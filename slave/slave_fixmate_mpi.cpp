#include <cstring>
#include <stdint.h>

#include "BamTools.h"

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

namespace {

const int kMinQuality = 15;

static inline uint16_t fm_le16(const uint8_t *p) {
    return (uint16_t)p[0] | ((uint16_t)p[1] << 8);
}

static inline uint32_t fm_le32(const uint8_t *p) {
    return (uint32_t)p[0] |
           ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) |
           ((uint32_t)p[3] << 24);
}

static inline void fm_store_le32(uint8_t *p, uint32_t value) {
    p[0] = (uint8_t)value;
    p[1] = (uint8_t)(value >> 8);
    p[2] = (uint8_t)(value >> 16);
    p[3] = (uint8_t)(value >> 24);
}

static const uint8_t *fm_aux_payload_end(const uint8_t *type,
                                         const uint8_t *end) {
    if (!type || type >= end) return nullptr;
    const uint8_t code = *type++;
    size_t bytes = 0;
    switch (code) {
    case 'A':
    case 'c':
    case 'C':
        bytes = 1;
        break;
    case 's':
    case 'S':
        bytes = 2;
        break;
    case 'i':
    case 'I':
    case 'f':
        bytes = 4;
        break;
    case 'd':
        bytes = 8;
        break;
    case 'Z':
    case 'H':
        while (type < end && *type != 0) ++type;
        return type < end ? type + 1 : nullptr;
    case 'B': {
        if ((size_t)(end - type) < 5) return nullptr;
        const uint8_t subtype = *type++;
        const uint32_t count = fm_le32(type);
        type += 4;
        size_t width = 0;
        switch (subtype) {
        case 'c':
        case 'C':
            width = 1;
            break;
        case 's':
        case 'S':
            width = 2;
            break;
        case 'i':
        case 'I':
        case 'f':
            width = 4;
            break;
        case 'd':
            width = 8;
            break;
        default:
            return nullptr;
        }
        if (count > SIZE_MAX / width) return nullptr;
        bytes = (size_t)count * width;
        break;
    }
    default:
        return nullptr;
    }
    return (size_t)(end - type) >= bytes ? type + bytes : nullptr;
}

static inline uint32_t fm_mate_score(const bam1_t *record) {
    const uint8_t *qual = bam_get_qual(record);
    uint32_t score = 0;
    for (int i = 0; i < record->core.l_qseq; ++i) {
        if (qual[i] >= kMinQuality) score += qual[i];
    }
    return score;
}

static inline uint32_t fm_decimal_digits(uint32_t value) {
    uint32_t digits = 1;
    while (value >= 10) {
        value /= 10;
        ++digits;
    }
    return digits;
}

static uint32_t fm_cigar_text_len(const bam1_t *record) {
    if (record->core.n_cigar == 0) return 1;
    const uint32_t *cigar = bam_get_cigar(record);
    uint32_t total = 0;
    for (uint32_t i = 0; i < record->core.n_cigar; ++i) {
        total += fm_decimal_digits(bam_cigar_oplen(cigar[i])) + 1;
    }
    return total;
}

static uint8_t *fm_write_decimal(uint8_t *dst, uint32_t value) {
    char buf[16];
    int n = 0;
    do {
        buf[n++] = (char)('0' + value % 10);
        value /= 10;
    } while (value != 0);
    while (n > 0) *dst++ = (uint8_t)buf[--n];
    return dst;
}

static uint8_t *fm_write_cigar(uint8_t *dst, const bam1_t *record) {
    if (record->core.n_cigar == 0) {
        *dst++ = '*';
        return dst;
    }
    const uint32_t *cigar = bam_get_cigar(record);
    for (uint32_t i = 0; i < record->core.n_cigar; ++i) {
        dst = fm_write_decimal(dst, bam_cigar_oplen(cigar[i]));
        *dst++ = (uint8_t)bam_cigar_opchr(cigar[i]);
    }
    return dst;
}

static inline void fm_sync_unmapped_pos(bam1_core_t *src,
                                        bam1_core_t *dest) {
    if ((dest->flag & BAM_FUNMAP) && !(src->flag & BAM_FUNMAP)) {
        dest->tid = src->tid;
        dest->pos = src->pos;
    }
}

static inline void fm_sync_mate_inner(const bam1_core_t *src,
                                      bam1_core_t *dest) {
    dest->mtid = src->tid;
    dest->mpos = src->pos;
    if (src->flag & BAM_FREVERSE) {
        dest->flag |= BAM_FMREVERSE;
    } else {
        dest->flag &= (uint16_t)~BAM_FMREVERSE;
    }
    if (src->flag & BAM_FUNMAP) dest->flag |= BAM_FMUNMAP;
}

static inline void fm_sync_mate(bam1_core_t *a, bam1_core_t *b) {
    fm_sync_unmapped_pos(a, b);
    fm_sync_unmapped_pos(b, a);
    fm_sync_mate_inner(a, b);
    fm_sync_mate_inner(b, a);
}

static inline int fm_plausibly_properly_paired(
        const bam1_core_t *a, const bam1_t *a_record,
        const bam1_core_t *b, const bam1_t *b_record) {
    if ((a->flag & BAM_FUNMAP) || (b->flag & BAM_FUNMAP)) return 0;
    if (a->tid != b->tid) return 0;
    int64_t a_pos = (a->flag & BAM_FREVERSE)
        ? bam_endpos(a_record) : a->pos;
    int64_t b_pos = (b->flag & BAM_FREVERSE)
        ? bam_endpos(b_record) : b->pos;
    const bam1_core_t *first = a;
    const bam1_core_t *second = b;
    if (a_pos > b_pos) {
        first = b;
        second = a;
    }
    return !(first->flag & BAM_FREVERSE) &&
           (second->flag & BAM_FREVERSE);
}

static int fm_plan_output_length(
        bam1_t *record, MpiFixmateRecordPlanShared *plan,
        bam1_t **records, int n_records) {
    const uint8_t *aux = bam_get_aux(record);
    const uint8_t *end = record->data + record->l_data;
    if (aux < record->data || aux > end) return -1;
    size_t output = (size_t)(aux - record->data);
    const uint8_t *p = aux;
    while (p < end) {
        if ((size_t)(end - p) < 3) return -1;
        const uint8_t *next = fm_aux_payload_end(p + 2, end);
        if (!next) return -1;
        const int drop =
            (plan->update_mq && p[0] == 'M' && p[1] == 'Q') ||
            (plan->update_mc && p[0] == 'M' && p[1] == 'C') ||
            (plan->update_ms && p[0] == 'm' && p[1] == 's');
        if (!drop) output += (size_t)(next - p);
        p = next;
    }
    if (plan->update_mq) output += 7;
    if (plan->update_mc) {
        if (plan->mate_index < 0 || plan->mate_index >= n_records) return -1;
        output += 4 + fm_cigar_text_len(records[plan->mate_index]);
    }
    if (plan->update_ms) output += 7;
    if (output > UINT32_MAX) return -1;
    plan->output_data_len = (uint32_t)output;
    return 0;
}

static void fm_pair_records(
        bam1_t **records, MpiFixmateRecordPlanShared *plans,
        int pre, int cur) {
    bam1_core_t *pre_core = &plans[pre].core;
    bam1_core_t *cur_core = &plans[cur].core;
    pre_core->flag |= BAM_FPAIRED;
    cur_core->flag |= BAM_FPAIRED;
    fm_sync_mate(pre_core, cur_core);

    if (pre_core->tid == cur_core->tid &&
        !(cur_core->flag & (BAM_FUNMAP | BAM_FMUNMAP)) &&
        !(pre_core->flag & (BAM_FUNMAP | BAM_FMUNMAP))) {
        const int64_t pre_end = bam_endpos(records[pre]);
        const int64_t cur_end = bam_endpos(records[cur]);
        const int64_t pre5 = (pre_core->flag & BAM_FREVERSE)
            ? pre_end : pre_core->pos;
        const int64_t cur5 = (cur_core->flag & BAM_FREVERSE)
            ? cur_end : cur_core->pos;
        cur_core->isize = pre5 - cur5;
        pre_core->isize = cur5 - pre5;
    } else {
        pre_core->isize = 0;
        cur_core->isize = 0;
    }

    if (!fm_plausibly_properly_paired(
            pre_core, records[pre], cur_core, records[cur])) {
        pre_core->flag &= (uint16_t)~BAM_FPROPER_PAIR;
        cur_core->flag &= (uint16_t)~BAM_FPROPER_PAIR;
    }

    plans[pre].mate_index = cur;
    plans[cur].mate_index = pre;
    plans[pre].rewrite_tags = 1;
    plans[cur].rewrite_tags = 1;
    plans[pre].update_mq =
        !(plans[cur].core.flag & BAM_FUNMAP);
    plans[cur].update_mq =
        !(plans[pre].core.flag & BAM_FUNMAP);
    plans[pre].update_mc =
        !(plans[cur].core.flag & BAM_FUNMAP) ||
        !(plans[pre].core.flag & BAM_FUNMAP);
    plans[cur].update_mc = plans[pre].update_mc;
    plans[pre].update_ms = 1;
    plans[cur].update_ms = 1;
}

static int fm_copy_aux_and_append(
        bam1_t *src, bam1_t *dst,
        const MpiFixmateRecordPlanShared *plan,
        bam1_t **records) {
    const uint8_t *src_aux = bam_get_aux(src);
    const uint8_t *src_end = src->data + src->l_data;
    const size_t prefix = (size_t)(src_aux - src->data);
    if (prefix > (size_t)src->l_data ||
        prefix > plan->output_data_len) {
        return -1;
    }
    memcpy(dst->data, src->data, prefix);
    uint8_t *out = dst->data + prefix;
    const uint8_t *p = src_aux;
    while (p < src_end) {
        if ((size_t)(src_end - p) < 3) return -1;
        const uint8_t *next = fm_aux_payload_end(p + 2, src_end);
        if (!next) return -1;
        const int drop =
            (plan->update_mq && p[0] == 'M' && p[1] == 'Q') ||
            (plan->update_mc && p[0] == 'M' && p[1] == 'C') ||
            (plan->update_ms && p[0] == 'm' && p[1] == 's');
        if (!drop) {
            memcpy(out, p, (size_t)(next - p));
            out += next - p;
        }
        p = next;
    }

    if (plan->mate_index < 0 &&
        (plan->update_mq || plan->update_mc || plan->update_ms)) {
        return -1;
    }
    bam1_t *mate =
        plan->mate_index >= 0 ? records[plan->mate_index] : nullptr;
    if (plan->update_mq) {
        *out++ = 'M';
        *out++ = 'Q';
        *out++ = 'i';
        fm_store_le32(out, mate->core.qual);
        out += 4;
    }
    if (plan->update_mc) {
        *out++ = 'M';
        *out++ = 'C';
        *out++ = 'Z';
        out = fm_write_cigar(out, mate);
        *out++ = 0;
    }
    if (plan->update_ms) {
        *out++ = 'm';
        *out++ = 's';
        *out++ = 'i';
        fm_store_le32(out, fm_mate_score(mate));
        out += 4;
    }
    return out == dst->data + plan->output_data_len ? 0 : -1;
}

} // namespace

extern "C" void slave_mpi_fixmate_plan(MpiFixmatePlanPara paras[64]) {
    const int id = _PEN;
    MpiFixmatePlanPara *para = &paras[id];
    if (para->status != 0 || !para->records || !para->plans ||
        !para->groups) {
        return;
    }

    para->paired_groups = 0;
    para->singleton_groups = 0;
    para->secondary_records = 0;
    para->supplementary_records = 0;
    para->mq_updates = 0;
    para->mc_updates = 0;
    para->ms_updates = 0;

    for (int g = para->group_begin; g < para->group_end; ++g) {
        const int begin = para->groups[g].begin;
        const int end = begin + para->groups[g].count;
        int pre = -1;
        int cur = -1;

        for (int i = begin; i < end; ++i) {
            bam1_t *record = para->records[i];
            if (!record) {
                para->status = -2;
                para->record_index = i;
                return;
            }
            MpiFixmateRecordPlanShared *plan = &para->plans[i];
            memset(plan, 0, sizeof(*plan));
            plan->core = record->core;
            plan->output_data_len = (uint32_t)record->l_data;
            plan->mate_index = -1;
            if (record->core.flag & BAM_FSECONDARY) {
                para->secondary_records++;
                continue;
            }
            if (record->core.flag & BAM_FSUPPLEMENTARY) {
                para->supplementary_records++;
                continue;
            }
            if (pre < 0) {
                pre = i;
                continue;
            }
            cur = i;
            fm_pair_records(para->records, para->plans, pre, cur);
        }

        if (cur < 0 && pre >= 0) {
            MpiFixmateRecordPlanShared *plan = &para->plans[pre];
            plan->core.mtid = -1;
            plan->core.mpos = -1;
            plan->core.isize = 0;
            plan->core.flag &=
                (uint16_t)~(BAM_FMREVERSE | BAM_FPROPER_PAIR);
            para->singleton_groups++;
        } else if (cur >= 0) {
            para->paired_groups++;
        }

        for (int i = begin; i < end; ++i) {
            if (fm_plan_output_length(
                    para->records[i], &para->plans[i],
                    para->records, end) != 0) {
                para->status = -3;
                para->record_index = i;
                return;
            }
            para->mq_updates += para->plans[i].update_mq;
            para->mc_updates += para->plans[i].update_mc;
            para->ms_updates += para->plans[i].update_ms;
        }
    }
    para->status = 0;
}

extern "C" void slave_mpi_fixmate_rewrite(
        MpiFixmateRewritePara paras[64]) {
    const int id = _PEN;
    MpiFixmateRewritePara *para = &paras[id];
    if (para->status != 0 || !para->records ||
        !para->output_records || !para->plans ||
        !para->output_offsets || !para->output_data) {
        return;
    }

    for (int i = para->record_begin; i < para->record_end; ++i) {
        const uint64_t offset = para->output_offsets[i];
        const uint64_t next = para->output_offsets[i + 1];
        const MpiFixmateRecordPlanShared *plan = &para->plans[i];
        if (next < offset ||
            next - offset != plan->output_data_len ||
            next > para->output_capacity) {
            para->status = -3;
            para->record_index = i;
            return;
        }
        bam1_t *src = para->records[i];
        bam1_t *dst = para->output_records + i;
        memset(dst, 0, sizeof(*dst));
        dst->core = plan->core;
        dst->id = src->id;
        dst->data = para->output_data + offset;
        dst->l_data = (int)plan->output_data_len;
        dst->m_data = plan->output_data_len;
        dst->mempolicy = BAM_USER_OWNS_DATA;
        if (fm_copy_aux_and_append(
                src, dst, plan, para->records) != 0) {
            para->status = -2;
            para->record_index = i;
            return;
        }
    }
    para->status = 0;
}
