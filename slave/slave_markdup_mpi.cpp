#include <cstring>
#include <stdint.h>
#include <stdlib.h>

#include "BamTools.h"

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

namespace {

const int kOrientationFF = 2;
const int kOrientationRR = 3;
const int kOrientationFR = 5;
const int kOrientationRF = 7;
const int kLeft = 11;
const int kRight = 13;
const int kMinQuality = 15;

static inline uint16_t md_le16(const uint8_t *p) {
    return (uint16_t)p[0] | ((uint16_t)p[1] << 8);
}

static inline uint32_t md_le32(const uint8_t *p) {
    return (uint32_t)p[0] |
           ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) |
           ((uint32_t)p[3] << 24);
}

static inline const uint8_t *md_cigar_bytes(const bam1_t *record) {
    return record->data + record->core.l_qname;
}

static inline uint32_t md_cigar_word(const bam1_t *record, uint32_t index) {
    return md_le32(md_cigar_bytes(record) + (size_t)index * 4);
}

static inline uint64_t md_le64(const uint8_t *p) {
    return (uint64_t)md_le32(p) | ((uint64_t)md_le32(p + 4) << 32);
}

static const uint8_t *md_aux_payload_end(const uint8_t *type,
                                         const uint8_t *end) {
    if (type >= end) return nullptr;
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
        const uint32_t count = md_le32(type);
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

static uint8_t *md_aux_find(bam1_t *record, char a, char b) {
    uint8_t *p = bam_get_aux(record);
    uint8_t *end = record->data + record->l_data;
    while ((size_t)(end - p) >= 3) {
        uint8_t *type = p + 2;
        const uint8_t *next = md_aux_payload_end(type, end);
        if (!next) return nullptr;
        if (p[0] == (uint8_t)a && p[1] == (uint8_t)b) return type;
        p = (uint8_t *)next;
    }
    return nullptr;
}

static int md_aux_to_i64(const uint8_t *type, int64_t *value) {
    if (!type || !value) return -1;
    switch (*type) {
    case 'c':
        *value = (int8_t)type[1];
        return 0;
    case 'C':
        *value = type[1];
        return 0;
    case 's':
        *value = (int16_t)md_le16(type + 1);
        return 0;
    case 'S':
        *value = md_le16(type + 1);
        return 0;
    case 'i':
        *value = (int32_t)md_le32(type + 1);
        return 0;
    case 'I':
        *value = md_le32(type + 1);
        return 0;
    default:
        return -1;
    }
}

static inline int64_t md_current_score(const bam1_t *record) {
    const uint8_t *qual = bam_get_qual(record);
    int64_t score = 0;
    for (int i = 0; i < record->core.l_qseq; ++i) {
        if (qual[i] >= kMinQuality) score += qual[i];
    }
    return score;
}

static inline int64_t md_unclipped_start(const bam1_t *record) {
    int64_t clipped = 0;
    for (uint32_t i = 0; i < record->core.n_cigar; ++i) {
        const uint32_t cigar = md_cigar_word(record, i);
        const int op = bam_cigar_op(cigar);
        if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            clipped += bam_cigar_oplen(cigar);
        } else if (op != BAM_CHARD_CLIP) {
            break;
        }
    }
    return (int64_t)record->core.pos - clipped + 1;
}

static inline int64_t md_unclipped_end(const bam1_t *record) {
    int64_t ref_len = 0;
    for (uint32_t i = 0; i < record->core.n_cigar; ++i) {
        const uint32_t cigar = md_cigar_word(record, i);
        const int op = bam_cigar_op(cigar);
        if (bam_cigar_type(op) & 2) ref_len += bam_cigar_oplen(cigar);
    }
    if (ref_len == 0) ref_len = 1;
    int64_t clipped = 0;
    for (int i = (int)record->core.n_cigar - 1; i >= 0; --i) {
        const uint32_t cigar = md_cigar_word(record, (uint32_t)i);
        const int op = bam_cigar_op(cigar);
        if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            clipped += bam_cigar_oplen(cigar);
        } else if (op != BAM_CHARD_CLIP) {
            break;
        }
    }
    return (int64_t)record->core.pos + ref_len + clipped;
}

static inline int md_cigar_digit(char c) {
    return c >= '0' && c <= '9';
}

static int64_t md_other_start(int64_t pos, const char *cigar) {
    const char *p = cigar;
    int64_t clipped = 0;
    while (*p && *p != '*') {
        long len = 0;
        while (md_cigar_digit(*p)) {
            len = len * 10 + (*p - '0');
            ++p;
        }
        if (len == 0) len = 1;
        if (*p == 'S' || *p == 'H') {
            clipped += len;
        } else if (*p != 'H') {
            break;
        }
        if (*p) ++p;
    }
    return pos - clipped + 1;
}

static int64_t md_other_end(int64_t pos, const char *cigar) {
    const char *p = cigar;
    int64_t ref_pos = 0;
    int skip = 1;
    while (*p && *p != '*') {
        long len = 0;
        while (md_cigar_digit(*p)) {
            len = len * 10 + (*p - '0');
            ++p;
        }
        if (len == 0) len = 1;
        switch (*p) {
        case 'M':
        case 'D':
        case 'N':
        case '=':
        case 'X':
            ref_pos += len;
            skip = 0;
            break;
        case 'S':
        case 'H':
            if (!skip) ref_pos += len;
            break;
        default:
            break;
        }
        if (*p) ++p;
    }
    return pos + ref_pos;
}

static inline int md_has_mate(const bam1_t *record) {
    return (record->core.flag & BAM_FPAIRED) &&
           !(record->core.flag & BAM_FMUNMAP) &&
           !(record->core.mtid == -1 && record->core.mpos == -1);
}

static void md_make_single_key(const bam1_t *record,
                               MpiMarkdupKeyShared *key) {
    memset(key, 0, sizeof(*key));
    key->single = 1;
    key->this_ref = record->core.tid + 1;
    if (record->core.flag & BAM_FREVERSE) {
        key->this_coord = md_unclipped_end(record);
        key->orientation = kOrientationRR;
    } else {
        key->this_coord = md_unclipped_start(record);
        key->orientation = kOrientationFF;
    }
}

static int md_make_pair_key(bam1_t *record,
                            const char *mate_cigar,
                            MpiMarkdupKeyShared *key) {
    const int32_t this_ref = record->core.tid + 1;
    const int32_t other_ref = record->core.mtid + 1;
    int64_t this_coord = md_unclipped_start(record);
    const int64_t this_end = md_unclipped_end(record);
    int64_t other_coord = md_other_start(record->core.mpos, mate_cigar);
    const int64_t other_end = md_other_end(record->core.mpos, mate_cigar);
    const int this_reverse = (record->core.flag & BAM_FREVERSE) != 0;
    const int mate_reverse = (record->core.flag & BAM_FMREVERSE) != 0;
    int leftmost = 0;

    if (this_ref != other_ref) {
        leftmost = this_ref < other_ref;
    } else if (this_reverse == mate_reverse) {
        leftmost = !this_reverse ? this_coord <= other_coord
                                 : this_end <= other_end;
    } else {
        leftmost = this_reverse ? this_end <= other_coord
                                : this_coord <= other_end;
    }

    int orientation = 0;
    if (leftmost) {
        if (this_reverse == mate_reverse) {
            other_coord = other_end;
            if (!this_reverse) {
                orientation = (record->core.flag & BAM_FREAD1)
                    ? kOrientationFF : kOrientationRR;
            } else {
                orientation = (record->core.flag & BAM_FREAD1)
                    ? kOrientationRR : kOrientationFF;
            }
        } else if (!this_reverse) {
            orientation = kOrientationFR;
            other_coord = other_end;
        } else {
            orientation = kOrientationRF;
            this_coord = this_end;
        }
    } else {
        if (this_reverse == mate_reverse) {
            this_coord = this_end;
            if (!this_reverse) {
                orientation = (record->core.flag & BAM_FREAD1)
                    ? kOrientationRR : kOrientationFF;
            } else {
                orientation = (record->core.flag & BAM_FREAD1)
                    ? kOrientationFF : kOrientationRR;
            }
        } else if (!this_reverse) {
            orientation = kOrientationRF;
            other_coord = other_end;
        } else {
            orientation = kOrientationFR;
            this_coord = this_end;
        }
    }

    memset(key, 0, sizeof(*key));
    key->single = 0;
    key->this_ref = this_ref;
    key->other_ref = other_ref;
    key->this_coord = this_coord;
    key->other_coord = other_coord;
    key->leftmost = leftmost ? kLeft : kRight;
    key->orientation = orientation;
    return 0;
}

static int md_append_candidate(MpiMarkdupExtractPara *para,
                               const MpiMarkdupCandidateShared &candidate) {
    if (para->n_candidates >= para->candidate_capacity) {
        para->status = -3;
        para->actual_value = para->n_candidates + 1;
        para->limit_value = para->candidate_capacity;
        para->limit_id = BOUNDS_LIMIT_MAX_RECORDS_PER_BLOCK;
        return -1;
    }
    para->candidates[para->n_candidates++] = candidate;
    return 0;
}

static int md_delete_aux_tag(bam1_t *record, char a, char b) {
    uint8_t *p = bam_get_aux(record);
    uint8_t *end = record->data + record->l_data;
    while ((size_t)(end - p) >= 3) {
        uint8_t *type = p + 2;
        const uint8_t *next = md_aux_payload_end(type, end);
        if (!next) return -1;
        if (p[0] == (uint8_t)a && p[1] == (uint8_t)b) {
            size_t remove_len = (size_t)(next - p);
            memmove(p, next, (size_t)(end - next));
            record->l_data -= (int)remove_len;
            return 1;
        }
        p = (uint8_t *)next;
    }
    return 0;
}

} // namespace

extern "C" void slave_mpi_markdup_extract(MpiMarkdupExtractPara paras[64]) {
    const int id = _PEN;
    MpiMarkdupExtractPara *para = &paras[id];
    if (para->status != 0 || !para->records || !para->candidates ||
        !para->qname_arena) {
        return;
    }

    para->n_candidates = 0;
    para->qname_used = 0;
    para->examined = 0;
    para->excluded = 0;
    para->pair_candidates = 0;
    para->single_candidates = 0;
    para->has_records = para->n_records > 0;
    para->first_tid = para->last_tid = -1;
    para->first_pos = para->last_pos = -1;

    int prev_tid = -1;
    int prev_pos = -1;
    for (int i = 0; i < para->n_records; ++i) {
        bam1_t *record = para->records[i];
        if (!record) {
            para->status = -2;
            para->record_index = i;
            return;
        }
        if (i == 0) {
            para->first_tid = record->core.tid;
            para->first_pos = record->core.pos;
        }
        para->last_tid = record->core.tid;
        para->last_pos = record->core.pos;
        if (record->core.tid >= 0 &&
            (record->core.tid < prev_tid ||
             (record->core.tid == prev_tid && record->core.pos < prev_pos))) {
            para->status = -20;
            para->record_index = i;
            return;
        }
        prev_tid = record->core.tid;
        prev_pos = record->core.pos;

        const uint16_t flag = record->core.flag;
        uint16_t excluded = BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP;
        if (!para->include_fails) excluded |= BAM_FQCFAIL;
        if (flag & excluded) {
            para->excluded++;
            continue;
        }
        para->examined++;

        const uint64_t ordinal = para->ordinal_base + (uint64_t)i;
        const uint64_t global_order =
            (para->global_block_index << 32) | (uint32_t)i;
        const int64_t current_score = md_current_score(record);

        if (md_has_mate(record)) {
            uint8_t *mc = md_aux_find(record, 'M', 'C');
            if (!mc) {
                para->status = -10;
                para->record_index = i;
                return;
            }
            if (*mc != 'Z') {
                para->status = -11;
                para->record_index = i;
                return;
            }
            uint8_t *ms = md_aux_find(record, 'm', 's');
            if (!ms) {
                para->status = -12;
                para->record_index = i;
                return;
            }
            int64_t mate_score = 0;
            if (md_aux_to_i64(ms, &mate_score) != 0) {
                para->status = -13;
                para->record_index = i;
                return;
            }

            const char *qname = bam_get_qname(record);
            size_t qname_len = strlen(qname);
            if (qname_len > UINT16_MAX ||
                qname_len > para->qname_capacity - para->qname_used) {
                para->status = -3;
                para->record_index = i;
                para->actual_value = (long long)qname_len;
                para->limit_value =
                    (long long)(para->qname_capacity - para->qname_used);
                para->limit_id = BOUNDS_LIMIT_INIT_DATA_SIZE;
                return;
            }

            MpiMarkdupCandidateShared pair_candidate;
            memset(&pair_candidate, 0, sizeof(pair_candidate));
            md_make_pair_key(record, (const char *)(mc + 1),
                             &pair_candidate.key);
            pair_candidate.ordinal = ordinal;
            pair_candidate.global_order = global_order;
            pair_candidate.score = current_score + mate_score;
            pair_candidate.qname_offset = (uint32_t)para->qname_used;
            pair_candidate.qname_len = (uint16_t)qname_len;
            pair_candidate.paired_marker = 1;
            pair_candidate.qc_fail = (flag & BAM_FQCFAIL) != 0;
            pair_candidate.source_rank = para->source_rank;
            memcpy(para->qname_arena + para->qname_used, qname, qname_len);
            para->qname_used += qname_len;
            if (md_append_candidate(para, pair_candidate) != 0) return;
            para->pair_candidates++;

            MpiMarkdupCandidateShared single_marker;
            memset(&single_marker, 0, sizeof(single_marker));
            md_make_single_key(record, &single_marker.key);
            single_marker.ordinal = ordinal;
            single_marker.global_order = global_order;
            single_marker.score = current_score;
            single_marker.paired_marker = 1;
            single_marker.qc_fail = (flag & BAM_FQCFAIL) != 0;
            single_marker.source_rank = para->source_rank;
            if (md_append_candidate(para, single_marker) != 0) return;
            para->single_candidates++;
        } else {
            MpiMarkdupCandidateShared single_candidate;
            memset(&single_candidate, 0, sizeof(single_candidate));
            md_make_single_key(record, &single_candidate.key);
            single_candidate.ordinal = ordinal;
            single_candidate.global_order = global_order;
            single_candidate.score = current_score;
            single_candidate.paired_marker = 0;
            single_candidate.qc_fail = (flag & BAM_FQCFAIL) != 0;
            single_candidate.source_rank = para->source_rank;
            if (md_append_candidate(para, single_candidate) != 0) return;
            para->single_candidates++;
        }
    }
    para->status = 0;
}

extern "C" void slave_mpi_markdup_rewrite(MpiMarkdupRewritePara paras[64]) {
    const int id = _PEN;
    MpiMarkdupRewritePara *para = &paras[id];
    if (para->status != 0 || !para->records || !para->bam_lens) return;

    int kept = 0;
    uint32_t kept_len = 0;
    para->marked_records = 0;
    para->cleared_records = 0;
    para->removed_records = 0;

    for (int i = 0; i < para->n_records; ++i) {
        bam1_t *record = para->records[i];
        if (!record) {
            para->status = -2;
            para->record_index = i;
            return;
        }
        const uint64_t ordinal = para->ordinal_base + (uint64_t)i;
        const size_t byte = (size_t)(ordinal >> 3);
        const uint8_t mask = (uint8_t)(1u << (ordinal & 7u));
        const int duplicate = byte < para->duplicate_bitmap_bytes &&
                              (para->duplicate_bitmap[byte] & mask);

        if (para->clear_old) {
            if (record->core.flag & BAM_FDUP) para->cleared_records++;
            record->core.flag &= (uint16_t)~BAM_FDUP;
            int ret = md_delete_aux_tag(record, 'd', 't');
            if (ret < 0) {
                para->status = -2;
                para->record_index = i;
                return;
            }
            ret = md_delete_aux_tag(record, 'd', 'o');
            if (ret < 0) {
                para->status = -2;
                para->record_index = i;
                return;
            }
        }
        if (duplicate) {
            record->core.flag |= BAM_FDUP;
            para->marked_records++;
        }
        if (para->remove_dups && (record->core.flag & BAM_FDUP)) {
            para->removed_records++;
            continue;
        }

        const uint32_t bam_len =
            (uint32_t)(record->l_data - record->core.l_extranul + 32);
        para->records[kept] = record;
        para->bam_lens[kept] = bam_len;
        kept_len += bam_len + 4;
        kept++;
    }

    para->n_kept_records = kept;
    para->kept_total_len = kept_len;
    para->status = 0;
}
