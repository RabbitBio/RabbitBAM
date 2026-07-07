#include "BamTools.h"
#include "swbam/cpe_bam_parser.h"
#include "swbam/cpe_codec.h"

#include <climits>
#include <cstring>

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

#define slave_cycle_now swbam_cpe_cycle_now

static inline int rb_slave_stats_unclipped_length(const bam1_t *record) {
    int read_len = record->core.l_qseq;
    const uint32_t *cigar = bam_get_cigar(record);
    for (int i = 0; i < record->core.n_cigar; ++i) {
        if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
            read_len += bam_cigar_oplen(cigar[i]);
        }
    }
    return read_len;
}

static inline long long rb_slave_stats_mapped_cigar_bases(const bam1_t *record) {
    long long total = 0;
    const uint32_t *cigar = bam_get_cigar(record);
    for (int i = 0; i < record->core.n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CEQUAL || op == BAM_CDIFF) {
            total += bam_cigar_oplen(cigar[i]);
        }
    }
    return total;
}

static inline uint32_t rb_slave_stats_read_order(const bam1_t *record) {
    const uint16_t flag = record->core.flag;
    if ((flag & BAM_FPAIRED) == 0) return 1;
    return ((flag & BAM_FREAD1) ? 1u : 0u) + ((flag & BAM_FREAD2) ? 2u : 0u);
}

static inline void rb_slave_stats_update_sort(const bam1_t *record,
                                              MpiStatsBlockSortState *sort_state) {
    const long long tid = record->core.tid;
    const long long pos = record->core.pos;
    if (!sort_state->has_coord) {
        sort_state->has_coord = 1;
        sort_state->first_tid = tid;
        sort_state->first_pos = pos;
        sort_state->last_tid = tid;
        sort_state->last_pos = pos;
        return;
    }
    if (tid < sort_state->last_tid ||
        (tid == sort_state->last_tid && pos < sort_state->last_pos)) {
        sort_state->sorted = 0;
    }
    sort_state->last_tid = tid;
    sort_state->last_pos = pos;
}

enum RbSlaveStatsOrientClass {
    RB_SLAVE_STATS_ORIENT_IN = 0,
    RB_SLAVE_STATS_ORIENT_OUT = 1,
    RB_SLAVE_STATS_ORIENT_OTHER = 2
};

static inline RbSlaveStatsOrientClass rb_slave_stats_reference_orient_class(const bam1_t *record) {
    const uint16_t flag = record->core.flag;
    const long long pos_fst = record->core.mpos - record->core.pos;
    const int is_fst = (flag & BAM_FREAD1) ? 1 : -1;
    const int is_fwd = (flag & BAM_FREVERSE) ? -1 : 1;
    const int is_mfwd = (flag & BAM_FMREVERSE) ? -1 : 1;

    if (is_fwd * is_mfwd > 0) return RB_SLAVE_STATS_ORIENT_OTHER;
    if (is_fst * pos_fst >= 0) {
        return (is_fst * is_fwd > 0) ? RB_SLAVE_STATS_ORIENT_IN : RB_SLAVE_STATS_ORIENT_OUT;
    }
    return (is_fst * is_fwd > 0) ? RB_SLAVE_STATS_ORIENT_OUT : RB_SLAVE_STATS_ORIENT_IN;
}

static inline RbSlaveStatsOrientClass rb_slave_stats_read_strand_orient_class(const bam1_t *record) {
    const uint16_t flag = record->core.flag;
    const int is_fst = (flag & BAM_FREAD1) ? 1 : -1;
    const int is_fwd = (flag & BAM_FREVERSE) ? -1 : 1;
    const int is_mfwd = (flag & BAM_FMREVERSE) ? -1 : 1;

    if (is_fwd * is_mfwd > 0) return RB_SLAVE_STATS_ORIENT_OTHER;
    return (is_fst * is_fwd > 0) ? RB_SLAVE_STATS_ORIENT_IN : RB_SLAVE_STATS_ORIENT_OUT;
}

static inline void rb_slave_stats_add_insert_diag(const bam1_t *record,
                                                  int isize,
                                                  long long signed_isize,
                                                  MpiStatsBasicCountSlice *counts) {
    const uint16_t flag = record->core.flag;
    const RbSlaveStatsOrientClass ref_class = rb_slave_stats_reference_orient_class(record);
    const RbSlaveStatsOrientClass doc_class = rb_slave_stats_read_strand_orient_class(record);

    if (ref_class == RB_SLAVE_STATS_ORIENT_IN) counts->orient_diag[RB_ORIENT_DIAG_REF_IN]++;
    else if (ref_class == RB_SLAVE_STATS_ORIENT_OUT) counts->orient_diag[RB_ORIENT_DIAG_REF_OUT]++;
    else counts->orient_diag[RB_ORIENT_DIAG_REF_OTHER]++;

    if (doc_class == RB_SLAVE_STATS_ORIENT_IN) counts->orient_diag[RB_ORIENT_DIAG_DOC_IN]++;
    else if (doc_class == RB_SLAVE_STATS_ORIENT_OUT) counts->orient_diag[RB_ORIENT_DIAG_DOC_OUT]++;
    else counts->orient_diag[RB_ORIENT_DIAG_DOC_OTHER]++;

    const int read_len = record->core.l_qseq > 0 ? record->core.l_qseq : 1;
    int size_bucket = 2;
    if (isize < read_len) size_bucket = 0;
    else if (isize < 2 * read_len) size_bucket = 1;

    if (ref_class == RB_SLAVE_STATS_ORIENT_IN) {
        counts->orient_diag[RB_ORIENT_DIAG_REF_IN_LT_READ + size_bucket]++;
    } else if (ref_class == RB_SLAVE_STATS_ORIENT_OUT) {
        counts->orient_diag[RB_ORIENT_DIAG_REF_OUT_LT_READ + size_bucket]++;
    }

    const long long pos_fst = record->core.mpos - record->core.pos;
    const int is_fst = (flag & BAM_FREAD1) ? 1 : -1;
    const int is_fwd = (flag & BAM_FREVERSE) ? -1 : 1;
    const long long pos_term = is_fst * pos_fst;
    const int strand_term = is_fst * is_fwd;
    int sign_base = RB_ORIENT_DIAG_POS_ZERO_STRAND_NEG;
    if (pos_term < 0) sign_base = RB_ORIENT_DIAG_POS_NEG_STRAND_NEG;
    else if (pos_term > 0) sign_base = RB_ORIENT_DIAG_POS_POS_STRAND_NEG;
    counts->orient_diag[sign_base + (strand_term > 0 ? 1 : 0)]++;

    if (flag & BAM_FREAD1) {
        if (pos_fst > 0) counts->orient_diag[RB_ORIENT_DIAG_READ1_LEFT]++;
        else if (pos_fst < 0) counts->orient_diag[RB_ORIENT_DIAG_READ1_RIGHT]++;
        else counts->orient_diag[RB_ORIENT_DIAG_READ1_SAME_POS]++;
    } else if (flag & BAM_FREAD2) {
        if (pos_fst > 0) counts->orient_diag[RB_ORIENT_DIAG_READ2_LEFT]++;
        else if (pos_fst < 0) counts->orient_diag[RB_ORIENT_DIAG_READ2_RIGHT]++;
        else counts->orient_diag[RB_ORIENT_DIAG_READ2_SAME_POS]++;
    }

    if (signed_isize < 0) counts->orient_diag[RB_ORIENT_DIAG_ISIZE_NEG]++;
    else if (signed_isize > 0) counts->orient_diag[RB_ORIENT_DIAG_ISIZE_POS]++;
    else counts->orient_diag[RB_ORIENT_DIAG_ISIZE_ZERO]++;
}

static inline void rb_slave_stats_add_insert(const bam1_t *record,
                                             MpiStatsBasicCountSlice *counts,
                                             int collect_diag) {
    const uint16_t flag = record->core.flag;
    if ((flag & BAM_FPAIRED) == 0 || (flag & BAM_FUNMAP) || (flag & BAM_FMUNMAP)) return;
    if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) return;

    const long long signed_isize = record->core.isize;
    long long isize_ll = signed_isize;
    if (isize_ll < 0) isize_ll = -isize_ll;
    if (isize_ll > RB_MPI_STATS_MAX_INSERT_SIZE) isize_ll = RB_MPI_STATS_MAX_INSERT_SIZE;
    const int isize = (int)isize_ll;
    if (isize == 0 && record->core.tid != record->core.mtid) return;

    if (collect_diag) rb_slave_stats_add_insert_diag(record, isize, signed_isize, counts);

    const RbSlaveStatsOrientClass orient_class = rb_slave_stats_reference_orient_class(record);
    if (orient_class == RB_SLAVE_STATS_ORIENT_OTHER) {
        counts->isize_other[isize]++;
    } else if (orient_class == RB_SLAVE_STATS_ORIENT_OUT) {
        counts->isize_outward[isize]++;
    } else {
        counts->isize_inward[isize]++;
    }
}

static inline int rb_slave_aux_type_size(int type) {
    switch (type) {
    case 'A':
    case 'c':
    case 'C':
        return 1;
    case 's':
    case 'S':
        return 2;
    case 'i':
    case 'I':
    case 'f':
        return 4;
    case 'd':
        return 8;
    default:
        return 0;
    }
}

static inline const uint8_t *rb_slave_aux_get_type_ptr(const bam1_t *record,
                                                       char tag0,
                                                       char tag1) {
    const int l_qseq = record->core.l_qseq;
    size_t off = (size_t)record->core.l_qname +
                 ((size_t)record->core.n_cigar << 2) +
                 (((size_t)l_qseq + 1u) >> 1) +
                 (size_t)l_qseq;
    const uint8_t *p = record->data + off;
    const uint8_t *end = record->data + record->l_data;
    while (p + 3 <= end) {
        if (p[0] == (uint8_t)tag0 && p[1] == (uint8_t)tag1) return p + 2;
        const int type = p[2];
        p += 3;
        if (type == 'Z' || type == 'H') {
            while (p < end && *p) p++;
            if (p >= end) return NULL;
            p++;
        } else if (type == 'B') {
            if (p + 5 > end) return NULL;
            const int subtype = p[0];
            uint32_t n = 0;
            memcpy(&n, p + 1, 4);
            const int elem_size = rb_slave_aux_type_size(subtype);
            if (elem_size <= 0) return NULL;
            const uint64_t skip = 5ull + (uint64_t)n * (uint64_t)elem_size;
            if ((uint64_t)(end - p) < skip) return NULL;
            p += skip;
        } else {
            const int elem_size = rb_slave_aux_type_size(type);
            if (elem_size <= 0 || p + elem_size > end) return NULL;
            p += elem_size;
        }
    }
    return NULL;
}

static inline long long rb_slave_aux2i(const uint8_t *type_ptr) {
    if (!type_ptr) return 0;
    const uint8_t *p = type_ptr + 1;
    switch (*type_ptr) {
    case 'c': {
        int8_t v = 0;
        memcpy(&v, p, 1);
        return v;
    }
    case 'C': {
        uint8_t v = 0;
        memcpy(&v, p, 1);
        return v;
    }
    case 's': {
        int16_t v = 0;
        memcpy(&v, p, 2);
        return v;
    }
    case 'S': {
        uint16_t v = 0;
        memcpy(&v, p, 2);
        return v;
    }
    case 'i': {
        int32_t v = 0;
        memcpy(&v, p, 4);
        return v;
    }
    case 'I': {
        uint32_t v = 0;
        memcpy(&v, p, 4);
        return v;
    }
    default:
        return 0;
    }
}

static inline void rb_slave_stats_add_record(bam1_t *record,
                                             MpiStatsBasicCountSlice *counts,
                                             MpiStatsBlockSortState *sort_state,
                                             int collect_diag) {
    const uint16_t flag = record->core.flag;
    const int is_secondary = (flag & BAM_FSECONDARY) != 0;
    const int is_supplementary = (flag & BAM_FSUPPLEMENTARY) != 0;
    const int is_original = !is_secondary && !is_supplementary;

    if (is_secondary) {
        counts->values[RB_STATS_NREADS_SECONDARY]++;
        return;
    }
    if (is_supplementary) {
        counts->values[RB_STATS_NREADS_SUPPLEMENTARY]++;
    }

    const int seq_len = record->core.l_qseq;
    if (seq_len == 0) return;

    if (flag & BAM_FDUP) {
        counts->values[RB_STATS_TOTAL_LEN_DUP] += seq_len;
        counts->values[RB_STATS_NREADS_DUP]++;
    }

    const uint32_t order = rb_slave_stats_read_order(record);
    const int read_len = rb_slave_stats_unclipped_length(record);
    if (read_len > counts->values[RB_STATS_MAX_LEN]) counts->values[RB_STATS_MAX_LEN] = read_len;
    if (order == 1 && read_len > counts->values[RB_STATS_MAX_LEN_1ST]) {
        counts->values[RB_STATS_MAX_LEN_1ST] = read_len;
    }
    if (order == 2 && read_len > counts->values[RB_STATS_MAX_LEN_2ND]) {
        counts->values[RB_STATS_MAX_LEN_2ND] = read_len;
    }

    if (is_original) {
        counts->values[RB_STATS_TOTAL_LEN] += seq_len;
        if (flag & BAM_FQCFAIL) counts->values[RB_STATS_NREADS_QCFAILED]++;
        if (flag & BAM_FPAIRED) counts->values[RB_STATS_NREADS_PAIRED_TECH]++;

        if (order == 1) {
            counts->values[RB_STATS_NREADS_1ST]++;
            counts->values[RB_STATS_TOTAL_LEN_1ST] += seq_len;
        } else if (order == 2) {
            counts->values[RB_STATS_NREADS_2ND]++;
            counts->values[RB_STATS_TOTAL_LEN_2ND] += seq_len;
        } else {
            counts->values[RB_STATS_NREADS_OTHER]++;
        }

        if (order == 1 || order == 2) {
            uint8_t *qual = bam_get_qual(record);
            for (int i = 0; i < seq_len; ++i) {
                counts->values[RB_STATS_SUM_QUAL] += qual[i];
            }
        }

        if (flag & BAM_FUNMAP) {
            counts->values[RB_STATS_NREADS_UNMAPPED]++;
        } else {
            counts->values[RB_STATS_NBASES_MAPPED] += seq_len;
            if (record->core.qual == 0) counts->values[RB_STATS_NREADS_MQ0]++;
            if ((flag & BAM_FPAIRED) && !(flag & BAM_FMUNMAP)) {
                counts->values[RB_STATS_NREADS_PAIRED_AND_MAPPED]++;
                if ((flag & (BAM_FPAIRED | BAM_FPROPER_PAIR)) == (BAM_FPAIRED | BAM_FPROPER_PAIR)) {
                    counts->values[RB_STATS_NREADS_PROPERLY_PAIRED]++;
                }
                if (record->core.tid != record->core.mtid) counts->values[RB_STATS_NREADS_ANOMALOUS]++;
            } else {
                counts->values[RB_STATS_NREADS_SINGLE_MAPPED]++;
            }
        }
    }

    if (flag & BAM_FUNMAP) return;

    rb_slave_stats_add_insert(record, counts, collect_diag);
    const uint8_t *nm = rb_slave_aux_get_type_ptr(record, 'N', 'M');
    if (nm) counts->values[RB_STATS_NMISMATCHES] += rb_slave_aux2i(nm);
    counts->values[RB_STATS_NBASES_MAPPED_CIGAR] += rb_slave_stats_mapped_cigar_bases(record);
    rb_slave_stats_update_sort(record, sort_state);
}

extern "C" void slave_mpi_stats_basic_count(MpiStatsBasicCountPara paras[64]) {
    int id = _PEN;
    MpiStatsBasicCountPara *para = &paras[id];

    para->decomp_alloc_cycles = 0;
    para->decomp_inflate_cycles = 0;
    para->decomp_crc_cycles = 0;
    para->decomp_parse_cycles = 0;
    para->decomp_total_cycles = 0;
    memset(&para->sort_state, 0, sizeof(para->sort_state));
    para->sort_state.sorted = 1;
    para->sort_state.first_tid = -1;
    para->sort_state.last_tid = -1;

    bam_block *comp = para->input_block;
    bam_block *un_comp = para->un_comp_block;
    if (comp == NULL) return;

    unsigned long total_t0 = slave_cycle_now();
    struct libdeflate_decompressor *z =
        swbam_cpe_get_decompressor(id, &para->decomp_alloc_cycles);
    if (!z || para->scratch_data == NULL || para->scratch_capacity == 0 || para->counts == NULL) {
        para->status = -2;
        para->decomp_total_cycles = (uint64_t)(slave_cycle_now() - total_t0);
        return;
    }

    if (swbam_cpe_decode_bgzf(comp, un_comp, z,
                                &para->decomp_inflate_cycles,
                                &para->decomp_crc_cycles) != 0) {
        para->status = -2;
        para->decomp_total_cycles = (uint64_t)(slave_cycle_now() - total_t0);
        return;
    }

    bam1_t record;
    memset(&record, 0, sizeof(record));
    record.mempolicy = BAM_USER_OWNS_DATA;

    int total_count = 0;
    int ret = -1;
    unsigned long parse_t0 = slave_cycle_now();
    while (1) {
        record.data = para->scratch_data;
        record.m_data = para->scratch_capacity > UINT32_MAX ?
            UINT32_MAX : (uint32_t)para->scratch_capacity;
        record.l_data = 0;
        record.mempolicy = BAM_USER_OWNS_DATA;
        ret = swbam_cpe_read_bam_record(un_comp, &record);
        if (ret < 0) break;
        rb_slave_stats_add_record(&record, para->counts, &para->sort_state, para->collect_diag);
        total_count++;
    }
    para->decomp_parse_cycles = (uint64_t)(slave_cycle_now() - parse_t0);
    para->n_total_records = total_count;

    if (ret == -5) {
        para->record_index = total_count;
        para->actual_value = record.l_data;
        para->limit_value = para->scratch_capacity;
        para->limit_id = BOUNDS_LIMIT_INIT_DATA_SIZE;
        para->status = -3;
        para->decomp_total_cycles = (uint64_t)(slave_cycle_now() - total_t0);
        return;
    }
    if (ret < -1) {
        para->record_index = total_count;
        para->actual_value = ret;
        para->limit_value = 0;
        para->limit_id = BOUNDS_LIMIT_GENERIC_RUNTIME_ERROR;
        para->status = -2;
        para->decomp_total_cycles = (uint64_t)(slave_cycle_now() - total_t0);
        return;
    }

    para->record_index = 0;
    para->actual_value = 0;
    para->limit_value = 0;
    para->limit_id = BOUNDS_LIMIT_NONE;
    para->decomp_total_cycles = (uint64_t)(slave_cycle_now() - total_t0);
    para->status = 0;
}


