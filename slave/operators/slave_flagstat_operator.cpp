#include "BamTools.h"
#include "swbam/cpe_bam_parser.h"
#include "swbam/cpe_codec.h"

#include <climits>
#include <cstring>

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

#define slave_cycle_now swbam_cpe_cycle_now

static inline void rb_slave_flagstat_add_record(const bam1_t *record,
                                                MpiFlagstatCountSlice *counts) {
    const uint16_t flag = record->core.flag;
    const int bucket = (flag & BAM_FQCFAIL) ? 1 : 0;
    const int is_secondary = (flag & BAM_FSECONDARY) != 0;
    const int is_supplementary = (flag & BAM_FSUPPLEMENTARY) != 0;
    const int is_primary = !is_secondary && !is_supplementary;
    const int is_paired = (flag & BAM_FPAIRED) != 0;
    const int is_mapped = (flag & BAM_FUNMAP) == 0;
    const int mate_mapped = (flag & BAM_FMUNMAP) == 0;

    counts->values[RB_FLAGSTAT_TOTAL][bucket]++;
    if (is_primary) counts->values[RB_FLAGSTAT_PRIMARY][bucket]++;
    if (is_secondary) counts->values[RB_FLAGSTAT_SECONDARY][bucket]++;
    if (is_supplementary) counts->values[RB_FLAGSTAT_SUPPLEMENTARY][bucket]++;
    if (flag & BAM_FDUP) counts->values[RB_FLAGSTAT_DUPLICATES][bucket]++;
    if (is_primary && (flag & BAM_FDUP)) counts->values[RB_FLAGSTAT_PRIMARY_DUPLICATES][bucket]++;
    if (is_mapped) counts->values[RB_FLAGSTAT_MAPPED][bucket]++;
    if (is_primary && is_mapped) counts->values[RB_FLAGSTAT_PRIMARY_MAPPED][bucket]++;

    if (is_primary && is_paired) {
        counts->values[RB_FLAGSTAT_PAIRED][bucket]++;
        if (flag & BAM_FREAD1) counts->values[RB_FLAGSTAT_READ1][bucket]++;
        if (flag & BAM_FREAD2) counts->values[RB_FLAGSTAT_READ2][bucket]++;
        if ((flag & BAM_FPROPER_PAIR) && is_mapped) {
            counts->values[RB_FLAGSTAT_PROPERLY_PAIRED][bucket]++;
        }
        if (is_mapped && mate_mapped) {
            counts->values[RB_FLAGSTAT_PAIR_MAPPED][bucket]++;
            if (record->core.tid != record->core.mtid) {
                counts->values[RB_FLAGSTAT_DIFF_CHR][bucket]++;
                if (record->core.qual >= 5) counts->values[RB_FLAGSTAT_DIFF_CHR_MAPQ5][bucket]++;
            }
        }
        if (is_mapped && !mate_mapped) {
            counts->values[RB_FLAGSTAT_SINGLETONS][bucket]++;
        }
    }
}

extern "C" void slave_mpi_flagstat_count(MpiFlagstatCountPara paras[64]) {
    int id = _PEN;
    MpiFlagstatCountPara *para = &paras[id];

    para->decomp_alloc_cycles = 0;
    para->decomp_inflate_cycles = 0;
    para->decomp_crc_cycles = 0;
    para->decomp_parse_cycles = 0;
    para->decomp_total_cycles = 0;

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
        rb_slave_flagstat_add_record(&record, para->counts);
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


