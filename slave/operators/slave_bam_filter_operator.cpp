#include "BamTools.h"
#include "swbam/cpe_bam_parser.h"
#include "swbam/cpe_codec.h"

#include <climits>

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

namespace {

void ResetDiagnostics(Bam2BamPara *para) {
    para->decomp_alloc_cycles = 0;
    para->decomp_inflate_cycles = 0;
    para->decomp_crc_cycles = 0;
    para->decomp_parse_cycles = 0;
    para->decomp_total_cycles = 0;
    para->data_arena_used = 0;
    para->record_index = 0;
    para->actual_value = 0;
    para->limit_value = 0;
    para->limit_id = BOUNDS_LIMIT_NONE;
}

int PrepareRecord(Bam2BamPara *para, bam1_t *record,
                  size_t arena_used, int record_index) {
    if (!para->data_arena || para->data_arena_capacity == 0) return 0;
    if (arena_used >= para->data_arena_capacity) {
        para->record_index = record_index;
        para->actual_value = 1;
        para->limit_value = 0;
        para->limit_id = BOUNDS_LIMIT_INIT_DATA_SIZE;
        return -1;
    }
    const size_t remaining = para->data_arena_capacity - arena_used;
    record->data = para->data_arena + arena_used;
    record->m_data = remaining > UINT32_MAX
        ? UINT32_MAX : (uint32_t)remaining;
    record->l_data = 0;
    record->mempolicy = BAM_USER_OWNS_DATA;
    return 0;
}

void SetArenaError(Bam2BamPara *para, int record_index,
                   size_t arena_used, const bam1_t *record) {
    para->record_index = record_index;
    para->actual_value = record ? (long long)record->l_data
                                : para->actual_value;
    para->limit_value = para->data_arena
        ? (long long)(para->data_arena_capacity - arena_used) : 0;
    para->limit_id = BOUNDS_LIMIT_INIT_DATA_SIZE;
    para->status = -3;
}

} // namespace

extern "C" void slave_mpi_decompress_filterfunc(Bam2BamPara paras[64]) {
    const int id = _PEN;
    Bam2BamPara *para = &paras[id];
    ResetDiagnostics(para);

    bam_block *compressed = para->input_block;
    bam_block *decoded = para->un_comp_block;
    if (!compressed) return;

    const unsigned long total_t0 = swbam_cpe_cycle_now();
    struct libdeflate_decompressor *decompressor =
        swbam_cpe_get_decompressor(id, &para->decomp_alloc_cycles);
    if (!decompressor) {
        para->status = -19;
        para->decomp_total_cycles =
            (uint64_t)(swbam_cpe_cycle_now() - total_t0);
        return;
    }
    if (swbam_cpe_decode_bgzf(
            compressed, decoded, decompressor,
            &para->decomp_inflate_cycles,
            &para->decomp_crc_cycles) != 0) {
        para->status = -20;
        para->decomp_total_cycles =
            (uint64_t)(swbam_cpe_cycle_now() - total_t0);
        return;
    }

    const int capacity = para->record_capacity > 0
        ? para->record_capacity : (int)MAX_RECORDS_PER_BLOCK;
    size_t arena_used = 0;
    int total_count = 0;
    int kept_count = 0;
    uint32_t kept_total_len = 0;
    int read_result = -1;
    const unsigned long parse_t0 = swbam_cpe_cycle_now();
    while (total_count < capacity) {
        bam1_t *record = para->record_base
            ? para->record_base + total_count
            : para->output_records[total_count];
        if (PrepareRecord(para, record, arena_used, total_count) != 0) {
            read_result = -5;
            break;
        }
        read_result = swbam_cpe_read_bam_record(decoded, record);
        if (read_result < 0) break;
        if (para->data_arena) {
            arena_used += ((size_t)record->l_data + 7u) & ~(size_t)7u;
            para->data_arena_used = arena_used;
        }

        total_count++;
        if (bam_filter_matches(record, para->filter)) {
            const uint32_t bam_len =
                (uint32_t)(record->l_data - record->core.l_extranul + 32);
            para->output_records[kept_count] = record;
            para->bam_lens[kept_count] = bam_len;
            kept_total_len += bam_len + 4;
            kept_count++;
        }
    }
    para->decomp_parse_cycles =
        (uint64_t)(swbam_cpe_cycle_now() - parse_t0);
    para->n_total_records = total_count;
    para->n_kept_records = kept_count;
    para->kept_total_len = kept_total_len;

    if (read_result == -5) {
        const bam1_t *record = para->record_base
            ? para->record_base + total_count : nullptr;
        SetArenaError(para, total_count, arena_used, record);
    } else if (read_result < -1) {
        para->status = -21;
    } else if (total_count == capacity && decoded->pos < decoded->length) {
        para->record_index = total_count;
        para->actual_value = total_count + 1;
        para->limit_value = capacity;
        para->limit_id = BOUNDS_LIMIT_MAX_RECORDS_PER_BLOCK;
        para->status = -3;
    } else {
        para->status = 0;
    }
    para->decomp_total_cycles =
        (uint64_t)(swbam_cpe_cycle_now() - total_t0);
}

extern "C" void slave_mpi_decompress_bam2bam_passthrough(
        Bam2BamPara paras[64]) {
    const int id = _PEN;
    Bam2BamPara *para = &paras[id];
    ResetDiagnostics(para);

    bam_block *compressed = para->input_block;
    bam_block *decoded = para->un_comp_block;
    if (!compressed) return;

    const unsigned long total_t0 = swbam_cpe_cycle_now();
    struct libdeflate_decompressor *decompressor =
        swbam_cpe_get_decompressor(id, &para->decomp_alloc_cycles);
    if (!decompressor) {
        para->status = -19;
        para->decomp_total_cycles =
            (uint64_t)(swbam_cpe_cycle_now() - total_t0);
        return;
    }
    if (swbam_cpe_decode_bgzf(
            compressed, decoded, decompressor,
            &para->decomp_inflate_cycles,
            &para->decomp_crc_cycles) != 0) {
        para->status = -20;
        para->decomp_total_cycles =
            (uint64_t)(swbam_cpe_cycle_now() - total_t0);
        return;
    }

    const int capacity = para->record_capacity > 0
        ? para->record_capacity : (int)MAX_RECORDS_PER_BLOCK;
    size_t arena_used = 0;
    int total_count = 0;
    uint32_t kept_total_len = 0;
    int read_result = -1;
    const unsigned long parse_t0 = swbam_cpe_cycle_now();
    while (total_count < capacity) {
        bam1_t *record = para->output_records[total_count];
        if (PrepareRecord(para, record, arena_used, total_count) != 0) {
            read_result = -5;
            break;
        }
        read_result = swbam_cpe_read_bam_record(decoded, record);
        if (read_result < 0) break;
        if (para->data_arena) {
            arena_used += ((size_t)record->l_data + 7u) & ~(size_t)7u;
            para->data_arena_used = arena_used;
        }

        const uint32_t bam_len =
            (uint32_t)(record->l_data - record->core.l_extranul + 32);
        para->bam_lens[total_count] = bam_len;
        kept_total_len += bam_len + 4;
        total_count++;
    }
    para->decomp_parse_cycles =
        (uint64_t)(swbam_cpe_cycle_now() - parse_t0);
    para->n_total_records = total_count;
    para->n_kept_records = total_count;
    para->kept_total_len = kept_total_len;

    if (read_result == -5) {
        const bam1_t *record = total_count < capacity
            ? para->output_records[total_count] : nullptr;
        SetArenaError(para, total_count, arena_used, record);
    } else if (read_result < -1) {
        para->status = -21;
    } else if (total_count == capacity && decoded->pos < decoded->length) {
        para->record_index = total_count;
        para->actual_value = total_count + 1;
        para->limit_value = capacity;
        para->limit_id = BOUNDS_LIMIT_MAX_RECORDS_PER_BLOCK;
        para->status = -3;
    } else {
        para->status = 0;
    }
    para->decomp_total_cycles =
        (uint64_t)(swbam_cpe_cycle_now() - total_t0);
}
