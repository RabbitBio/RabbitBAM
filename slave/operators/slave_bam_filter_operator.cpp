#include "BamTools.h"
#include "swbam/cpe_bam_read_steps.h"

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

template <bool Filter>
class BamTransformAction {
public:
    explicit BamTransformAction(Bam2BamPara *para)
        : para_(para), arena_used_(0), kept_(0), kept_len_(0) {}

    bam1_t *Prepare(int index) {
        bam1_t *record = para_->record_base
            ? para_->record_base + index : para_->output_records[index];
        return PrepareRecord(para_, record, arena_used_, index) == 0
            ? record : nullptr;
    }

    void Process(bam1_t *record, int) {
        if (para_->data_arena) {
            arena_used_ += ((size_t)record->l_data + 7u) & ~(size_t)7u;
            para_->data_arena_used = arena_used_;
        }
        if (!Filter || bam_filter_matches(record, para_->filter)) {
            const uint32_t bam_len =
                (uint32_t)(record->l_data - record->core.l_extranul + 32);
            para_->output_records[kept_] = record;
            para_->bam_lens[kept_] = bam_len;
            kept_len_ += bam_len + 4;
            ++kept_;
        }
    }

    size_t arena_used() const { return arena_used_; }
    int kept() const { return kept_; }
    uint32_t kept_len() const { return kept_len_; }

private:
    Bam2BamPara *para_;
    size_t arena_used_;
    int kept_;
    uint32_t kept_len_;
};

template <bool Filter>
void RunBamTransform(Bam2BamPara *para, int id) {
    ResetDiagnostics(para);
    bam_block *compressed = para->input_block;
    bam_block *decoded = para->un_comp_block;
    if (!compressed) return;
    const unsigned long total_t0 = swbam_cpe_cycle_now();
    const int decode = swbam_cpe_decode_bam_block(
        compressed, decoded, id, &para->decomp_alloc_cycles,
        &para->decomp_inflate_cycles, &para->decomp_crc_cycles);
    if (decode != 0) {
        para->status = decode == -1 ? -19 : -20;
        para->decomp_total_cycles =
            (uint64_t)(swbam_cpe_cycle_now() - total_t0);
        return;
    }
    const int capacity = para->record_capacity > 0
        ? para->record_capacity : (int)MAX_RECORDS_PER_BLOCK;
    BamTransformAction<Filter> action(para);
    const unsigned long parse_t0 = swbam_cpe_cycle_now();
    const SwbamCpeWalkResult result =
        swbam_cpe_walk_bam_records(decoded, capacity, action);
    para->decomp_parse_cycles =
        (uint64_t)(swbam_cpe_cycle_now() - parse_t0);
    para->n_total_records = result.count;
    para->n_kept_records = action.kept();
    para->kept_total_len = action.kept_len();
    if (result.read_result == -5) {
        SetArenaError(para, result.count, action.arena_used(),
                      result.last_record);
    } else if (result.read_result == -6) {
        para->record_index = result.count;
        para->actual_value = result.count + 1;
        para->limit_value = capacity;
        para->limit_id = BOUNDS_LIMIT_MAX_RECORDS_PER_BLOCK;
        para->status = -3;
    } else if (result.read_result < -1) {
        para->status = -21;
    } else {
        para->status = 0;
    }
    para->decomp_total_cycles =
        (uint64_t)(swbam_cpe_cycle_now() - total_t0);
}

} // namespace

extern "C" void slave_mpi_decompress_filterfunc(Bam2BamPara paras[64]) {
    RunBamTransform<true>(&paras[_PEN], _PEN);
}

extern "C" void slave_mpi_decompress_bam2bam_passthrough(
        Bam2BamPara paras[64]) {
    RunBamTransform<false>(&paras[_PEN], _PEN);
}
