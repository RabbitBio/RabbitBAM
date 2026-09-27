#include "swbam/operators/stats_basic.h"
#include "swbam/operators/bam_read_batch.h"

#include "swbam/bam_types.h"

#include <cstdio>
#include <cstring>
#include <stdint.h>

extern "C" {
    void slave_mpi_stats_basic_count();
}

namespace swbam {
namespace operators {

namespace {

const int kStatsBatchBlocks = 64;

#define SWBAM_STATS_ID_MATCH(public_id, cpe_id)                         \
    static_assert(static_cast<int>(public_id) == static_cast<int>(cpe_id), \
                  "stats ID mismatch")

SWBAM_STATS_ID_MATCH(kStatsFirstFragments, RB_STATS_NREADS_1ST);
SWBAM_STATS_ID_MATCH(kStatsLastFragments, RB_STATS_NREADS_2ND);
SWBAM_STATS_ID_MATCH(kStatsOtherFragments, RB_STATS_NREADS_OTHER);
SWBAM_STATS_ID_MATCH(kStatsFilteredSequences, RB_STATS_NREADS_FILTERED);
SWBAM_STATS_ID_MATCH(kStatsDuplicateReads, RB_STATS_NREADS_DUP);
SWBAM_STATS_ID_MATCH(kStatsUnmappedReads, RB_STATS_NREADS_UNMAPPED);
SWBAM_STATS_ID_MATCH(kStatsSingleMappedReads, RB_STATS_NREADS_SINGLE_MAPPED);
SWBAM_STATS_ID_MATCH(kStatsPairedAndMappedReads,
                     RB_STATS_NREADS_PAIRED_AND_MAPPED);
SWBAM_STATS_ID_MATCH(kStatsProperlyPairedReads,
                     RB_STATS_NREADS_PROPERLY_PAIRED);
SWBAM_STATS_ID_MATCH(kStatsPairedReads, RB_STATS_NREADS_PAIRED_TECH);
SWBAM_STATS_ID_MATCH(kStatsAnomalousReads, RB_STATS_NREADS_ANOMALOUS);
SWBAM_STATS_ID_MATCH(kStatsMq0Reads, RB_STATS_NREADS_MQ0);
SWBAM_STATS_ID_MATCH(kStatsQcFailedReads, RB_STATS_NREADS_QCFAILED);
SWBAM_STATS_ID_MATCH(kStatsSecondaryReads, RB_STATS_NREADS_SECONDARY);
SWBAM_STATS_ID_MATCH(kStatsSupplementaryReads,
                     RB_STATS_NREADS_SUPPLEMENTARY);
SWBAM_STATS_ID_MATCH(kStatsTotalLength, RB_STATS_TOTAL_LEN);
SWBAM_STATS_ID_MATCH(kStatsFirstTotalLength, RB_STATS_TOTAL_LEN_1ST);
SWBAM_STATS_ID_MATCH(kStatsLastTotalLength, RB_STATS_TOTAL_LEN_2ND);
SWBAM_STATS_ID_MATCH(kStatsDuplicateTotalLength, RB_STATS_TOTAL_LEN_DUP);
SWBAM_STATS_ID_MATCH(kStatsMappedBases, RB_STATS_NBASES_MAPPED);
SWBAM_STATS_ID_MATCH(kStatsMappedCigarBases, RB_STATS_NBASES_MAPPED_CIGAR);
SWBAM_STATS_ID_MATCH(kStatsTrimmedBases, RB_STATS_NBASES_TRIMMED);
SWBAM_STATS_ID_MATCH(kStatsMismatches, RB_STATS_NMISMATCHES);
SWBAM_STATS_ID_MATCH(kStatsMaxLength, RB_STATS_MAX_LEN);
SWBAM_STATS_ID_MATCH(kStatsFirstMaxLength, RB_STATS_MAX_LEN_1ST);
SWBAM_STATS_ID_MATCH(kStatsLastMaxLength, RB_STATS_MAX_LEN_2ND);
SWBAM_STATS_ID_MATCH(kStatsQualitySum, RB_STATS_SUM_QUAL);
SWBAM_STATS_ID_MATCH(kStatsValueCount, RB_STATS_LONG_COUNT);
SWBAM_STATS_ID_MATCH(kStatsOrientationDiagnosticCount,
                     RB_ORIENT_DIAG_COUNT);
static_assert(kStatsMaxInsertSize == RB_MPI_STATS_MAX_INSERT_SIZE,
              "stats insert-size limit mismatch");

#undef SWBAM_STATS_ID_MATCH

void AccumulateDecompDetail(const MpiStatsBasicCountPara *paras,
                            int active_blocks,
                            double kernel_wall,
                            StatsBasicMetrics *metrics) {
    if (!metrics || active_blocks <= 0 || kernel_wall <= 0.0) return;

    const MpiStatsBasicCountPara *critical = nullptr;
    uint64_t critical_total = 0;
    for (int b = 0; b < active_blocks; ++b) {
        if (paras[b].decomp_total_cycles >= critical_total) {
            critical_total = paras[b].decomp_total_cycles;
            critical = &paras[b];
        }
    }
    if (!critical || critical_total == 0) {
        metrics->decomp_other += kernel_wall;
        return;
    }

    const double scale = kernel_wall / static_cast<double>(critical_total);
    const double alloc = scale * critical->decomp_alloc_cycles;
    const double inflate = scale * critical->decomp_inflate_cycles;
    const double crc = scale * critical->decomp_crc_cycles;
    const double parse = scale * critical->decomp_parse_cycles;
    double other = kernel_wall - alloc - inflate - crc - parse;
    if (other < 0.0) other = 0.0;

    metrics->decomp_alloc += alloc;
    metrics->decomp_inflate += inflate;
    metrics->decomp_crc += crc;
    metrics->decomp_parse += parse;
    metrics->decomp_other += other;
}

void MergeSlice(const MpiStatsBasicCountSlice &slice,
                StatsBasicCounts *counts) {
    for (int i = 0; i < kStatsValueCount; ++i) {
        if (i == kStatsMaxLength || i == kStatsFirstMaxLength ||
            i == kStatsLastMaxLength) {
            if (slice.values[i] > counts->values[i]) {
                counts->values[i] = slice.values[i];
            }
            continue;
        }
        counts->values[i] += slice.values[i];
    }
    for (int i = 0; i < kStatsInsertBins; ++i) {
        counts->isize_inward[i] += slice.isize_inward[i];
        counts->isize_outward[i] += slice.isize_outward[i];
        counts->isize_other[i] += slice.isize_other[i];
    }
    for (int i = 0; i < kStatsOrientationDiagnosticCount; ++i) {
        counts->orientation_diagnostics[i] += slice.orient_diag[i];
    }
}

void MergeBlockSort(const MpiStatsBlockSortState &block_sort,
                    StatsSortState *rank_sort) {
    if (!block_sort.has_coord) return;
    if (!rank_sort->has_coord) {
        rank_sort->has_coord = 1;
        rank_sort->first_tid = block_sort.first_tid;
        rank_sort->first_pos = block_sort.first_pos;
        rank_sort->last_tid = block_sort.last_tid;
        rank_sort->last_pos = block_sort.last_pos;
        if (!block_sort.sorted) rank_sort->sorted = 0;
        return;
    }
    if (block_sort.first_tid < rank_sort->last_tid ||
        (block_sort.first_tid == rank_sort->last_tid &&
         block_sort.first_pos < rank_sort->last_pos)) {
        rank_sort->sorted = 0;
    }
    if (!block_sort.sorted) rank_sort->sorted = 0;
    rank_sort->last_tid = block_sort.last_tid;
    rank_sort->last_pos = block_sort.last_pos;
}

} // namespace

StatsBasicCounts::StatsBasicCounts() {
    std::memset(this, 0, sizeof(*this));
}

StatsSortState::StatsSortState()
    : has_coord(0), sorted(1), first_tid(-1), first_pos(0),
      last_tid(-1), last_pos(0) {}

StatsBasicMetrics::StatsBasicMetrics()
    : input_blocks(0), batch_count(0), total_records(0),
      read(0.0), kernel(0.0), decomp_alloc(0.0),
      decomp_inflate(0.0), decomp_crc(0.0), decomp_parse(0.0),
      decomp_other(0.0), post_process(0.0), total(0.0) {}

class StatsBasicOperator::Impl {
public:
    Impl(StatsBasicCounts *counts_in,
         StatsSortState *sort_state_in,
         StatsBasicMetrics *metrics_in,
         bool collect_orientation_diagnostics_in)
        : counts(counts_in), sort_state(sort_state_in), metrics(metrics_in),
          collect_orientation_diagnostics(
              collect_orientation_diagnostics_in ? 1 : 0),
          count_slices(nullptr), scratch_data(nullptr) {
        std::memset(paras, 0, sizeof(paras));
    }

    StatsBasicCounts *counts;
    StatsSortState *sort_state;
    StatsBasicMetrics *metrics;
    int collect_orientation_diagnostics;
    MpiStatsBasicCountPara paras[kStatsBatchBlocks];
    MpiStatsBasicCountSlice *count_slices;
    unsigned char *scratch_data;
};

StatsBasicOperator::StatsBasicOperator(
        StatsBasicCounts *counts,
        StatsSortState *sort_state,
        StatsBasicMetrics *metrics,
        bool collect_orientation_diagnostics)
    : BamReadBatchOperator({"stats-basic",
          reinterpret_cast<void *>(slave_mpi_stats_basic_count),
          kStatsBatchBlocks}),
      impl_(new Impl(counts, sort_state, metrics,
                     collect_orientation_diagnostics)) {}

StatsBasicOperator::~StatsBasicOperator() {
    Shutdown();
    delete impl_;
}

int StatsBasicOperator::Initialize() {
    if (!impl_ || !impl_->counts || !impl_->sort_state) return -1;
    Shutdown();
    *impl_->counts = StatsBasicCounts();
    *impl_->sort_state = StatsSortState();
    if (impl_->metrics) *impl_->metrics = StatsBasicMetrics();
    impl_->count_slices =
        reinterpret_cast<MpiStatsBasicCountSlice *>(aligned_alloc_custom(
            64, static_cast<size_t>(kStatsBatchBlocks) *
                    sizeof(MpiStatsBasicCountSlice)));
    impl_->scratch_data = aligned_alloc_custom(
        64, static_cast<size_t>(kStatsBatchBlocks) *
                MPI_BAM_BLOCK_ARENA_SIZE);
    if (!impl_->count_slices || !impl_->scratch_data) return -1;
    std::memset(impl_->count_slices, 0,
                static_cast<size_t>(kStatsBatchBlocks) *
                    sizeof(MpiStatsBasicCountSlice));
    return 0;
}

void StatsBasicOperator::Shutdown() {
    if (!impl_) return;
    if (impl_->count_slices) {
        aligned_free_custom(
            reinterpret_cast<unsigned char *>(impl_->count_slices));
        impl_->count_slices = nullptr;
    }
    if (impl_->scratch_data) {
        aligned_free_custom(impl_->scratch_data);
        impl_->scratch_data = nullptr;
    }
}

int StatsBasicOperator::Prepare(const BgzfBlockBatch &compressed,
                                BgzfBlockBatch *decoded,
                                size_t active_blocks) {
    if (!impl_ || !impl_->count_slices || !impl_->scratch_data ||
        !decoded || active_blocks > kStatsBatchBlocks) {
        return -1;
    }
    BindBamReadBatch(impl_->paras, kStatsBatchBlocks,
                     compressed, decoded, active_blocks);
    for (int b = 0; b < kStatsBatchBlocks; ++b) {
        MpiStatsBasicCountPara &para = impl_->paras[b];
        para.block_id = b;
        para.scratch_data = impl_->scratch_data +
            static_cast<size_t>(b) * MPI_BAM_BLOCK_ARENA_SIZE;
        para.scratch_capacity = MPI_BAM_BLOCK_ARENA_SIZE;
        para.counts = &impl_->count_slices[b];
        para.collect_diag = impl_->collect_orientation_diagnostics;
        para.n_total_records = 0;
        para.record_index = 0;
        para.actual_value = 0;
        para.limit_value = 0;
        para.limit_id = BOUNDS_LIMIT_NONE;
        std::memset(&para.sort_state, 0, sizeof(para.sort_state));
        para.sort_state.sorted = 1;
        para.sort_state.first_tid = -1;
        para.sort_state.last_tid = -1;
        para.decomp_alloc_cycles = 0;
        para.decomp_inflate_cycles = 0;
        para.decomp_crc_cycles = 0;
        para.decomp_parse_cycles = 0;
        para.decomp_total_cycles = 0;
    }
    return 0;
}

void *StatsBasicOperator::kernel_arguments() {
    return impl_ ? impl_->paras : nullptr;
}

void StatsBasicOperator::ObserveKernel(double wall_seconds,
                                       size_t active_blocks) {
    if (!impl_) return;
    AccumulateDecompDetail(impl_->paras, static_cast<int>(active_blocks),
                           wall_seconds, impl_->metrics);
}

int StatsBasicOperator::Validate(size_t active_blocks) const {
    if (!impl_ || active_blocks > kStatsBatchBlocks) return -1;
    for (size_t b = 0; b < active_blocks; ++b) {
        const MpiStatsBasicCountPara &para = impl_->paras[b];
        if (para.status == 0) continue;
        if (para.status == -3) {
            std::fprintf(
                stderr,
                "ERROR: MPI stats capacity exceeded on input block %zu. "
                "limit_id=%d limit=%lld actual=%lld record=%d.\n",
                b, para.limit_id, para.limit_value,
                para.actual_value, para.record_index);
        } else {
            std::fprintf(
                stderr,
                "ERROR: MPI stats count failed on input block %zu with "
                "status %d. record=%d actual=%lld limit_id=%d.\n",
                b, para.status, para.record_index,
                para.actual_value, para.limit_id);
        }
        return -1;
    }
    return 0;
}

int StatsBasicOperator::PostProcessBatch(size_t active_blocks,
                                         long long *records_processed) {
    if (!impl_ || !impl_->sort_state ||
        active_blocks > kStatsBatchBlocks) {
        return -1;
    }
    long long records = 0;
    for (size_t b = 0; b < active_blocks; ++b) {
        records += impl_->paras[b].n_total_records;
        MergeBlockSort(impl_->paras[b].sort_state, impl_->sort_state);
    }
    if (records_processed) *records_processed = records;
    return 0;
}

int StatsBasicOperator::Finish() {
    if (!impl_ || !impl_->counts || !impl_->count_slices) return -1;
    for (int b = 0; b < kStatsBatchBlocks; ++b) {
        MergeSlice(impl_->count_slices[b], impl_->counts);
    }
    return 0;
}

int RunStatsBasicPipeline(
        const BamInputBackend &input,
        const BgzfBlockSpan *spans,
        size_t span_count,
        StatsBasicCounts *counts,
        StatsSortState *sort_state,
        StatsBasicMetrics *metrics,
        bool collect_orientation_diagnostics,
        const cpe::CpeReadPipelineOptions &options) {
    if (!counts || !sort_state) return -1;

    StatsBasicOperator op(counts, sort_state, metrics,
                          collect_orientation_diagnostics);
    cpe::CpeReadPipelineTiming timing;
    const int ret = cpe::RunCpeReadPipeline(
        input, spans, span_count, &op, &timing, options);
    if (metrics) {
        metrics->input_blocks = timing.input_blocks;
        metrics->batch_count = timing.batch_count;
        metrics->total_records = timing.total_records;
        metrics->read = timing.read;
        metrics->kernel = timing.kernel;
        metrics->post_process = timing.post_process;
        metrics->total = timing.total;
    }
    return ret;
}

} // namespace operators
} // namespace swbam
