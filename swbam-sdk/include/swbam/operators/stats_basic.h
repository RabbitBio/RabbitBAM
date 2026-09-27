#ifndef SWBAM_OPERATORS_STATS_BASIC_H
#define SWBAM_OPERATORS_STATS_BASIC_H

#include "swbam/operators/bam_read_batch.h"

#include <cstddef>

namespace swbam {
namespace operators {

static const int kStatsMaxInsertSize = 8000;
static const int kStatsInsertBins = kStatsMaxInsertSize + 1;

enum StatsBasicValueId {
    kStatsFirstFragments = 0,
    kStatsLastFragments,
    kStatsOtherFragments,
    kStatsFilteredSequences,
    kStatsDuplicateReads,
    kStatsUnmappedReads,
    kStatsSingleMappedReads,
    kStatsPairedAndMappedReads,
    kStatsProperlyPairedReads,
    kStatsPairedReads,
    kStatsAnomalousReads,
    kStatsMq0Reads,
    kStatsQcFailedReads,
    kStatsSecondaryReads,
    kStatsSupplementaryReads,
    kStatsTotalLength,
    kStatsFirstTotalLength,
    kStatsLastTotalLength,
    kStatsDuplicateTotalLength,
    kStatsMappedBases,
    kStatsMappedCigarBases,
    kStatsTrimmedBases,
    kStatsMismatches,
    kStatsMaxLength,
    kStatsFirstMaxLength,
    kStatsLastMaxLength,
    kStatsQualitySum,
    kStatsValueCount
};

enum StatsOrientationDiagnosticId {
    kStatsOrientRefIn = 0,
    kStatsOrientRefOut,
    kStatsOrientRefOther,
    kStatsOrientDocIn,
    kStatsOrientDocOut,
    kStatsOrientDocOther,
    kStatsOrientRefInLtRead,
    kStatsOrientRefInLt2Read,
    kStatsOrientRefInGe2Read,
    kStatsOrientRefOutLtRead,
    kStatsOrientRefOutLt2Read,
    kStatsOrientRefOutGe2Read,
    kStatsOrientPosNegStrandNeg,
    kStatsOrientPosNegStrandPos,
    kStatsOrientPosZeroStrandNeg,
    kStatsOrientPosZeroStrandPos,
    kStatsOrientPosPosStrandNeg,
    kStatsOrientPosPosStrandPos,
    kStatsOrientRead1Left,
    kStatsOrientRead1Right,
    kStatsOrientRead1SamePos,
    kStatsOrientRead2Left,
    kStatsOrientRead2Right,
    kStatsOrientRead2SamePos,
    kStatsOrientIsizeNeg,
    kStatsOrientIsizeZero,
    kStatsOrientIsizePos,
    kStatsOrientationDiagnosticCount
};

struct StatsBasicCounts {
    long long values[kStatsValueCount];
    long long isize_inward[kStatsInsertBins];
    long long isize_outward[kStatsInsertBins];
    long long isize_other[kStatsInsertBins];
    long long orientation_diagnostics[kStatsOrientationDiagnosticCount];

    StatsBasicCounts();
};

struct StatsSortState {
    long long has_coord;
    long long sorted;
    long long first_tid;
    long long first_pos;
    long long last_tid;
    long long last_pos;

    StatsSortState();
};

struct StatsBasicMetrics {
    long long input_blocks;
    long long batch_count;
    long long total_records;
    double read;
    double kernel;
    double decomp_alloc;
    double decomp_inflate;
    double decomp_crc;
    double decomp_parse;
    double decomp_other;
    double post_process;
    double total;

    StatsBasicMetrics();
};

// Rank-local fused decode + parse + stats operator. Count histograms remain
// private to each CPE slot across batches and are merged once in Finish().
class StatsBasicOperator : public BamReadBatchOperator {
public:
    StatsBasicOperator(StatsBasicCounts *counts,
                       StatsSortState *sort_state,
                       StatsBasicMetrics *metrics,
                       bool collect_orientation_diagnostics);
    ~StatsBasicOperator();

    int Initialize();
    void Shutdown();
    int Prepare(const BgzfBlockBatch &compressed,
                BgzfBlockBatch *decoded,
                size_t active_blocks);
    void *kernel_arguments();
    void ObserveKernel(double wall_seconds, size_t active_blocks);
    int Validate(size_t active_blocks) const;
    int PostProcessBatch(size_t active_blocks,
                         long long *records_processed);
    int Finish();

private:
    class Impl;
    Impl *impl_;

    StatsBasicOperator(const StatsBasicOperator &);
    StatsBasicOperator &operator=(const StatsBasicOperator &);
};

int RunStatsBasicPipeline(
    const BamInputBackend &input,
    const BgzfBlockSpan *spans,
    size_t span_count,
    StatsBasicCounts *counts,
    StatsSortState *sort_state,
    StatsBasicMetrics *metrics,
    bool collect_orientation_diagnostics,
    const cpe::CpeReadPipelineOptions &options =
        cpe::CpeReadPipelineOptions());

} // namespace operators
} // namespace swbam

#endif
