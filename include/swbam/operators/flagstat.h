#ifndef SWBAM_OPERATORS_FLAGSTAT_H
#define SWBAM_OPERATORS_FLAGSTAT_H

#include "swbam/cpe_pipeline.h"

#include <cstddef>

namespace swbam {
namespace operators {

enum FlagstatCounterId {
    kFlagstatTotal = 0,
    kFlagstatPrimary,
    kFlagstatSecondary,
    kFlagstatSupplementary,
    kFlagstatDuplicates,
    kFlagstatPrimaryDuplicates,
    kFlagstatMapped,
    kFlagstatPrimaryMapped,
    kFlagstatPaired,
    kFlagstatRead1,
    kFlagstatRead2,
    kFlagstatProperlyPaired,
    kFlagstatPairMapped,
    kFlagstatSingletons,
    kFlagstatDiffChr,
    kFlagstatDiffChrMapq5,
    kFlagstatCounterCount
};

struct FlagstatCounts {
    long long values[kFlagstatCounterCount][2];

    FlagstatCounts();
};

struct FlagstatMetrics {
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

    FlagstatMetrics();
};

// Rank-local fused decode + parse + flagstat operator. The common read
// pipeline owns BGZF batches and scheduling; this object owns only operator
// parameters, per-CPE scratch space and rank-local results.
class FlagstatOperator : public cpe::CpeBatchOperator {
public:
    FlagstatOperator(FlagstatCounts *counts, FlagstatMetrics *metrics);
    ~FlagstatOperator();

    const char *name() const;
    size_t batch_capacity() const;
    int Initialize();
    void Shutdown();
    int Prepare(const BgzfBlockBatch &compressed,
                BgzfBlockBatch *decoded,
                size_t active_blocks);
    void *kernel_entry() const;
    void *kernel_arguments();
    void ObserveKernel(double wall_seconds, size_t active_blocks);
    int Validate(size_t active_blocks) const;
    int PostProcessBatch(size_t active_blocks,
                         long long *records_processed);
    int Finish();

private:
    class Impl;
    Impl *impl_;

    FlagstatOperator(const FlagstatOperator &);
    FlagstatOperator &operator=(const FlagstatOperator &);
};

int RunFlagstatPipeline(
    const BamInputBackend &input,
    const BgzfBlockSpan *spans,
    size_t span_count,
    FlagstatCounts *counts,
    FlagstatMetrics *metrics,
    const cpe::CpeReadPipelineOptions &options =
        cpe::CpeReadPipelineOptions());

} // namespace operators
} // namespace swbam

#endif
