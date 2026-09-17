#ifndef SWBAM_BAM1_H
#define SWBAM_BAM1_H

#include "swbam/raw_bam.h"

#include <cstddef>

#include <htslib/sam.h>

namespace swbam {
namespace cpe {

struct ComposableBam1Metrics {
    long long records;
    long long encoded_bytes;
    long long data_bytes;
    size_t peak_batch_records;
    double materialize;
    double post_process;

    ComposableBam1Metrics();
};

// MPE-side post-processing extension point for materialized bam1_t records.
class Bam1BatchPostProcessor {
public:
    virtual ~Bam1BatchPostProcessor() {}

    // Records are owned by the pipeline and remain valid only during this
    // call. Use bam_dup1() when a record must outlive the callback.
    virtual int PostProcessBam1Batch(const bam1_t *const *records,
                                     size_t count) = 0;
};

int RunComposableBam1Pipeline(
    const BamInputBackend &input,
    const BgzfBlockSpan *spans,
    size_t span_count,
    Bam1BatchPostProcessor *post_processor,
    CpeReadPipelineTiming *timing,
    ComposableDecodeMetrics *decode_metrics,
    ComposableBam1Metrics *bam1_metrics,
    const CpeReadPipelineOptions &options = CpeReadPipelineOptions());

} // namespace cpe
} // namespace swbam

#endif
