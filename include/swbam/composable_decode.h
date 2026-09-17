#ifndef SWBAM_COMPOSABLE_DECODE_H
#define SWBAM_COMPOSABLE_DECODE_H

#include "swbam/cpe_pipeline.h"

namespace swbam {
namespace cpe {

struct ComposableDecodeMetrics {
    long long decoded_blocks;
    long long decoded_bytes;
    double alloc;
    double inflate;
    double crc;
    double other;

    ComposableDecodeMetrics();
};

// MPE-side post-processing extension point invoked after a decoded batch
// completes on the CPEs.
class DecodedBgzfBatchPostProcessor {
public:
    virtual ~DecodedBgzfBatchPostProcessor() {}

    // The batch remains valid only until this call returns.
    virtual int PostProcessDecodedBatch(const bam_block *blocks, size_t count) = 0;
};

int RunComposableDecodePipeline(
    const BamInputBackend &input,
    const BgzfBlockSpan *spans,
    size_t span_count,
    DecodedBgzfBatchPostProcessor *post_processor,
    CpeReadPipelineTiming *timing,
    ComposableDecodeMetrics *metrics,
    const CpeReadPipelineOptions &options = CpeReadPipelineOptions());

} // namespace cpe
} // namespace swbam

#endif
