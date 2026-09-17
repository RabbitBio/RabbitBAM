#ifndef SWBAM_RAW_BAM_H
#define SWBAM_RAW_BAM_H

#include "swbam/composable_decode.h"

#include <cstddef>
#include <cstdint>

namespace swbam {
namespace cpe {

struct RawBamRecordView {
    const unsigned char *encoded;
    uint32_t encoded_size;
    uint32_t block_size;
    uint32_t block_offset;
    uint32_t block_index;
};

struct ComposableRawBamMetrics {
    long long records;
    long long encoded_bytes;

    ComposableRawBamMetrics() : records(0), encoded_bytes(0) {}
};

// MPE-side post-processing extension point for zero-copy raw BAM views.
class RawBamBatchPostProcessor {
public:
    virtual ~RawBamBatchPostProcessor() {}

    // Views point into the decoded batch and remain valid only during this call.
    virtual int PostProcessRawBatch(const RawBamRecordView *records,
                                    size_t count) = 0;
};

int RunComposableRawBamPipeline(
    const BamInputBackend &input,
    const BgzfBlockSpan *spans,
    size_t span_count,
    RawBamBatchPostProcessor *post_processor,
    CpeReadPipelineTiming *timing,
    ComposableDecodeMetrics *decode_metrics,
    ComposableRawBamMetrics *raw_metrics,
    const CpeReadPipelineOptions &options = CpeReadPipelineOptions());

} // namespace cpe
} // namespace swbam

#endif
