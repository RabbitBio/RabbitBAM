#ifndef SWBAM_RAW_BAM_H
#define SWBAM_RAW_BAM_H

#include "swbam/generic_decode.h"

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

struct GenericRawBamMetrics {
    long long records;
    long long encoded_bytes;

    GenericRawBamMetrics() : records(0), encoded_bytes(0) {}
};

class RawBamRecordConsumer {
public:
    virtual ~RawBamRecordConsumer() {}

    // Views point into the decoded batch and remain valid only during this call.
    virtual int ConsumeRaw(const RawBamRecordView *records,
                           size_t count) = 0;
};

int RunGenericRawBamPipeline(
    const BamInputBackend &input,
    const BgzfBlockSpan *spans,
    size_t span_count,
    RawBamRecordConsumer *consumer,
    CpeReadPipelineTiming *timing,
    GenericDecodeMetrics *decode_metrics,
    GenericRawBamMetrics *raw_metrics,
    const CpeReadPipelineOptions &options = CpeReadPipelineOptions());

} // namespace cpe
} // namespace swbam

#endif
