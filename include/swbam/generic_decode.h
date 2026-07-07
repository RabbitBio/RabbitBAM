#ifndef SWBAM_GENERIC_DECODE_H
#define SWBAM_GENERIC_DECODE_H

#include "swbam/cpe_pipeline.h"

namespace swbam {
namespace cpe {

struct GenericDecodeMetrics {
    long long decoded_blocks;
    long long decoded_bytes;
    double alloc;
    double inflate;
    double crc;
    double other;

    GenericDecodeMetrics();
};

class DecodedBgzfConsumer {
public:
    virtual ~DecodedBgzfConsumer() {}

    // The view remains valid only until this call returns.
    virtual int ConsumeDecoded(const bam_block *blocks, size_t count) = 0;
};

int RunGenericDecodePipeline(
    const BamInputBackend &input,
    const BgzfBlockSpan *spans,
    size_t span_count,
    DecodedBgzfConsumer *consumer,
    CpeReadPipelineTiming *timing,
    GenericDecodeMetrics *metrics,
    const CpeReadPipelineOptions &options = CpeReadPipelineOptions());

} // namespace cpe
} // namespace swbam

#endif
