#ifndef SWBAM_GENERIC_COMPRESS_H
#define SWBAM_GENERIC_COMPRESS_H

#include "swbam/cpe_write_pipeline.h"

namespace swbam {
namespace cpe {

struct GenericCompressMetrics {
    double alloc;
    double deflate;
    double footer;
    double other;

    GenericCompressMetrics();
};

int RunGenericCompressPipeline(
    UncompressedBgzfSource *source,
    CompressedBgzfConsumer *consumer,
    int compression_level,
    CpeWritePipelineTiming *timing,
    GenericCompressMetrics *metrics,
    const CpeWritePipelineOptions &options = CpeWritePipelineOptions());

} // namespace cpe
} // namespace swbam

#endif
