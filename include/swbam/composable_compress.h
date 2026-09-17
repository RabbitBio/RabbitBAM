#ifndef SWBAM_COMPOSABLE_COMPRESS_H
#define SWBAM_COMPOSABLE_COMPRESS_H

#include "swbam/cpe_write_pipeline.h"

namespace swbam {
namespace cpe {

struct ComposableCompressMetrics {
    double alloc;
    double deflate;
    double footer;
    double other;

    ComposableCompressMetrics();
};

int RunComposableCompressPipeline(
    UncompressedBgzfSource *source,
    CompressedBgzfBatchPostProcessor *post_processor,
    int compression_level,
    CpeWritePipelineTiming *timing,
    ComposableCompressMetrics *metrics,
    const CpeWritePipelineOptions &options = CpeWritePipelineOptions());

} // namespace cpe
} // namespace swbam

#endif
