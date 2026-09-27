#ifndef SWBAM_OPERATORS_BGZF_COMPRESS_H
#define SWBAM_OPERATORS_BGZF_COMPRESS_H

#include "swbam/cpe_write_pipeline.h"

namespace swbam {
namespace cpe {

struct BgzfCompressMetrics {
    double alloc;
    double deflate;
    double footer;
    double other;

    BgzfCompressMetrics();
};

int RunBgzfCompressPipeline(
    UncompressedBgzfSource *source,
    CompressedBgzfBatchPostProcessor *post_processor,
    int compression_level,
    CpeWritePipelineTiming *timing,
    BgzfCompressMetrics *metrics,
    const CpeWritePipelineOptions &options = CpeWritePipelineOptions());

} // namespace cpe
} // namespace swbam

#endif
