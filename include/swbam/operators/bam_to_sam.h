#ifndef SWBAM_OPERATORS_BAM_TO_SAM_H
#define SWBAM_OPERATORS_BAM_TO_SAM_H

#include "swbam/cpe_pipeline.h"
#include "swbam/cpe_sam_write_pipeline.h"

namespace swbam {
namespace operators {

struct BamToSamMetrics {
    long long input_blocks;
    long long group_count;
    long long total_records;
    long long format_tiles;
    double t_decomp;
    double t_decomp_alloc;
    double t_decomp_inflate;
    double t_decomp_crc;
    double t_decomp_parse;
    double t_decomp_other;
    double t_format;
    double t_collect;
    double t_read;
    double t_write;
    double t_gather;
    double t_optimized_total;
    double t_alloc_init;
    double t_free_workspace;
};

int RunBamToSamPipeline(MemReader *reader, RankBodySink *sink,
                        const sam_hdr_t *header, BamToSamMetrics *metrics);

int RunBamToSamPipeline(const BamInputBackend &input,
                        const BgzfBlockSpan *spans, size_t span_count,
                        RankBodySink *sink, const sam_hdr_t *header,
                        BamToSamMetrics *metrics);

} // namespace operators
} // namespace swbam

#endif
