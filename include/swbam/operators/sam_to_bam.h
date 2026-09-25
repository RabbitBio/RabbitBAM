#ifndef SWBAM_OPERATORS_SAM_TO_BAM_H
#define SWBAM_OPERATORS_SAM_TO_BAM_H

#include "swbam/io.h"

namespace swbam {
namespace operators {

struct SamToBamMetrics {
    long long input_chunks;
    long long chunk_groups;
    long long total_records;
    long long compress_groups;
    long long bgzf_blocks;
    long long parse_fast_records;
    long long parse_fallback_records;
    double t_split;
    double t_copy_count;
    double t_parse;
    double t_parse_core;
    double t_parse_aux;
    double t_parse_cg;
    double t_parse_fallback;
    double t_parse_other;
    double t_pack;
    double t_compress;
    double t_compress_serialize;
    double t_compress_alloc;
    double t_compress_deflate;
    double t_compress_footer;
    double t_compress_other;
    double t_write;
    double t_gather;
    double t_optimized_total;
    double t_alloc_init;
    double t_setup_reset;
    double t_compress_setup;
    double t_status_check;
    double t_free_workspace;
};

int RunSamToBamPipeline(MemReader *reader, RankBodySink *sink,
                        const sam_hdr_t *header, int compression_level,
                        SamToBamMetrics *metrics);

} // namespace operators
} // namespace swbam

#endif
