#ifndef SWBAM_OPERATORS_BAM_TRANSFORM_H
#define SWBAM_OPERATORS_BAM_TRANSFORM_H

#include "swbam/cpe_pipeline.h"
#include "swbam/cpe_record_write_adapter.h"

namespace swbam {
namespace operators {

struct BamTransformMetrics {
    long long input_blocks;
    long long group_count;
    long long total_records;
    long long kept_records;
    long long dropped_records;
    long long bgzf_blocks;
    long long pack_records;
    double t_decomp_filter;
    double t_decomp_alloc;
    double t_decomp_inflate;
    double t_decomp_crc;
    double t_decomp_parse;
    double t_decomp_other;
    double t_pack;
    double t_compress;
    double t_compress_serialize;
    double t_compress_alloc;
    double t_compress_deflate;
    double t_compress_footer;
    double t_compress_other;
    double t_read;
    double t_write;
    double t_mpi_write;
    double t_optimized_total;
    double t_alloc_init;
    double t_initial_read;
    double t_prepare_decomp;
    double t_decomp_check;
    double t_compress_prepare;
    double t_compress_check;
    double t_empty_group_path;
    double t_final_flush;
    double t_free_workspace;
};

int RunBamTransformPipeline(
    const BamInputBackend &input,
    const BgzfBlockSpan *spans, size_t span_count,
    RankBodySink *sink, const BamFilterOptions &filter,
    int compression_level, BamTransformMetrics *metrics);

// The memory overload borrows the caller's compressed rank body in place.
int RunBamTransformPipeline(
    MemReader *reader, RankBodySink *sink,
    const BamFilterOptions &filter, int compression_level,
    BamTransformMetrics *metrics);

} // namespace operators
} // namespace swbam

#endif
