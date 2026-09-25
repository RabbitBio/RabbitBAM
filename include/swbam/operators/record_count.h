#ifndef SWBAM_OPERATORS_RECORD_COUNT_H
#define SWBAM_OPERATORS_RECORD_COUNT_H

#include "swbam/cpe_pipeline.h"

namespace swbam {
namespace operators {

struct RecordCountPara {
    bam_block *input_block;
    bam_block *un_comp_block;
    unsigned char *scratch_data;
    size_t scratch_capacity;
    long long records;
    int status;
};

struct RecordCountMetrics {
    long long records;
    long long blocks;
    double read;
    double kernel;
    double total;
};

// Parses each block-aligned BAM record into a reusable bam1_t before counting.
int RunRecordCountPipeline(const BamInputBackend &input,
                           const BgzfBlockSpan *spans, size_t span_count,
                           RecordCountMetrics *metrics);

} // namespace operators
} // namespace swbam

#endif
