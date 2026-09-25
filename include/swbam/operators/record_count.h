#ifndef SWBAM_OPERATORS_RECORD_COUNT_H
#define SWBAM_OPERATORS_RECORD_COUNT_H

#include "swbam/cpe_pipeline.h"

namespace swbam {
namespace operators {

struct RecordCountPara {
    bam_block *compressed;
    bam_block *decoded;
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

// Counts block-aligned BAM records without materializing bam1_t objects.
int RunRecordCountPipeline(const BamInputBackend &input,
                           const BgzfBlockSpan *spans, size_t span_count,
                           RecordCountMetrics *metrics);

} // namespace operators
} // namespace swbam

#endif
