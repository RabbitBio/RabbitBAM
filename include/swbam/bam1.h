#ifndef SWBAM_BAM1_H
#define SWBAM_BAM1_H

#include "swbam/raw_bam.h"

#include <cstddef>

#include <htslib/sam.h>

namespace swbam {
namespace cpe {

struct GenericBam1Metrics {
    long long records;
    long long encoded_bytes;
    long long data_bytes;
    size_t peak_batch_records;
    double materialize;
    double consume;

    GenericBam1Metrics();
};

class Bam1RecordConsumer {
public:
    virtual ~Bam1RecordConsumer() {}

    // Records are owned by the pipeline and remain valid only during this
    // call. Use bam_dup1() when a record must outlive the callback.
    virtual int ConsumeBam1(const bam1_t *const *records,
                            size_t count) = 0;
};

int RunGenericBam1Pipeline(
    const BamInputBackend &input,
    const BgzfBlockSpan *spans,
    size_t span_count,
    Bam1RecordConsumer *consumer,
    CpeReadPipelineTiming *timing,
    GenericDecodeMetrics *decode_metrics,
    GenericBam1Metrics *bam1_metrics,
    const CpeReadPipelineOptions &options = CpeReadPipelineOptions());

} // namespace cpe
} // namespace swbam

#endif
