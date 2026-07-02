#ifndef SWBAM_CPE_PIPELINE_H
#define SWBAM_CPE_PIPELINE_H

#include "swbam/io.h"

#include <cstddef>

namespace swbam {
namespace cpe {

struct CpeReadPipelineTiming {
    long long input_blocks;
    long long batch_count;
    long long total_records;
    double read;
    double kernel;
    double consume;
    double total;

    CpeReadPipelineTiming();
};

struct CpeReadPipelineOptions {
    bool overlap_input;

    CpeReadPipelineOptions() : overlap_input(true) {}
};

// A batch-level extension point. Operators keep their fused CPE parameter
// layout and kernel while the runtime owns I/O, double buffering and overlap.
class CpeBatchOperator {
public:
    virtual ~CpeBatchOperator() {}

    virtual const char *name() const = 0;
    virtual size_t batch_capacity() const = 0;
    virtual int Initialize() = 0;
    virtual void Shutdown() = 0;

    virtual int Prepare(const BgzfBlockBatch &compressed,
                        BgzfBlockBatch *decoded,
                        size_t active_blocks) = 0;
    virtual void *kernel_entry() const = 0;
    virtual void *kernel_arguments() = 0;
    virtual void ObserveKernel(double wall_seconds,
                               size_t active_blocks) = 0;
    virtual int Validate(size_t active_blocks) const = 0;
    virtual int Consume(size_t active_blocks,
                        long long *records_processed) = 0;
    virtual int Finish() = 0;
};

int RunCpeReadPipeline(const BamInputBackend &input,
                       const BgzfBlockSpan *spans,
                       size_t span_count,
                       CpeBatchOperator *op,
                       CpeReadPipelineTiming *timing,
                       const CpeReadPipelineOptions &options =
                           CpeReadPipelineOptions());

} // namespace cpe
} // namespace swbam

#endif
