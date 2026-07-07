#ifndef SWBAM_CPE_WRITE_PIPELINE_H
#define SWBAM_CPE_WRITE_PIPELINE_H

#include "swbam/io.h"

#include <cstddef>

namespace swbam {
namespace cpe {

struct CpeWritePipelineTiming {
    long long input_blocks;
    long long output_blocks;
    long long output_bytes;
    long long batch_count;
    double source;
    double kernel;
    double consume;
    double total;

    CpeWritePipelineTiming();
};

struct CpeWritePipelineOptions {
    bool overlap_output;

    CpeWritePipelineOptions() : overlap_output(true) {}
};

class UncompressedBgzfSource {
public:
    virtual ~UncompressedBgzfSource() {}
    virtual int Fill(BgzfBlockBatch *batch, size_t *count) = 0;
};

class CompressedBgzfConsumer {
public:
    virtual ~CompressedBgzfConsumer() {}
    virtual int ConsumeCompressed(const bam_block *blocks,
                                  size_t count) = 0;
};

class CpeWriteBatchOperator {
public:
    virtual ~CpeWriteBatchOperator() {}

    virtual const char *name() const = 0;
    virtual size_t batch_capacity() const = 0;
    virtual int Initialize() = 0;
    virtual void Shutdown() = 0;
    virtual int Prepare(const BgzfBlockBatch &input,
                        BgzfBlockBatch *output,
                        size_t active_blocks) = 0;
    virtual void *kernel_entry() const = 0;
    virtual void *kernel_arguments() = 0;
    virtual void ObserveKernel(double wall_seconds,
                               size_t active_blocks) = 0;
    virtual int Validate(size_t active_blocks) const = 0;
    virtual long long OutputBytes(size_t active_blocks) const = 0;
    virtual int Finish() = 0;
};

int RunCpeWritePipeline(
    UncompressedBgzfSource *source,
    CompressedBgzfConsumer *consumer,
    CpeWriteBatchOperator *op,
    CpeWritePipelineTiming *timing,
    const CpeWritePipelineOptions &options = CpeWritePipelineOptions());

} // namespace cpe
} // namespace swbam

#endif
