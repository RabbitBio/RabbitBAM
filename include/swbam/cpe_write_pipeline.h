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
    double post_process;
    double total;

    CpeWritePipelineTiming();
};

struct CpeWritePipelineOptions {
    bool overlap_output;
    // Pull sources may fill the next input buffer during current compression.
    bool overlap_source;

    CpeWritePipelineOptions() : overlap_output(true), overlap_source(false) {}
};

// MPE work only: the callback must not launch another CPE kernel or touch
// the current input/output batch. It completes before Submit returns.
typedef int (*CpeWriteBatchPrefetch)(void *context);

class UncompressedBgzfSource {
public:
    virtual ~UncompressedBgzfSource() {}
    virtual int Fill(BgzfBlockBatch *batch, size_t *count) = 0;
};

// MPE-side post-processing extension point for a compressed output batch.
class CompressedBgzfBatchPostProcessor {
public:
    virtual ~CompressedBgzfBatchPostProcessor() {}
    virtual int PostProcessCompressedBatch(const bam_block *blocks,
                                           size_t count) = 0;
};

struct BamRecordPackBlock {
    bam1_t **records;
    int n_records;
    uint32_t total_len;
};

// Inputs are borrowed until Submit returns, or until Complete after Start.
struct CpeWriteBatchInput {
    enum Kind { kUncompressedBgzf, kBamRecords };

    CpeWriteBatchInput(const BgzfBlockBatch *blocks, size_t n)
        : kind(kUncompressedBgzf), uncompressed(blocks), records(nullptr), count(n) {}
    CpeWriteBatchInput(const BamRecordPackBlock *blocks, size_t n)
        : kind(kBamRecords), uncompressed(nullptr), records(blocks), count(n) {}

    Kind kind;
    const BgzfBlockBatch *uncompressed;
    const BamRecordPackBlock *records;
    size_t count;
};

class CpeWriteBatchOperator {
public:
    virtual ~CpeWriteBatchOperator() {}

    virtual const char *name() const = 0;
    virtual size_t batch_capacity() const = 0;
    virtual int Initialize() = 0;
    virtual void Shutdown() = 0;
    virtual bool needs_scratch() const { return false; }
    virtual int Prepare(const CpeWriteBatchInput &input,
                        BgzfBlockBatch *scratch,
                        BgzfBlockBatch *output) = 0;
    virtual void *kernel_entry() const = 0;
    virtual void *kernel_arguments() = 0;
    virtual void ObserveKernel(double wall_seconds,
                               size_t active_blocks) = 0;
    virtual int Validate(size_t active_blocks) const = 0;
    virtual long long OutputBytes(size_t active_blocks) const = 0;
    virtual int Finish() = 0;
};

struct CpeWriteKernelSpec {
    const char *name;
    void *entry;
    size_t batch_capacity;
    bool needs_scratch;
};

// Common batch-level dispatch for block compression and record serialization
// plus compression. Operator-specific Prepare/Validate stay in the derived type.
class CpeWriteKernelOperator : public CpeWriteBatchOperator {
public:
    explicit CpeWriteKernelOperator(const CpeWriteKernelSpec &spec)
        : spec_(spec) {}

    const char *name() const { return spec_.name; }
    size_t batch_capacity() const { return spec_.batch_capacity; }
    bool needs_scratch() const { return spec_.needs_scratch; }
    void *kernel_entry() const { return spec_.entry; }

private:
    CpeWriteKernelSpec spec_;
};

// Shared double-buffered scheduler for pull-style block sources and push-style
// record batches. The post-processor and operator must outlive this session.
// The previous compressed batch is flushed during CPE work.
class CpeWritePipelineSession {
public:
    CpeWritePipelineSession();
    ~CpeWritePipelineSession();

    int Initialize(CompressedBgzfBatchPostProcessor *post_processor,
                   CpeWriteBatchOperator *op,
                   const CpeWritePipelineOptions &options = CpeWritePipelineOptions());
    int Submit(const CpeWriteBatchInput &input);
    int SubmitWithPrefetch(const CpeWriteBatchInput &input,
                           CpeWriteBatchPrefetch prefetch, void *context);
    // Split submission for stateful apps: do MPE-only work between these calls.
    // Input bytes and operator arguments must remain valid and unchanged until
    // Complete. No other CPE kernel may be launched during this interval.
    // Empty Start is a no-op; Complete/Flush/Finish require no other submission
    // in flight. Destruction joins unfinished work but does not flush output.
    int Start(const CpeWriteBatchInput &input);
    int Complete();
    int Flush();
    int Finish();
    const CpeWritePipelineTiming &timing() const;

private:
    class Impl;
    Impl *impl_;

    CpeWritePipelineSession(const CpeWritePipelineSession &);
    CpeWritePipelineSession &operator=(const CpeWritePipelineSession &);
};

int RunCpeWritePipeline(
    UncompressedBgzfSource *source,
    CompressedBgzfBatchPostProcessor *post_processor,
    CpeWriteBatchOperator *op,
    CpeWritePipelineTiming *timing,
    const CpeWritePipelineOptions &options = CpeWritePipelineOptions());

} // namespace cpe
} // namespace swbam

#endif
