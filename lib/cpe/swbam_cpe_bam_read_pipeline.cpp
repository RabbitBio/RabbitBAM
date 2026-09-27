#include "swbam/cpe_pipeline.h"

#include <cstdio>
#include <cstring>
#include <utility>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

namespace swbam {
namespace cpe {

namespace {
class SpanBatchSource : public CpeReadBatchSource {
public:
    SpanBatchSource(const BamInputBackend &input,
                    const BgzfBlockSpan *spans, size_t count)
        : reader_(&input, spans, count) {}

    int ReadNext(BgzfBlockBatch *batch, size_t *count) {
        return reader_.ReadNext(batch, count);
    }

private:
    BgzfSpanBatchReader reader_;
};

struct PrefetchContext {
    CpeReadBatchSource *source;
    BgzfBlockBatch *batch;
    size_t *count;
    size_t capacity;
    CpeReadPipelineTiming *timing;
    bool called;
};

int PrefetchNext(void *opaque) {
    PrefetchContext *context = static_cast<PrefetchContext *>(opaque);
    if (context->called) return -1;
    context->called = true;
    const double t0 = GetTime();
    const int status = context->source->ReadNext(context->batch, context->count);
    context->timing->read += GetTime() - t0;
    return status == 0 && *context->count <= context->capacity ? 0 : -1;
}
} // namespace

CpeReadPipelineTiming::CpeReadPipelineTiming()
    : input_blocks(0), batch_count(0), total_records(0),
      read(0.0), kernel(0.0), post_process(0.0), total(0.0) {}

int RunCpeReadPipeline(CpeReadBatchSource *source,
                       CpeBatchOperator *op,
                       CpeReadPipelineTiming *timing,
                       const CpeReadPipelineOptions &options) {
    if (!source || !op ||
        op->batch_capacity() == 0) {
        return -1;
    }

    CpeReadPipelineTiming local_timing;
    const double total_t0 = GetTime();
    int ret = -1;
    bool initialized = false;
    BgzfBlockBatch compressed_a;
    BgzfBlockBatch compressed_b;
    BgzfBlockBatch decoded_a;
    BgzfBlockBatch decoded_b;
    BgzfBlockBatch *current_compressed = &compressed_a;
    BgzfBlockBatch *next_compressed = &compressed_b;
    BgzfBlockBatch *current_decoded = &decoded_a;
    BgzfBlockBatch *next_decoded = &decoded_b;
    size_t active_blocks = 0;

    const int initialize_ret = op->Initialize();
    initialized = true;
    if (initialize_ret != 0) {
        fprintf(stderr, "ERROR: failed to initialize CPE operator %s.\n",
                op->name());
        goto cleanup;
    }

    if (compressed_a.Allocate(op->batch_capacity()) != 0 ||
        compressed_b.Allocate(op->batch_capacity()) != 0 ||
        decoded_a.Allocate(op->batch_capacity()) != 0 ||
        decoded_b.Allocate(op->batch_capacity()) != 0) {
        fprintf(stderr, "ERROR: failed to allocate CPE pipeline buffers for %s.\n",
                op->name());
        goto cleanup;
    }

    {
        const double read_t0 = GetTime();
        const int read_ret = source->ReadNext(current_compressed, &active_blocks);
        local_timing.read += GetTime() - read_t0;
        if (read_ret != 0 || active_blocks > op->batch_capacity()) goto cleanup;
    }

    while (active_blocks > 0) {
        local_timing.input_blocks += static_cast<long long>(active_blocks);
        local_timing.batch_count++;
        if (op->Prepare(*current_compressed, current_decoded,
                        active_blocks) != 0) {
            fprintf(stderr, "ERROR: failed to prepare CPE operator %s.\n",
                    op->name());
            goto cleanup;
        }

        const double kernel_t0 = GetTime();
#ifdef PLATFORM_SUNWAY
        __real_athread_spawn(op->kernel_entry(), op->kernel_arguments(), 1);
#else
        fprintf(stderr, "ERROR: CPE pipeline requires the Sunway backend.\n");
        goto cleanup;
#endif

        size_t next_active_blocks = 0;
        if (options.overlap_input && !options.prefetch_during_post_process) {
            const double read_t0 = GetTime();
            const int read_ret = source->ReadNext(next_compressed,
                                                 &next_active_blocks);
            local_timing.read += GetTime() - read_t0;
            if (read_ret != 0 || next_active_blocks > op->batch_capacity()) {
                athread_join();
                goto cleanup;
            }
        }
        const int overlap_status = op->DuringKernel();
#ifdef PLATFORM_SUNWAY
        athread_join();
#endif
        if (overlap_status != 0) goto cleanup;
        const double kernel_wall = GetTime() - kernel_t0;
        local_timing.kernel += kernel_wall;
        op->ObserveKernel(kernel_wall, active_blocks);

        if (!options.overlap_input && !options.prefetch_during_post_process) {
            const double read_t0 = GetTime();
            const int read_ret = source->ReadNext(next_compressed,
                                                 &next_active_blocks);
            local_timing.read += GetTime() - read_t0;
            if (read_ret != 0 || next_active_blocks > op->batch_capacity()) goto cleanup;
        }

        if (op->Validate(active_blocks) != 0) goto cleanup;

        const double post_process_t0 = GetTime();
        long long batch_records = 0;
        if (options.prefetch_during_post_process) {
            PrefetchContext context = {source, next_compressed,
                                       &next_active_blocks, op->batch_capacity(),
                                       &local_timing, false};
            if (op->PostProcessBatchWithPrefetch(
                    active_blocks, &batch_records, PrefetchNext, &context) != 0)
                goto cleanup;
            if (!context.called && PrefetchNext(&context) != 0) goto cleanup;
        } else if (op->PostProcessBatch(active_blocks, &batch_records) != 0) {
            goto cleanup;
        }
        local_timing.post_process += GetTime() - post_process_t0;
        local_timing.total_records += batch_records;

        std::swap(current_compressed, next_compressed);
        std::swap(current_decoded, next_decoded);
        active_blocks = next_active_blocks;
    }

    {
        const double post_process_t0 = GetTime();
        if (op->Finish() != 0) goto cleanup;
        local_timing.post_process += GetTime() - post_process_t0;
    }
    ret = 0;

cleanup:
    local_timing.total = GetTime() - total_t0;
    if (initialized) op->Shutdown();
    if (timing) *timing = local_timing;
    return ret;
}

int RunCpeReadPipeline(const BamInputBackend &input,
                       const BgzfBlockSpan *spans,
                       size_t span_count,
                       CpeBatchOperator *op,
                       CpeReadPipelineTiming *timing,
                       const CpeReadPipelineOptions &options) {
    if (!spans && span_count != 0) return -1;
    SpanBatchSource source(input, spans, span_count);
    return RunCpeReadPipeline(&source, op, timing, options);
}

} // namespace cpe
} // namespace swbam
