#include "swbam/cpe_pipeline.h"

#include <cstdio>
#include <cstring>
#include <utility>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

namespace swbam {
namespace cpe {

CpeReadPipelineTiming::CpeReadPipelineTiming()
    : input_blocks(0), batch_count(0), total_records(0),
      read(0.0), kernel(0.0), consume(0.0), total(0.0) {}

int RunCpeReadPipeline(const BamInputBackend &input,
                       const BgzfBlockSpan *spans,
                       size_t span_count,
                       CpeBatchOperator *op,
                       CpeReadPipelineTiming *timing,
                       const CpeReadPipelineOptions &options) {
    if (!op || (!spans && span_count != 0) ||
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
    BgzfSpanBatchReader reader(&input, spans, span_count);
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
        const int read_ret = reader.ReadNext(current_compressed, &active_blocks);
        local_timing.read += GetTime() - read_t0;
        if (read_ret != 0) goto cleanup;
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
        if (options.overlap_input) {
            const double read_t0 = GetTime();
            const int read_ret = reader.ReadNext(next_compressed,
                                                 &next_active_blocks);
            local_timing.read += GetTime() - read_t0;
            if (read_ret != 0) {
                athread_join();
                goto cleanup;
            }
        }
#ifdef PLATFORM_SUNWAY
        athread_join();
#endif
        const double kernel_wall = GetTime() - kernel_t0;
        local_timing.kernel += kernel_wall;
        op->ObserveKernel(kernel_wall, active_blocks);

        if (!options.overlap_input) {
            const double read_t0 = GetTime();
            const int read_ret = reader.ReadNext(next_compressed,
                                                 &next_active_blocks);
            local_timing.read += GetTime() - read_t0;
            if (read_ret != 0) goto cleanup;
        }

        if (op->Validate(active_blocks) != 0) goto cleanup;

        const double consume_t0 = GetTime();
        long long batch_records = 0;
        if (op->Consume(active_blocks, &batch_records) != 0) goto cleanup;
        local_timing.consume += GetTime() - consume_t0;
        local_timing.total_records += batch_records;

        std::swap(current_compressed, next_compressed);
        std::swap(current_decoded, next_decoded);
        active_blocks = next_active_blocks;
    }

    {
        const double consume_t0 = GetTime();
        if (op->Finish() != 0) goto cleanup;
        local_timing.consume += GetTime() - consume_t0;
    }
    ret = 0;

cleanup:
    local_timing.total = GetTime() - total_t0;
    if (initialized) op->Shutdown();
    if (timing) *timing = local_timing;
    return ret;
}

} // namespace cpe
} // namespace swbam
