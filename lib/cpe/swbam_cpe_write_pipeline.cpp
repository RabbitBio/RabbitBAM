#include "swbam/cpe_write_pipeline.h"

#include <cstdio>
#include <utility>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

namespace swbam {
namespace cpe {

CpeWritePipelineTiming::CpeWritePipelineTiming()
    : input_blocks(0), output_blocks(0), output_bytes(0),
      batch_count(0), source(0.0), kernel(0.0), consume(0.0),
      total(0.0) {}

int RunCpeWritePipeline(
        UncompressedBgzfSource *source,
        CompressedBgzfConsumer *consumer,
        CpeWriteBatchOperator *op,
        CpeWritePipelineTiming *timing,
        const CpeWritePipelineOptions &options) {
    if (!source || !consumer || !op || op->batch_capacity() == 0) return -1;

    CpeWritePipelineTiming local;
    const double total_t0 = GetTime();
    int ret = -1;
    bool initialized = false;
    BgzfBlockBatch input_a;
    BgzfBlockBatch input_b;
    BgzfBlockBatch output_a;
    BgzfBlockBatch output_b;
    BgzfBlockBatch *current_input = &input_a;
    BgzfBlockBatch *next_input = &input_b;
    BgzfBlockBatch *current_output = &output_a;
    BgzfBlockBatch *next_output = &output_b;
    const bam_block *previous_output = nullptr;
    size_t previous_count = 0;
    size_t active_blocks = 0;

    const int initialize_ret = op->Initialize();
    initialized = true;
    if (initialize_ret != 0) goto cleanup;
    if (input_a.Allocate(op->batch_capacity()) != 0 ||
        input_b.Allocate(op->batch_capacity()) != 0 ||
        output_a.Allocate(op->batch_capacity()) != 0 ||
        output_b.Allocate(op->batch_capacity()) != 0) {
        fprintf(stderr, "ERROR: failed to allocate write pipeline buffers for %s.\n",
                op->name());
        goto cleanup;
    }

    {
        const double source_t0 = GetTime();
        const int source_ret = source->Fill(current_input, &active_blocks);
        local.source += GetTime() - source_t0;
        if (source_ret != 0 || active_blocks > op->batch_capacity()) goto cleanup;
    }

    while (active_blocks > 0) {
        local.input_blocks += (long long)active_blocks;
        local.batch_count++;
        if (op->Prepare(*current_input, current_output, active_blocks) != 0) {
            goto cleanup;
        }

        if (previous_output && !options.overlap_output) {
            const double consume_t0 = GetTime();
            if (consumer->ConsumeCompressed(
                    previous_output, previous_count) != 0) goto cleanup;
            local.consume += GetTime() - consume_t0;
            previous_output = nullptr;
            previous_count = 0;
        }

        const double kernel_t0 = GetTime();
#ifdef PLATFORM_SUNWAY
        __real_athread_spawn(op->kernel_entry(), op->kernel_arguments(), 1);
#else
        goto cleanup;
#endif

        int overlap_consume_ret = 0;
        if (previous_output && options.overlap_output) {
            const double consume_t0 = GetTime();
            overlap_consume_ret = consumer->ConsumeCompressed(
                previous_output, previous_count);
            local.consume += GetTime() - consume_t0;
        }
#ifdef PLATFORM_SUNWAY
        athread_join();
#endif
        const double kernel_wall = GetTime() - kernel_t0;
        local.kernel += kernel_wall;
        op->ObserveKernel(kernel_wall, active_blocks);
        if (overlap_consume_ret != 0 || op->Validate(active_blocks) != 0) {
            goto cleanup;
        }

        local.output_blocks += (long long)active_blocks;
        local.output_bytes += op->OutputBytes(active_blocks);
        previous_output = current_output->blocks();
        previous_count = active_blocks;

        size_t next_blocks = 0;
        const double source_t0 = GetTime();
        const int source_ret = source->Fill(next_input, &next_blocks);
        local.source += GetTime() - source_t0;
        if (source_ret != 0 || next_blocks > op->batch_capacity()) goto cleanup;

        std::swap(current_input, next_input);
        std::swap(current_output, next_output);
        active_blocks = next_blocks;
    }

    if (previous_output) {
        const double consume_t0 = GetTime();
        if (consumer->ConsumeCompressed(
                previous_output, previous_count) != 0) goto cleanup;
        local.consume += GetTime() - consume_t0;
    }
    if (op->Finish() != 0) goto cleanup;
    ret = 0;

cleanup:
    local.total = GetTime() - total_t0;
    if (initialized) op->Shutdown();
    if (timing) *timing = local;
    return ret;
}

} // namespace cpe
} // namespace swbam
