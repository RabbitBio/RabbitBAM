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
      batch_count(0), source(0.0), kernel(0.0), post_process(0.0),
      total(0.0) {}

class CpeWritePipelineSession::Impl {
public:
    Impl()
        : post_processor(nullptr), op(nullptr), current_output(&output_a),
          pending_output(nullptr), current_scratch(&scratch_a),
          pending_count(0), started_at(0.0), op_initialized(false),
          initialized(false), finished(false), failed(false) {}

    ~Impl() {
        if (op_initialized) op->Shutdown();
    }

    int FlushPending() {
        if (!pending_output) return 0;
        const double t0 = GetTime();
        const int ret = post_processor->PostProcessCompressedBatch(
            pending_output->blocks(), pending_count);
        timing.post_process += GetTime() - t0;
        if (ret != 0) {
            failed = true;
            return -1;
        }
        pending_output = nullptr;
        pending_count = 0;
        return 0;
    }

    CompressedBgzfBatchPostProcessor *post_processor;
    CpeWriteBatchOperator *op;
    CpeWritePipelineOptions options;
    CpeWritePipelineTiming timing;
    BgzfBlockBatch output_a;
    BgzfBlockBatch output_b;
    BgzfBlockBatch scratch_a;
    BgzfBlockBatch scratch_b;
    BgzfBlockBatch *current_output;
    BgzfBlockBatch *pending_output;
    BgzfBlockBatch *current_scratch;
    size_t pending_count;
    double started_at;
    bool op_initialized;
    bool initialized;
    bool finished;
    bool failed;
};

CpeWritePipelineSession::CpeWritePipelineSession() : impl_(new Impl) {}
CpeWritePipelineSession::~CpeWritePipelineSession() { delete impl_; }

int CpeWritePipelineSession::Initialize(
        CompressedBgzfBatchPostProcessor *post_processor,
        CpeWriteBatchOperator *op,
        const CpeWritePipelineOptions &options) {
    if (!impl_ || !post_processor || !op || op->batch_capacity() == 0 ||
        impl_->initialized || impl_->op_initialized) return -1;
    impl_->started_at = GetTime();
    impl_->post_processor = post_processor;
    impl_->op = op;
    impl_->options = options;
    impl_->op_initialized = true;
    if (op->Initialize() != 0 ||
        impl_->output_a.Allocate(op->batch_capacity()) != 0 ||
        impl_->output_b.Allocate(op->batch_capacity()) != 0 ||
        (op->needs_scratch() &&
         (impl_->scratch_a.Allocate(op->batch_capacity()) != 0 ||
          impl_->scratch_b.Allocate(op->batch_capacity()) != 0))) {
        fprintf(stderr, "ERROR: failed to initialize write pipeline for %s.\n",
                op->name());
        impl_->failed = true;
        return -1;
    }
    impl_->initialized = true;
    return 0;
}

int CpeWritePipelineSession::Submit(const CpeWriteBatchInput &input) {
    if (!impl_ || !impl_->initialized || impl_->failed || impl_->finished ||
        input.count > impl_->op->batch_capacity()) return -1;
    if (input.count == 0) return 0;
    if ((input.kind == CpeWriteBatchInput::kUncompressedBgzf &&
         (!input.uncompressed || input.count > input.uncompressed->capacity())) ||
        (input.kind == CpeWriteBatchInput::kBamRecords && !input.records)) {
        return -1;
    }
    BgzfBlockBatch *scratch = impl_->op->needs_scratch()
        ? impl_->current_scratch : nullptr;
    if (impl_->op->Prepare(input, scratch, impl_->current_output) != 0) {
        impl_->failed = true;
        return -1;
    }
    if (!impl_->options.overlap_output && impl_->FlushPending() != 0) {
        return -1;
    }

    const double kernel_t0 = GetTime();
    int flush_ret = 0;
#ifdef PLATFORM_SUNWAY
    __real_athread_spawn(impl_->op->kernel_entry(),
                         impl_->op->kernel_arguments(), 1);
    flush_ret = impl_->options.overlap_output
        ? impl_->FlushPending() : 0;
    athread_join();
#else
    impl_->failed = true;
    return -1;
#endif
    const double kernel_wall = GetTime() - kernel_t0;
    impl_->timing.kernel += kernel_wall;
    impl_->op->ObserveKernel(kernel_wall, input.count);
    if (flush_ret != 0 || impl_->op->Validate(input.count) != 0) {
        impl_->failed = true;
        return -1;
    }
    impl_->timing.input_blocks += static_cast<long long>(input.count);
    impl_->timing.output_blocks += static_cast<long long>(input.count);
    impl_->timing.output_bytes += impl_->op->OutputBytes(input.count);
    impl_->timing.batch_count++;
    impl_->pending_output = impl_->current_output;
    impl_->pending_count = input.count;
    impl_->current_output = impl_->current_output == &impl_->output_a
        ? &impl_->output_b : &impl_->output_a;
    impl_->current_scratch = impl_->current_scratch == &impl_->scratch_a
        ? &impl_->scratch_b : &impl_->scratch_a;
    return 0;
}

int CpeWritePipelineSession::Flush() {
    if (!impl_ || !impl_->initialized || impl_->failed || impl_->finished) {
        return -1;
    }
    return impl_->FlushPending();
}

int CpeWritePipelineSession::Finish() {
    if (!impl_ || !impl_->initialized || impl_->failed || impl_->finished ||
        impl_->FlushPending() != 0) return -1;
    const int ret = impl_->op->Finish();
    impl_->timing.total = GetTime() - impl_->started_at;
    impl_->op->Shutdown();
    impl_->op_initialized = false;
    impl_->finished = true;
    return ret;
}

const CpeWritePipelineTiming &CpeWritePipelineSession::timing() const {
    return impl_->timing;
}

int RunCpeWritePipeline(
        UncompressedBgzfSource *source,
        CompressedBgzfBatchPostProcessor *post_processor,
        CpeWriteBatchOperator *op,
        CpeWritePipelineTiming *timing,
        const CpeWritePipelineOptions &options) {
    if (!source || !post_processor || !op || op->batch_capacity() == 0) return -1;

    const double total_t0 = GetTime();
    CpeWritePipelineTiming local;
    CpeWritePipelineSession session;
    BgzfBlockBatch input_a;
    BgzfBlockBatch input_b;
    BgzfBlockBatch *current_input = &input_a;
    BgzfBlockBatch *next_input = &input_b;
    size_t active_blocks = 0;
    int ret = -1;

    if (session.Initialize(post_processor, op, options) != 0 ||
        input_a.Allocate(op->batch_capacity()) != 0 ||
        input_b.Allocate(op->batch_capacity()) != 0) goto cleanup;

    {
        const double t0 = GetTime();
        const int source_ret = source->Fill(current_input, &active_blocks);
        local.source += GetTime() - t0;
        if (source_ret != 0 || active_blocks > op->batch_capacity()) goto cleanup;
    }

    while (active_blocks > 0) {
        if (session.Submit(CpeWriteBatchInput(current_input, active_blocks)) != 0) {
            goto cleanup;
        }
        size_t next_blocks = 0;
        const double t0 = GetTime();
        const int source_ret = source->Fill(next_input, &next_blocks);
        local.source += GetTime() - t0;
        if (source_ret != 0 || next_blocks > op->batch_capacity()) goto cleanup;
        std::swap(current_input, next_input);
        active_blocks = next_blocks;
    }
    ret = session.Finish();

cleanup:
    {
        const double source_time = local.source;
        local = session.timing();
        local.source = source_time;
        local.total = GetTime() - total_t0;
    }
    if (timing) *timing = local;
    return ret;
}

} // namespace cpe
} // namespace swbam
