#include "swbam/cpe_record_write_adapter.h"

#include <cstdio>
#include <cstring>

extern "C" void slave_mpi_compressfunc();

namespace swbam {
namespace cpe {

namespace {
const size_t kBatchBlocks = 64;

void ResetPara(Comp_Para *para, size_t block_id) {
    memset(para, 0, sizeof(*para));
    para->block_id = static_cast<int>(block_id);
    para->status = -1;
    para->compress_level = 1;
}

void AccumulateDetail(const Comp_Para *paras, size_t count,
                      double wall, CpeRecordWriteMetrics *metrics) {
    if (count == 0 || wall <= 0.0) return;
    const Comp_Para *critical = nullptr;
    uint64_t total = 0;
    for (size_t i = 0; i < count; ++i) {
        if (paras[i].compress_total_cycles >= total) {
            total = paras[i].compress_total_cycles;
            critical = &paras[i];
        }
    }
    if (!critical || total == 0) {
        metrics->other += wall;
        return;
    }
    const double scale = wall / static_cast<double>(total);
    const double serialize = scale * critical->compress_serialize_cycles;
    const double alloc = scale * critical->compress_alloc_cycles;
    const double deflate = scale * critical->compress_deflate_cycles;
    const double footer = scale * critical->compress_footer_cycles;
    double other = wall - serialize - alloc - deflate - footer;
    if (other < 0.0) other = 0.0;
    metrics->serialize += serialize;
    metrics->alloc += alloc;
    metrics->deflate += deflate;
    metrics->footer += footer;
    metrics->other += other;
}

class RecordCompressOperator : public CpeWriteBatchOperator {
public:
    explicit RecordCompressOperator(CpeRecordWriteMetrics *metrics)
        : level_(1), metrics_(metrics) {
        for (size_t i = 0; i < kBatchBlocks; ++i) ResetPara(&paras_[i], i);
    }

    void set_level(int level) { level_ = level; }
    const char *name() const { return "bam-record-compress"; }
    size_t batch_capacity() const { return kBatchBlocks; }
    bool needs_scratch() const { return true; }
    int Initialize() {
        return (level_ == 0 || level_ == 1 || level_ == 6) ? 0 : -1;
    }
    void Shutdown() {}

    int Prepare(const CpeWriteBatchInput &input,
                BgzfBlockBatch *scratch,
                BgzfBlockBatch *output) {
        if (input.kind != CpeWriteBatchInput::kBamRecords ||
            !input.records || !scratch || !output) return -1;
        for (size_t i = 0; i < kBatchBlocks; ++i) {
            ResetPara(&paras_[i], i);
            if (i >= input.count) continue;
            const BamRecordPackBlock &block = input.records[i];
            if (!block.records || block.n_records <= 0 ||
                block.total_len > BGZF_BLOCK_SIZE) return -1;
            Comp_Para &para = paras_[i];
            para.input_records = block.records;
            para.n_records = block.n_records;
            para.un_comp_block = &scratch->blocks()[i];
            para.un_comp_size = static_cast<int>(block.total_len);
            para.output_block = &output->blocks()[i];
            para.compress_level = level_;
            para.status = 0;
        }
        return 0;
    }

    void *kernel_entry() const { return (void *)slave_mpi_compressfunc; }
    void *kernel_arguments() { return paras_; }
    void ObserveKernel(double wall, size_t count) {
        AccumulateDetail(paras_, count, wall, metrics_);
    }
    int Validate(size_t count) const {
        for (size_t i = 0; i < count; ++i) {
            const Comp_Para &para = paras_[i];
            if (para.status != 0 || !para.output_block ||
                para.output_block->length == 0 ||
                para.output_block->length > BGZF_MAX_BLOCK_SIZE) {
                fprintf(stderr,
                        "ERROR: BAM record compression failed on block %zu with status %d.\n",
                        i, para.status);
                return -1;
            }
        }
        return 0;
    }
    long long OutputBytes(size_t count) const {
        long long total = 0;
        for (size_t i = 0; i < count; ++i) total += paras_[i].output_size;
        return total;
    }
    int Finish() { return 0; }

private:
    int level_;
    CpeRecordWriteMetrics *metrics_;
    Comp_Para paras_[kBatchBlocks];
};

class SinkPostProcessor : public CompressedBgzfBatchPostProcessor {
public:
    SinkPostProcessor() : sink_(nullptr) {}
    void set_sink(RankBodySink *sink) { sink_ = sink; }

    int PostProcessCompressedBatch(const bam_block *blocks, size_t count) {
        for (size_t i = 0; i < count; ++i) {
            if (!sink_ || blocks[i].length == 0 ||
                blocks[i].length > BGZF_MAX_BLOCK_SIZE ||
                sink_->Append(blocks[i].data, blocks[i].length) != 0) {
                fprintf(stderr,
                        "ERROR: BAM record writer failed to append BGZF block.\n");
                return -1;
            }
        }
        return 0;
    }

private:
    RankBodySink *sink_;
};
} // namespace

CpeRecordWriteMetrics::CpeRecordWriteMetrics()
    : bgzf_blocks(0), compress(0.0), serialize(0.0), alloc(0.0),
      deflate(0.0), footer(0.0), other(0.0), write(0.0) {}

class CpeRecordWriteSession::Impl {
public:
    Impl() : op(&metrics), initialized(false) {}

    CpeRecordWriteMetrics metrics;
    RecordCompressOperator op;
    SinkPostProcessor post_processor;
    CpeWritePipelineSession pipeline;
    bool initialized;
};

CpeRecordWriteSession::CpeRecordWriteSession() : impl_(new Impl) {}
CpeRecordWriteSession::~CpeRecordWriteSession() { delete impl_; }

int CpeRecordWriteSession::Initialize(RankBodySink *sink, int level) {
    if (!impl_ || !sink || impl_->initialized ||
        (level != 0 && level != 1 && level != 6)) return -1;
    impl_->metrics = CpeRecordWriteMetrics();
    impl_->op.set_level(level);
    impl_->post_processor.set_sink(sink);
    if (impl_->pipeline.Initialize(&impl_->post_processor, &impl_->op) != 0) {
        return -1;
    }
    impl_->initialized = true;
    return 0;
}

int CpeRecordWriteSession::Submit(const BamRecordPackBlock *blocks,
                                  size_t count) {
    if (!impl_ || !impl_->initialized || count > kBatchBlocks ||
        (count && !blocks)) return -1;
    if (count == 0) return impl_->pipeline.Flush();
    return impl_->pipeline.Submit(CpeWriteBatchInput(blocks, count));
}

int CpeRecordWriteSession::Finish() {
    return impl_ && impl_->initialized ? impl_->pipeline.Finish() : -1;
}

const CpeRecordWriteMetrics &CpeRecordWriteSession::metrics() const {
    const CpeWritePipelineTiming &timing = impl_->pipeline.timing();
    impl_->metrics.bgzf_blocks = timing.output_blocks;
    impl_->metrics.compress = timing.kernel;
    impl_->metrics.write = timing.post_process;
    return impl_->metrics;
}

} // namespace cpe
} // namespace swbam
