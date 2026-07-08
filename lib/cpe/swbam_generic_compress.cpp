#include "swbam/generic_compress.h"
#include "swbam/cpe_codec.h"

#include <cstdio>
#include <cstring>

extern "C" void slave_swbam_compress_bgzf();

namespace swbam {
namespace cpe {
namespace {

const size_t kGenericCompressBatch = 64;

class GenericCompressOperator : public CpeWriteBatchOperator {
public:
    GenericCompressOperator(int level, GenericCompressMetrics *metrics)
        : level_(level), metrics_(metrics) {
        memset(paras_, 0, sizeof(paras_));
    }

    const char *name() const { return "generic-compress"; }
    size_t batch_capacity() const { return kGenericCompressBatch; }
    int Initialize() {
        if (level_ != 0 && level_ != 1 && level_ != 6) return -1;
        if (metrics_) *metrics_ = GenericCompressMetrics();
        return 0;
    }
    void Shutdown() {}

    int Prepare(const BgzfBlockBatch &input,
                BgzfBlockBatch *output,
                size_t active_blocks) {
        for (size_t i = 0; i < kGenericCompressBatch; ++i) {
            SwbamCpeCompressPara &para = paras_[i];
            para.alloc_cycles = 0;
            para.deflate_cycles = 0;
            para.footer_cycles = 0;
            para.total_cycles = 0;
            para.output_size = 0;
            para.level = level_;
            if (i < active_blocks) {
                para.uncompressed = const_cast<bam_block *>(
                    &input.blocks()[i]);
                para.compressed = &output->blocks()[i];
                para.status = 0;
            } else {
                para.uncompressed = nullptr;
                para.compressed = nullptr;
                para.status = -1;
            }
        }
        return 0;
    }

    void *kernel_entry() const { return (void *)slave_swbam_compress_bgzf; }
    void *kernel_arguments() { return paras_; }

    void ObserveKernel(double wall_seconds, size_t active_blocks) {
        if (!metrics_ || active_blocks == 0 || wall_seconds <= 0.0) return;
        const SwbamCpeCompressPara *critical = nullptr;
        uint64_t critical_cycles = 0;
        for (size_t i = 0; i < active_blocks; ++i) {
            if (paras_[i].total_cycles >= critical_cycles) {
                critical_cycles = paras_[i].total_cycles;
                critical = &paras_[i];
            }
        }
        if (!critical || critical_cycles == 0) {
            metrics_->other += wall_seconds;
            return;
        }
        const double scale = wall_seconds / (double)critical_cycles;
        const double alloc = scale * critical->alloc_cycles;
        const double deflate = scale * critical->deflate_cycles;
        const double footer = scale * critical->footer_cycles;
        double other = wall_seconds - alloc - deflate - footer;
        if (other < 0.0) other = 0.0;
        metrics_->alloc += alloc;
        metrics_->deflate += deflate;
        metrics_->footer += footer;
        metrics_->other += other;
    }

    int Validate(size_t active_blocks) const {
        for (size_t i = 0; i < active_blocks; ++i) {
            if (paras_[i].status != 0 || paras_[i].output_size <= 0) {
                fprintf(stderr,
                        "ERROR: generic BGZF compression failed on block %zu with status %d.\n",
                        i, paras_[i].status);
                return -1;
            }
        }
        return 0;
    }

    long long OutputBytes(size_t active_blocks) const {
        long long total = 0;
        for (size_t i = 0; i < active_blocks; ++i) {
            total += paras_[i].output_size;
        }
        return total;
    }
    int Finish() { return 0; }

private:
    int level_;
    GenericCompressMetrics *metrics_;
    SwbamCpeCompressPara paras_[kGenericCompressBatch];
};

} // namespace

GenericCompressMetrics::GenericCompressMetrics()
    : alloc(0.0), deflate(0.0), footer(0.0), other(0.0) {}

int RunGenericCompressPipeline(
        UncompressedBgzfSource *source,
        CompressedBgzfConsumer *consumer,
        int compression_level,
        CpeWritePipelineTiming *timing,
        GenericCompressMetrics *metrics,
        const CpeWritePipelineOptions &options) {
    GenericCompressOperator op(compression_level, metrics);
    return RunCpeWritePipeline(source, consumer, &op, timing, options);
}

} // namespace cpe
} // namespace swbam
