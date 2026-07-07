#include "swbam/generic_decode.h"
#include "swbam/cpe_codec.h"

#include <cstdio>
#include <cstring>

extern "C" void slave_swbam_decode_bgzf();

namespace swbam {
namespace cpe {
namespace {

const size_t kGenericDecodeBatch = 64;

class GenericDecodeOperator : public CpeBatchOperator {
public:
    GenericDecodeOperator(DecodedBgzfConsumer *consumer,
                          GenericDecodeMetrics *metrics)
        : consumer_(consumer), metrics_(metrics), decoded_(nullptr) {
        memset(paras_, 0, sizeof(paras_));
    }

    const char *name() const { return "generic-decode"; }
    size_t batch_capacity() const { return kGenericDecodeBatch; }

    int Initialize() {
        if (!consumer_) return -1;
        if (metrics_) *metrics_ = GenericDecodeMetrics();
        return 0;
    }
    void Shutdown() { decoded_ = nullptr; }

    int Prepare(const BgzfBlockBatch &compressed,
                BgzfBlockBatch *decoded,
                size_t active_blocks) {
        decoded_ = decoded;
        for (size_t i = 0; i < kGenericDecodeBatch; ++i) {
            SwbamCpeDecodePara &para = paras_[i];
            para.alloc_cycles = 0;
            para.inflate_cycles = 0;
            para.crc_cycles = 0;
            para.total_cycles = 0;
            if (i < active_blocks) {
                para.compressed = const_cast<bam_block *>(
                    &compressed.blocks()[i]);
                para.decoded = &decoded->blocks()[i];
                para.status = 0;
            } else {
                para.compressed = nullptr;
                para.decoded = nullptr;
                para.status = -1;
            }
        }
        return 0;
    }

    void *kernel_entry() const {
        return (void *)slave_swbam_decode_bgzf;
    }
    void *kernel_arguments() { return paras_; }

    void ObserveKernel(double wall_seconds, size_t active_blocks) {
        if (!metrics_ || active_blocks == 0 || wall_seconds <= 0.0) return;
        const SwbamCpeDecodePara *critical = nullptr;
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
        const double inflate = scale * critical->inflate_cycles;
        const double crc = scale * critical->crc_cycles;
        double other = wall_seconds - alloc - inflate - crc;
        if (other < 0.0) other = 0.0;
        metrics_->alloc += alloc;
        metrics_->inflate += inflate;
        metrics_->crc += crc;
        metrics_->other += other;
    }

    int Validate(size_t active_blocks) const {
        for (size_t i = 0; i < active_blocks; ++i) {
            if (paras_[i].status != 0) {
                fprintf(stderr,
                        "ERROR: generic BGZF decode failed on block %zu with status %d.\n",
                        i, paras_[i].status);
                return -1;
            }
        }
        return 0;
    }

    int Consume(size_t active_blocks, long long *records_processed) {
        if (!decoded_) return -1;
        if (metrics_) {
            metrics_->decoded_blocks += (long long)active_blocks;
            for (size_t i = 0; i < active_blocks; ++i) {
                metrics_->decoded_bytes += decoded_->blocks()[i].length;
            }
        }
        if (records_processed) *records_processed = 0;
        return consumer_->ConsumeDecoded(decoded_->blocks(), active_blocks);
    }

    int Finish() { return 0; }

private:
    DecodedBgzfConsumer *consumer_;
    GenericDecodeMetrics *metrics_;
    BgzfBlockBatch *decoded_;
    SwbamCpeDecodePara paras_[kGenericDecodeBatch];
};

} // namespace

GenericDecodeMetrics::GenericDecodeMetrics()
    : decoded_blocks(0), decoded_bytes(0), alloc(0.0),
      inflate(0.0), crc(0.0), other(0.0) {}

int RunGenericDecodePipeline(
        const BamInputBackend &input,
        const BgzfBlockSpan *spans,
        size_t span_count,
        DecodedBgzfConsumer *consumer,
        CpeReadPipelineTiming *timing,
        GenericDecodeMetrics *metrics,
        const CpeReadPipelineOptions &options) {
    GenericDecodeOperator op(consumer, metrics);
    return RunCpeReadPipeline(
        input, spans, span_count, &op, timing, options);
}

} // namespace cpe
} // namespace swbam
