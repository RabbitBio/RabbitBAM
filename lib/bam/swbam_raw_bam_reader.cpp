#include "swbam/raw_bam.h"

#include <cstring>
#include <vector>

namespace swbam {
namespace cpe {
namespace {

class RawBamBatchPostProcessorAdapter : public DecodedBgzfBatchPostProcessor {
public:
    RawBamBatchPostProcessorAdapter(RawBamBatchPostProcessor *post_processor,
                       ComposableRawBamMetrics *metrics)
        : post_processor_(post_processor), metrics_(metrics) {
        records_.reserve(64 * 1024);
    }

    int PostProcessDecodedBatch(const bam_block *blocks, size_t count) {
        if (!post_processor_ || (!blocks && count != 0)) return -1;
        records_.clear();

        long long encoded_bytes = 0;
        for (size_t block_index = 0; block_index < count; ++block_index) {
            const bam_block &block = blocks[block_index];
            uint32_t offset = 0;
            while (offset < block.length) {
                const uint32_t remaining = block.length - offset;
                if (remaining < sizeof(int32_t)) return -1;

                int32_t block_size = 0;
                memcpy(&block_size, block.data + offset, sizeof(block_size));
                if (block_size < 32) return -1;
                const uint64_t encoded_size =
                    (uint64_t)(uint32_t)block_size + sizeof(int32_t);
                if (encoded_size > remaining || encoded_size > UINT32_MAX) {
                    return -1;
                }

                RawBamRecordView view;
                view.encoded = block.data + offset;
                view.encoded_size = (uint32_t)encoded_size;
                view.block_size = (uint32_t)block_size;
                view.block_offset = offset;
                view.block_index = (uint32_t)block_index;
                records_.push_back(view);
                encoded_bytes += (long long)encoded_size;
                offset += (uint32_t)encoded_size;
            }
        }

        if (post_processor_->PostProcessRawBatch(
                records_.empty() ? nullptr : records_.data(),
                records_.size()) != 0) {
            return -1;
        }
        if (metrics_) {
            metrics_->records += (long long)records_.size();
            metrics_->encoded_bytes += encoded_bytes;
        }
        return 0;
    }

private:
    RawBamBatchPostProcessor *post_processor_;
    ComposableRawBamMetrics *metrics_;
    std::vector<RawBamRecordView> records_;
};

} // namespace

int RunComposableRawBamPipeline(
        const BamInputBackend &input,
        const BgzfBlockSpan *spans,
        size_t span_count,
        RawBamBatchPostProcessor *post_processor,
        CpeReadPipelineTiming *timing,
        ComposableDecodeMetrics *decode_metrics,
        ComposableRawBamMetrics *raw_metrics,
        const CpeReadPipelineOptions &options) {
    if (!post_processor) return -1;
    if (raw_metrics) *raw_metrics = ComposableRawBamMetrics();
    RawBamBatchPostProcessorAdapter adapter(post_processor, raw_metrics);
    return RunComposableDecodePipeline(
        input, spans, span_count, &adapter,
        timing, decode_metrics, options);
}

} // namespace cpe
} // namespace swbam
