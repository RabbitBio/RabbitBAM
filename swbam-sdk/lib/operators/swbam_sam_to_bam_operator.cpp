#include "swbam/operators/sam_to_bam.h"

#include "swbam/cpe_record_write_adapter.h"
#include "swbam/cpe_sam_read_pipeline.h"

#include <cstdio>
#include <cstring>

namespace swbam {
namespace operators {
namespace {

class SamPackWorkspace {
public:
    SamPackWorkspace()
        : records_(nullptr), plans_(nullptr), active_blocks_(0),
          total_records_(0), current_begin_(nullptr), current_records_(0),
          current_len_(0) {}
    ~SamPackWorkspace() {
        if (records_) aligned_free_custom(reinterpret_cast<unsigned char *>(records_));
        if (plans_) aligned_free_custom(reinterpret_cast<unsigned char *>(plans_));
    }

    int Initialize() {
        records_ = reinterpret_cast<bam1_t **>(aligned_alloc_custom(
            64, static_cast<size_t>(FUSED_SAM2BAM_BAM_POOL_SIZE) * sizeof(bam1_t *)));
        plans_ = reinterpret_cast<cpe::BamRecordPackBlock *>(aligned_alloc_custom(
            64, 64 * sizeof(cpe::BamRecordPackBlock)));
        if (!records_ || !plans_) return -1;
        Reset();
        return 0;
    }

    void Reset() {
        active_blocks_ = total_records_ = current_records_ = 0;
        current_begin_ = records_;
        current_len_ = 0;
    }

    int Seal() {
        if (!current_records_) return 0;
        if (active_blocks_ >= 64) return -3;
        cpe::BamRecordPackBlock &plan = plans_[active_blocks_++];
        plan.records = current_begin_;
        plan.n_records = current_records_;
        plan.total_len = current_len_;
        current_begin_ = records_ + total_records_;
        current_records_ = 0;
        current_len_ = 0;
        return 0;
    }

    int Append(bam1_t *record, uint32_t bam_len) {
        const uint32_t packed_len = bam_len + 4;
        if (packed_len > BGZF_BLOCK_SIZE) return -1;
        if (current_records_ > 0 && current_len_ + packed_len > BGZF_BLOCK_SIZE)
            return 1;
        if (total_records_ >= FUSED_SAM2BAM_BAM_POOL_SIZE) return -2;
        if (!current_records_) current_begin_ = records_ + total_records_;
        records_[total_records_++] = record;
        ++current_records_;
        current_len_ += packed_len;
        return 0;
    }

    cpe::BamRecordPackBlock *plans() { return plans_; }
    int active_blocks() const { return active_blocks_; }

private:
    bam1_t **records_;
    cpe::BamRecordPackBlock *plans_;
    int active_blocks_;
    int total_records_;
    bam1_t **current_begin_;
    int current_records_;
    uint32_t current_len_;
};

class SamToBamPostProcessor : public cpe::SamParsedBatchPostProcessor {
public:
    explicit SamToBamPostProcessor(SamToBamMetrics *metrics)
        : metrics_(metrics) {}

    int Initialize(RankBodySink *sink, int level) {
        return pack_.Initialize() == 0 && writer_.Initialize(sink, level) == 0
            ? 0 : -1;
    }

    int FlushPendingOutput() { return writer_.Submit(nullptr, 0); }

    int PostProcessParsedBatch(const MpiSamParseChunk *chunks, size_t count) {
        double pack_t0 = GetTime();
        for (size_t i = 0; i < count; ++i) {
            const MpiSamParseChunk &chunk = chunks[i];
            for (int j = 0; j < chunk.count; ++j) {
                bam1_t *record = chunk.bams + j;
                const uint32_t length = chunk.bam_lens[j];
                ++metrics_->total_records;
                for (;;) {
                    const int ret = pack_.Append(record, length);
                    if (ret == 0) break;
                    if (ret == -1) {
                        fprintf(stderr, "ERROR: MPI sam2bam pack encountered an oversized BAM record.\n");
                        return -1;
                    }
                    if (ret == -2) {
                        fprintf(stderr, "ERROR: MPI sam2bam pack workspace capacity exceeded.\n");
                        return -1;
                    }
                    if (pack_.Seal() != 0) return -1;
                    if (pack_.active_blocks() == 64) {
                        metrics_->t_pack += GetTime() - pack_t0;
                        if (SubmitPlans() != 0) return -1;
                        pack_t0 = GetTime();
                    }
                }
            }
        }
        if (pack_.Seal() != 0) {
            fprintf(stderr, "ERROR: MPI sam2bam final output block plan capacity exceeded.\n");
            return -1;
        }
        metrics_->t_pack += GetTime() - pack_t0;
        return pack_.active_blocks() ? SubmitPlans() : 0;
    }

    int Finish() {
        if (writer_.Finish() != 0) return -1;
        const cpe::CpeRecordWriteMetrics &written = writer_.metrics();
        metrics_->bgzf_blocks = written.bgzf_blocks;
        metrics_->t_compress = written.compress;
        metrics_->t_compress_serialize = written.serialize;
        metrics_->t_compress_alloc = written.alloc;
        metrics_->t_compress_deflate = written.deflate;
        metrics_->t_compress_footer = written.footer;
        metrics_->t_compress_other = written.other;
        metrics_->t_write = written.write;
        return 0;
    }

private:
    int SubmitPlans() {
        if (writer_.Submit(pack_.plans(), static_cast<size_t>(pack_.active_blocks())) != 0)
            return -1;
        ++metrics_->compress_groups;
        pack_.Reset();
        return 0;
    }

    SamToBamMetrics *metrics_;
    SamPackWorkspace pack_;
    cpe::CpeRecordWriteSession writer_;
};

} // namespace

int RunSamToBamPipeline(MemReader *reader, RankBodySink *sink,
                        const sam_hdr_t *header, int compression_level,
                        SamToBamMetrics *metrics) {
    if (!reader || !sink || !header || !metrics) return -1;
    *metrics = SamToBamMetrics();
    const double t0 = GetTime();
    int ret = -1;
    {
        SamToBamPostProcessor post(metrics);
        cpe::CpeSamReadTiming read;
        if (post.Initialize(sink, compression_level) == 0 &&
            cpe::RunCpeSamReadPipeline(reader, header, &post, &read) == 0 &&
            post.Finish() == 0) ret = 0;
        metrics->input_chunks = read.input_chunks;
        metrics->chunk_groups = read.chunk_groups;
        metrics->parse_fast_records = read.parse_fast_records;
        metrics->parse_fallback_records = read.parse_fallback_records;
        metrics->t_split = read.split;
        metrics->t_copy_count = read.copy_count;
        metrics->t_parse = read.parse;
        metrics->t_parse_core = read.parse_core;
        metrics->t_parse_aux = read.parse_aux;
        metrics->t_parse_cg = read.parse_cg;
        metrics->t_parse_fallback = read.parse_fallback;
        metrics->t_parse_other = read.parse_other;
    }
    metrics->t_optimized_total = GetTime() - t0;
    return ret;
}

} // namespace operators
} // namespace swbam
