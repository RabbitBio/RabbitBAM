#include "swbam/operators/bam_transform.h"

#include <cstdio>
#include <cstring>
#include <vector>

extern "C" {
void slave_mpi_decompress_filterfunc();
void slave_mpi_decompress_bam2bam_passthrough();
}

namespace swbam {
namespace operators {
namespace {

const int kBatchBlocks = 64;
const int kRecordsPerBlock = static_cast<int>(MPI_RECORDS_PER_BLOCK);

class RecordWorkspace {
public:
    RecordWorkspace()
        : records_(nullptr), data_(nullptr), ptrs_(nullptr), lengths_(nullptr) {}
    ~RecordWorkspace() { Release(); }

    int Allocate() {
        const size_t count = static_cast<size_t>(kBatchBlocks) * kRecordsPerBlock;
        records_ = reinterpret_cast<bam1_t *>(
            aligned_alloc_custom(64, count * sizeof(bam1_t)));
        data_ = aligned_alloc_custom(
            64, static_cast<size_t>(kBatchBlocks) * MPI_BAM_BLOCK_ARENA_SIZE);
        ptrs_ = reinterpret_cast<bam1_t **>(
            aligned_alloc_custom(64, count * sizeof(bam1_t *)));
        lengths_ = reinterpret_cast<uint32_t *>(
            aligned_alloc_custom(64, count * sizeof(uint32_t)));
        if (!records_ || !data_ || !ptrs_ || !lengths_) return -1;
        memset(records_, 0, count * sizeof(bam1_t));
        memset(data_, 0, static_cast<size_t>(kBatchBlocks) *
                         MPI_BAM_BLOCK_ARENA_SIZE);
        memset(lengths_, 0, count * sizeof(uint32_t));
        for (size_t i = 0; i < count; ++i) {
            records_[i].mempolicy = BAM_USER_OWNS_DATA;
            ptrs_[i] = &records_[i];
        }
        return 0;
    }

    void Release() {
        if (records_) aligned_free_custom(reinterpret_cast<unsigned char *>(records_));
        if (data_) aligned_free_custom(data_);
        if (ptrs_) aligned_free_custom(reinterpret_cast<unsigned char *>(ptrs_));
        if (lengths_) aligned_free_custom(reinterpret_cast<unsigned char *>(lengths_));
        records_ = nullptr;
        data_ = nullptr;
        ptrs_ = nullptr;
        lengths_ = nullptr;
    }

    bam1_t *records(int b) { return records_ + static_cast<size_t>(b) * kRecordsPerBlock; }
    unsigned char *data(int b) {
        return data_ + static_cast<size_t>(b) * MPI_BAM_BLOCK_ARENA_SIZE;
    }
    bam1_t **ptrs(int b) { return ptrs_ + static_cast<size_t>(b) * kRecordsPerBlock; }
    uint32_t *lengths(int b) {
        return lengths_ + static_cast<size_t>(b) * kRecordsPerBlock;
    }

private:
    bam1_t *records_;
    unsigned char *data_;
    bam1_t **ptrs_;
    uint32_t *lengths_;
};

class PackWorkspace {
public:
    PackWorkspace() : active_blocks(0), total_records(0),
        current_begin(nullptr), current_records(0), current_len(0) {
        records.resize(static_cast<size_t>(kBatchBlocks) * kRecordsPerBlock);
        Reset();
    }

    void Reset() {
        active_blocks = 0;
        total_records = 0;
        current_begin = records.data();
        current_records = 0;
        current_len = 0;
    }

    int Seal() {
        if (!current_records) return 0;
        if (active_blocks >= kBatchBlocks) return -3;
        cpe::BamRecordPackBlock &plan = plans[active_blocks++];
        plan.records = current_begin;
        plan.n_records = current_records;
        plan.total_len = current_len;
        current_begin = records.data() + total_records;
        current_records = 0;
        current_len = 0;
        return 0;
    }

    int AppendRange(bam1_t **input, const uint32_t *lengths,
                    int count, uint32_t total_len) {
        if (count <= 0) return 0;
        if (total_len > 0 && total_len <= BGZF_BLOCK_SIZE &&
            (!current_records || current_len + total_len <= BGZF_BLOCK_SIZE)) {
            if (total_records + count > static_cast<int>(records.size())) return -2;
            if (!current_records) current_begin = records.data() + total_records;
            memcpy(records.data() + total_records, input,
                   static_cast<size_t>(count) * sizeof(bam1_t *));
            total_records += count;
            current_records += count;
            current_len += total_len;
            return 0;
        }
        int offset = 0;
        while (offset < count) {
            if (!current_records) current_begin = records.data() + total_records;
            const uint32_t room = BGZF_BLOCK_SIZE - current_len;
            int take = 0;
            uint32_t take_len = 0;
            while (offset + take < count) {
                const uint32_t packed_len = lengths[offset + take] + 4;
                if (packed_len > BGZF_BLOCK_SIZE) return -1;
                if (take_len + packed_len > room) break;
                take_len += packed_len;
                ++take;
            }
            if (!take) {
                if (Seal() != 0) return -3;
                continue;
            }
            if (total_records + take > static_cast<int>(records.size())) return -2;
            memcpy(records.data() + total_records, input + offset,
                   static_cast<size_t>(take) * sizeof(bam1_t *));
            total_records += take;
            current_records += take;
            current_len += take_len;
            offset += take;
        }
        return 0;
    }

    int BuildPassthrough(const Bam2BamPara *paras,
                         const BgzfBlockBatch &decoded, size_t count) {
        for (size_t b = 0; b < count; ++b) {
            const int n = paras[b].n_kept_records;
            if (n <= 0) continue;
            uint32_t total_len = paras[b].kept_total_len;
            if (!total_len) total_len = decoded.blocks()[b].length;
            if (total_len > BGZF_BLOCK_SIZE) {
                const int ret = AppendRange(paras[b].output_records,
                                            paras[b].bam_lens, n, total_len);
                if (ret != 0) return ret;
                continue;
            }
            if (current_records && Seal() != 0) return -3;
            if (active_blocks >= kBatchBlocks) return -3;
            cpe::BamRecordPackBlock &plan = plans[active_blocks++];
            plan.records = paras[b].output_records;
            plan.n_records = n;
            plan.total_len = total_len;
            total_records += n;
        }
        return Seal();
    }

    std::vector<bam1_t *> records;
    cpe::BamRecordPackBlock plans[kBatchBlocks];
    int active_blocks;
    int total_records;
    bam1_t **current_begin;
    int current_records;
    uint32_t current_len;
};

void AccumulateDecodeDetail(const Bam2BamPara *paras, size_t count,
                            double wall, BamTransformMetrics *metrics) {
    if (!metrics || !count || wall <= 0.0) return;
    const Bam2BamPara *critical = nullptr;
    uint64_t total = 0;
    for (size_t b = 0; b < count; ++b) {
        if (paras[b].decomp_total_cycles >= total) {
            total = paras[b].decomp_total_cycles;
            critical = &paras[b];
        }
    }
    if (!critical || !total) {
        metrics->t_decomp_other += wall;
        return;
    }
    const double scale = wall / static_cast<double>(total);
    const double alloc = scale * critical->decomp_alloc_cycles;
    const double inflate = scale * critical->decomp_inflate_cycles;
    const double crc = scale * critical->decomp_crc_cycles;
    const double parse = scale * critical->decomp_parse_cycles;
    double other = wall - alloc - inflate - crc - parse;
    if (other < 0.0) other = 0.0;
    metrics->t_decomp_alloc += alloc;
    metrics->t_decomp_inflate += inflate;
    metrics->t_decomp_crc += crc;
    metrics->t_decomp_parse += parse;
    metrics->t_decomp_other += other;
}

class MemoryBatchSource : public cpe::CpeReadBatchSource {
public:
    explicit MemoryBatchSource(MemReader *reader) : reader_(reader) {}

    int ReadNext(BgzfBlockBatch *batch, size_t *count) {
        if (!reader_ || !batch || !count || reader_->pos > reader_->size ||
            (reader_->size && !reader_->base)) return -1;
        *count = 0;
        while (*count < batch->capacity() && reader_->pos < reader_->size) {
            const size_t pos = reader_->pos;
            if (reader_->size - pos < BLOCK_HEADER_LENGTH) return -1;
            const unsigned char *src =
                reinterpret_cast<const unsigned char *>(reader_->base + pos);
            const size_t len = static_cast<size_t>(src[16] | (src[17] << 8)) + 1;
            if (len > BGZF_MAX_BLOCK_SIZE || len > reader_->size - pos ||
                len < BLOCK_HEADER_LENGTH) return -1;
            reader_->pos += len;
            if (len == 28) break; // BGZF EOF is not a data block.
            bam_block &block = batch->blocks()[*count];
            memcpy(block.data, src, len);
            block.length = static_cast<unsigned int>(len);
            block.pos = 0;
            block.errcode = 0;
            block.block_address = static_cast<int64_t>(pos);
            block.block_id = static_cast<int>(*count);
            ++*count;
        }
        return 0;
    }

private:
    MemReader *reader_;
};

class BamTransformOperator : public cpe::CpeBatchOperator {
public:
    BamTransformOperator(RankBodySink *sink, const BamFilterOptions &filter,
                         int level, BamTransformMetrics *metrics)
        : sink_(sink), filter_(filter), level_(level), metrics_(metrics),
          no_filter_(bam_filter_is_noop(filter)), decoded_(nullptr) {
        memset(paras_, 0, sizeof(paras_));
    }

    const char *name() const { return "bam-transform"; }
    size_t batch_capacity() const { return kBatchBlocks; }

    int Initialize() {
        if (!sink_ || !metrics_ || workspace_.Allocate() != 0 ||
            writer_.Initialize(sink_, level_) != 0) {
            fprintf(stderr, "ERROR: failed to allocate BAM transform workspace.\n");
            return -1;
        }
        return 0;
    }

    void Shutdown() { workspace_.Release(); }

    int Prepare(const BgzfBlockBatch &compressed,
                BgzfBlockBatch *decoded, size_t count) {
        if (!decoded || count > kBatchBlocks) return -1;
        decoded_ = decoded;
        for (int b = 0; b < kBatchBlocks; ++b) {
            Bam2BamPara &para = paras_[b];
            para.filter = filter_;
            para.block_id = b;
            para.output_records = workspace_.ptrs(b);
            para.record_base = workspace_.records(b);
            para.bam_lens = workspace_.lengths(b);
            para.record_capacity = kRecordsPerBlock;
            para.data_arena = workspace_.data(b);
            para.data_arena_capacity = MPI_BAM_BLOCK_ARENA_SIZE;
            para.data_arena_used = 0;
            para.record_index = 0;
            para.actual_value = 0;
            para.limit_value = 0;
            para.limit_id = BOUNDS_LIMIT_NONE;
            para.n_total_records = 0;
            para.n_kept_records = 0;
            para.kept_total_len = 0;
            para.decomp_alloc_cycles = 0;
            para.decomp_inflate_cycles = 0;
            para.decomp_crc_cycles = 0;
            para.decomp_parse_cycles = 0;
            para.decomp_total_cycles = 0;
            if (b < static_cast<int>(count)) {
                para.input_block = const_cast<bam_block *>(&compressed.blocks()[b]);
                para.un_comp_block = &decoded->blocks()[b];
                para.status = 0;
            } else {
                para.input_block = nullptr;
                para.un_comp_block = nullptr;
                para.status = -1;
            }
        }
        return 0;
    }

    void *kernel_entry() const {
        return reinterpret_cast<void *>(
            no_filter_ ? slave_mpi_decompress_bam2bam_passthrough
                       : slave_mpi_decompress_filterfunc);
    }
    void *kernel_arguments() { return paras_; }

    void ObserveKernel(double wall, size_t count) {
        metrics_->t_decomp_filter += wall;
        AccumulateDecodeDetail(paras_, count, wall, metrics_);
    }

    int Validate(size_t count) const {
        for (size_t b = 0; b < count; ++b) {
            if (paras_[b].status == 0) continue;
            if (paras_[b].status == -3) {
                fprintf(stderr,
                        "ERROR: MPI bam2bam capacity exceeded on input block %zu. "
                        "limit_id=%d limit=%lld actual=%lld record=%d.\n",
                        b, paras_[b].limit_id, paras_[b].limit_value,
                        paras_[b].actual_value, paras_[b].record_index);
            } else {
                fprintf(stderr,
                        "ERROR: MPI bam2bam decompress/filter failed on input block %zu "
                        "with status %d.\n", b, paras_[b].status);
            }
            return -1;
        }
        return 0;
    }

    int PostProcessBatch(size_t count, long long *records_processed) {
        pack_.Reset();
        const double pack_t0 = GetTime();
        long long batch_records = 0;
        for (size_t b = 0; b < count; ++b) {
            const Bam2BamPara &para = paras_[b];
            metrics_->total_records += para.n_total_records;
            metrics_->kept_records += para.n_kept_records;
            metrics_->dropped_records += para.n_total_records - para.n_kept_records;
            metrics_->pack_records += para.n_kept_records;
            batch_records += para.n_total_records;
        }
        int ret = 0;
        if (no_filter_) {
            ret = pack_.BuildPassthrough(paras_, *decoded_, count);
        } else {
            for (size_t b = 0; b < count; ++b) {
                ret = pack_.AppendRange(paras_[b].output_records,
                                        paras_[b].bam_lens,
                                        paras_[b].n_kept_records,
                                        paras_[b].kept_total_len);
                if (ret != 0) break;
            }
            if (ret == 0) ret = pack_.Seal();
        }
        if (ret != 0) {
            if (ret == -1) {
                fprintf(stderr, "ERROR: MPI bam2bam pack encountered an oversized BAM record.\n");
            } else if (ret == -2) {
                fprintf(stderr, "ERROR: MPI bam2bam pack workspace capacity exceeded.\n");
            } else {
                fprintf(stderr, "ERROR: MPI bam2bam produced more than %d output blocks from one input group.\n",
                        kBatchBlocks);
            }
            return -1;
        }
        metrics_->t_pack += GetTime() - pack_t0;
        metrics_->group_count++;
        if (records_processed) *records_processed = batch_records;
        return writer_.Submit(pack_.plans, static_cast<size_t>(pack_.active_blocks));
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
    RankBodySink *sink_;
    BamFilterOptions filter_;
    int level_;
    BamTransformMetrics *metrics_;
    bool no_filter_;
    BgzfBlockBatch *decoded_;
    Bam2BamPara paras_[kBatchBlocks];
    RecordWorkspace workspace_;
    PackWorkspace pack_;
    cpe::CpeRecordWriteSession writer_;
};

int RunCore(cpe::CpeReadBatchSource *source,
            RankBodySink *sink, const BamFilterOptions &filter,
            int level, BamTransformMetrics *metrics) {
    if (!source || !sink || !metrics) return -1;
    *metrics = BamTransformMetrics();
    const double t0 = GetTime();
    BamTransformOperator op(sink, filter, level, metrics);
    cpe::CpeReadPipelineTiming read_timing;
    const int ret = cpe::RunCpeReadPipeline(source, &op, &read_timing);
    metrics->input_blocks = read_timing.input_blocks;
    metrics->t_read = read_timing.read;
    metrics->t_optimized_total = GetTime() - t0;
    return ret;
}

class SpanSource : public cpe::CpeReadBatchSource {
public:
    SpanSource(const BamInputBackend &input,
               const BgzfBlockSpan *spans, size_t count)
        : reader_(&input, spans, count) {}
    int ReadNext(BgzfBlockBatch *batch, size_t *count) {
        return reader_.ReadNext(batch, count);
    }
private:
    BgzfSpanBatchReader reader_;
};

} // namespace

int RunBamTransformPipeline(const BamInputBackend &input,
                            const BgzfBlockSpan *spans, size_t count,
                            RankBodySink *sink, const BamFilterOptions &filter,
                            int level, BamTransformMetrics *metrics) {
    if (count && !spans) return -1;
    SpanSource source(input, spans, count);
    return RunCore(&source, sink, filter, level, metrics);
}

int RunBamTransformPipeline(MemReader *reader, RankBodySink *sink,
                            const BamFilterOptions &filter, int level,
                            BamTransformMetrics *metrics) {
    if (!reader) return -1;
    MemoryBatchSource source(reader);
    return RunCore(&source, sink, filter, level, metrics);
}

} // namespace operators
} // namespace swbam
