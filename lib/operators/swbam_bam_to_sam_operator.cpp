#include "swbam/operators/bam_to_sam.h"

#include <cstdio>
#include <cstring>

extern "C" void slave_mpi_decompress_bam2bam_passthrough();

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
        const size_t n = static_cast<size_t>(kBatchBlocks) * kRecordsPerBlock;
        records_ = reinterpret_cast<bam1_t *>(
            aligned_alloc_custom(64, n * sizeof(bam1_t)));
        data_ = aligned_alloc_custom(
            64, static_cast<size_t>(kBatchBlocks) * MPI_BAM_BLOCK_ARENA_SIZE);
        ptrs_ = reinterpret_cast<bam1_t **>(
            aligned_alloc_custom(64, n * sizeof(bam1_t *)));
        lengths_ = reinterpret_cast<uint32_t *>(
            aligned_alloc_custom(64, n * sizeof(uint32_t)));
        if (!records_ || !data_ || !ptrs_ || !lengths_) return -1;
        memset(records_, 0, n * sizeof(bam1_t));
        memset(data_, 0, static_cast<size_t>(kBatchBlocks) *
                         MPI_BAM_BLOCK_ARENA_SIZE);
        memset(lengths_, 0, n * sizeof(uint32_t));
        for (size_t i = 0; i < n; ++i) {
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
            const unsigned char *src = reinterpret_cast<const unsigned char *>(
                reader_->base + pos);
            const size_t len = static_cast<size_t>(src[16] | (src[17] << 8)) + 1;
            if (len < BLOCK_HEADER_LENGTH || len > BGZF_MAX_BLOCK_SIZE ||
                len > reader_->size - pos) return -1;
            reader_->pos += len;
            if (len == 28) break;
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

void AccumulateDecodeDetail(const Bam2BamPara *paras, size_t count,
                            double wall, BamToSamMetrics *metrics) {
    if (!count || wall <= 0.0) return;
    const Bam2BamPara *critical = nullptr;
    uint64_t cycles = 0;
    for (size_t b = 0; b < count; ++b) {
        if (paras[b].decomp_total_cycles >= cycles) {
            cycles = paras[b].decomp_total_cycles;
            critical = &paras[b];
        }
    }
    if (!critical || !cycles) {
        metrics->t_decomp_other += wall;
        return;
    }
    const double scale = wall / static_cast<double>(cycles);
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

class BamToSamOperator : public cpe::CpeBatchOperator {
public:
    BamToSamOperator(RankBodySink *sink, const sam_hdr_t *header,
                     BamToSamMetrics *metrics)
        : sink_(sink), header_(header), metrics_(metrics) {
        memset(paras_, 0, sizeof(paras_));
    }

    const char *name() const { return "bam-to-sam"; }
    size_t batch_capacity() const { return kBatchBlocks; }

    int Initialize() {
        const double t0 = GetTime();
        const int status = workspace_.Allocate() == 0 &&
                           writer_.Initialize(sink_, header_) == 0 ? 0 : -1;
        metrics_->t_alloc_init += GetTime() - t0;
        return status;
    }

    void Shutdown() {
        const double t0 = GetTime();
        writer_.Close();
        workspace_.Release();
        metrics_->t_free_workspace += GetTime() - t0;
    }

    int Prepare(const BgzfBlockBatch &compressed,
                BgzfBlockBatch *decoded, size_t count) {
        if (!decoded || count > kBatchBlocks) return -1;
        for (int b = 0; b < kBatchBlocks; ++b) {
            Bam2BamPara &para = paras_[b];
            para.block_id = b;
            para.output_records = workspace_.ptrs(b);
            para.record_base = workspace_.records(b);
            para.bam_lens = workspace_.lengths(b);
            para.record_capacity = kRecordsPerBlock;
            para.data_arena = workspace_.data(b);
            para.data_arena_capacity = MPI_BAM_BLOCK_ARENA_SIZE;
            para.data_arena_used = 0;
            para.n_total_records = 0;
            para.n_kept_records = 0;
            para.record_index = 0;
            para.actual_value = 0;
            para.limit_value = 0;
            para.limit_id = BOUNDS_LIMIT_NONE;
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
        return reinterpret_cast<void *>(slave_mpi_decompress_bam2bam_passthrough);
    }
    void *kernel_arguments() { return paras_; }

    int DuringKernel() { return writer_.FlushPrevious(); }

    void ObserveKernel(double wall, size_t count) {
        metrics_->t_decomp += wall;
        AccumulateDecodeDetail(paras_, count, wall, metrics_);
    }

    int Validate(size_t count) const {
        for (size_t b = 0; b < count; ++b) {
            if (paras_[b].status == 0) continue;
            if (paras_[b].status == -3) {
                fprintf(stderr,
                        "ERROR: MPI bam2sam capacity exceeded on input block %zu. "
                        "limit_id=%d limit=%lld actual=%lld record=%d.\n",
                        b, paras_[b].limit_id, paras_[b].limit_value,
                        paras_[b].actual_value, paras_[b].record_index);
            } else {
                fprintf(stderr,
                        "ERROR: MPI bam2sam decompress/parse failed on input block %zu "
                        "with status %d.\n", b, paras_[b].status);
            }
            return -1;
        }
        return 0;
    }

    int PostProcessBatch(size_t count, long long *processed) {
        return CollectAndSubmit(count, processed, nullptr, nullptr);
    }

    int PostProcessBatchWithPrefetch(
            size_t count, long long *processed,
            cpe::CpeBatchPrefetch prefetch, void *context) {
        return CollectAndSubmit(count, processed, prefetch, context);
    }

    int CollectAndSubmit(size_t count, long long *processed,
                         cpe::CpeBatchPrefetch prefetch, void *context) {
        const double collect_t0 = GetTime();
        size_t n = 0;
        bam1_t **slots = writer_.RecordSlots();
        for (size_t b = 0; b < count; ++b) {
            const int block_records = paras_[b].n_total_records;
            if (block_records < 0 || n + block_records > writer_.capacity()) {
                fprintf(stderr, "ERROR: MPI bam2sam format batch capacity exceeded.\n");
                return -1;
            }
            for (int r = 0; r < block_records; ++r) {
                slots[n++] = paras_[b].output_records[r];
            }
        }
        metrics_->t_collect += GetTime() - collect_t0;
        if (writer_.SubmitWithPrefetch(n, prefetch, context) != 0) return -1;
        *processed = static_cast<long long>(n);
        metrics_->group_count++;
        return 0;
    }

    int Finish() { return writer_.Finish(); }
    const cpe::CpeSamWriteTiming &write_timing() const { return writer_.timing(); }

private:
    RankBodySink *sink_;
    const sam_hdr_t *header_;
    BamToSamMetrics *metrics_;
    RecordWorkspace workspace_;
    cpe::CpeSamWriteSession writer_;
    Bam2BamPara paras_[kBatchBlocks];
};

int RunWithSource(cpe::CpeReadBatchSource *source, RankBodySink *sink,
                  const sam_hdr_t *header, BamToSamMetrics *metrics) {
    if (!source || !sink || !header) return -1;
    BamToSamMetrics local_metrics = {};
    if (!metrics) metrics = &local_metrics;
    const double t0 = GetTime();
    BamToSamOperator op(sink, header, metrics);
    cpe::CpeReadPipelineTiming read_timing;
    cpe::CpeReadPipelineOptions options;
    options.prefetch_during_post_process = true;
    const int status = cpe::RunCpeReadPipeline(source, &op, &read_timing,
                                                 options);
    metrics->input_blocks += read_timing.input_blocks;
    metrics->total_records += read_timing.total_records;
    metrics->t_read += read_timing.read;
    metrics->format_tiles += op.write_timing().batches;
    metrics->t_format += op.write_timing().format;
    metrics->t_write += op.write_timing().write;
    metrics->t_optimized_total += GetTime() - t0;
    return status;
}

class SpanBatchSource : public cpe::CpeReadBatchSource {
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

} // namespace

int RunBamToSamPipeline(MemReader *reader, RankBodySink *sink,
                        const sam_hdr_t *header, BamToSamMetrics *metrics) {
    if (!reader) return -1;
    MemoryBatchSource source(reader);
    return RunWithSource(&source, sink, header, metrics);
}

int RunBamToSamPipeline(const BamInputBackend &input,
                        const BgzfBlockSpan *spans, size_t span_count,
                        RankBodySink *sink, const sam_hdr_t *header,
                        BamToSamMetrics *metrics) {
    if (!spans && span_count) return -1;
    SpanBatchSource source(input, spans, span_count);
    return RunWithSource(&source, sink, header, metrics);
}

} // namespace operators
} // namespace swbam
