#include "swbam_mpi.h"
#include "swbam/mpi_runtime.h"
#include "swbam/operators/record_count.h"
#include "swbam/operators/bgzf_compress.h"

#include <cstdio>
#include <cstring>
#include <vector>
#include "libdeflate.h"
#include <mpi.h>

namespace {

unsigned char SyntheticByte(int block_id, size_t offset) {
    return (unsigned char)(((size_t)block_id * 17u + offset) % 251u);
}

class SyntheticBgzfSource : public swbam::cpe::UncompressedBgzfSource {
public:
    explicit SyntheticBgzfSource(int total_blocks)
        : total_blocks_(total_blocks), next_block_(0) {}

    int Fill(swbam::BgzfBlockBatch *batch, size_t *count) {
        if (!batch || !count) return -1;
        size_t produced = 0;
        while (produced < batch->capacity() && next_block_ < total_blocks_) {
            bam_block &block = batch->blocks()[produced];
            const size_t length = 16384u +
                (size_t)(next_block_ % 8) * 1024u;
            for (size_t i = 0; i < length; ++i) {
                block.data[i] = SyntheticByte(next_block_, i);
            }
            block.pos = 0;
            block.length = (unsigned int)length;
            block.errcode = 0;
            block.block_id = next_block_;
            produced++;
            next_block_++;
        }
        *count = produced;
        return 0;
    }

private:
    int total_blocks_;
    int next_block_;
};

uint32_t ReadLe32(const unsigned char *data) {
    return (uint32_t)data[0] |
           ((uint32_t)data[1] << 8) |
           ((uint32_t)data[2] << 16) |
           ((uint32_t)data[3] << 24);
}

class VerifyingBgzfBatchPostProcessor : public swbam::cpe::CompressedBgzfBatchPostProcessor {
public:
    VerifyingBgzfBatchPostProcessor()
        : decompressor_(libdeflate_alloc_decompressor()), blocks_(0),
          bytes_(0), scratch_(BGZF_MAX_BLOCK_SIZE) {}
    ~VerifyingBgzfBatchPostProcessor() {
        if (decompressor_) libdeflate_free_decompressor(decompressor_);
    }

    int PostProcessCompressedBatch(const bam_block *blocks, size_t count) {
        if (!decompressor_ || (!blocks && count != 0)) return -1;
        for (size_t i = 0; i < count; ++i) {
            const bam_block &block = blocks[i];
            if (block.length < BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH ||
                block.data[0] != 0x1f || block.data[1] != 0x8b ||
                block.data[12] != 'B' || block.data[13] != 'C') {
                return -1;
            }
            const uint32_t bsize =
                (uint32_t)block.data[16] |
                ((uint32_t)block.data[17] << 8);
            if (bsize + 1 != block.length) return -1;

            const uint32_t expected_crc =
                ReadLe32(block.data + block.length - 8);
            const uint32_t expected_size =
                ReadLe32(block.data + block.length - 4);
            size_t decoded_size = scratch_.size();
            const int ret = libdeflate_deflate_decompress(
                decompressor_, block.data + BLOCK_HEADER_LENGTH,
                block.length - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH,
                scratch_.data(), scratch_.size(), &decoded_size);
            if (ret != 0 || decoded_size != expected_size ||
                libdeflate_crc32(0, scratch_.data(), decoded_size) !=
                    expected_crc) {
                return -1;
            }
            for (size_t j = 0; j < decoded_size; ++j) {
                if (scratch_[j] != SyntheticByte(block.block_id, j)) {
                    return -1;
                }
            }
            blocks_++;
            bytes_ += block.length;
        }
        return 0;
    }

    long long blocks() const { return blocks_; }
    long long bytes() const { return bytes_; }

private:
    struct libdeflate_decompressor *decompressor_;
    long long blocks_;
    long long bytes_;
    std::vector<unsigned char> scratch_;
};

int SameBamHeader(const sam_hdr_t *lhs, const sam_hdr_t *rhs) {
    if (!lhs || !rhs || lhs->n_targets != rhs->n_targets) return 0;
    sam_hdr_t *mutable_lhs = const_cast<sam_hdr_t *>(lhs);
    sam_hdr_t *mutable_rhs = const_cast<sam_hdr_t *>(rhs);
    const size_t lhs_text_size = sam_hdr_length(mutable_lhs);
    const size_t rhs_text_size = sam_hdr_length(mutable_rhs);
    const char *lhs_text = sam_hdr_str(mutable_lhs);
    const char *rhs_text = sam_hdr_str(mutable_rhs);
    if (lhs_text_size != rhs_text_size ||
        (lhs_text_size > 0 &&
         (!lhs_text || !rhs_text ||
          memcmp(lhs_text, rhs_text, lhs_text_size) != 0))) {
        return 0;
    }
    for (int32_t i = 0; i < lhs->n_targets; ++i) {
        if (!lhs->target_name || !rhs->target_name ||
            !lhs->target_name[i] || !rhs->target_name[i] ||
            strcmp(lhs->target_name[i], rhs->target_name[i]) != 0 ||
            !lhs->target_len || !rhs->target_len ||
            lhs->target_len[i] != rhs->target_len[i]) {
            return 0;
        }
    }
    return 1;
}

struct ReferenceCounts {
    long long records = 0;
    long long bytes = 0;
    long long kept = 0;
    long long kept_bytes = 0;
};

// Independent MPE validation, not part of the performance measurement. Check
// CRC, record boundaries and MAPQ against the new CPE count/filter pipelines.
int ScanReference(const swbam::MemoryBamInput &input,
                  const swbam::BgzfBlockSpan *spans, size_t count,
                  ReferenceCounts *counts) {
    struct libdeflate_decompressor *decoder = libdeflate_alloc_decompressor();
    if (!decoder) return -1;
    std::vector<unsigned char> raw(BGZF_MAX_BLOCK_SIZE);
    int status = 0;
    for (size_t i = 0; i < count && status == 0; ++i) {
        if (spans[i].offset > input.size() ||
            spans[i].compressed_size > input.size() - spans[i].offset ||
            spans[i].compressed_size < BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH) {
            status = -1;
            break;
        }
        const unsigned char *block = reinterpret_cast<const unsigned char *>(
            input.data()) + spans[i].offset;
        const size_t length = spans[i].compressed_size;
        size_t size = 0;
        if (libdeflate_deflate_decompress(decoder, block + BLOCK_HEADER_LENGTH,
                length - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH,
                raw.data(), raw.size(), &size) != 0 ||
            size != ReadLe32(block + length - 4) ||
            libdeflate_crc32(0, raw.data(), size) != ReadLe32(block + length - 8)) {
            status = -1;
            break;
        }
        for (size_t offset = 0; offset < size;) {
            if (size - offset < 4) { status = -1; break; }
            const uint64_t length = (uint64_t)ReadLe32(raw.data() + offset) + 4;
            if (length < 36 || length > size - offset) { status = -1; break; }
            ++counts->records;
            counts->bytes += length;
            if (raw[offset + 13] >= 30) {
                ++counts->kept;
                counts->kept_bytes += length;
            }
            offset += length;
        }
    }
    libdeflate_free_decompressor(decoder);
    return status;
}

int CheckRoundtrip(const swbam::MemoryBamInput &input,
                   const swbam::mpi::MpiBamInputPlan &plan,
                   const ReferenceCounts &expected, bool filter, int rank) {
    swbam::AdaptiveRankBodySink body;
    size_t reserve = 1024 * 1024;
    for (size_t i = plan.rank_begin; i < plan.rank_end; ++i)
        reserve += plan.blocks[i].compressed_size;
    int ok = body.Open("memory", 0, reserve, "io-check-body", "") == 0;
    swbam::operators::BamTransformMetrics metrics = {};
    swbam::BamFilterOptions options = swbam::DefaultBamFilterOptions();
    if (filter) options.min_mapq = 30;
    if (ok) ok = swbam::operators::RunBamTransformPipeline(
        input, plan.rank_spans(), plan.rank_block_count(), &body,
        options, 1, &metrics) == 0;
    const long long expected_records = filter ? expected.kept : expected.records;
    const long long expected_bytes = filter ? expected.kept_bytes : expected.bytes;
    ok = ok && metrics.total_records == expected.records &&
         metrics.kept_records == expected_records &&
         metrics.dropped_records == expected.records - expected_records;

    swbam::MemoryBamOutput output;
    char *header = nullptr;
    size_t header_size = 0;
    if (ok) ok = MpiCommonBuildBamHeaderMemory(
        const_cast<sam_hdr_t *>(input.header()), 1, &header, &header_size) == 0;
    if (ok) ok = output.Write(header, header_size) == 0;
    free(header);
    std::vector<unsigned char> buffer(8u * 1024u * 1024u);
    for (uint64_t offset = 0; ok && offset < body.size();) {
        const size_t n = (size_t)std::min<uint64_t>(buffer.size(), body.size() - offset);
        ok = body.ReadAt(offset, buffer.data(), n) == 0 &&
             output.Write(buffer.data(), n) == 0;
        offset += n;
    }
    if (ok) ok = swbam::WriteBgzfEof(&output) == 0;
    swbam::MemoryBamInput reread;
    std::vector<swbam::BgzfBlockSpan> spans;
    ReferenceCounts actual;
    if (ok) ok = reread.OpenMemoryCopy(output.data(), output.size()) == 0 &&
        SameBamHeader(reread.header(), input.header()) &&
        reread.ScanBlocks(&spans) == 0 &&
        ScanReference(reread, spans.data(), spans.size(), &actual) == 0;
    if (ok) {
        const uint64_t end = spans.empty() ? reread.body_offset()
            : spans.back().offset + spans.back().compressed_size;
        ok = end + 28 == output.size() && actual.records == expected_records &&
             actual.bytes == expected_bytes && (!filter || actual.kept == actual.records);
    }
    if (!ok) fprintf(stderr, "[rank %d] io-check %s roundtrip FAILED.\n",
                     rank, filter ? "filter" : "BAM");
    if (!swbam::mpi::AllRanksOk(ok)) return -1;
    long long local[3] = {metrics.total_records, metrics.kept_records,
                          (long long)output.size()};
    long long global[3] = {};
    MPI_Reduce(local, global, 3, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    const double cost = swbam::mpi::ReduceMaxCost(metrics.t_optimized_total);
    if (rank == 0) printf("SWBAM BAM read -> BAM write %s check PASS. "
        "total=%lld kept=%lld output_bytes=%lld pipeline_max=%.6f\n",
        filter ? "MAPQ>=30 filter" : "roundtrip", global[0], global[1], global[2], cost);
    return 0;
}

} // namespace

int ProcessIoCheckMPI(CmdInfo *cmd_info) {
    if (!cmd_info) return 1;
    int rank = 0, ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    swbam::MemoryBamInput input;
    swbam::mpi::MpiBamInputPlan plan;
    int ok = input.Open(cmd_info->in_file_name_) == 0 && input.format() == bam;
    if (!swbam::mpi::AllRanksOk(ok)) return 1;
    if (swbam::mpi::PrepareMpiBamInputPlan(input, &plan) != 0) return 1;
    ReferenceCounts expected;
    ok = ScanReference(input, plan.rank_spans(), plan.rank_block_count(), &expected) == 0;
    if (!swbam::mpi::AllRanksOk(ok)) return 1;

    swbam::PosixBamInput posix_input;
    swbam::mpi::MpiIoBamInput mpiio_input;
    swbam::BamInputBackend *backends[] = {&input, &posix_input, &mpiio_input};
    const char *names[] = {"memory", "posix", "mpiio"};
    for (size_t b = 0; b < 3; ++b) {
        swbam::BamInputBackend &backend = *backends[b];
        ok = b == 0 || backend.Open(cmd_info->in_file_name_) == 0;
        if (!swbam::mpi::AllRanksOk(ok)) return 1;
        ok = backend.format() == bam && backend.size() == input.size() &&
             backend.body_offset() == input.body_offset() &&
             SameBamHeader(backend.header(), input.header());
        if (!swbam::mpi::AllRanksOk(ok)) return 1;
        swbam::mpi::MpiBamInputPlan current;
        if (swbam::mpi::PrepareMpiBamInputPlan(backend, &current) != 0) return 1;
        ok = current.blocks.size() == plan.blocks.size() &&
             current.rank_begin == plan.rank_begin && current.rank_end == plan.rank_end;
        for (size_t i = 0; ok && i < plan.blocks.size(); ++i) {
            ok = current.blocks[i].offset == plan.blocks[i].offset &&
                 current.blocks[i].compressed_size == plan.blocks[i].compressed_size;
        }
        if (!swbam::mpi::AllRanksOk(ok)) return 1;
        swbam::operators::RecordCountMetrics metrics = {};
        ok = swbam::operators::RunRecordCountPipeline(
            backend, current.rank_spans(), current.rank_block_count(), &metrics) == 0 &&
            metrics.records == expected.records &&
            metrics.blocks == (long long)current.rank_block_count();
        if (!ok) fprintf(stderr, "[rank %d] io-check %s count FAILED.\n", rank, names[b]);
        if (!swbam::mpi::AllRanksOk(ok)) return 1;
        long long local[3] = {metrics.blocks, metrics.records, expected.bytes};
        long long global[3] = {};
        MPI_Reduce(local, global, 3, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        const double total = swbam::mpi::ReduceMaxCost(metrics.total);
        const double kernel = swbam::mpi::ReduceMaxCost(metrics.kernel);
        if (rank == 0) printf("SWBAM BAM read/count %s check PASS. ranks=%d "
            "blocks=%lld records=%lld decoded_bytes=%lld pipeline_max=%.6f kernel_max=%.6f\n",
            names[b], ranks, global[0], global[1], global[2], total, kernel);
    }

    // Two full batches exercise output flushing while CPE compression runs.
    SyntheticBgzfSource source(128);
    VerifyingBgzfBatchPostProcessor verifier;
    swbam::cpe::CpeWritePipelineTiming timing;
    swbam::cpe::BgzfCompressMetrics compress;
    ok = swbam::cpe::RunBgzfCompressPipeline(&source, &verifier, 1,
        &timing, &compress) == 0 && verifier.blocks() == 128 &&
        timing.input_blocks == 128 && timing.output_blocks == 128 &&
        timing.output_bytes == verifier.bytes();
    if (!swbam::mpi::AllRanksOk(ok)) return 1;
    if (rank == 0) printf("SWBAM BAM write/compress check PASS. blocks_per_rank=128\n");
    if (CheckRoundtrip(input, plan, expected, false, rank) != 0 ||
        CheckRoundtrip(input, plan, expected, true, rank) != 0) return 1;
    if (rank == 0) printf("io-check: PASS (BAM read/write, count/filter, memory/POSIX/MPI-IO)\n");
    return 0;
}
