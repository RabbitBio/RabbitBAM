#include "swbam_mpi.h"
#include "swbam/generic_compress.h"
#include "swbam/io.h"
#include "swbam/mpi_runtime.h"
#include "swbam/raw_bam.h"
#include "swbam/raw_bam_filter.h"
#include "swbam/raw_bam_writer.h"

#include <cstdio>
#include <cstring>
#include <vector>

#include "libdeflate.h"

#include <mpi.h>

namespace {

class DecodeCheckConsumer : public swbam::cpe::RawBamRecordConsumer {
public:
    DecodeCheckConsumer() : records_(0), encoded_bytes_(0) {}

    int ConsumeRaw(const swbam::cpe::RawBamRecordView *records,
                   size_t count) {
        if (!records && count != 0) return -1;
        for (size_t i = 0; i < count; ++i) {
            if (!records[i].encoded || records[i].block_size < 32 ||
                records[i].encoded_size != records[i].block_size + 4) {
                return -1;
            }
            records_++;
            encoded_bytes_ += records[i].encoded_size;
        }
        return 0;
    }

    long long records() const { return records_; }
    long long encoded_bytes() const { return encoded_bytes_; }

private:
    long long records_;
    long long encoded_bytes_;
};

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

class VerifyingBgzfConsumer : public swbam::cpe::CompressedBgzfConsumer {
public:
    VerifyingBgzfConsumer()
        : decompressor_(libdeflate_alloc_decompressor()), blocks_(0),
          bytes_(0), scratch_(BGZF_MAX_BLOCK_SIZE) {}
    ~VerifyingBgzfConsumer() {
        if (decompressor_) libdeflate_free_decompressor(decompressor_);
    }

    int ConsumeCompressed(const bam_block *blocks, size_t count) {
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

int VerifyMemoryBam(const swbam::MemoryBamOutput &output,
                    const sam_hdr_t *expected_header,
                    long long expected_records,
                    long long expected_raw_bytes) {
    if (!output.data() || output.size() < 28) return -1;
    swbam::MemoryBamInput verified_input;
    std::vector<swbam::BgzfBlockSpan> spans;
    if (verified_input.OpenMemoryCopy(output.data(), output.size()) != 0 ||
        verified_input.format() != bam ||
        !SameBamHeader(verified_input.header(), expected_header) ||
        verified_input.ScanBlocks(&spans) != 0) {
        return -1;
    }
    const size_t body_end = spans.empty()
        ? (size_t)verified_input.body_offset()
        : (size_t)(spans.back().offset + spans.back().compressed_size);
    if (body_end + 28 != output.size()) return -1;

    struct libdeflate_decompressor *decompressor =
        libdeflate_alloc_decompressor();
    if (!decompressor) return -1;
    std::vector<unsigned char> decoded(BGZF_MAX_BLOCK_SIZE);
    long long records = 0;
    long long raw_bytes = 0;
    int result = 0;
    for (size_t i = 0; i < spans.size(); ++i) {
        const unsigned char *compressed =
            output.data() + (size_t)spans[i].offset;
        const size_t compressed_size = spans[i].compressed_size;
        if (compressed_size < BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH) {
            result = -1;
            break;
        }
        size_t decoded_size = decoded.size();
        const int ret = libdeflate_deflate_decompress(
            decompressor, compressed + BLOCK_HEADER_LENGTH,
            compressed_size - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH,
            decoded.data(), decoded.size(), &decoded_size);
        if (ret != 0 ||
            decoded_size != ReadLe32(compressed + compressed_size - 4) ||
            libdeflate_crc32(0, decoded.data(), decoded_size) !=
                ReadLe32(compressed + compressed_size - 8)) {
            result = -1;
            break;
        }
        size_t offset = 0;
        while (offset < decoded_size) {
            if (decoded_size - offset < 4) {
                result = -1;
                break;
            }
            const uint32_t block_size = ReadLe32(decoded.data() + offset);
            const uint64_t encoded_size = (uint64_t)block_size + 4;
            if (block_size < 32 || encoded_size > decoded_size - offset) {
                result = -1;
                break;
            }
            records++;
            raw_bytes += (long long)encoded_size;
            offset += (size_t)encoded_size;
        }
        if (result != 0) break;
    }
    libdeflate_free_decompressor(decompressor);
    if (result != 0 || records != expected_records ||
        raw_bytes != expected_raw_bytes) {
        return -1;
    }
    return 0;
}

} // namespace

int ProcessIoCheckMPI(CmdInfo *cmd_info) {
    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int exit_code = 1;
    int local_ok = 1;
    swbam::MemoryBamInput input;
    swbam::mpi::MpiBamInputPlan plan;
    DecodeCheckConsumer consumer;
    swbam::cpe::CpeReadPipelineTiming timing;
    swbam::cpe::GenericDecodeMetrics metrics;
    swbam::cpe::GenericRawBamMetrics raw_metrics;
    swbam::cpe::CpeWritePipelineTiming write_timing;
    swbam::cpe::GenericCompressMetrics compress_metrics;

    if (!cmd_info || input.Load(cmd_info->in_file_name_) != 0 ||
        input.ParseHeader() != 0 || input.format() != bam) {
        local_ok = 0;
    }
    if (!swbam::mpi::AllRanksOk(local_ok)) goto cleanup;

    if (swbam::mpi::PrepareMpiBamInputPlan(input, &plan) != 0) {
        local_ok = 0;
    }
    if (!swbam::mpi::AllRanksOk(local_ok)) goto cleanup;

    if (swbam::cpe::RunGenericRawBamPipeline(
            input, plan.rank_spans(), plan.rank_block_count(),
            &consumer, &timing, &metrics, &raw_metrics) != 0 ||
        metrics.decoded_blocks != (long long)plan.rank_block_count() ||
        raw_metrics.records != consumer.records() ||
        raw_metrics.encoded_bytes != consumer.encoded_bytes() ||
        raw_metrics.encoded_bytes != metrics.decoded_bytes) {
        local_ok = 0;
    }
    if (!swbam::mpi::AllRanksOk(local_ok)) goto cleanup;

    {
        long long local_values[3] = {
            metrics.decoded_blocks, metrics.decoded_bytes,
            consumer.records()
        };
        long long global_values[3] = {};
        MPI_Reduce(local_values, global_values, 3, MPI_LONG_LONG,
                   MPI_SUM, 0, MPI_COMM_WORLD);
        const double max_total = swbam::mpi::ReduceMaxCost(timing.total);
        const double max_kernel = swbam::mpi::ReduceMaxCost(timing.kernel);
        const double max_consume = swbam::mpi::ReduceMaxCost(timing.consume);
        if (rank == 0) {
            const long long expected_blocks =
                (long long)plan.blocks.size();
            printf("SWBAM generic raw BAM check finished. ranks=%d blocks=%lld expected=%lld decoded_bytes=%lld records=%lld\n",
                   comm_size, global_values[0], expected_blocks,
                   global_values[1], global_values[2]);
            printf("  pipeline_max=%.6f kernel_max=%.6f consume_max=%.6f\n",
                   max_total, max_kernel, max_consume);
            if (global_values[0] != expected_blocks) local_ok = 0;
        }
        MPI_Bcast(&local_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (!local_ok) goto cleanup;
    }

    {
        swbam::PosixBamInput posix_input;
        swbam::mpi::MpiBamInputPlan posix_plan;
        DecodeCheckConsumer posix_consumer;
        swbam::cpe::CpeReadPipelineTiming posix_timing;
        swbam::cpe::GenericDecodeMetrics posix_decode_metrics;
        swbam::cpe::GenericRawBamMetrics posix_raw_metrics;

        if (posix_input.Open(cmd_info->in_file_name_) != 0 ||
            posix_input.format() != bam ||
            posix_input.size() != input.size() ||
            posix_input.body_offset() != input.body_offset() ||
            !SameBamHeader(posix_input.header(), input.header()) ||
            swbam::mpi::PrepareMpiBamInputPlan(
                posix_input, &posix_plan) != 0 ||
            posix_plan.blocks.size() != plan.blocks.size()) {
            local_ok = 0;
        }
        if (local_ok) {
            for (size_t i = 0; i < plan.blocks.size(); ++i) {
                if (posix_plan.blocks[i].offset != plan.blocks[i].offset ||
                    posix_plan.blocks[i].compressed_size !=
                        plan.blocks[i].compressed_size) {
                    local_ok = 0;
                    break;
                }
            }
        }
        if (!swbam::mpi::AllRanksOk(local_ok)) goto cleanup;

        if (swbam::cpe::RunGenericRawBamPipeline(
                posix_input, posix_plan.rank_spans(),
                posix_plan.rank_block_count(), &posix_consumer,
                &posix_timing, &posix_decode_metrics,
                &posix_raw_metrics) != 0 ||
            posix_decode_metrics.decoded_blocks !=
                metrics.decoded_blocks ||
            posix_decode_metrics.decoded_bytes != metrics.decoded_bytes ||
            posix_raw_metrics.records != raw_metrics.records ||
            posix_raw_metrics.encoded_bytes != raw_metrics.encoded_bytes ||
            posix_consumer.records() != consumer.records() ||
            posix_consumer.encoded_bytes() != consumer.encoded_bytes()) {
            local_ok = 0;
        }
        if (!swbam::mpi::AllRanksOk(local_ok)) goto cleanup;

        const double posix_total_max =
            swbam::mpi::ReduceMaxCost(posix_timing.total);
        const double posix_read_max =
            swbam::mpi::ReduceMaxCost(posix_timing.read);
        long long local_posix[3] = {
            posix_decode_metrics.decoded_blocks,
            posix_decode_metrics.decoded_bytes,
            posix_raw_metrics.records
        };
        long long global_posix[3] = {};
        MPI_Reduce(local_posix, global_posix, 3, MPI_LONG_LONG,
                   MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            printf("SWBAM POSIX streaming input check finished. blocks=%lld decoded_bytes=%lld records=%lld\n",
                   global_posix[0], global_posix[1], global_posix[2]);
            printf("  pipeline_max=%.6f pread_max=%.6f\n",
                   posix_total_max, posix_read_max);
        }
    }

    {
        swbam::mpi::MpiIoBamInput mpiio_input;
        swbam::mpi::MpiBamInputPlan mpiio_plan;
        DecodeCheckConsumer mpiio_consumer;
        swbam::cpe::CpeReadPipelineTiming mpiio_timing;
        swbam::cpe::GenericDecodeMetrics mpiio_decode_metrics;
        swbam::cpe::GenericRawBamMetrics mpiio_raw_metrics;

        if (mpiio_input.Open(cmd_info->in_file_name_) != 0 ||
            mpiio_input.format() != bam ||
            mpiio_input.size() != input.size() ||
            mpiio_input.body_offset() != input.body_offset() ||
            !SameBamHeader(mpiio_input.header(), input.header()) ||
            swbam::mpi::PrepareMpiBamInputPlan(
                mpiio_input, &mpiio_plan) != 0 ||
            mpiio_plan.blocks.size() != plan.blocks.size()) {
            local_ok = 0;
        }
        if (local_ok) {
            for (size_t i = 0; i < plan.blocks.size(); ++i) {
                if (mpiio_plan.blocks[i].offset != plan.blocks[i].offset ||
                    mpiio_plan.blocks[i].compressed_size !=
                        plan.blocks[i].compressed_size) {
                    local_ok = 0;
                    break;
                }
            }
        }
        if (!swbam::mpi::AllRanksOk(local_ok)) goto cleanup;

        if (swbam::cpe::RunGenericRawBamPipeline(
                mpiio_input, mpiio_plan.rank_spans(),
                mpiio_plan.rank_block_count(), &mpiio_consumer,
                &mpiio_timing, &mpiio_decode_metrics,
                &mpiio_raw_metrics) != 0 ||
            mpiio_decode_metrics.decoded_blocks != metrics.decoded_blocks ||
            mpiio_decode_metrics.decoded_bytes != metrics.decoded_bytes ||
            mpiio_raw_metrics.records != raw_metrics.records ||
            mpiio_raw_metrics.encoded_bytes != raw_metrics.encoded_bytes ||
            mpiio_consumer.records() != consumer.records() ||
            mpiio_consumer.encoded_bytes() != consumer.encoded_bytes()) {
            local_ok = 0;
        }
        if (!swbam::mpi::AllRanksOk(local_ok)) goto cleanup;

        const double mpiio_total_max =
            swbam::mpi::ReduceMaxCost(mpiio_timing.total);
        const double mpiio_read_max =
            swbam::mpi::ReduceMaxCost(mpiio_timing.read);
        long long local_mpiio[3] = {
            mpiio_decode_metrics.decoded_blocks,
            mpiio_decode_metrics.decoded_bytes,
            mpiio_raw_metrics.records
        };
        long long global_mpiio[3] = {};
        MPI_Reduce(local_mpiio, global_mpiio, 3, MPI_LONG_LONG,
                   MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            printf("SWBAM MPI-IO streaming input check finished. blocks=%lld decoded_bytes=%lld records=%lld\n",
                   global_mpiio[0], global_mpiio[1], global_mpiio[2]);
            printf("  pipeline_max=%.6f mpi_read_max=%.6f\n",
                   mpiio_total_max, mpiio_read_max);
        }
    }

    {
        const int synthetic_blocks = 128;
        SyntheticBgzfSource source(synthetic_blocks);
        VerifyingBgzfConsumer write_consumer;
        if (swbam::cpe::RunGenericCompressPipeline(
                &source, &write_consumer, 1,
                &write_timing, &compress_metrics) != 0 ||
            write_consumer.blocks() != synthetic_blocks ||
            write_timing.input_blocks != synthetic_blocks ||
            write_timing.output_blocks != synthetic_blocks ||
            write_timing.output_bytes != write_consumer.bytes()) {
            local_ok = 0;
        }
        if (!swbam::mpi::AllRanksOk(local_ok)) goto cleanup;

        const double write_max =
            swbam::mpi::ReduceMaxCost(write_timing.total);
        const double compress_max =
            swbam::mpi::ReduceMaxCost(write_timing.kernel);
        const double output_max =
            swbam::mpi::ReduceMaxCost(write_timing.consume);
        if (rank == 0) {
            printf("SWBAM generic write check finished. blocks_per_rank=%d compression=1\n",
                   synthetic_blocks);
            printf("  pipeline_max=%.6f compress_max=%.6f output_verify_max=%.6f\n",
                   write_max, compress_max, output_max);
        }
    }

    {
        swbam::MemoryBamOutput roundtrip_output;
        swbam::cpe::RawBamMemoryWriter writer;
        swbam::cpe::CpeReadPipelineTiming roundtrip_read_timing;
        swbam::cpe::GenericDecodeMetrics roundtrip_decode_metrics;
        swbam::cpe::GenericRawBamMetrics roundtrip_raw_metrics;
        const size_t reserve_size = input.size() /
            (size_t)(comm_size > 0 ? comm_size : 1) + 1024 * 1024;
        if (roundtrip_output.Reserve(reserve_size) != 0 ||
            writer.InitializeBam(
                &roundtrip_output, input.header(), 1, 256) != 0 ||
            swbam::cpe::RunGenericRawBamPipeline(
                input, plan.rank_spans(), plan.rank_block_count(),
                &writer, &roundtrip_read_timing,
                &roundtrip_decode_metrics, &roundtrip_raw_metrics) != 0 ||
            writer.Finish(true) != 0 ||
            writer.metrics().records != roundtrip_raw_metrics.records ||
            writer.metrics().raw_bytes != roundtrip_raw_metrics.encoded_bytes ||
            VerifyMemoryBam(
                roundtrip_output, input.header(),
                roundtrip_raw_metrics.records,
                roundtrip_raw_metrics.encoded_bytes) != 0) {
            local_ok = 0;
        }
        if (!swbam::mpi::AllRanksOk(local_ok)) goto cleanup;

        const double roundtrip_max =
            swbam::mpi::ReduceMaxCost(roundtrip_read_timing.total);
        const double writer_max =
            swbam::mpi::ReduceMaxCost(writer.metrics().total);
        long long local_output[3] = {
            writer.metrics().records,
            writer.metrics().packed_blocks,
            (long long)roundtrip_output.size()
        };
        long long global_output[3] = {};
        MPI_Reduce(local_output, global_output, 3, MPI_LONG_LONG,
                   MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            printf("SWBAM complete BAM memory roundtrip finished. records=%lld bgzf_blocks=%lld output_bytes=%lld\n",
                   global_output[0], global_output[1], global_output[2]);
            printf("  read_pack_write_max=%.6f writer_flush_max=%.6f\n",
                   roundtrip_max, writer_max);
        }
    }

    {
        BamFilterOptions filter;
        filter.min_mapq = 30;
        filter.max_mapq = -1;
        filter.require_flag = 0;
        filter.exclude_flag = 0;
        filter.ref_tid = -2;
        filter.min_read_len = -1;
        filter.max_read_len = -1;

        swbam::MemoryBamOutput filtered_output;
        swbam::cpe::RawBamMemoryWriter writer;
        swbam::cpe::RawBamFilterConsumer filter_consumer(filter, &writer);
        swbam::cpe::CpeReadPipelineTiming filter_timing;
        swbam::cpe::GenericDecodeMetrics filter_decode_metrics;
        swbam::cpe::GenericRawBamMetrics filter_raw_metrics;
        const size_t reserve_size = input.size() /
            (size_t)(comm_size > 0 ? comm_size : 1) + 1024 * 1024;
        if (filtered_output.Reserve(reserve_size) != 0 ||
            writer.InitializeBam(
                &filtered_output, input.header(), 1, 256) != 0 ||
            swbam::cpe::RunGenericRawBamPipeline(
                input, plan.rank_spans(), plan.rank_block_count(),
                &filter_consumer, &filter_timing,
                &filter_decode_metrics, &filter_raw_metrics) != 0 ||
            writer.Finish(true) != 0 ||
            filter_consumer.metrics().total_records !=
                filter_raw_metrics.records ||
            filter_consumer.metrics().kept_records !=
                writer.metrics().records ||
            filter_consumer.metrics().kept_bytes !=
                writer.metrics().raw_bytes ||
            VerifyMemoryBam(
                filtered_output, input.header(),
                writer.metrics().records,
                writer.metrics().raw_bytes) != 0) {
            local_ok = 0;
        }
        if (!swbam::mpi::AllRanksOk(local_ok)) goto cleanup;

        long long local_filter[4] = {
            filter_consumer.metrics().total_records,
            filter_consumer.metrics().kept_records,
            filter_consumer.metrics().dropped_records,
            (long long)filtered_output.size()
        };
        long long global_filter[4] = {};
        MPI_Reduce(local_filter, global_filter, 4, MPI_LONG_LONG,
                   MPI_SUM, 0, MPI_COMM_WORLD);
        const double filter_max = swbam::mpi::ReduceMaxCost(
            filter_consumer.metrics().filter);
        const double pipeline_max =
            swbam::mpi::ReduceMaxCost(filter_timing.total);
        if (rank == 0) {
            printf("SWBAM generic BAM filter check finished. min_mapq=30 total=%lld kept=%lld dropped=%lld output_bytes=%lld\n",
                   global_filter[0], global_filter[1],
                   global_filter[2], global_filter[3]);
            printf("  pipeline_max=%.6f filter_max=%.6f\n",
                   pipeline_max, filter_max);
        }
    }

    exit_code = 0;

cleanup:
    input.Close();
    return exit_code;
}
