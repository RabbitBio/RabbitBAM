#include "swbam/raw_bam_writer.h"
#include "swbam/generic_compress.h"

#include <algorithm>
#include <cstring>
#include <vector>

namespace swbam {
namespace cpe {
namespace {

void AppendLe32(std::vector<unsigned char> *output, uint32_t value) {
    output->push_back((unsigned char)(value & 0xff));
    output->push_back((unsigned char)((value >> 8) & 0xff));
    output->push_back((unsigned char)((value >> 16) & 0xff));
    output->push_back((unsigned char)((value >> 24) & 0xff));
}

int SerializeBamHeader(const sam_hdr_t *header,
                       std::vector<unsigned char> *output) {
    if (!header || !output || header->n_targets < 0) return -1;
    sam_hdr_t *mutable_header = const_cast<sam_hdr_t *>(header);
    const char *text = sam_hdr_str(mutable_header);
    const size_t text_size = sam_hdr_length(mutable_header);
    if ((text_size > 0 && !text) || text_size > UINT32_MAX) return -1;

    size_t reserve_size = 12;
    if (text_size > output->max_size() - reserve_size) return -1;
    reserve_size += text_size;
    if ((size_t)header->n_targets >
            (output->max_size() - reserve_size) / 32) {
        return -1;
    }
    reserve_size += (size_t)header->n_targets * 32;

    try {
        output->clear();
        output->reserve(reserve_size);
        const unsigned char magic[4] = {'B', 'A', 'M', 1};
        output->insert(output->end(), magic, magic + 4);
        AppendLe32(output, (uint32_t)text_size);
        if (text_size > 0) output->insert(output->end(), text, text + text_size);
        AppendLe32(output, (uint32_t)header->n_targets);
        for (int32_t i = 0; i < header->n_targets; ++i) {
            if (!header->target_name || !header->target_name[i] ||
                !header->target_len) {
                return -1;
            }
            const size_t name_size = strlen(header->target_name[i]) + 1;
            if (name_size > UINT32_MAX ||
                output->size() > output->max_size() - name_size - 8) {
                return -1;
            }
            AppendLe32(output, (uint32_t)name_size);
            output->insert(output->end(), header->target_name[i],
                           header->target_name[i] + name_size);
            AppendLe32(output, header->target_len[i]);
        }
    } catch (...) {
        output->clear();
        return -1;
    }
    return 0;
}

class PackedChunkSource : public UncompressedBgzfSource {
public:
    PackedChunkSource(const BgzfBlockBatch *blocks, size_t count)
        : blocks_(blocks), count_(count), position_(0) {}

    int Fill(BgzfBlockBatch *batch, size_t *count) {
        if (!batch || !count || !blocks_) return -1;
        size_t n = count_ - position_;
        if (n > batch->capacity()) n = batch->capacity();
        for (size_t i = 0; i < n; ++i) {
            const bam_block &src = blocks_->blocks()[position_ + i];
            bam_block &dst = batch->blocks()[i];
            if (src.length == 0 || src.length > BGZF_BLOCK_SIZE) return -1;
            memcpy(dst.data, src.data, src.length);
            dst.pos = 0;
            dst.length = src.length;
            dst.errcode = 0;
            dst.block_id = src.block_id;
        }
        position_ += n;
        *count = n;
        return 0;
    }

private:
    const BgzfBlockBatch *blocks_;
    size_t count_;
    size_t position_;
};

class BackendOutputConsumer : public CompressedBgzfConsumer {
public:
    explicit BackendOutputConsumer(BamOutputBackend *output)
        : output_(output) {}

    int ConsumeCompressed(const bam_block *blocks, size_t count) {
        return WriteBgzfBlocks(output_, blocks, count);
    }

private:
    BamOutputBackend *output_;
};

} // namespace

RawBamWriterMetrics::RawBamWriterMetrics()
    : records(0), raw_bytes(0), header_raw_bytes(0), header_blocks(0),
      packed_blocks(0), compressed_bytes(0),
      pack(0.0), source(0.0), kernel(0.0), consume(0.0), total(0.0) {}

class RawBamWriter::Impl {
public:
    Impl()
        : output(nullptr), compression_level(1), chunk_blocks(0),
          finalized_blocks(0), current_size(0), next_block_id(0),
          initialized(false), finished(false) {}

    int Initialize(BamOutputBackend *new_output, int level,
                   size_t new_chunk_blocks) {
        if (!new_output || new_chunk_blocks == 0 ||
            (level != 0 && level != 1 && level != 6)) {
            return -1;
        }
        output = new_output;
        compression_level = level;
        chunk_blocks = new_chunk_blocks;
        finalized_blocks = 0;
        current_size = 0;
        next_block_id = 0;
        metrics = RawBamWriterMetrics();
        finished = false;
        if (packed.Allocate(chunk_blocks) != 0) return -1;
        initialized = true;
        return 0;
    }

    int InitializeBam(BamOutputBackend *new_output,
                      const sam_hdr_t *header, int level,
                      size_t new_chunk_blocks) {
        std::vector<unsigned char> raw_header;
        if (Initialize(new_output, level, new_chunk_blocks) != 0 ||
            SerializeBamHeader(header, &raw_header) != 0) {
            return -1;
        }
        if (AppendBytes(raw_header.data(), raw_header.size()) != 0 ||
            FinalizeCurrent() != 0) {
            return -1;
        }
        metrics.header_raw_bytes = (long long)raw_header.size();
        metrics.header_blocks = metrics.packed_blocks;
        return 0;
    }

    int AppendBytes(const unsigned char *data, size_t size) {
        if ((!data && size != 0) || !initialized || finished) return -1;
        double pack_t0 = GetTime();
        size_t offset = 0;
        while (offset < size) {
            if (finalized_blocks == chunk_blocks) {
                metrics.pack += GetTime() - pack_t0;
                if (Flush() != 0) return -1;
                pack_t0 = GetTime();
            }
            const size_t available = BGZF_BLOCK_SIZE - current_size;
            const size_t copy_size = std::min(available, size - offset);
            bam_block &block = packed.blocks()[finalized_blocks];
            memcpy(block.data + current_size, data + offset, copy_size);
            current_size += copy_size;
            offset += copy_size;
            if (current_size == BGZF_BLOCK_SIZE &&
                FinalizeCurrent() != 0) {
                return -1;
            }
        }
        metrics.pack += GetTime() - pack_t0;
        return 0;
    }

    int FinalizeCurrent() {
        if (current_size == 0) return 0;
        if (finalized_blocks >= chunk_blocks) return -1;
        bam_block &block = packed.blocks()[finalized_blocks];
        block.pos = 0;
        block.length = (unsigned int)current_size;
        block.errcode = 0;
        block.block_id = next_block_id++;
        finalized_blocks++;
        current_size = 0;
        metrics.packed_blocks++;
        return 0;
    }

    int Flush() {
        if (finalized_blocks == 0) return 0;
        PackedChunkSource source(&packed, finalized_blocks);
        BackendOutputConsumer consumer(output);
        CpeWritePipelineTiming timing;
        GenericCompressMetrics compress_metrics;
        const uint64_t output_before = output->bytes_written();
        const double flush_t0 = GetTime();
        if (RunGenericCompressPipeline(
                &source, &consumer, compression_level,
                &timing, &compress_metrics) != 0) {
            return -1;
        }
        metrics.source += timing.source;
        metrics.kernel += timing.kernel;
        metrics.consume += timing.consume;
        metrics.total += GetTime() - flush_t0;
        metrics.compressed_bytes +=
            (long long)(output->bytes_written() - output_before);
        finalized_blocks = 0;
        return 0;
    }

    int Consume(const RawBamRecordView *records, size_t count) {
        if (!initialized || finished || (!records && count != 0)) return -1;
        double pack_t0 = GetTime();
        for (size_t i = 0; i < count; ++i) {
            const RawBamRecordView &record = records[i];
            if (!record.encoded || record.encoded_size == 0 ||
                record.encoded_size > BGZF_BLOCK_SIZE) {
                return -1;
            }
            if (current_size + record.encoded_size > BGZF_BLOCK_SIZE &&
                FinalizeCurrent() != 0) {
                return -1;
            }
            if (finalized_blocks == chunk_blocks) {
                metrics.pack += GetTime() - pack_t0;
                if (Flush() != 0) return -1;
                pack_t0 = GetTime();
            }
            bam_block &block = packed.blocks()[finalized_blocks];
            memcpy(block.data + current_size,
                   record.encoded, record.encoded_size);
            current_size += record.encoded_size;
            metrics.records++;
            metrics.raw_bytes += record.encoded_size;
            if (current_size == BGZF_BLOCK_SIZE && FinalizeCurrent() != 0) {
                return -1;
            }
            if (finalized_blocks == chunk_blocks) {
                metrics.pack += GetTime() - pack_t0;
                if (Flush() != 0) return -1;
                pack_t0 = GetTime();
            }
        }
        metrics.pack += GetTime() - pack_t0;
        return 0;
    }

    int Finish(bool append_eof) {
        if (!initialized || finished) return -1;
        const double pack_t0 = GetTime();
        if (FinalizeCurrent() != 0) return -1;
        metrics.pack += GetTime() - pack_t0;
        if (Flush() != 0) return -1;
        if (append_eof && WriteBgzfEof(output) != 0) return -1;
        if (output->Flush() != 0) return -1;
        finished = true;
        return 0;
    }

    BamOutputBackend *output;
    int compression_level;
    size_t chunk_blocks;
    BgzfBlockBatch packed;
    size_t finalized_blocks;
    size_t current_size;
    int next_block_id;
    bool initialized;
    bool finished;
    RawBamWriterMetrics metrics;
};

RawBamWriter::RawBamWriter() : impl_(new Impl()) {}
RawBamWriter::~RawBamWriter() { delete impl_; }

int RawBamWriter::Initialize(BamOutputBackend *output,
                             int compression_level,
                             size_t chunk_blocks) {
    return impl_->Initialize(output, compression_level, chunk_blocks);
}

int RawBamWriter::InitializeBam(
        BamOutputBackend *output, const sam_hdr_t *header,
        int compression_level, size_t chunk_blocks) {
    return impl_->InitializeBam(
        output, header, compression_level, chunk_blocks);
}

int RawBamWriter::ConsumeRaw(
        const RawBamRecordView *records, size_t count) {
    return impl_->Consume(records, count);
}

int RawBamWriter::Finish(bool append_bgzf_eof) {
    return impl_->Finish(append_bgzf_eof);
}

const RawBamWriterMetrics &RawBamWriter::metrics() const {
    return impl_->metrics;
}

} // namespace cpe
} // namespace swbam
