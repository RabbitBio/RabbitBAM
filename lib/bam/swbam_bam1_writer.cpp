#include "swbam/bam1_writer.h"

#include <climits>
#include <cstring>
#include <vector>

namespace swbam {
namespace cpe {
namespace {

void WriteLe32(unsigned char *output, uint32_t value) {
    output[0] = (unsigned char)(value & 0xff);
    output[1] = (unsigned char)((value >> 8) & 0xff);
    output[2] = (unsigned char)((value >> 16) & 0xff);
    output[3] = (unsigned char)((value >> 24) & 0xff);
}

uint32_t ReadHost32(const unsigned char *data) {
    uint32_t value = 0;
    memcpy(&value, data, sizeof(value));
    return value;
}

int ValidateCigarAndBin(const bam1_t *record, uint16_t *encoded_bin) {
    if (!record || !encoded_bin) return -1;
    const bam1_core_t &core = record->core;
    const uint32_t *cigar = bam_get_cigar(record);
    uint64_t query_length = 0;
    uint64_t reference_length = 0;
    for (uint32_t i = 0; i < core.n_cigar; ++i) {
        const uint32_t operation = bam_cigar_op(cigar[i]);
        if (operation > BAM_CBACK) return -1;
        const uint64_t length = bam_cigar_oplen(cigar[i]);
        const int type = bam_cigar_type(operation);
        if ((type & 1) != 0) query_length += length;
        if ((type & 2) != 0) reference_length += length;
    }
    if (!(core.flag & BAM_FUNMAP) && core.l_qseq > 0 &&
        query_length != (uint64_t)core.l_qseq) {
        return -1;
    }
    if ((core.flag & BAM_FUNMAP) || reference_length == 0) {
        reference_length = 1;
    }
    if (core.tid >= 0 && core.pos >= 0) {
        if (reference_length > (uint64_t)INT64_MAX - (uint64_t)core.pos) {
            return -1;
        }
        const hts_pos_t end = core.pos + (hts_pos_t)reference_length;
        *encoded_bin = (uint16_t)hts_reg2bin(core.pos, end, 14, 5);
    } else {
        *encoded_bin = core.bin;
    }
    return 0;
}

int EncodedRecordSize(const bam1_t *record, uint32_t *encoded_size,
                      uint16_t *encoded_bin) {
    if (!record || !encoded_size || record->l_data < 0 ||
        (record->l_data > 0 && !record->data)) {
        return -1;
    }
    const bam1_core_t &core = record->core;
    if (core.l_extranul > 3 || core.l_qname <= core.l_extranul ||
        core.l_qname > record->l_data || (core.l_qname & 3) != 0 ||
        core.l_qseq < 0 || core.n_cigar > 0xffff ||
        core.pos < INT_MIN || core.pos > INT_MAX ||
        core.mpos < INT_MIN || core.mpos > INT_MAX ||
        core.isize < INT_MIN || core.isize > INT_MAX) {
        return -1;
    }

    const uint32_t raw_qname_length =
        (uint32_t)core.l_qname - core.l_extranul;
    if (raw_qname_length < 1 || raw_qname_length > 255 ||
        record->data[raw_qname_length - 1] != '\0') {
        return -1;
    }
    const uint64_t minimum_data =
        (uint64_t)core.l_qname + (uint64_t)core.n_cigar * 4 +
        (((uint64_t)(uint32_t)core.l_qseq + 1) >> 1) +
        (uint64_t)(uint32_t)core.l_qseq;
    if (minimum_data > (uint64_t)record->l_data) return -1;

    const uint64_t block_size = 32 + (uint64_t)raw_qname_length +
        (uint64_t)((size_t)record->l_data - core.l_qname);
    const uint64_t total_size = block_size + 4;
    if (block_size > INT_MAX || total_size > BGZF_BLOCK_SIZE ||
        total_size > UINT32_MAX) {
        return -1;
    }
    if (ValidateCigarAndBin(record, encoded_bin) != 0) return -1;
    *encoded_size = (uint32_t)total_size;
    return 0;
}

int EncodeRecord(const bam1_t *record, uint32_t encoded_size,
                 uint16_t encoded_bin,
                 unsigned char *output) {
    if (!record || !output || encoded_size < 36) return -1;
    const bam1_core_t &core = record->core;
    const uint32_t raw_qname_length =
        (uint32_t)core.l_qname - core.l_extranul;
    const uint32_t block_size = encoded_size - 4;

    WriteLe32(output, block_size);
    WriteLe32(output + 4, (uint32_t)core.tid);
    WriteLe32(output + 8, (uint32_t)core.pos);
    WriteLe32(output + 12,
              (uint32_t)encoded_bin << 16 |
              (uint32_t)core.qual << 8 |
              raw_qname_length);
    WriteLe32(output + 16,
              (uint32_t)core.flag << 16 | core.n_cigar);
    WriteLe32(output + 20, (uint32_t)core.l_qseq);
    WriteLe32(output + 24, (uint32_t)core.mtid);
    WriteLe32(output + 28, (uint32_t)core.mpos);
    WriteLe32(output + 32, (uint32_t)core.isize);

    unsigned char *payload = output + 36;
    memcpy(payload, record->data, raw_qname_length);
    const unsigned char *cigar = record->data + core.l_qname;
    unsigned char *encoded_cigar = payload + raw_qname_length;
    for (uint32_t i = 0; i < core.n_cigar; ++i) {
        WriteLe32(encoded_cigar + (size_t)i * 4,
                  ReadHost32(cigar + (size_t)i * 4));
    }
    const size_t cigar_bytes = (size_t)core.n_cigar * 4;
    const size_t tail_offset = (size_t)core.l_qname + cigar_bytes;
    const size_t tail_size = (size_t)record->l_data - tail_offset;
    memcpy(encoded_cigar + cigar_bytes,
           record->data + tail_offset, tail_size);
    return 0;
}

} // namespace

Bam1WriterMetrics::Bam1WriterMetrics()
    : records(0), bam1_data_bytes(0), encoded_bytes(0), encode(0.0) {}

class Bam1Writer::Impl {
public:
    Impl() : initialized(false), finished(false) {}

    int Initialize(BamOutputBackend *output, int compression_level,
                   size_t chunk_blocks, const sam_hdr_t *header) {
        metrics = Bam1WriterMetrics();
        arena.clear();
        encoded_sizes.clear();
        encoded_bins.clear();
        views.clear();
        finished = false;
        const int status = header
            ? raw_writer.InitializeBam(
                  output, header, compression_level, chunk_blocks)
            : raw_writer.Initialize(
                  output, compression_level, chunk_blocks);
        initialized = status == 0;
        return status;
    }

    int Consume(const bam1_t *const *records, size_t count) {
        if (!initialized || finished || (!records && count != 0) ||
            count > UINT32_MAX) {
            return -1;
        }
        const double encode_t0 = GetTime();
        size_t total_size = 0;
        long long data_bytes = 0;
        try {
            encoded_sizes.resize(count);
            encoded_bins.resize(count);
            views.resize(count);
        } catch (...) {
            return -1;
        }
        for (size_t i = 0; i < count; ++i) {
            uint32_t encoded_size = 0;
            uint16_t encoded_bin = 0;
            if (EncodedRecordSize(records[i], &encoded_size,
                                  &encoded_bin) != 0 ||
                encoded_size > arena.max_size() - total_size) {
                return -1;
            }
            encoded_sizes[i] = encoded_size;
            encoded_bins[i] = encoded_bin;
            total_size += encoded_size;
            data_bytes += records[i]->l_data;
        }
        try {
            arena.resize(total_size);
        } catch (...) {
            return -1;
        }

        size_t offset = 0;
        for (size_t i = 0; i < count; ++i) {
            if (EncodeRecord(records[i], encoded_sizes[i], encoded_bins[i],
                             arena.data() + offset) != 0) {
                return -1;
            }
            RawBamRecordView &view = views[i];
            view.encoded = arena.data() + offset;
            view.encoded_size = encoded_sizes[i];
            view.block_size = encoded_sizes[i] - 4;
            view.block_offset = 0;
            view.block_index = (uint32_t)i;
            offset += encoded_sizes[i];
        }
        const double encode = GetTime() - encode_t0;
        if (raw_writer.ConsumeRaw(
                count == 0 ? nullptr : views.data(), count) != 0) {
            return -1;
        }
        metrics.records += (long long)count;
        metrics.bam1_data_bytes += data_bytes;
        metrics.encoded_bytes += (long long)total_size;
        metrics.encode += encode;
        return 0;
    }

    RawBamWriter raw_writer;
    Bam1WriterMetrics metrics;
    std::vector<unsigned char> arena;
    std::vector<uint32_t> encoded_sizes;
    std::vector<uint16_t> encoded_bins;
    std::vector<RawBamRecordView> views;
    bool initialized;
    bool finished;
};

Bam1Writer::Bam1Writer() : impl_(new Impl()) {}
Bam1Writer::~Bam1Writer() { delete impl_; }

int Bam1Writer::Initialize(BamOutputBackend *output,
                           int compression_level,
                           size_t chunk_blocks) {
    return impl_->Initialize(
        output, compression_level, chunk_blocks, nullptr);
}

int Bam1Writer::InitializeBam(BamOutputBackend *output,
                              const sam_hdr_t *header,
                              int compression_level,
                              size_t chunk_blocks) {
    if (!header) return -1;
    return impl_->Initialize(
        output, compression_level, chunk_blocks, header);
}

int Bam1Writer::ConsumeBam1(
        const bam1_t *const *records, size_t count) {
    return impl_->Consume(records, count);
}

int Bam1Writer::Finish(bool append_bgzf_eof) {
    if (!impl_->initialized || impl_->finished) return -1;
    if (impl_->raw_writer.Finish(append_bgzf_eof) != 0) return -1;
    impl_->finished = true;
    return 0;
}

const Bam1WriterMetrics &Bam1Writer::metrics() const {
    return impl_->metrics;
}

const RawBamWriterMetrics &Bam1Writer::raw_metrics() const {
    return impl_->raw_writer.metrics();
}

} // namespace cpe
} // namespace swbam
