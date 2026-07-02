#include "swbam/io.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

namespace swbam {
namespace {

const unsigned char kBgzfEofBlock[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
};

int NormalizeFormat(int format) {
    if (format == bam || format == binary_format) return bam;
    if (format == sam || format == text_format) return sam;
    return format;
}

int BgzfBlockLengthAt(const char *base, size_t size,
                      size_t pos, size_t *block_len) {
    if (!base || !block_len || pos > size ||
        size - pos < BLOCK_HEADER_LENGTH) {
        return -1;
    }
    const unsigned char *block =
        reinterpret_cast<const unsigned char *>(base + pos);
    if (block[0] != 0x1f || block[1] != 0x8b ||
        block[2] != 0x08 || block[3] != 0x04 ||
        block[10] != 0x06 || block[12] != 'B' ||
        block[13] != 'C' || block[14] != 0x02 ||
        block[15] != 0x00) {
        return -1;
    }
    size_t len = static_cast<size_t>(block[16]) |
                 (static_cast<size_t>(block[17]) << 8);
    ++len;
    if (len == 0 || len > BGZF_MAX_BLOCK_SIZE || len > size - pos) {
        return -1;
    }
    *block_len = len;
    return 0;
}

} // namespace

BgzfBlockBatch::BgzfBlockBatch()
    : blocks_(nullptr), data_(nullptr), capacity_(0) {}

BgzfBlockBatch::~BgzfBlockBatch() {
    Release();
}

int BgzfBlockBatch::Allocate(size_t capacity) {
    Release();
    if (capacity == 0 ||
        capacity > std::numeric_limits<size_t>::max() / sizeof(bam_block) ||
        capacity > std::numeric_limits<size_t>::max() / BGZF_MAX_BLOCK_SIZE) {
        return -1;
    }
    blocks_ = reinterpret_cast<bam_block *>(
        aligned_alloc_custom(64, capacity * sizeof(bam_block)));
    data_ = aligned_alloc_custom(64, capacity * BGZF_MAX_BLOCK_SIZE);
    if (!blocks_ || !data_) {
        Release();
        return -1;
    }
    memset(blocks_, 0, capacity * sizeof(bam_block));
    for (size_t i = 0; i < capacity; ++i) {
        blocks_[i].data = data_ + i * BGZF_MAX_BLOCK_SIZE;
        blocks_[i].block_id = static_cast<int>(i);
    }
    capacity_ = capacity;
    return 0;
}

void BgzfBlockBatch::Release() {
    if (blocks_) {
        aligned_free_custom(reinterpret_cast<unsigned char *>(blocks_));
    }
    if (data_) aligned_free_custom(data_);
    blocks_ = nullptr;
    data_ = nullptr;
    capacity_ = 0;
}

int LoadFileToMemory(const std::string &path, char **data, size_t *size) {
    if (!data || !size) return -1;
    *data = nullptr;
    *size = 0;

    FILE *fp = fopen(path.c_str(), "rb");
    if (!fp) return -1;
    if (fseeko(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        return -1;
    }
    off_t end = ftello(fp);
    if (end < 0 || fseeko(fp, 0, SEEK_SET) != 0 ||
        static_cast<unsigned long long>(end) >
            static_cast<unsigned long long>(SIZE_MAX)) {
        fclose(fp);
        return -1;
    }

    *size = static_cast<size_t>(end);
    *data = *size ? static_cast<char *>(malloc(*size)) : nullptr;
    if (*size > 0 && !*data) {
        fclose(fp);
        return -1;
    }
    if (*size > 0 && fread(*data, 1, *size, fp) != *size) {
        fclose(fp);
        free(*data);
        *data = nullptr;
        *size = 0;
        return -1;
    }
    fclose(fp);
    return 0;
}

int ScanBgzfBlocksInMemory(const char *base, size_t size, long long body_start,
                           std::vector<long long> *offsets,
                           std::vector<long long> *lengths) {
    if (!offsets || !lengths || !base || body_start < 0 ||
        static_cast<unsigned long long>(body_start) >
            static_cast<unsigned long long>(size)) {
        return -1;
    }

    long long pos = body_start;
    while (static_cast<unsigned long long>(pos) <
           static_cast<unsigned long long>(size)) {
        size_t block_len = 0;
        if (BgzfBlockLengthAt(base, size, static_cast<size_t>(pos),
                              &block_len) != 0) {
            return -1;
        }
        bool is_eof = block_len == sizeof(kBgzfEofBlock) &&
                      memcmp(base + pos, kBgzfEofBlock,
                             sizeof(kBgzfEofBlock)) == 0;
        if (is_eof) break;
        offsets->push_back(pos);
        lengths->push_back(static_cast<long long>(block_len));
        pos += static_cast<long long>(block_len);
    }
    return 0;
}

int SelectBlockRangeFromMemory(char *base, size_t input_size,
                               const std::vector<long long> &offsets,
                               const std::vector<long long> &lengths,
                               long long begin, long long end,
                               char **data, size_t *size) {
    if (!data || !size || begin < 0 || end < begin ||
        static_cast<size_t>(end) > offsets.size() ||
        offsets.size() != lengths.size()) {
        return -1;
    }
    *data = nullptr;
    *size = 0;
    if (begin >= end) return 0;

    long long start = offsets[static_cast<size_t>(begin)];
    long long stop = offsets[static_cast<size_t>(end - 1)] +
                     lengths[static_cast<size_t>(end - 1)];
    if (!base || start < 0 || stop < start ||
        static_cast<unsigned long long>(stop) >
            static_cast<unsigned long long>(input_size)) {
        return -1;
    }
    unsigned long long total = static_cast<unsigned long long>(stop - start);
    if (total > static_cast<unsigned long long>(SIZE_MAX)) return -1;
    *data = base + start;
    *size = static_cast<size_t>(total);
    return 0;
}

MemoryBamInput::MemoryBamInput()
    : data_(nullptr), size_(0), input_(nullptr), pending_hfile_(nullptr),
      header_(nullptr), format_(-1), body_offset_(0) {}

MemoryBamInput::~MemoryBamInput() {
    Close();
}

int MemoryBamInput::Open(const std::string &path) {
    if (Load(path) != 0) return -1;
    if (ParseHeader() != 0) {
        Close();
        return -1;
    }
    return 0;
}

int MemoryBamInput::Load(const std::string &path) {
    Close();
    if (LoadFileToMemory(path, &data_, &size_) != 0) return -1;
    return 0;
}

int MemoryBamInput::ParseHeader() {
    if (!data_ || size_ == 0 || input_ || pending_hfile_ || header_) {
        return -1;
    }
    pending_hfile_ = hopen("mem:", "rb:", data_, size_);
    if (!pending_hfile_) {
        free(data_);
        data_ = nullptr;
        size_ = 0;
        return -1;
    }
    input_ = reinterpret_cast<samFile *>(
        hts_hopen(pending_hfile_, "data", "rb"));
    if (!input_) {
        int close_ret = hclose(pending_hfile_);
        (void)close_ret;
        pending_hfile_ = nullptr;
        data_ = nullptr;
        size_ = 0;
        return -1;
    }
    pending_hfile_ = nullptr;

    header_ = sam_hdr_read(input_);
    if (!header_) {
        Close();
        return -1;
    }
    format_ = NormalizeFormat(input_->format.format);
    if (!input_->fp.bgzf || input_->fp.bgzf->block_address < 0 ||
        static_cast<unsigned long long>(input_->fp.bgzf->block_address) >
            static_cast<unsigned long long>(size_)) {
        Close();
        return -1;
    }
    body_offset_ = static_cast<uint64_t>(input_->fp.bgzf->block_address);
    return 0;
}

void MemoryBamInput::Close() {
    if (header_) sam_hdr_destroy(header_);
    header_ = nullptr;
    if (input_) {
        hts_close(input_);
        input_ = nullptr;
        data_ = nullptr;
    } else if (pending_hfile_) {
        int close_ret = hclose(pending_hfile_);
        (void)close_ret;
        pending_hfile_ = nullptr;
        data_ = nullptr;
    } else if (data_) {
        free(data_);
        data_ = nullptr;
    }
    size_ = 0;
    format_ = -1;
    body_offset_ = 0;
}

int MemoryBamInput::ScanBlocks(
        std::vector<BgzfBlockSpan> *spans) const {
    if (!spans || !data_) return -1;
    std::vector<long long> offsets;
    std::vector<long long> lengths;
    if (ScanBgzfBlocksInMemory(
            data_, size_, static_cast<long long>(body_offset_),
            &offsets, &lengths) != 0) {
        return -1;
    }
    spans->clear();
    spans->reserve(offsets.size());
    for (size_t i = 0; i < offsets.size(); ++i) {
        if (offsets[i] < 0 || lengths[i] <= 0 ||
            static_cast<unsigned long long>(lengths[i]) > UINT32_MAX) {
            return -1;
        }
        BgzfBlockSpan span;
        span.offset = static_cast<uint64_t>(offsets[i]);
        span.compressed_size = static_cast<uint32_t>(lengths[i]);
        spans->push_back(span);
    }
    return 0;
}

int MemoryBamInput::ReadBatch(const BgzfBlockSpan *spans, size_t count,
                              BgzfBlockBatch *batch) const {
    if (!batch || count > batch->capacity() || (count > 0 && !spans)) {
        return -1;
    }
    bam_block *blocks = batch->blocks();
    if (count > 0 && (!data_ || !blocks)) return -1;
    for (size_t i = 0; i < count; ++i) {
        uint64_t offset = spans[i].offset;
        size_t length = spans[i].compressed_size;
        if (length == 0 || length > BGZF_MAX_BLOCK_SIZE ||
            offset > size_ || length > size_ - static_cast<size_t>(offset)) {
            return -1;
        }
        memcpy(blocks[i].data, data_ + static_cast<size_t>(offset), length);
        blocks[i].length = static_cast<unsigned int>(length);
        blocks[i].pos = 0;
        blocks[i].errcode = 0;
        blocks[i].block_address = static_cast<int64_t>(offset);
        blocks[i].block_id = static_cast<int>(i);
    }
    return 0;
}

BgzfSpanBatchReader::BgzfSpanBatchReader(
        const BamInputBackend *input,
        const BgzfBlockSpan *spans, size_t count)
    : input_(input), spans_(spans), span_count_(count), position_(0) {}

int BgzfSpanBatchReader::ReadNext(BgzfBlockBatch *batch, size_t *count) {
    if (!batch || !count || !input_ ||
        (span_count_ > 0 && !spans_) || position_ > span_count_) {
        return -1;
    }
    size_t n = span_count_ - position_;
    if (n > batch->capacity()) n = batch->capacity();
    if (n > 0 && input_->ReadBatch(spans_ + position_, n, batch) != 0) {
        return -1;
    }
    position_ += n;
    *count = n;
    return 0;
}

} // namespace swbam
