#include "swbam/io.h"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <sys/stat.h>
#include <sys/uio.h>
#include <unistd.h>

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

int ReadAt(int fd, uint64_t offset, void *data, size_t size) {
    if (fd < 0 || (!data && size != 0) || offset > INT64_MAX) return -1;
    unsigned char *bytes = static_cast<unsigned char *>(data);
    size_t done = 0;
    while (done < size) {
        if (offset + done > (uint64_t)INT64_MAX) return -1;
        const ssize_t result = pread(
            fd, bytes + done, size - done, (off_t)(offset + done));
        if (result < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        if (result == 0) return -1;
        done += (size_t)result;
    }
    return 0;
}

int ReadVectorAt(int fd, uint64_t offset,
                 struct iovec *vectors, int count) {
    if (fd < 0 || !vectors || count <= 0 || offset > INT64_MAX) return -1;
    int first = 0;
    uint64_t done = 0;
    while (first < count) {
        if (offset + done > (uint64_t)INT64_MAX) return -1;
        const ssize_t result = preadv(
            fd, vectors + first, count - first,
            (off_t)(offset + done));
        if (result < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        if (result == 0) return -1;
        done += (uint64_t)result;
        size_t consumed = (size_t)result;
        while (first < count && consumed >= vectors[first].iov_len) {
            consumed -= vectors[first].iov_len;
            first++;
        }
        if (consumed > 0 && first < count) {
            vectors[first].iov_base =
                static_cast<unsigned char *>(vectors[first].iov_base) +
                consumed;
            vectors[first].iov_len -= consumed;
        }
    }
    return 0;
}

int BgzfBlockLengthFromHeader(const unsigned char *block,
                              size_t *block_len) {
    if (!block || !block_len ||
        block[0] != 0x1f || block[1] != 0x8b ||
        block[2] != 0x08 || block[3] != 0x04 ||
        block[10] != 0x06 || block[12] != 'B' ||
        block[13] != 'C' || block[14] != 0x02 ||
        block[15] != 0x00) {
        return -1;
    }
    const size_t length =
        (size_t)block[16] | ((size_t)block[17] << 8);
    const size_t total = length + 1;
    if (total < BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH ||
        total > BGZF_MAX_BLOCK_SIZE) {
        return -1;
    }
    *block_len = total;
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

int FindSamBodyStartInMemory(const char *base, size_t size,
                             long long *body_start) {
    if ((!base && size > 0) || !body_start ||
        size > static_cast<size_t>(
            std::numeric_limits<long long>::max())) {
        return -1;
    }
    size_t pos = 0;
    while (pos < size) {
        const size_t line_start = pos;
        while (pos < size && base[pos] != '\n') pos++;
        size_t length = pos - line_start;
        if (length > 0 && base[line_start + length - 1] == '\r') length--;
        if (length > 0 && base[line_start] != '@') {
            *body_start = static_cast<long long>(line_start);
            return 0;
        }
        if (pos < size) pos++;
    }
    *body_start = static_cast<long long>(size);
    return 0;
}

int SplitSamRangesInMemory(const char *base, size_t size,
                           long long body_start, int range_count,
                           std::vector<long long> *offsets,
                           std::vector<long long> *lengths) {
    if ((!base && size > 0) || body_start < 0 || !offsets || !lengths ||
        static_cast<unsigned long long>(body_start) >
            static_cast<unsigned long long>(size) ||
        size > static_cast<size_t>(
            std::numeric_limits<long long>::max()) || range_count <= 0) {
        return -1;
    }

    std::vector<long long> boundaries(
        static_cast<size_t>(range_count) + 1, body_start);
    const long long file_end = static_cast<long long>(size);
    const long long body_len = file_end - body_start;
    boundaries[0] = body_start;
    boundaries[static_cast<size_t>(range_count)] = file_end;
    for (int r = 1; r < range_count; ++r) {
        long long pos = body_start + body_len * r / range_count;
        while (pos > body_start && pos < file_end &&
               base[static_cast<size_t>(pos - 1)] != '\n') {
            pos++;
        }
        if (pos < boundaries[static_cast<size_t>(r - 1)]) {
            pos = boundaries[static_cast<size_t>(r - 1)];
        }
        boundaries[static_cast<size_t>(r)] = pos;
    }

    offsets->resize(static_cast<size_t>(range_count));
    lengths->resize(static_cast<size_t>(range_count));
    for (int r = 0; r < range_count; ++r) {
        const long long start = boundaries[static_cast<size_t>(r)];
        long long stop = boundaries[static_cast<size_t>(r + 1)];
        if (stop < start) stop = start;
        (*offsets)[static_cast<size_t>(r)] = start;
        (*lengths)[static_cast<size_t>(r)] = stop - start;
    }
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

int MemoryBamInput::OpenMemoryCopy(const void *data, size_t size) {
    if (LoadMemoryCopy(data, size) != 0) return -1;
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

int MemoryBamInput::LoadMemoryCopy(const void *data, size_t size) {
    if (!data || size == 0) return -1;
    Close();
    data_ = static_cast<char *>(malloc(size));
    if (!data_) return -1;
    memcpy(data_, data, size);
    size_ = size;
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

class PosixBamInput::Impl {
public:
    Impl()
        : fd(-1), file_size(0), header(nullptr), input_format(-1),
          bam_body_offset(0) {}

    int fd;
    size_t file_size;
    sam_hdr_t *header;
    int input_format;
    uint64_t bam_body_offset;
    mutable std::vector<struct iovec> vectors;
};

PosixBamInput::PosixBamInput() : impl_(new Impl()) {}

PosixBamInput::~PosixBamInput() {
    Close();
    delete impl_;
}

int PosixBamInput::Open(const std::string &path) {
    if (path.empty()) return -1;
    Close();

    struct stat file_stat;
    if (stat(path.c_str(), &file_stat) != 0 || file_stat.st_size <= 0 ||
        (unsigned long long)file_stat.st_size >
            (unsigned long long)SIZE_MAX) {
        return -1;
    }

    samFile *input = sam_open(path.c_str(), "rb");
    if (!input) return -1;
    sam_hdr_t *parsed_header = sam_hdr_read(input);
    const int parsed_format = NormalizeFormat(input->format.format);
    uint64_t parsed_body_offset = 0;
    int valid = parsed_header && parsed_format == bam && input->fp.bgzf &&
        input->fp.bgzf->block_address >= 0;
    if (valid) {
        parsed_body_offset =
            (uint64_t)input->fp.bgzf->block_address;
        valid = parsed_body_offset <= (uint64_t)file_stat.st_size;
    }
    if (sam_close(input) != 0) valid = 0;
    if (!valid) {
        if (parsed_header) sam_hdr_destroy(parsed_header);
        return -1;
    }

    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        sam_hdr_destroy(parsed_header);
        return -1;
    }
    impl_->fd = fd;
    impl_->file_size = (size_t)file_stat.st_size;
    impl_->header = parsed_header;
    impl_->input_format = parsed_format;
    impl_->bam_body_offset = parsed_body_offset;
    return 0;
}

void PosixBamInput::Close() {
    if (impl_->fd >= 0) {
        close(impl_->fd);
        impl_->fd = -1;
    }
    if (impl_->header) {
        sam_hdr_destroy(impl_->header);
        impl_->header = nullptr;
    }
    impl_->file_size = 0;
    impl_->input_format = -1;
    impl_->bam_body_offset = 0;
}

const sam_hdr_t *PosixBamInput::header() const {
    return impl_->header;
}

int PosixBamInput::format() const {
    return impl_->input_format;
}

size_t PosixBamInput::size() const {
    return impl_->file_size;
}

uint64_t PosixBamInput::body_offset() const {
    return impl_->bam_body_offset;
}

int PosixBamInput::ScanBlocks(
        std::vector<BgzfBlockSpan> *spans) const {
    if (!spans || impl_->fd < 0 || !impl_->header ||
        impl_->bam_body_offset > impl_->file_size) {
        return -1;
    }
    spans->clear();
    uint64_t offset = impl_->bam_body_offset;
    unsigned char header_bytes[BLOCK_HEADER_LENGTH];
    unsigned char eof_bytes[sizeof(kBgzfEofBlock)];
    while (offset < impl_->file_size) {
        if (impl_->file_size - (size_t)offset < BLOCK_HEADER_LENGTH ||
            ReadAt(impl_->fd, offset, header_bytes,
                   sizeof(header_bytes)) != 0) {
            return -1;
        }
        size_t block_size = 0;
        if (BgzfBlockLengthFromHeader(
                header_bytes, &block_size) != 0 ||
            block_size > impl_->file_size - (size_t)offset) {
            return -1;
        }
        if (block_size == sizeof(kBgzfEofBlock)) {
            if (ReadAt(impl_->fd, offset, eof_bytes,
                       sizeof(eof_bytes)) != 0) {
                return -1;
            }
            if (memcmp(eof_bytes, kBgzfEofBlock,
                       sizeof(kBgzfEofBlock)) == 0) {
                break;
            }
        }
        BgzfBlockSpan span;
        span.offset = offset;
        span.compressed_size = (uint32_t)block_size;
        spans->push_back(span);
        offset += block_size;
    }
    return 0;
}

int PosixBamInput::ReadBatch(
        const BgzfBlockSpan *spans, size_t count,
        BgzfBlockBatch *batch) const {
    if (!batch || count > batch->capacity() ||
        (count > 0 && !spans) || impl_->fd < 0) {
        return -1;
    }
    bam_block *blocks = batch->blocks();
    if (count > 0 && !blocks) return -1;

    bool contiguous = count > 0;
    uint64_t range_offset = count > 0 ? spans[0].offset : 0;
    size_t range_size = 0;
    for (size_t i = 0; i < count; ++i) {
        const uint64_t offset = spans[i].offset;
        const size_t length = spans[i].compressed_size;
        if (length < BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH ||
            length > BGZF_MAX_BLOCK_SIZE || offset > impl_->file_size ||
            length > impl_->file_size - (size_t)offset ||
            length > SIZE_MAX - range_size) {
            return -1;
        }
        if (i > 0 && offset != range_offset + range_size) contiguous = false;
        range_size += length;
    }

    if (contiguous) {
        impl_->vectors.resize(count);
        for (size_t i = 0; i < count; ++i) {
            impl_->vectors[i].iov_base = blocks[i].data;
            impl_->vectors[i].iov_len = spans[i].compressed_size;
        }
        if (ReadVectorAt(impl_->fd, range_offset,
                         impl_->vectors.data(), (int)count) != 0) {
            return -1;
        }
    }

    for (size_t i = 0; i < count; ++i) {
        const uint64_t offset = spans[i].offset;
        const size_t length = spans[i].compressed_size;
        if (!contiguous && ReadAt(
                       impl_->fd, offset,
                       blocks[i].data, length) != 0) {
            return -1;
        }
        blocks[i].length = (unsigned int)length;
        blocks[i].pos = 0;
        blocks[i].errcode = 0;
        blocks[i].block_address = (int64_t)offset;
        blocks[i].block_id = (int)i;
    }
    return 0;
}

class PosixSamInput::Impl {
public:
    Impl() : fd(-1), file_size(0), header(nullptr), sam_body_offset(0) {}

    int fd;
    size_t file_size;
    sam_hdr_t *header;
    uint64_t sam_body_offset;
};

static int FindSamBodyOffsetAt(int fd, size_t file_size,
                               uint64_t *body_offset) {
    if (fd < 0 || !body_offset) return -1;
    const size_t chunk_size = 64 * 1024;
    std::vector<char> buffer(chunk_size);
    uint64_t file_offset = 0;
    uint64_t line_start = 0;
    size_t line_length = 0;
    char first = 0;
    char last = 0;
    while (file_offset < file_size) {
        const size_t length = std::min(
            chunk_size, file_size - static_cast<size_t>(file_offset));
        if (ReadAt(fd, file_offset, buffer.data(), length) != 0) return -1;
        for (size_t i = 0; i < length; ++i) {
            const char value = buffer[i];
            if (value == '\n') {
                const size_t effective =
                    line_length > 0 && last == '\r'
                        ? line_length - 1 : line_length;
                if (effective > 0 && first != '@') {
                    *body_offset = line_start;
                    return 0;
                }
                line_start = file_offset + i + 1;
                line_length = 0;
                first = 0;
                last = 0;
            } else {
                if (line_length == 0) first = value;
                last = value;
                line_length++;
            }
        }
        file_offset += length;
    }
    const size_t effective =
        line_length > 0 && last == '\r' ? line_length - 1 : line_length;
    *body_offset = effective > 0 && first != '@'
        ? line_start : static_cast<uint64_t>(file_size);
    return 0;
}

static int FindSamLineBoundaryAt(int fd, size_t file_size,
                                 uint64_t target,
                                 uint64_t *boundary) {
    if (fd < 0 || !boundary || target > file_size) return -1;
    if (target == 0 || target == file_size) {
        *boundary = target;
        return 0;
    }
    char previous = 0;
    if (ReadAt(fd, target - 1, &previous, 1) != 0) return -1;
    if (previous == '\n') {
        *boundary = target;
        return 0;
    }

    const size_t chunk_size = 64 * 1024;
    std::vector<char> buffer(chunk_size);
    uint64_t offset = target;
    while (offset < file_size) {
        const size_t length = std::min(
            chunk_size, file_size - static_cast<size_t>(offset));
        if (ReadAt(fd, offset, buffer.data(), length) != 0) return -1;
        for (size_t i = 0; i < length; ++i) {
            if (buffer[i] == '\n') {
                *boundary = offset + i + 1;
                return 0;
            }
        }
        offset += length;
    }
    *boundary = static_cast<uint64_t>(file_size);
    return 0;
}

PosixSamInput::PosixSamInput() : impl_(new Impl()) {}

PosixSamInput::~PosixSamInput() {
    Close();
    delete impl_;
}

int PosixSamInput::Open(const std::string &path) {
    if (path.empty()) return -1;
    Close();
    struct stat file_stat;
    if (stat(path.c_str(), &file_stat) != 0 || file_stat.st_size <= 0 ||
        static_cast<unsigned long long>(file_stat.st_size) >
            static_cast<unsigned long long>(SIZE_MAX)) {
        return -1;
    }

    samFile *input = sam_open(path.c_str(), "r");
    if (!input) return -1;
    sam_hdr_t *header = sam_hdr_read(input);
    const int format = NormalizeFormat(input->format.format);
    const int close_result = sam_close(input);
    if (!header || format != sam || close_result != 0) {
        if (header) sam_hdr_destroy(header);
        return -1;
    }

    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        sam_hdr_destroy(header);
        return -1;
    }
    uint64_t body_offset = 0;
    if (FindSamBodyOffsetAt(
            fd, static_cast<size_t>(file_stat.st_size),
            &body_offset) != 0) {
        close(fd);
        sam_hdr_destroy(header);
        return -1;
    }
    impl_->fd = fd;
    impl_->file_size = static_cast<size_t>(file_stat.st_size);
    impl_->header = header;
    impl_->sam_body_offset = body_offset;
    return 0;
}

void PosixSamInput::Close() {
    if (impl_->fd >= 0) {
        close(impl_->fd);
        impl_->fd = -1;
    }
    if (impl_->header) {
        sam_hdr_destroy(impl_->header);
        impl_->header = nullptr;
    }
    impl_->file_size = 0;
    impl_->sam_body_offset = 0;
}

const sam_hdr_t *PosixSamInput::header() const { return impl_->header; }
size_t PosixSamInput::size() const { return impl_->file_size; }
uint64_t PosixSamInput::body_offset() const { return impl_->sam_body_offset; }

int PosixSamInput::SplitRanges(
        int range_count, std::vector<long long> *offsets,
        std::vector<long long> *lengths) const {
    if (impl_->fd < 0 || !impl_->header || range_count <= 0 ||
        !offsets || !lengths ||
        impl_->file_size > static_cast<size_t>(
            std::numeric_limits<long long>::max())) {
        return -1;
    }
    std::vector<uint64_t> boundaries(
        static_cast<size_t>(range_count) + 1, impl_->sam_body_offset);
    boundaries.back() = static_cast<uint64_t>(impl_->file_size);
    const uint64_t body_size =
        static_cast<uint64_t>(impl_->file_size) - impl_->sam_body_offset;
    const uint64_t quotient = body_size / static_cast<uint64_t>(range_count);
    const uint64_t remainder = body_size % static_cast<uint64_t>(range_count);
    for (int rank = 1; rank < range_count; ++rank) {
        const uint64_t r = static_cast<uint64_t>(rank);
        const uint64_t target = impl_->sam_body_offset + quotient * r +
            remainder * r / static_cast<uint64_t>(range_count);
        if (FindSamLineBoundaryAt(
                impl_->fd, impl_->file_size, target,
                &boundaries[static_cast<size_t>(rank)]) != 0) {
            return -1;
        }
    }
    offsets->resize(static_cast<size_t>(range_count));
    lengths->resize(static_cast<size_t>(range_count));
    for (int rank = 0; rank < range_count; ++rank) {
        const uint64_t begin = boundaries[static_cast<size_t>(rank)];
        const uint64_t end = boundaries[static_cast<size_t>(rank + 1)];
        if (end < begin || end > static_cast<uint64_t>(
                std::numeric_limits<long long>::max())) {
            return -1;
        }
        (*offsets)[static_cast<size_t>(rank)] = static_cast<long long>(begin);
        (*lengths)[static_cast<size_t>(rank)] =
            static_cast<long long>(end - begin);
    }
    return 0;
}

int PosixSamInput::ReadRange(uint64_t offset, size_t length,
                             std::vector<char> *data) const {
    if (impl_->fd < 0 || !data || offset > impl_->file_size ||
        length > impl_->file_size - static_cast<size_t>(offset)) {
        return -1;
    }
    try {
        data->resize(length);
    } catch (...) {
        return -1;
    }
    return length > 0
        ? ReadAt(impl_->fd, offset, data->data(), length) : 0;
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

MemoryBamOutput::MemoryBamOutput() {}

void MemoryBamOutput::Clear() {
    data_.clear();
}

int MemoryBamOutput::Reserve(size_t capacity) {
    try {
        data_.reserve(capacity);
    } catch (...) {
        return -1;
    }
    return 0;
}

int MemoryBamOutput::Write(const void *data, size_t size) {
    if (!data && size != 0) return -1;
    if (size == 0) return 0;
    if (size > data_.max_size() - data_.size()) return -1;
    try {
        const unsigned char *bytes =
            static_cast<const unsigned char *>(data);
        data_.insert(data_.end(), bytes, bytes + size);
    } catch (...) {
        return -1;
    }
    return 0;
}

int MemoryBamOutput::Append(const void *data, size_t size) {
    return Write(data, size);
}

int WriteBgzfBlocks(BamOutputBackend *output,
                    const bam_block *blocks, size_t count) {
    if (!output || (!blocks && count != 0)) return -1;
    for (size_t i = 0; i < count; ++i) {
        if (!blocks[i].data || blocks[i].length == 0 ||
            blocks[i].length > BGZF_MAX_BLOCK_SIZE ||
            output->Write(blocks[i].data, blocks[i].length) != 0) {
            return -1;
        }
    }
    return 0;
}

int WriteBgzfEof(BamOutputBackend *output) {
    return output
        ? output->Write(kBgzfEofBlock, sizeof(kBgzfEofBlock)) : -1;
}

int MemoryBamOutput::AppendBlocks(const bam_block *blocks, size_t count) {
    return WriteBgzfBlocks(this, blocks, count);
}

int MemoryBamOutput::AppendBgzfEof() {
    return WriteBgzfEof(this);
}

class PosixBamOutput::Impl {
public:
    Impl() : file(nullptr), written(0) {}

    FILE *file;
    uint64_t written;
};

PosixBamOutput::PosixBamOutput() : impl_(new Impl()) {}

PosixBamOutput::~PosixBamOutput() {
    (void)Close();
    delete impl_;
}

int PosixBamOutput::Open(const std::string &path) {
    if (path.empty()) return -1;
    if (Close() != 0) return -1;
    impl_->file = fopen(path.c_str(), "wb");
    impl_->written = 0;
    return impl_->file ? 0 : -1;
}

int PosixBamOutput::Close() {
    if (impl_->file) {
        const int result = fclose(impl_->file);
        impl_->file = nullptr;
        return result == 0 ? 0 : -1;
    }
    return 0;
}

int PosixBamOutput::Write(const void *data, size_t size) {
    if (!impl_->file || (!data && size != 0) ||
        size > UINT64_MAX - impl_->written) {
        return -1;
    }
    if (size == 0) return 0;
    const unsigned char *bytes =
        static_cast<const unsigned char *>(data);
    size_t offset = 0;
    while (offset < size) {
        const size_t written =
            fwrite(bytes + offset, 1, size - offset, impl_->file);
        if (written == 0) return -1;
        offset += written;
    }
    impl_->written += size;
    return 0;
}

int PosixBamOutput::Flush() {
    return impl_->file && fflush(impl_->file) == 0 ? 0 : -1;
}

uint64_t PosixBamOutput::bytes_written() const {
    return impl_->written;
}

} // namespace swbam
