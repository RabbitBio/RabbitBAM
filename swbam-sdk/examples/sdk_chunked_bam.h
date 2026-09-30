#ifndef SWBAM_SDK_CHUNKED_BAM_H
#define SWBAM_SDK_CHUNKED_BAM_H

#include "swbam/io.h"

#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <mpi.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

namespace sdk_bench {

inline size_t ChunkLimit(const char *argument) {
    if (!argument) return 0;
    if (strcmp(argument, "--chunked-memory") == 0 ||
        strcmp(argument, "--chunked-memory-progress") == 0)
        return size_t(2) << 30;
    const char prefix[] = "--chunked-memory=";
    if (strncmp(argument, prefix, sizeof(prefix) - 1) != 0) return 0;
    const char *value = argument + sizeof(prefix) - 1;
    if (!*value || *value == '-') return 0;
    char *end = nullptr;
    errno = 0;
    const unsigned long long mib = strtoull(value, &end, 10);
    if (errno || *end || mib < 1 || mib > 2048) return 0;
    return static_cast<size_t>(mib) << 20;
}

class ChunkedBamInput : public swbam::BamInputBackend {
public:
    ChunkedBamInput() : fd_(-1), base_(0), length_(0), index_ready_(false) {}
    ~ChunkedBamInput() { Close(); }

    int Open(const std::string &path) {
        Close();
        if (metadata_.Open(path) != 0) return -1;
        fd_ = open(path.c_str(), O_RDONLY);
        if (fd_ < 0) {
            metadata_.Close();
            return -1;
        }
        path_ = path;
        return 0;
    }

    void Close() {
        if (fd_ >= 0) close(fd_);
        fd_ = -1;
        data_.clear();
        index_.clear();
        index_ready_ = false;
        path_.clear();
        metadata_.Close();
        base_ = 0;
        length_ = 0;
    }

    const sam_hdr_t *header() const { return metadata_.header(); }
    int format() const { return metadata_.format(); }
    size_t size() const { return metadata_.size(); }
    uint64_t body_offset() const { return metadata_.body_offset(); }
    int ScanBlocks(std::vector<swbam::BgzfBlockSpan> *spans) const {
        if (!spans) return -1;
        if (index_ready_) {
            *spans = index_;
            return 0;
        }
        return metadata_.ScanBlocks(spans);
    }

    int PrepareIndex() {
        if (fd_ < 0 || index_ready_) return -1;
        const char *cache = std::getenv("SWBAM_SDK_INDEX_CACHE");
        if (!cache || !*cache || LoadIndex(cache) != 0) {
            if (metadata_.ScanBlocks(&index_) != 0) {
                fprintf(stderr, "ERROR: BGZF index scan failed for %s.\n",
                        path_.c_str());
                return -1;
            }
            if (cache && *cache && SaveIndex(cache) != 0) {
                fprintf(stderr, "ERROR: cannot save BGZF index to %s: %s.\n",
                        cache, strerror(errno));
                return -1;
            }
        }
        index_ready_ = true;
        return 0;
    }

    int Stage(const swbam::BgzfBlockSpan *spans, size_t count,
              size_t max_bytes,
              double *read_seconds, double *prepare_seconds) {
        if (fd_ < 0 || !spans || count == 0) return -1;
        const uint64_t begin = spans[0].offset;
        const uint64_t end = spans[count - 1].offset +
                             spans[count - 1].compressed_size;
        if (end < begin || end > size() || end - begin > SIZE_MAX) return -1;
        const size_t length = static_cast<size_t>(end - begin);
        if (length == 0 || length > max_bytes) return -1;
        const char *progress_env = std::getenv("SWBAM_SDK_PROGRESS");
        const bool progress = progress_env && strcmp(progress_env, "1") == 0;
        int rank = 0;
        if (progress) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const double prepare_start = MPI_Wtime();
        try {
            if (data_.capacity() < max_bytes) data_.reserve(max_bytes);
            data_.resize(length);
        } catch (...) {
            fprintf(stderr, "ERROR: cannot allocate %zu bytes for BAM chunk.\n",
                    length);
            return -1;
        }
        if (prepare_seconds) *prepare_seconds += MPI_Wtime() - prepare_start;
        if (progress) {
            fprintf(stderr, "[sdk-stage rank=%d] prepared=%zu elapsed=%.3f\n",
                    rank, length, MPI_Wtime() - prepare_start);
            fflush(stderr);
        }
        const double read_start = MPI_Wtime();
        size_t done = 0;
        size_t next_report = size_t(256) << 20;
        while (done < length) {
            const size_t remaining = length - done;
            const size_t amount = remaining < 8u * 1024u * 1024u
                                  ? remaining : 8u * 1024u * 1024u;
            const ssize_t got = pread(fd_, data_.data() + done, amount,
                                      static_cast<off_t>(begin + done));
            if (got < 0 && errno == EINTR) continue;
            if (got <= 0) {
                fprintf(stderr, "ERROR: BAM chunk pread failed at offset %llu "
                        "size %zu: %s.\n",
                        static_cast<unsigned long long>(begin + done), amount,
                        got < 0 ? strerror(errno) : "unexpected EOF");
                return -1;
            }
            done += static_cast<size_t>(got);
            if (progress && done >= next_report) {
                fprintf(stderr,
                        "[sdk-stage rank=%d] read=%zu/%zu elapsed=%.3f\n",
                        rank, done, length, MPI_Wtime() - read_start);
                fflush(stderr);
                next_report += size_t(256) << 20;
            }
        }
        if (read_seconds) *read_seconds += MPI_Wtime() - read_start;
        base_ = begin;
        length_ = length;
        return 0;
    }

    int ReadBatch(const swbam::BgzfBlockSpan *spans, size_t count,
                  swbam::BgzfBlockBatch *batch) const {
        if (!batch || count > batch->capacity() || (count && !spans)) return -1;
        bam_block *blocks = batch->blocks();
        if (count && !blocks) return -1;
        for (size_t i = 0; i < count; ++i) {
            const uint64_t offset = spans[i].offset;
            const size_t length = spans[i].compressed_size;
            if (!length || length > BGZF_MAX_BLOCK_SIZE || offset < base_ ||
                offset - base_ > length_ ||
                length > length_ - static_cast<size_t>(offset - base_)) {
                return -1;
            }
            memcpy(blocks[i].data, data_.data() +
                   static_cast<size_t>(offset - base_), length);
            blocks[i].length = static_cast<unsigned int>(length);
            blocks[i].pos = 0;
            blocks[i].errcode = 0;
            blocks[i].block_address = static_cast<int64_t>(offset);
            blocks[i].block_id = static_cast<int>(i);
        }
        return 0;
    }

    size_t peak_capacity() const { return data_.capacity(); }

private:
    struct IndexHeader {
        uint64_t magic;
        uint64_t file_size;
        uint64_t body_offset;
        uint64_t mtime;
        uint64_t block_count;
    };
    static const uint64_t kIndexMagic = UINT64_C(0x535742414d424731);

    IndexHeader IndexIdentity() const {
        struct stat details = {};
        IndexHeader header = {};
        if (stat(path_.c_str(), &details) == 0) {
            header.magic = kIndexMagic;
            header.file_size = size();
            header.body_offset = body_offset();
            header.mtime = static_cast<uint64_t>(details.st_mtime);
        }
        return header;
    }

    int LoadIndex(const char *path) {
        FILE *file = fopen(path, "rb");
        if (!file) return -1;
        IndexHeader saved = {}, current = IndexIdentity();
        int ok = fread(&saved, sizeof(saved), 1, file) == 1 &&
                 saved.magic == current.magic &&
                 saved.file_size == current.file_size &&
                 saved.body_offset == current.body_offset &&
                 saved.mtime == current.mtime &&
                 saved.block_count <= static_cast<uint64_t>(INT_MAX);
        if (ok) {
            try {
                index_.resize(static_cast<size_t>(saved.block_count));
            } catch (...) {
                ok = 0;
            }
        }
        if (ok && !index_.empty()) {
            ok = fread(index_.data(), sizeof(index_[0]), index_.size(), file)
                 == index_.size();
        }
        if (ok) {
            for (size_t i = 0; i < index_.size(); ++i) {
                const swbam::BgzfBlockSpan &span = index_[i];
                if (span.compressed_size == 0 ||
                    span.compressed_size > BGZF_MAX_BLOCK_SIZE ||
                    span.offset < body_offset() ||
                    span.offset > size() ||
                    span.compressed_size > size() - span.offset ||
                    (i && index_[i - 1].offset +
                          index_[i - 1].compressed_size != span.offset)) {
                    ok = 0;
                    break;
                }
            }
        }
        if (fclose(file) != 0) ok = 0;
        if (!ok) index_.clear();
        return ok ? 0 : -1;
    }

    int SaveIndex(const char *path) const {
        const std::string temporary = std::string(path) + "." +
                                      std::to_string(static_cast<long>(getpid()));
        FILE *file = fopen(temporary.c_str(), "wb");
        if (!file) return -1;
        IndexHeader header = IndexIdentity();
        header.block_count = index_.size();
        int ok = header.magic == kIndexMagic &&
                 fwrite(&header, sizeof(header), 1, file) == 1;
        if (ok && !index_.empty()) {
            ok = fwrite(index_.data(), sizeof(index_[0]), index_.size(), file)
                 == index_.size();
        }
        if (fclose(file) != 0) ok = 0;
        if (ok && rename(temporary.c_str(), path) != 0) ok = 0;
        if (!ok) unlink(temporary.c_str());
        return ok ? 0 : -1;
    }

    ChunkedBamInput(const ChunkedBamInput &);
    ChunkedBamInput &operator=(const ChunkedBamInput &);

    swbam::PosixBamInput metadata_;
    int fd_;
    std::string path_;
    std::vector<unsigned char> data_;
    std::vector<swbam::BgzfBlockSpan> index_;
    uint64_t base_;
    size_t length_;
    bool index_ready_;
};

inline size_t NextChunkEnd(const swbam::BgzfBlockSpan *spans,
                           size_t count, size_t begin, size_t max_bytes) {
    if (begin >= count) return begin;
    size_t end = begin + 1;
    const uint64_t first = spans[begin].offset;
    while (end < count && spans[end].offset - first +
           spans[end].compressed_size <= max_bytes) ++end;
    return end;
}

class CircularBodySink : public swbam::RankBodySink {
public:
    CircularBodySink() : size_(0), position_(0) {}
    int Open(size_t capacity) {
        if (!capacity) return -1;
        try {
            data_.resize(capacity);
        } catch (...) {
            return -1;
        }
        size_ = 0;
        position_ = 0;
        return 0;
    }
    int Append(const void *source, size_t length) {
        if ((!source && length) || data_.empty() ||
            length > std::numeric_limits<uint64_t>::max() - size_) return -1;
        const unsigned char *bytes = static_cast<const unsigned char *>(source);
        size_t remaining = length;
        while (remaining) {
            const size_t amount = remaining < data_.size() - position_
                                  ? remaining : data_.size() - position_;
            memcpy(data_.data() + position_, bytes, amount);
            bytes += amount;
            remaining -= amount;
            position_ += amount;
            if (position_ == data_.size()) position_ = 0;
        }
        size_ += length;
        return 0;
    }
    int Flush() { return 0; }
    uint64_t size() const { return size_; }
    int ReadAt(uint64_t, void *, size_t) const { return -1; }
    bool is_memory() const { return true; }
    size_t capacity() const { return data_.capacity(); }

private:
    std::vector<unsigned char> data_;
    uint64_t size_;
    size_t position_;
};

} // namespace sdk_bench
#endif
