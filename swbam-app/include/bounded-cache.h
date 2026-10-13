#ifndef SWBAM_APP_BOUNDED_CACHE_H
#define SWBAM_APP_BOUNDED_CACHE_H

#ifdef SWBAM_ENABLE_RAM_STORAGE
#include "storage.h"
#else
#include "storage-native.h"
#endif
#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <unistd.h>

namespace swbam_usage {

// A bounded resident cache over the selected temporary storage. The slow path
// is lazy: ordinary boundary groups do not create temporary files.
class CachedTemporary {
public:
    explicit CachedTemporary(size_t limit = 0)
        : fd_(-1), size_(0), capacity_(0), limit_(limit) {}
    ~CachedTemporary() { if (fd_ >= 0) CloseTemporary(fd_); free(data_); }
    uint64_t size() const { return size_; }
    size_t capacity() const { return capacity_; }
    int WriteAt(uint64_t offset, const void *data, size_t bytes) {
        if (offset > size_ || bytes > UINT64_MAX - offset ||
            offset + bytes > (uint64_t)INT64_MAX || (!data && bytes)) return -1;
        const uint64_t end = offset + bytes;
        if (fd_ < 0 && end <= limit_) {
            if (end > capacity_) {
                const size_t capacity = std::min<uint64_t>(limit_,
                    std::max<uint64_t>(end, std::max<size_t>(4096, capacity_ * 2)));
                void *next = realloc(data_, capacity);
                if (!next) return -1;
                data_ = static_cast<char *>(next); capacity_ = capacity;
            }
            if (bytes) memcpy(data_ + offset, data, bytes);
        } else {
            if (Spill() || Transfer(offset, const_cast<void *>(data), bytes, true)) return -1;
        }
        size_ = std::max(size_, end);
        return 0;
    }
    int Append(const void *data, size_t bytes) { return WriteAt(size_, data, bytes); }
    int ReadAt(uint64_t offset, void *data, size_t bytes) const {
        if (offset > size_ || bytes > size_ - offset || (!data && bytes)) return -1;
        if (fd_ >= 0) return Transfer(offset, data, bytes, false);
        if (bytes) memcpy(data, data_ + offset, bytes);
        return 0;
    }
    int Spill() {
        if (fd_ >= 0) return 0;
        static unsigned long serial = 0;
        const char *tmp = getenv("TMPDIR");
        const std::string path = std::string(tmp ? tmp : "/tmp") +
            "/swbam-cache." + std::to_string((long long)getpid()) + "." +
            std::to_string(++serial);
        fd_ = OpenTemporary(path);
        if (fd_ < 0) return -1;
        if (Transfer(0, data_, size_, true)) { CloseTemporary(fd_); fd_ = -1; return -1; }
        free(data_); data_ = nullptr; capacity_ = 0;
        return 0;
    }
private:
    int Transfer(uint64_t offset, void *data, size_t bytes, bool write) const {
        char *p = static_cast<char *>(data);
        while (bytes) {
            const size_t count = std::min<size_t>(bytes, 8u << 20);
            const ssize_t n = write ? WriteTemporary(fd_, p, count, offset)
                                    : ReadTemporary(fd_, p, count, offset);
            if (n < 0 && errno == EINTR) continue;
            if (n <= 0) return -1;
            bytes -= n; offset += n; p += n;
        }
        return 0;
    }
    CachedTemporary(const CachedTemporary &);
    CachedTemporary &operator=(const CachedTemporary &);
    int fd_;
    uint64_t size_;
    char *data_ = nullptr;
    size_t capacity_, limit_;
};
}
#endif
