#include "ram_file.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <limits>

namespace swbam_usage {
const size_t RamPool::kPage;
RamPool::RamPool() : data_(nullptr), used_(0), peak_(0), cursor_(0) {}
RamPool::~RamPool() { free(data_); }
int RamPool::Initialize(size_t capacity) {
    if (data_ || capacity < kPage) { errno = EINVAL; return -1; }
    const size_t count = capacity / kPage;
    pages_.assign(count, 0);
    data_ = static_cast<char *>(malloc(count * kPage));
    if (!data_) return -1;
    // -xmalloc routes ordinary malloc to the cross-CG heap before MPI_Init.
    for (size_t pos = 0; pos < count * kPage; pos += 4096)
        data_[pos] = 0;
    return 0;
}
int RamPool::Allocate(size_t count, size_t *first) {
    if (!first || !count || count > pages_.size() - used_) {
        errno = ENOSPC; return -1;
    }
    size_t run = 0;
    for (size_t step = 0; step < pages_.size(); ++step) {
        const size_t i = (cursor_ + step) % pages_.size();
        if (i == 0) run = 0;
        run = pages_[i] ? 0 : run + 1;
        if (run == count) {
            *first = i + 1 - count;
            for (size_t p = *first; p <= i; ++p) pages_[p] = 1;
            used_ += count;
            peak_ = std::max(peak_, used_);
            cursor_ = (i + 1) % pages_.size();
            return 0;
        }
    }
    // A contiguous request can straddle the search origin, not the pool end.
    run = 0;
    for (size_t i = 0; i < pages_.size(); ++i) {
        run = pages_[i] ? 0 : run + 1;
        if (run != count) continue;
        *first = i + 1 - count;
        for (size_t p = *first; p <= i; ++p) pages_[p] = 1;
        used_ += count;
        peak_ = std::max(peak_, used_);
        cursor_ = (i + 1) % pages_.size();
        return 0;
    }
    errno = ENOSPC; return -1;
}
void RamPool::Release(size_t first, size_t count) {
    for (size_t i = first; i < first + count; ++i) {
        if (pages_[i]) { pages_[i] = 0; --used_; }
    }
}
char *RamPool::Address(size_t page) const { return data_ + page * kPage; }
RamFile::RamFile(RamPool *pool) : pool_(pool), size_(0) {}
RamFile::~RamFile() {
    for (size_t i = released_pages_; i < pages_.size(); ++i)
        pool_->Release(pages_[i], 1);
}
int RamFile::ReleaseBefore(uint64_t offset) {
    if (offset > size_) { errno = EINVAL; return -1; }
    const size_t end = offset == size_ ? pages_.size() : offset / RamPool::kPage;
    while (released_pages_ < end) pool_->Release(pages_[released_pages_++], 1);
    return 0;
}
int RamFile::ReserveContiguous(size_t bytes) {
    if (!pages_.empty() || bytes > SIZE_MAX - RamPool::kPage + 1) {
        errno = EINVAL; return -1;
    }
    const size_t count = (bytes + RamPool::kPage - 1) / RamPool::kPage;
    if (!count) return 0;
    size_t first = 0;
    if (pool_->Allocate(count, &first)) return -1;
    try {
        pages_.reserve(count);
        for (size_t i = 0; i < count; ++i) pages_.push_back(first + i);
    } catch (...) { pool_->Release(first, count); pages_.clear(); errno = ENOMEM; return -1; }
    return 0;
}
int RamFile::Grow(uint64_t size) {
    if (size > SIZE_MAX - RamPool::kPage + 1) { errno = EFBIG; return -1; }
    const size_t needed = (size + RamPool::kPage - 1) / RamPool::kPage;
    const size_t old = pages_.size();
    if (needed <= old) return 0;
    try { pages_.reserve(needed); } catch (...) { errno = ENOMEM; return -1; }
    while (pages_.size() < needed) {
        size_t page = 0;
        if (pool_->Allocate(1, &page)) {
            while (pages_.size() > old) {
                pool_->Release(pages_.back(), 1); pages_.pop_back();
            }
            return -1;
        }
        pages_.push_back(page);
    }
    return 0;
}
void RamFile::Copy(uint64_t offset, void *buffer, size_t bytes, bool write) const {
    char *data = static_cast<char *>(buffer);
    while (bytes) {
        const size_t in_page = offset % RamPool::kPage;
        const size_t count = std::min(bytes, RamPool::kPage - in_page);
        char *address = pool_->Address(pages_[offset / RamPool::kPage]) + in_page;
        if (write) memcpy(address, data, count); else memcpy(data, address, count);
        offset += count; data += count; bytes -= count;
    }
}
int RamFile::Resize(uint64_t size) {
    if (released_pages_) { errno = EINVAL; return -1; }
    if (Grow(size)) return -1;
    if (size > size_) {
        uint64_t pos = size_;
        while (pos < size) {
            const size_t in_page = pos % RamPool::kPage;
            const size_t count = std::min<uint64_t>(size - pos, RamPool::kPage - in_page);
            memset(pool_->Address(pages_[pos / RamPool::kPage]) + in_page, 0, count);
            pos += count;
        }
    }
    size_ = size;
    const size_t needed = (size + RamPool::kPage - 1) / RamPool::kPage;
    while (pages_.size() > needed) {
        pool_->Release(pages_.back(), 1); pages_.pop_back();
    }
    return 0;
}
int RamFile::WriteAt(uint64_t offset, const void *data, size_t bytes) {
    if (released_pages_) { errno = EINVAL; return -1; }
    if ((!data && bytes) || offset > UINT64_MAX - bytes) { errno = EINVAL; return -1; }
    if (!bytes) return 0;
    if (Grow(offset + bytes)) return -1;
    // Only holes need zeroing; ordinary appends do not make a second full copy.
    uint64_t hole = size_;
    while (hole < offset) {
        const size_t in_page = hole % RamPool::kPage;
        const size_t count = std::min<uint64_t>(offset - hole, RamPool::kPage - in_page);
        memset(pool_->Address(pages_[hole / RamPool::kPage]) + in_page, 0, count);
        hole += count;
    }
    Copy(offset, const_cast<void *>(data), bytes, true);
    size_ = std::max(size_, offset + bytes);
    return 0;
}
int RamFile::ReadAt(uint64_t offset, void *data, size_t bytes) const {
    if (offset < released_pages_ * RamPool::kPage) { errno = EINVAL; return -1; }
    if ((!data && bytes) || offset > size_ || bytes > size_ - offset) {
        errno = EIO; return -1;
    }
    if (bytes) Copy(offset, data, bytes, false);
    return 0;
}
char *RamFile::Range(uint64_t offset, size_t bytes) const {
    if (offset < released_pages_ * RamPool::kPage) return nullptr;
    if (offset > size_ || bytes > size_ - offset || pages_.empty()) return nullptr;
    if (bytes) {
        const size_t first = offset / RamPool::kPage;
        const size_t last = (offset + bytes - 1) / RamPool::kPage;
        for (size_t i = first + 1; i <= last; ++i)
            if (pages_[i] != pages_[i - 1] + 1) return nullptr;
    }
    if (offset == size_ && offset / RamPool::kPage == pages_.size())
        return pool_->Address(pages_.back()) + RamPool::kPage;
    return pool_->Address(pages_[offset / RamPool::kPage]) + offset % RamPool::kPage;
}
}
