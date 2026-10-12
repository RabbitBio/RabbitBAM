#ifndef SWBAM_USAGE_RAM_FILE_H
#define SWBAM_USAGE_RAM_FILE_H

#include <cstddef>
#include <cstdint>
#include <vector>

namespace swbam_usage {

class RamPool {
public:
    static const size_t kPage = 4u * 1024u * 1024u;
    RamPool();
    ~RamPool();
    int Initialize(size_t capacity);
    int Allocate(size_t pages, size_t *first);
    void Release(size_t first, size_t pages);
    char *Address(size_t page) const;
    size_t capacity() const { return pages_.size() * kPage; }
    size_t used() const { return used_ * kPage; }
    size_t peak() const { return peak_ * kPage; }
private:
    RamPool(const RamPool &);
    RamPool &operator=(const RamPool &);
    char *data_;
    std::vector<unsigned char> pages_;
    size_t used_, peak_, cursor_;
};

class RamFile {
public:
    explicit RamFile(RamPool *pool);
    ~RamFile();
    int ReserveContiguous(size_t bytes);
    int WriteAt(uint64_t offset, const void *data, size_t bytes);
    int ReadAt(uint64_t offset, void *data, size_t bytes) const;
    int Resize(uint64_t size);
    char *Range(uint64_t offset, size_t bytes) const;
    uint64_t size() const { return size_; }
private:
    RamFile(const RamFile &);
    RamFile &operator=(const RamFile &);
    int Grow(uint64_t size);
    void Copy(uint64_t offset, void *data, size_t bytes, bool write) const;
    RamPool *pool_;
    std::vector<size_t> pages_;
    uint64_t size_;
};
}
#endif
