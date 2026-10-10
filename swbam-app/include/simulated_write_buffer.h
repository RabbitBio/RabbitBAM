#ifndef SWBAM_APP_SIMULATED_WRITE_BUFFER_H
#define SWBAM_APP_SIMULATED_WRITE_BUFFER_H

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <limits>

namespace swbam_app {

class SimulatedWriteBuffer {
public:
    static size_t WindowBytes() { return 64ull * 1024 * 1024; }

    SimulatedWriteBuffer()
        : data_(nullptr), body_capacity_(0), size_(0), bounded_(false) {}
    ~SimulatedWriteBuffer() { std::free(data_); }

    void Reset() {
        std::free(data_);
        data_ = nullptr;
        body_capacity_ = size_ = 0;
        bounded_ = false;
    }

    int Allocate(size_t body_size, size_t envelope_size,
                 void *(*allocate)(size_t) = std::malloc) {
        Reset();
        body_capacity_ = body_size;
        if (envelope_size > std::numeric_limits<size_t>::max() - body_size) {
            errno = EOVERFLOW;
            return -1;
        }
        size_ = body_size + envelope_size;
        if (size_ == 0) return 0;
        data_ = static_cast<char *>(allocate(size_));
        // Preserve the full-copy fast path; only the extra simulation copy
        // falls back to a window. The actual rank body remains unchanged.
        if (!data_ && body_size > WindowBytes()) {
            body_capacity_ = WindowBytes();
            size_ = body_capacity_ + envelope_size;
            data_ = static_cast<char *>(allocate(size_));
            bounded_ = true;
        }
        return data_ ? 0 : -1;
    }

    template <typename Source>
    int CopyBodyWindowed(const Source &source, size_t prefix_size,
                         volatile unsigned long long *guard) {
        if (source.size() == 0) return 0;
        if (!data_ || !body_capacity_ || prefix_size > size_ ||
            body_capacity_ > size_ - prefix_size || !guard) {
            errno = EINVAL;
            return -1;
        }
        uint64_t copied = 0;
        char *destination = data_ + prefix_size;
        while (copied < source.size()) {
            const size_t count = static_cast<size_t>(std::min<uint64_t>(
                source.size() - copied, body_capacity_));
            if (source.ReadAt(copied, destination, count) != 0) return -1;
            *guard += static_cast<unsigned char>(destination[0]);
            *guard += static_cast<unsigned char>(destination[count - 1]);
            copied += count;
        }
        return 0;
    }

    char *data() const { return data_; }
    size_t body_capacity() const { return body_capacity_; }
    size_t size() const { return size_; }
    bool bounded() const { return bounded_; }

private:
    SimulatedWriteBuffer(const SimulatedWriteBuffer &);
    SimulatedWriteBuffer &operator=(const SimulatedWriteBuffer &);

    char *data_;
    size_t body_capacity_;
    size_t size_;
    bool bounded_;
};

} // namespace swbam_app

#endif
