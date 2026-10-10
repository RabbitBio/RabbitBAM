// Host-only regression of the app's output simulation, not an x86 app build.
#include "../swbam-app/include/simulated_write_buffer.h"

#include <cstdio>
#include <cstring>

static void Check(bool ok, const char *message) {
    if (ok) return;
    fprintf(stderr, "FAIL simulated write buffer: %s\n", message);
    std::exit(1);
}

static size_t allocation_limit;
static size_t allocation_calls;

static void *LimitedMalloc(size_t size) {
    ++allocation_calls;
    if (size > allocation_limit) {
        errno = ENOMEM;
        return nullptr;
    }
    return std::malloc(size);
}

struct PatternSource {
    uint64_t bytes;
    mutable uint64_t copied;
    mutable size_t calls;
    size_t fail_call;

    explicit PatternSource(uint64_t size, size_t fail = 0)
        : bytes(size), copied(0), calls(0), fail_call(fail) {}

    uint64_t size() const { return bytes; }
    int ReadAt(uint64_t offset, void *data, size_t count) const {
        ++calls;
        Check(offset == copied && offset <= bytes && count <= bytes - offset,
              "contiguous reads and exact bounds");
        if (calls == fail_call) return -1;
        memset(data, static_cast<unsigned char>(calls), count);
        copied += count;
        return 0;
    }
};

int main() {
    typedef swbam_app::SimulatedWriteBuffer Buffer;
    Buffer buffer;
    volatile unsigned long long guard = 0;

    Check(buffer.Allocate(0, 0) == 0 && !buffer.data(), "empty allocation");
    PatternSource empty(0);
    Check(buffer.CopyBodyWindowed(empty, 0, &guard) == 0 && !empty.calls,
          "empty body does not read");
    Check(buffer.Allocate(0, 28) == 0 && buffer.size() == 28,
          "envelope-only allocation");

    allocation_calls = 0;
    allocation_limit = 1024;
    Check(buffer.Allocate(128, 17, LimitedMalloc) == 0 && !buffer.bounded() &&
          allocation_calls == 1 && buffer.body_capacity() == 128 &&
          buffer.size() == 145, "full-copy fast allocation");
    memset(buffer.data(), 0x5a, buffer.size());
    PatternSource small(128);
    Check(buffer.CopyBodyWindowed(small, 8, &guard) == 0 && small.calls == 1 &&
          small.copied == 128, "full-copy body");
    Check(static_cast<unsigned char>(buffer.data()[7]) == 0x5a &&
          buffer.data()[8] == 1 && buffer.data()[135] == 1 &&
          static_cast<unsigned char>(buffer.data()[136]) == 0x5a,
          "prefix and suffix preserved");

    const size_t window = Buffer::WindowBytes();
    allocation_calls = 0;
    allocation_limit = std::numeric_limits<size_t>::max();
    Check(buffer.Allocate(window + 1, 17, LimitedMalloc) == 0 &&
          !buffer.bounded() && allocation_calls == 1 &&
          buffer.body_capacity() == window + 1,
          "large successful allocation keeps full-copy path");
    allocation_calls = 0;
    allocation_limit = window + 17;
    const size_t large_size = 2 * window + 37;
    Check(buffer.Allocate(large_size, 17, LimitedMalloc) == 0 &&
          buffer.bounded() && allocation_calls == 2 &&
          buffer.body_capacity() == window && buffer.size() == window + 17,
          "failed full allocation falls back to bounded window");
    memset(buffer.data(), 0x5a, buffer.size());
    PatternSource large(large_size);
    guard = 0;
    Check(buffer.CopyBodyWindowed(large, 8, &guard) == 0 && large.calls == 3 &&
          large.copied == large_size && guard == 12,
          "every byte including short tail is copied and observed");
    Check(buffer.data()[8] == 3 && buffer.data()[44] == 3 &&
          buffer.data()[45] == 2 &&
          static_cast<unsigned char>(buffer.data()[7]) == 0x5a &&
          static_cast<unsigned char>(buffer.data()[8 + window]) == 0x5a,
          "tail and envelope bounds");
    PatternSource bad_prefix(1);
    Check(buffer.CopyBodyWindowed(bad_prefix, buffer.size(), &guard) != 0 &&
          !bad_prefix.calls, "invalid window rejected before read");
    PatternSource failed(large_size, 2);
    Check(buffer.CopyBodyWindowed(failed, 8, &guard) != 0 &&
          failed.calls == 2 && failed.copied == window,
          "read failure stops without skipping bytes");

    allocation_calls = 0;
    allocation_limit = 0;
    Check(buffer.Allocate(large_size, 17, LimitedMalloc) != 0 && !buffer.data() &&
          allocation_calls == 2, "both allocation failures propagate");
    allocation_calls = 0;
    Check(buffer.Allocate(std::numeric_limits<size_t>::max(), 1,
                          LimitedMalloc) != 0 && errno == EOVERFLOW &&
          !allocation_calls, "size overflow rejected");
    buffer.Reset();
    Check(!buffer.data() && !buffer.size() && !buffer.body_capacity() &&
          !buffer.bounded(), "reset releases buffer and clears state");
    puts("PASS simulated write buffer: full, bounded, empty, tail, failure, overflow");
    return 0;
}
