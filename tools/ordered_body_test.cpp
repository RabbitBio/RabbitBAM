#include "ordered-body.h"
#include <cassert>
#include <dirent.h>
#include <map>
#include <string>

#ifdef SWBAM_ENABLE_RAM_STORAGE
// Exercise the real RAM file implementation without launching an MPI session.
namespace swbam_usage {
static RamPool test_pool;
static std::map<int, std::unique_ptr<RamFile> > test_files;
static int next_handle = 1;
bool IsRam() { return true; }
int OpenTemporary(const std::string &) {
    const int handle = next_handle++;
    test_files[handle].reset(new RamFile(&test_pool));
    return handle;
}
int CloseTemporary(int handle) { return test_files.erase(handle) == 1 ? 0 : -1; }
ssize_t ReadTemporary(int handle, void *data, size_t bytes, off_t offset) {
    auto it = test_files.find(handle);
    return it != test_files.end() && offset >= 0 &&
        it->second->ReadAt(offset, data, bytes) == 0 ? (ssize_t)bytes : -1;
}
ssize_t WriteTemporary(int handle, const void *data, size_t bytes, off_t offset) {
    auto it = test_files.find(handle);
    return it != test_files.end() && offset >= 0 &&
        it->second->WriteAt(offset, data, bytes) == 0 ? (ssize_t)bytes : -1;
}
}
#endif

class Prefix : public swbam::RankBodySink {
public:
    std::string data;
    mutable uint64_t released = 0;
    int flushes = 0;
    int Append(const void *src, size_t bytes) {
        if (bytes) data.append(static_cast<const char *>(src), bytes);
        return 0;
    }
    int Flush() { ++flushes; return 0; }
    uint64_t size() const { return data.size(); }
    bool is_memory() const { return true; }
    int ReadAt(uint64_t offset, void *dst, size_t bytes) const {
        if (offset < released || offset > data.size() || bytes > data.size() - offset) return -1;
        if (bytes) memcpy(dst, data.data() + offset, bytes);
        return 0;
    }
    int ReleaseBefore(uint64_t offset) const {
        if (offset > data.size()) return -1;
        released = std::max(released, offset);
        return 0;
    }
};

static size_t OpenFiles() {
#ifdef SWBAM_ENABLE_RAM_STORAGE
    return swbam_usage::test_files.size();
#else
    DIR *dir = opendir("/proc/self/fd");
    assert(dir);
    size_t count = 0;
    while (dirent *entry = readdir(dir)) {
        if (strcmp(entry->d_name, ".") && strcmp(entry->d_name, "..")) ++count;
    }
    assert(closedir(dir) == 0);
    return count;
#endif
}

static void TestPrefix() {
    Prefix prefix;
    swbam_usage::OrderedBody body(&prefix, 0);
    assert(body.AppendOrdered(4, "abc", 3) == 0);
    assert(body.AppendOrdered(5, "de", 2) == 0);
    assert(body.size() == 5 && prefix.data == "abcde");
    assert(body.memory_bytes() == 0 && body.is_memory());
    assert(body.Flush() == 0 && body.Flush() == 0 && prefix.flushes == 1);
    char actual[5];
    assert(body.ReadAt(0, actual, sizeof(actual)) == 0);
    assert(memcmp(actual, "abcde", 5) == 0);
    assert(body.ReleaseBefore(2) == 0 && prefix.released == 2);
    assert(body.ReadAt(2, actual, 3) == 0 && memcmp(actual, "cde", 3) == 0);
    assert(body.ReadAt(body.size(), nullptr, 0) == 0);
    assert(body.ReadAt(0, actual, 1) != 0);
    assert(body.AppendOrdered(6, "x", 1) != 0);
    assert(body.Activate() != 0);
    assert(body.ReleaseBefore(5) == 0 && body.size() == 5);
}

static void TestOrder() {
    const size_t before = OpenFiles();
    Prefix prefix;
    swbam_usage::OrderedBody body(&prefix, 1u << 20);
    assert(body.AppendOrdered(0, "prefix", 6) == 0);
    assert(body.Activate() == 0 && body.Activate() == 0);
    assert(body.AppendOrdered(7, "D", 1) == 0);
    assert(body.AppendOrdered(4, "BB", 2) == 0);
    assert(body.AppendOrdered(6, "CCC", 3) == 0);
    assert(body.AppendOrdered(2, "AAAA", 4) == 0);
    assert(prefix.data == "prefix" && body.size() == 16);
    assert(body.memory_bytes() > 0 && body.memory_bytes() <= (1u << 20));
    assert(body.Flush() == 0 && OpenFiles() == before + 1);
    const std::string expected = "prefixAAAABBCCCD";
    for (size_t offset = 0; offset <= expected.size(); ++offset) {
        std::vector<char> actual(expected.size() - offset);
        assert(body.ReadAt(offset, actual.data(), actual.size()) == 0);
        assert(std::string(actual.begin(), actual.end()) == expected.substr(offset));
    }
    const size_t metadata = body.memory_bytes();
    assert(body.ReleaseBefore(8) == 0 && prefix.released == 6);
    assert(OpenFiles() == before + 1);
    char tail[8];
    assert(body.ReadAt(8, tail, sizeof(tail)) == 0);
    assert(memcmp(tail, expected.data() + 8, sizeof(tail)) == 0);
    assert(body.ReleaseBefore(body.size()) == 0 && OpenFiles() == before);
    assert(body.memory_bytes() < metadata && body.size() == expected.size());
    assert(body.ReleaseBefore(8) == 0);
}

static void TestSegmentRelease() {
    const size_t before = OpenFiles();
    const size_t half = 32u << 20;
    Prefix prefix;
    swbam_usage::OrderedBody body(&prefix, 1u << 20);
    assert(body.Activate() == 0);
    std::vector<char> payload(half, 'A');
    assert(body.AppendOrdered(12, payload.data(), payload.size()) == 0);
    std::fill(payload.begin(), payload.end(), 'B');
    assert(body.AppendOrdered(4, payload.data(), payload.size()) == 0);
    std::fill(payload.begin(), payload.end(), 'C');
    assert(body.AppendOrdered(8, payload.data(), payload.size()) == 0);
    assert(body.Flush() == 0 && OpenFiles() == before + 2);
    // B and A share segment 0, but C occupies segment 1. Releasing logical B
    // must retain A, and finishing C must close segment 1 before segment 0.
    std::vector<char> chunk(1u << 20);
    const char expected[] = {'B', 'C', 'A'};
    for (size_t part = 0; part < 3; ++part) {
        for (size_t pos = 0; pos < half; pos += chunk.size()) {
            const uint64_t offset = part * half + pos;
            assert(body.ReadAt(offset, chunk.data(), chunk.size()) == 0);
            for (char c : chunk) assert(c == expected[part]);
            assert(body.ReleaseBefore(offset + chunk.size()) == 0);
        }
        assert(OpenFiles() == before + (part == 0 ? 2 : part == 1 ? 1 : 0));
    }
    assert(body.size() == 3ull * half);
#ifdef SWBAM_ENABLE_RAM_STORAGE
    assert(swbam_usage::test_pool.used() == 0);
#endif
}

static void TestSplitAndCopy() {
    const size_t before = OpenFiles();
    Prefix prefix;
    swbam_usage::OrderedBody body(&prefix, 1u << 20);
    assert(body.Activate() == 0);
    std::vector<char> payload((64u << 20) + 17);
    for (size_t i = 0; i < payload.size(); ++i) payload[i] = char(i * 13);
    assert(body.AppendOrdered(4, payload.data(), payload.size()) == 0);
    assert(body.AppendOrdered(2, "before", 6) == 0);
    assert(body.AppendOrdered(5, "after", 5) == 0);
    assert(body.Flush() == 0 && OpenFiles() == before + 2);
    swbam_usage::CachedTemporary destination(0);
    std::vector<char> chunk(1u << 20);
    for (uint64_t offset = 0; offset < body.size();) {
        const size_t count = (size_t)std::min<uint64_t>(chunk.size(), body.size() - offset);
        assert(body.ReadAt(offset, chunk.data(), count) == 0);
        assert(body.ReleaseBefore(offset + count) == 0);
        assert(destination.Append(chunk.data(), count) == 0);
        offset += count;
    }
    assert(OpenFiles() == before + 1 && destination.size() == payload.size() + 11);
    char envelope[6];
    assert(destination.ReadAt(0, envelope, 6) == 0 && memcmp(envelope, "before", 6) == 0);
    for (size_t pos = 0; pos < payload.size();) {
        const size_t count = std::min(chunk.size(), payload.size() - pos);
        assert(destination.ReadAt(pos + 6, chunk.data(), count) == 0);
        assert(memcmp(chunk.data(), payload.data() + pos, count) == 0);
        pos += count;
    }
    assert(destination.ReadAt(payload.size() + 6, envelope, 5) == 0);
    assert(memcmp(envelope, "after", 5) == 0);
}

static void TestErrorsAndLifetime() {
    Prefix prefix;
    {
        swbam_usage::OrderedBody body(&prefix, 2048);
        assert(body.Append("unkeyed", 7) != 0);
        assert(body.Activate() == 0);
        size_t i = 0;
        while (i < 1000 && body.AppendOrdered(i * 2, "x", 1) == 0) {
            assert(body.memory_bytes() <= 2048);
            ++i;
        }
        assert(i > 0 && i < 1000 && errno == ENOBUFS);
        assert(body.memory_bytes() <= 2048 && body.Flush() != 0);
    }
    {
        swbam_usage::OrderedBody body(&prefix, 1u << 20);
        assert(body.Activate() == 0);
        assert(body.AppendOrdered(2, "a", 1) == 0);
        assert(body.AppendOrdered(2, "b", 1) == 0);
        assert(body.Flush() != 0 && errno == EINVAL);
    }
    std::unique_ptr<swbam_usage::OrderedBody> body;
    {
        Prefix short_lived;
        body.reset(new swbam_usage::OrderedBody(&short_lived, 1u << 20));
        assert(body->AppendOrdered(0, "prefix", 6) == 0);
        assert(body->Activate() == 0);
        assert(body->AppendOrdered(2, "suffix", 6) == 0);
    }
    body.reset();
}

int main() {
#ifdef SWBAM_ENABLE_RAM_STORAGE
    assert(swbam_usage::test_pool.Initialize(160u << 20) == 0);
#endif
    const size_t before = OpenFiles();
    TestPrefix();
    TestOrder();
    TestSegmentRelease();
    TestSplitAndCopy();
    TestErrorsAndLifetime();
    assert(OpenFiles() == before);
#ifdef SWBAM_ENABLE_RAM_STORAGE
    assert(swbam_usage::test_pool.used() == 0);
#endif
    puts("Ordered body tests PASS");
}
