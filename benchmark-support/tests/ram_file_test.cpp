#include "ram_file.h"
#include <cassert>
#include <cerrno>
#include <cstring>
#include <cstdio>
#include <memory>
#include <vector>

int main() {
    using namespace swbam_usage;
    RamPool pool;
    assert(pool.Initialize(32u * 1024u * 1024u) == 0);
    {
        RamFile file(&pool);
        const char value[] = "retained";
        const uint64_t offset = 5u * 1024u * 1024u + 17;
        assert(file.WriteAt(offset, value, sizeof(value)) == 0);
        char actual[sizeof(value)] = {};
        assert(file.ReadAt(offset, actual, sizeof(actual)) == 0);
        assert(memcmp(actual, value, sizeof(value)) == 0);
        std::vector<char> hole(offset);
        assert(file.ReadAt(0, hole.data(), hole.size()) == 0);
        for (char c : hole) assert(c == 0);
        assert(file.WriteAt(offset + 2, "X", 1) == 0);
        assert(file.ReadAt(offset + 2, actual, 1) == 0 && actual[0] == 'X');
        const uint64_t before = file.size();
        const size_t used = pool.used();
        assert(file.WriteAt(1ULL << 32, value, sizeof(value)) == -1 && errno == ENOSPC);
        assert(file.size() == before && pool.used() == used);
        assert(file.ReadAt(before, actual, 1) == -1);
        assert(file.ReadAt(before, nullptr, 0) == 0);
        assert(file.Resize(0) == 0 && pool.used() == 0);
        assert(file.ReserveContiguous(8u * 1024u * 1024u) == 0);
        std::vector<char> data(8u * 1024u * 1024u, 42);
        assert(file.WriteAt(0, data.data(), data.size()) == 0);
        assert(file.Range(0, data.size()) && file.Range(0, data.size())[data.size() - 1] == 42);
        assert(file.Resize(17) == 0);
        assert(file.Resize(4096) == 0);
        assert(file.Range(17, 1)[0] == 0);
    }
    assert(pool.used() == 0);
    {
        std::vector<std::unique_ptr<RamFile> > occupied;
        for (int i = 0; i < 8; ++i) {
            occupied.emplace_back(new RamFile(&pool));
            assert(occupied.back()->WriteAt(0, "x", 1) == 0);
        }
        occupied[2].reset();
        occupied[4].reset();
        RamFile fragmented(&pool);
        std::vector<char> data(2 * RamPool::kPage, 37);
        assert(fragmented.WriteAt(0, data.data(), RamPool::kPage) == 0);
        char *first = fragmented.Range(0, 1);
        assert(fragmented.WriteAt(RamPool::kPage, data.data(), RamPool::kPage) == 0);
        assert(fragmented.Range(0, 1) == first);
        assert(fragmented.Range(0, data.size()) == nullptr);
        std::vector<char> actual(data.size());
        assert(fragmented.ReadAt(0, actual.data(), actual.size()) == 0);
        assert(actual == data);
    }
    assert(pool.used() == 0);
    puts("RAM file tests PASS");
}
