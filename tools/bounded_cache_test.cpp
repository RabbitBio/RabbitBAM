#include "bounded-cache.h"
#include <cassert>
#include <cstdio>
#include <vector>

int main() {
    swbam_usage::CachedTemporary cache(64u << 10);
    std::vector<char> expected(200000), actual(expected.size());
    for (size_t i = 0; i < expected.size(); ++i) expected[i] = char(i * 17);
    assert(cache.Append(expected.data(), 32000) == 0);
    assert(cache.capacity() <= 65536);
    assert(cache.Append(expected.data() + 32000, expected.size() - 32000) == 0);
    assert(cache.capacity() == 0);
    expected[317] = 42;
    assert(cache.WriteAt(317, expected.data() + 317, 1) == 0);
    assert(cache.ReadAt(0, actual.data(), actual.size()) == 0);
    assert(actual == expected);
    assert(cache.ReadAt(cache.size(), nullptr, 0) == 0);
    assert(cache.ReadAt(cache.size(), actual.data(), 1) != 0);
    assert(cache.WriteAt(cache.size() + 1, actual.data(), 1) != 0);
    puts("Bounded temporary cache tests PASS");
}
