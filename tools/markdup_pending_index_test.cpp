#include "../swbam-app/algorithms/markdup/markdup-pending-index.h"

#include <cassert>
#include <cstdio>
#include <list>
#include <vector>

struct Batch {
    uint64_t ordinal_base;
    uint64_t ordinal_end;
};

int main() {
    swbam_markdup::PendingIndex<Batch> index;
    assert(!index.Find(0) && !index.Find(UINT64_MAX));
    std::list<Batch> batches;
    for (unsigned i = 0; i < 400; ++i) {
        batches.push_back({i * 100u, i * 100u + 73});
        assert(index.Insert(&batches.back()));
    }
    Batch empty = {100, 100};
    Batch invalid = {10, 9};
    Batch overlap = {70, 120};
    assert(index.Insert(&empty));
    index.Erase(&empty);
    assert(!index.Insert(&invalid) && !index.Insert(&overlap));
    assert(!index.Insert(&batches.front()));
    assert(index.size() == 400);
    uint32_t random = 1;
    for (unsigned i = 0; i < 100000; ++i) {
        random = random * 1664525u + 1013904223u;
        uint64_t ordinal = (random >> 8) % 41000;
        Batch *expected = nullptr;
        for (auto &batch : batches) {
            if (ordinal >= batch.ordinal_base && ordinal < batch.ordinal_end) {
                expected = &batch;
                break;
            }
        }
        assert(index.Find(ordinal) == expected);
    }
    for (auto it = batches.begin(); it != batches.end();) {
        // Populate the lookup cache immediately before erasing its owner.
        assert(index.Find(it->ordinal_base) == &*it);
        const uint64_t ordinal = it->ordinal_base;
        index.Erase(&*it);
        it = batches.erase(it);
        assert(!index.Find(ordinal));
        if (it != batches.end()) ++it;
    }
    assert(index.size() == 200);
    assert(index.memory_bytes() == 200 * index.EntryBytes());
    index.Clear();
    batches.clear();
    assert(index.size() == 0 && !index.Find(0));
    Batch reused = {0, 10};
    assert(index.Insert(&reused) && index.Find(9) == &reused);
    index.Erase(&reused);
    assert(!index.Find(9));
    puts("markdup pending index range/lifetime/budget PASS");
}
