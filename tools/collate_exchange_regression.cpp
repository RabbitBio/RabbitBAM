// Test the production internal exchange functions. Link in place of main.cpp
// and collate.cpp; this is not an application or SDK translation unit.
#include "../swbam-app/algorithms/collate/swbam_collate_mpi.cpp"

static void Check(bool ok, const char *message) {
    if (ok) return;
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    fprintf(stderr, "[rank %d] collate exchange regression FAILED: %s\n",
            rank, message);
    MPI_Abort(MPI_COMM_WORLD, 2);
    abort();
}

static size_t Count(int test, int src, int dst, int ranks, size_t limit) {
    if (test == 0) return 0;
    if (test == 1) {
        if (src == 0 && dst == 1) return 2 * limit + 1;
        if (src == 1 && dst == 2) return 1;
        if (src == ranks - 1 && dst == 0) return 3;
        return src == 0 && dst == 0 ? 2 : 0;
    }
    if (test == 2) return (src * 7 + dst * 3) % 5 * 3;
    if (src == 0 && dst == ranks - 1) return limit + 3;
    return dst == (src + 2) % ranks ? 2 : 0;
}

static unsigned char Pattern(int src, size_t index) {
    return (unsigned char)(src * 17 + index * 13);
}

static void TestExchange(int test, bool external, int rank, int ranks,
                         const char *temp_dir) {
    const size_t raw_len = test == 3 ? 1 : external ? 1024 : 4096;
    const size_t raw_limit = external ? 4096 : kCollateExchangeRaw;
    const size_t meta_limit = external ? 8 : kCollateExchangeMeta / sizeof(CollateMeta);
    const size_t limit = test == 3 ? meta_limit : raw_limit / raw_len;
    std::vector<CollateMeta> records;
    CollateRawArena raw;
    size_t count = 0;
    for (int dst = 0; dst < ranks; ++dst)
        count += Count(test, rank, dst, ranks, limit);
    raw.size = raw.capacity = count * raw_len;
    if (count) {
        raw.data = aligned_alloc_custom(64, raw.capacity);
        Check(raw.data != nullptr, "source allocation");
    }
    for (int dst = 0; dst < ranks; ++dst) {
        for (size_t i = 0; i < Count(test, rank, dst, ranks, limit); ++i) {
            CollateMeta meta = {};
            meta.hash = rank;
            meta.bin = dst;
            meta.raw_offset = records.size() * raw_len;
            meta.raw_len = raw_len;
            meta.global_order = i;
            memset(raw.data + meta.raw_offset, Pattern(rank, i), raw_len);
            records.push_back(meta);
        }
    }

    MpiCollateStats stats = {};
    std::vector<CollateMemorySegment> memory_segments;
    std::vector<CollateSegment> file_segments;
    CollateTempStore store;
    if (!external) {
        Check(CollateExchangeMemory(records, raw.data, raw.size,
              rank, ranks, ranks, &memory_segments, &stats) == 0,
              "memory exchange");
    } else {
        Check(CollateOpenTemp(temp_dir, rank, &store) == 0, "temp files");
        CollateRunArena arena;
        arena.capacity = raw.size + records.size() * sizeof(CollateMeta);
        arena.data = aligned_alloc_custom(64, std::max<size_t>(64, arena.capacity));
        Check(arena.data != nullptr, "run allocation");
        arena.raw_used = arena.meta_begin = raw.size;
        arena.record_count = records.size();
        std::vector<CollateRun> runs;
        if (count) {
            memcpy(arena.data, raw.data, raw.size);
            memcpy(arena.records(), records.data(), records.size() * sizeof(CollateMeta));
            runs.resize(1);
            runs[0].resident = true;
        }
        unsigned char scratch[4096];
        Check(CollateExchangeRuns(runs, &arena, &store, rank, ranks, ranks,
              raw_limit, meta_limit * sizeof(CollateMeta), scratch, sizeof(scratch),
              &file_segments, &stats) == 0, "external exchange");
        aligned_free_custom(arena.data);
    }

    // The merge input must remain valid after freeing the extraction buffers.
    std::vector<CollateMeta>().swap(records);
    raw.Reset();
    Check(!raw.data && !raw.size && !raw.capacity, "arena reset");
    std::vector<size_t> seen(ranks, 0);
    auto verify = [&](const CollateMeta *meta, size_t n,
                      const unsigned char *data, size_t bytes) {
        size_t offset = 0;
        for (size_t i = 0; i < n; ++i) {
            const CollateMeta &m = meta[i];
            Check(m.hash < (unsigned)ranks && m.bin == (unsigned)rank, "owner/source");
            Check(m.raw_len == raw_len && m.raw_offset == offset &&
                  offset + raw_len <= bytes, "raw offsets");
            Check(m.global_order == seen[m.hash], "peer record order");
            const unsigned char expected = Pattern(m.hash, seen[m.hash]++);
            for (size_t j = 0; j < raw_len; ++j)
                Check(data[offset + j] == expected, "payload bytes");
            offset += raw_len;
        }
        Check(offset == bytes, "raw total");
    };
    for (const auto &segment : memory_segments)
        verify(segment.records.data(), segment.records.size(),
               segment.raw_data, segment.raw_size);
    for (const auto &segment : file_segments) {
        for (const auto &extent : segment.extents) {
            std::vector<CollateMeta> meta(extent.record_count);
            std::vector<unsigned char> data(extent.raw_size);
            Check(CollatePreadAll(extent.fd, meta.data(), meta.size() * sizeof(CollateMeta),
                  extent.meta_offset, &stats) == 0, "read metadata");
            Check(CollatePreadAll(extent.fd, data.data(), data.size(),
                  extent.raw_offset, &stats) == 0, "read payload");
            verify(meta.data(), meta.size(), data.data(), data.size());
        }
    }
    if (external) CollateCloseTemp(&store);
    for (int src = 0; src < ranks; ++src)
        Check(seen[src] == Count(test, src, rank, ranks, limit), "peer totals");
    MPI_Barrier(MPI_COMM_WORLD);
    if (!rank) {
        printf("PASS exchange mode=%s ranks=%d case=%d\n",
               external ? "external" : "memory", ranks, test);
        fflush(stdout);
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    Check(ranks >= 3 && argc == 2, "usage: >=3 ranks, existing temp directory");
    for (int mode = 0; mode < 2; ++mode)
        for (int test = 0; test < 4; ++test)
            TestExchange(test, mode != 0, rank, ranks, argv[1]);
    MPI_Finalize();
    return 0;
}
