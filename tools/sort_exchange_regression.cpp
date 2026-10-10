// Exercise production byte and external-run exchange, including CPE packing.
#include "../swbam-app/algorithms/sort/swbam_sort_optimized_mpi.cpp"

static void Check(bool ok, const char *message) {
    if (ok) return;
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    fprintf(stderr, "[rank %d] sort exchange regression FAILED: %s\n",
            rank, message);
    MPI_Abort(MPI_COMM_WORLD, 2);
    abort();
}

static void TestRunArena(int rank) {
    const size_t meta_bytes = sizeof(MpiSortRecordMeta);
    MpiSortStrictRunArena arena;
    std::vector<unsigned char> data(8 * meta_bytes + 512, 0x5a);
    arena.data = data.data();
    arena.capacity = data.size();
    arena.reset();
    Check(arena.can_append(arena.meta_begin - meta_bytes, 1),
          "arena exact fit");
    Check(!arena.can_append(arena.meta_begin - meta_bytes + 1, 1),
          "arena payload overflow");
    arena.raw_used = arena.meta_begin - meta_bytes + 1;
    Check(!arena.can_append(1, 1), "arena metadata overlap");
    const std::vector<unsigned char> before = data;
    MpiSortRecordMeta record = {};
    record.raw_len = 16;
    const unsigned char raw[16] = {};
    Check(MpiSortStrictAppendBlock(&arena, raw, sizeof(raw), 0,
                                  &record, 1) == 1,
          "arena spill before metadata overlap");
    Check(data == before && arena.record_count == 0,
          "arena rejected append preserves data");
    if (!rank) {
        printf("PASS sort external arena capacity\n");
        fflush(stdout);
    }
}

static size_t Bytes(int test, int src, int dst, int ranks) {
    const size_t chunk = kSortExchangeChunk;
    if (test == 0) return 0;
    if (test == 1) {
        if (src == 0 && dst == 1) return 1;
        if (src == 1 && dst == 2) return 7;
        return 0;
    }
    if (test == 2) {
        if (dst != (src + 1) % ranks) return 0;
        if (src == 0) return 2 * chunk + 17;
        if (src == 1) return 1;
        if (src == 2) return chunk + 3;
        return 13 + src;
    }
    if (src == 0 && dst == ranks - 1) return chunk + 1;
    return (src * 7 + dst * 3) % 5;
}

static unsigned char Pattern(int src, int dst, size_t index, int pass) {
    return (unsigned char)(src * 17 + dst * 11 + index * 13 + pass);
}

static void TestExchange(int test, int rank, int ranks) {
    for (int step = 1; step < ranks; ++step) {
        const int send_to = (rank + step) % ranks;
        const int recv_from = (rank - step + ranks) % ranks;
        const size_t send_size = Bytes(test, rank, send_to, ranks);
        const size_t recv_size = Bytes(test, recv_from, rank, ranks);
        std::vector<unsigned char> send(send_size);
        std::vector<unsigned char> recv(recv_size);
        // Reuse the production tags across metadata/raw calls and rounds.
        for (int pass = 0; pass < 2; ++pass) {
            for (size_t i = 0; i < send_size; ++i)
                send[i] = Pattern(rank, send_to, i, pass);
            std::fill(recv.begin(), recv.end(), 0xff);
            Check(MpiSortSendRecvBytes(send_to,
                      send_size ? (const char *)send.data() : nullptr, send_size,
                      recv_from, recv_size ? (char *)recv.data() : nullptr, recv_size,
                      pass ? kSortTagRaw : kSortTagMeta) == 0, "byte exchange");
            for (size_t i = 0; i < recv_size; ++i)
                Check(recv[i] == Pattern(recv_from, rank, i, pass), "payload bytes");
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        printf("PASS sort byte exchange ranks=%d case=%d\n", ranks, test);
        fflush(stdout);
    }
}

static size_t RunCount(int test, int src) {
    if (test != 3) return 1;
    return src == 0 ? 3 : src == 1 ? 1 : src % 2 ? 2 : 0;
}

static size_t RecordCount(int test, int src, int dst, int ranks, int round) {
    if ((size_t)round >= RunCount(test, src) || test == 0) return 0;
    if (test == 1) {
        if (src == 0 && dst == 1) return 5;
        if (src == 1 && dst == 2) return 1;
        if (src == ranks - 1 && dst == 0) return 3;
        return src == 0 && dst == 0 ? 2 : 0;
    }
    if (test == 2) return (src * 7 + dst * 3) % 5;
    if (src == 0 && dst == ranks - 1) return 5 + round;
    return dst == (src + 2) % ranks ? 2 : 0;
}

static void TestExternalExchange(int test, int rank, int ranks,
                                 const char *temp_dir) {
    const size_t raw_len = 64;
    const uint64_t order_stride = 100000;
    MpiSortMemoryTracker tracker(2 * 1024 * 1024);
    MpiSortStrictBuffer arena_memory;
    Check(MpiSortStrictAlloc(&tracker, 64 * 1024, &arena_memory) == 0,
          "external arena allocation");
    MpiSortStrictRunArena arena;
    arena.data = arena_memory.data;
    arena.capacity = arena_memory.size;
    arena.reset();
    MpiSortStrictExchangeWorkspace ws;
    Check(MpiSortStrictAllocateExchange(&tracker, 1024 * 1024, &ws) == 0,
          "external exchange allocation");
    // Force unequal chunk counts while leaving the underlying buffers aligned.
    ws.meta_capacity = ws.prefix_capacity = 4;
    ws.raw_capacity = 2 * raw_len;
    unsigned char scratch[4096];
    MpiSortStats stats = {};
    MpiSortStrictTempStore store;
    Check(MpiSortStrictOpenTemp(temp_dir, rank, &store) == 0, "external temp files");
    std::vector<MpiSortStrictRun> runs;
    for (size_t round = 0; round < RunCount(test, rank); ++round) {
        std::vector<MpiSortRecordMeta> records;
        std::vector<unsigned char> raw;
        for (int dst = 0; dst < ranks; ++dst) {
            for (size_t i = 0; i < RecordCount(test, rank, dst, ranks, round); ++i) {
                MpiSortRecordMeta meta = {};
                meta.tid = dst;
                meta.pos = (int32_t)(round * 100 + i);
                meta.raw_offset = raw.size();
                meta.raw_len = raw_len;
                meta.global_order = (rank * 3 + round) * order_stride + i;
                raw.resize(raw.size() + raw_len);
                for (size_t j = 0; j < raw_len; ++j)
                    raw[meta.raw_offset + j] = Pattern(rank, dst, j, meta.pos);
                uint32_t size = raw_len - 4;
                for (int j = 0; j < 4; ++j)
                    raw[meta.raw_offset + j] = (unsigned char)(size >> (8 * j));
                records.push_back(meta);
            }
        }
        MpiSortStrictRun run;
        Check(MpiSortStrictAppendExtent(store.runs_fd, &store.runs_end,
                  records.data(), records.size(), raw.data(), raw.size(),
                  scratch, sizeof(scratch), &stats, &run.extent) == 0,
              "store external run");
        runs.push_back(run);
    }
    std::vector<MpiSortKey> splitters;
    for (int dst = 0; dst < ranks - 1; ++dst) {
        MpiSortKey key = {};
        key.tid = dst;
        key.pos = INT_MAX;
        key.flag = UINT16_MAX;
        key.global_order = UINT64_MAX;
        splitters.push_back(key);
    }
    std::vector<MpiSortStrictSegment> segments;
    Check(MpiSortStrictExchangeRuns(runs, splitters, &arena, &ws, &store,
              scratch, sizeof(scratch), rank, ranks, &segments, &stats) == 0,
          "external run exchange");
    std::vector<size_t> seen((size_t)ranks * 3, 0);
    for (const auto &segment : segments) {
        for (const auto &extent : segment.extents) {
            std::vector<MpiSortRecordMeta> records(extent.record_count);
            std::vector<unsigned char> raw(extent.raw_size);
            Check(MpiSortStrictPreadAll(extent.fd, records.data(),
                      records.size() * sizeof(MpiSortRecordMeta),
                      extent.meta_offset, &stats) == 0, "read external metadata");
            Check(MpiSortStrictPreadAll(extent.fd, raw.data(), raw.size(),
                      extent.raw_offset, &stats) == 0, "read external payload");
            size_t pos = 0;
            for (const auto &meta : records) {
                const size_t source_round = meta.global_order / order_stride;
                Check(source_round < seen.size() && meta.tid == rank,
                      "external owner/source");
                Check(meta.global_order % order_stride == seen[source_round]++,
                      "external peer record order");
                Check(meta.raw_offset == pos && meta.raw_len == raw_len &&
                      pos + raw_len <= raw.size(), "external raw offsets");
                Check(MpiSortReadLe32(raw.data() + pos) == raw_len - 4,
                      "external record size");
                for (size_t j = 4; j < raw_len; ++j)
                    Check(raw[pos + j] == Pattern(source_round / 3, rank, j, meta.pos),
                          "external payload bytes");
                pos += raw_len;
            }
            Check(pos == raw.size(), "external raw total");
        }
    }
    for (int src = 0; src < ranks; ++src)
        for (int round = 0; round < 3; ++round)
            Check(seen[(size_t)src * 3 + round] ==
                      RecordCount(test, src, rank, ranks, round),
                  "external peer totals");
    MpiSortStrictCloseTemp(&store);
    MpiSortStrictFreeExchange(&tracker, &ws);
    MpiSortStrictFree(&tracker, &arena_memory);
    Check(tracker.current == 0, "external allocation cleanup");
    MPI_Barrier(MPI_COMM_WORLD);
    if (!rank) {
        printf("PASS sort external exchange ranks=%d case=%d\n", ranks, test);
        fflush(stdout);
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    Check(ranks >= 3 && argc <= 2, "usage: >=3 ranks, optional existing temp directory");
    TestRunArena(rank);
    for (int test = 0; test < 4; ++test) TestExchange(test, rank, ranks);
    if (argc == 2) {
        athread_init();
        for (int test = 0; test < 4; ++test)
            TestExternalExchange(test, rank, ranks, argv[1]);
        athread_halt();
    }
    MPI_Finalize();
    return 0;
}
