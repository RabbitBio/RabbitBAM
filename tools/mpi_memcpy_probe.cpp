#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#if defined(__GNUC__)
#define RB_NO_INLINE __attribute__((noinline))
#else
#define RB_NO_INLINE
#endif

static volatile uint64_t g_probe_sink = 0;

static inline void compiler_barrier(const void *ptr) {
#if defined(__GNUC__)
    __asm__ __volatile__("" : : "r"(ptr) : "memory");
#else
    (void)ptr;
#endif
}

static RB_NO_INLINE void copy_by_memcpy(unsigned char *dst, const unsigned char *src, size_t n) {
    std::memcpy(dst, src, n);
    compiler_barrier(dst);
}

static RB_NO_INLINE void copy_by_for_u8(unsigned char *dst, const unsigned char *src, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        dst[i] = src[i];
    }
    compiler_barrier(dst);
}

static RB_NO_INLINE void copy_by_for_u64(unsigned char *dst, const unsigned char *src, size_t n) {
    size_t qwords = n / sizeof(uint64_t);
    uint64_t *dst64 = reinterpret_cast<uint64_t *>(dst);
    const uint64_t *src64 = reinterpret_cast<const uint64_t *>(src);

    for (size_t i = 0; i < qwords; ++i) {
        dst64[i] = src64[i];
    }
    for (size_t i = qwords * sizeof(uint64_t); i < n; ++i) {
        dst[i] = src[i];
    }
    compiler_barrier(dst);
}

static uint64_t sample_checksum(const unsigned char *buf, size_t n) {
    if (n == 0) return 0;

    uint64_t h = 1469598103934665603ULL;
    size_t stride = n / 4096;
    if (stride == 0) stride = 1;

    for (size_t i = 0; i < n; i += stride) {
        h ^= static_cast<uint64_t>(buf[i]);
        h *= 1099511628211ULL;
    }
    h ^= static_cast<uint64_t>(buf[n - 1]);
    return h;
}

static void fill_source(unsigned char *buf, size_t n, int rank) {
    uint32_t x = 0x9e3779b9u ^ static_cast<uint32_t>(rank * 2654435761u);
    for (size_t i = 0; i < n; ++i) {
        x = x * 1664525u + 1013904223u;
        buf[i] = static_cast<unsigned char>((x >> 24) ^ (i & 0xffu));
    }
}

typedef void (*CopyFn)(unsigned char *, const unsigned char *, size_t);

struct ProbeResult {
    double seconds;
    double gbps;
    uint64_t checksum;
};

static ProbeResult run_probe(const char *name,
                             CopyFn fn,
                             unsigned char *dst,
                             const unsigned char *src,
                             size_t bytes,
                             int iterations,
                             int warmup,
                             MPI_Comm comm) {
    for (int i = 0; i < warmup; ++i) {
        fn(dst, src, bytes);
    }
    MPI_Barrier(comm);

    double t0 = MPI_Wtime();
    for (int i = 0; i < iterations; ++i) {
        fn(dst, src, bytes);
    }
    double t1 = MPI_Wtime();

    ProbeResult result;
    result.seconds = t1 - t0;
    result.gbps = (static_cast<double>(bytes) * static_cast<double>(iterations)) /
                  std::max(result.seconds, 1e-12) / 1e9;
    result.checksum = sample_checksum(dst, bytes);
    g_probe_sink ^= result.checksum;

    (void)name;
    return result;
}

static unsigned long long parse_ull_arg(const char *s, unsigned long long fallback) {
    if (s == NULL || *s == '\0') return fallback;
    char *end = NULL;
    unsigned long long v = std::strtoull(s, &end, 10);
    if (end == s || v == 0) return fallback;
    return v;
}

static void print_rank_results(const char *name,
                               const ProbeResult &local,
                               int rank,
                               int world_size,
                               MPI_Comm comm) {
    double max_seconds = 0.0;
    double min_gbps = 0.0;
    double max_gbps = 0.0;
    MPI_Reduce(&local.seconds, &max_seconds, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&local.gbps, &min_gbps, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
    MPI_Reduce(&local.gbps, &max_gbps, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    for (int r = 0; r < world_size; ++r) {
        MPI_Barrier(comm);
        if (rank == r) {
            std::printf("[rank %d] %-12s time=%.6f s  bandwidth=%.3f GB/s  checksum=%llu\n",
                        rank,
                        name,
                        local.seconds,
                        local.gbps,
                        static_cast<unsigned long long>(local.checksum));
            std::fflush(stdout);
        }
    }
    MPI_Barrier(comm);

    if (rank == 0) {
        std::printf("[summary] %-12s max_time=%.6f s  rank_bandwidth_min=%.3f GB/s  rank_bandwidth_max=%.3f GB/s\n",
                    name,
                    max_seconds,
                    min_gbps,
                    max_gbps);
        std::fflush(stdout);
    }
    MPI_Barrier(comm);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    unsigned long long mib = parse_ull_arg(argc > 1 ? argv[1] : NULL, 256ULL);
    unsigned long long iterations_ull = parse_ull_arg(argc > 2 ? argv[2] : NULL, 20ULL);
    unsigned long long warmup_ull = parse_ull_arg(argc > 3 ? argv[3] : NULL, 3ULL);

    size_t bytes = static_cast<size_t>(mib) * 1024ULL * 1024ULL;
    int iterations = static_cast<int>(iterations_ull);
    int warmup = static_cast<int>(warmup_ull);

    unsigned char *src = NULL;
    unsigned char *dst = NULL;
    int ok = 1;
    if (posix_memalign(reinterpret_cast<void **>(&src), 64, bytes) != 0) ok = 0;
    if (posix_memalign(reinterpret_cast<void **>(&dst), 64, bytes) != 0) ok = 0;

    int global_ok = 0;
    MPI_Allreduce(&ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!global_ok) {
        if (rank == 0) {
            std::fprintf(stderr, "Failed to allocate %llu MiB source/destination buffers per rank.\n", mib);
        }
        std::free(src);
        std::free(dst);
        MPI_Finalize();
        return 1;
    }

    fill_source(src, bytes, rank);
    std::memset(dst, 0, bytes);

    if (rank == 0) {
        std::printf("MPI memcpy probe: ranks=%d buffer=%llu MiB iterations=%d warmup=%d\n",
                    world_size,
                    mib,
                    iterations,
                    warmup);
        std::fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    ProbeResult memcpy_result = run_probe("memcpy", copy_by_memcpy, dst, src, bytes, iterations, warmup, MPI_COMM_WORLD);
    print_rank_results("memcpy", memcpy_result, rank, world_size, MPI_COMM_WORLD);

    ProbeResult for_u8_result = run_probe("for_u8", copy_by_for_u8, dst, src, bytes, iterations, warmup, MPI_COMM_WORLD);
    print_rank_results("for_u8", for_u8_result, rank, world_size, MPI_COMM_WORLD);

    ProbeResult for_u64_result = run_probe("for_u64", copy_by_for_u64, dst, src, bytes, iterations, warmup, MPI_COMM_WORLD);
    print_rank_results("for_u64", for_u64_result, rank, world_size, MPI_COMM_WORLD);

    if (rank == 0) {
        std::printf("probe_sink=%llu\n", static_cast<unsigned long long>(g_probe_sink));
        std::fflush(stdout);
    }

    std::free(src);
    std::free(dst);
    MPI_Finalize();
    return 0;
}
