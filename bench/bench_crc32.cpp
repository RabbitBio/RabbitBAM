#include "crc32_bench.h"
#include "Globals.h"

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

#include <cerrno>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {

const char *kVariantNames[CRC32_BENCH_VARIANTS] = {
    "slice8",
    "slice8x2",
    "slice8_ldm",
    "slice16_ldm",
    "simd_folding"
};

uint32_t g_host_crc_table[256];

void InitHostCrcTable() {
    for (uint32_t i = 0; i < 256; i++) {
        uint32_t crc = i;
        for (int bit = 0; bit < 8; bit++) {
            crc = (crc >> 1) ^ ((crc & 1) ? 0xedb88320U : 0U);
        }
        g_host_crc_table[i] = crc;
    }
}

uint32_t HostCrc32Update(uint32_t crc, const uint8_t *data, size_t len) {
    crc = ~crc;
    for (size_t i = 0; i < len; i++) {
        crc = (crc >> 8) ^ g_host_crc_table[(crc ^ data[i]) & 0xffU];
    }
    return ~crc;
}

size_t AlignUp(size_t value, size_t align) {
    return (value + align - 1) & ~(align - 1);
}

bool ParseSize(const char *text, size_t *out) {
    if (text == nullptr || *text == '\0') return false;
    char *end = nullptr;
    errno = 0;
    unsigned long long value = strtoull(text, &end, 10);
    if (errno != 0 || end == text) return false;
    size_t scale = 1;
    if (*end != '\0') {
        if ((end[1] != '\0')) return false;
        switch (*end) {
        case 'k':
        case 'K':
            scale = 1024ULL;
            break;
        case 'm':
        case 'M':
            scale = 1024ULL * 1024ULL;
            break;
        case 'g':
        case 'G':
            scale = 1024ULL * 1024ULL * 1024ULL;
            break;
        default:
            return false;
        }
    }
    *out = (size_t)value * scale;
    return true;
}

void PrintUsage(const char *argv0) {
    std::fprintf(stderr,
                 "Usage: %s [--len 64K] [--iters 200] [--warmup 2] "
                 "[--cpes 64] [--misalign 0]\n",
                 argv0);
}

const char *StatusName(int status) {
    switch (status) {
    case CRC32_BENCH_OK:
        return "ok";
    case CRC32_BENCH_MISMATCH:
        return "mismatch";
    case CRC32_BENCH_UNSUPPORTED:
        return "unsupported";
    case CRC32_BENCH_NO_LDM:
        return "no_ldm";
    default:
        return "unknown";
    }
}

}  // namespace

int main(int argc, char **argv) {
#ifndef PLATFORM_SUNWAY
    std::fprintf(stderr, "RabbitBAM-CRC32-Bench is only supported on Sunway.\n");
    return 1;
#else
    size_t len = 64 * 1024;
    int iterations = 200;
    int warmup = 2;
    int active_cpes = 64;
    int misalign = 0;

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--len") == 0 && i + 1 < argc) {
            if (!ParseSize(argv[++i], &len)) {
                PrintUsage(argv[0]);
                return 1;
            }
        } else if (std::strcmp(argv[i], "--iters") == 0 && i + 1 < argc) {
            iterations = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--warmup") == 0 && i + 1 < argc) {
            warmup = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--cpes") == 0 && i + 1 < argc) {
            active_cpes = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--misalign") == 0 && i + 1 < argc) {
            misalign = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--help") == 0 ||
                   std::strcmp(argv[i], "-h") == 0) {
            PrintUsage(argv[0]);
            return 0;
        } else {
            PrintUsage(argv[0]);
            return 1;
        }
    }

    if (len == 0 || iterations <= 0 || warmup < 0 ||
        active_cpes <= 0 || active_cpes > CRC32_BENCH_MAX_CPES ||
        misalign < 0 || misalign >= 64) {
        PrintUsage(argv[0]);
        return 1;
    }

    InitHostCrcTable();

    const size_t stride = AlignUp(len + (size_t)misalign + 64, 64);
    uint8_t *input = nullptr;
    Crc32BenchArgs *args = nullptr;
    if (posix_memalign(reinterpret_cast<void **>(&input), 64,
                       stride * CRC32_BENCH_MAX_CPES) != 0) {
        std::fprintf(stderr, "failed to allocate input buffer\n");
        return 1;
    }
    if (posix_memalign(reinterpret_cast<void **>(&args), 64,
                       sizeof(Crc32BenchArgs)) != 0) {
        std::fprintf(stderr, "failed to allocate args\n");
        std::free(input);
        return 1;
    }
    std::memset(args, 0, sizeof(*args));

    for (int cpe = 0; cpe < CRC32_BENCH_MAX_CPES; cpe++) {
        uint8_t *data = input + (size_t)cpe * stride + misalign;
        for (size_t i = 0; i < len; i++) {
            data[i] = static_cast<uint8_t>((i * 131U + (i >> 3) + 17U) & 0xffU);
        }
    }

    uint32_t expected = 0;
    const uint8_t *reference_data = input + misalign;
    for (int i = 0; i < iterations; i++) {
        expected = HostCrc32Update(expected, reference_data, len);
    }

    args->input = input;
    args->len = len;
    args->stride = stride;
    args->iterations = iterations;
    args->warmup_iterations = warmup;
    args->active_cpes = active_cpes;
    args->misalign = misalign;
    args->expected_crc = expected;

    std::printf("RabbitBAM CRC32 CPE benchmark: len=%zu iters=%d warmup=%d cpes=%d misalign=%d expected=%08x\n",
                len, iterations, warmup, active_cpes, misalign, expected);

    athread_init();
    double t0 = GetTime();
    __real_athread_spawn(reinterpret_cast<void *>(slave_crc32_bench), args, 1);
    athread_join();
    double elapsed = GetTime() - t0;
    athread_halt();

    double baseline_avg_cycles = 0.0;
    for (int variant = 0; variant < CRC32_BENCH_VARIANTS; variant++) {
        int ok = 0;
        int bad = 0;
        uint64_t min_cycles = UINT64_MAX;
        uint64_t max_cycles = 0;
        long double sum_cycles = 0.0;
        uint64_t bytes_per_cpe = 0;
        int first_bad_status = CRC32_BENCH_OK;

        for (int cpe = 0; cpe < active_cpes; cpe++) {
            int status = args->status[cpe][variant];
            if (status == CRC32_BENCH_OK) {
                uint64_t cycles = args->cycles[cpe][variant];
                ok++;
                if (cycles < min_cycles) min_cycles = cycles;
                if (cycles > max_cycles) max_cycles = cycles;
                sum_cycles += (long double)cycles;
                bytes_per_cpe = args->bytes[cpe][variant];
            } else {
                bad++;
                if (first_bad_status == CRC32_BENCH_OK) first_bad_status = status;
            }
        }

        if (ok == 0) {
            std::printf("%-14s status=%s ok=%d bad=%d\n",
                        kVariantNames[variant], StatusName(first_bad_status), ok, bad);
            continue;
        }

        double avg_cycles = static_cast<double>(sum_cycles / ok);
        if (variant == CRC32_BENCH_SLICE8) baseline_avg_cycles = avg_cycles;
        double cycles_per_byte = avg_cycles / static_cast<double>(bytes_per_cpe);
        double speedup = baseline_avg_cycles > 0.0 ?
            baseline_avg_cycles / avg_cycles : 1.0;
        double gib_per_cycle = (static_cast<double>(bytes_per_cpe) / avg_cycles) /
            (1024.0 * 1024.0 * 1024.0);

        std::printf("%-14s ok=%d bad=%d avg_cycles=%.0f min=%" PRIu64
                    " max=%" PRIu64 " cycles_per_byte=%.6f speedup_vs_slice8=%.4f"
                    " bytes_per_cpe=%" PRIu64 " gib_per_cycle=%.9f\n",
                    kVariantNames[variant], ok, bad, avg_cycles,
                    min_cycles, max_cycles, cycles_per_byte, speedup,
                    bytes_per_cpe, gib_per_cycle);
    }

    std::printf("wall_time=%.6f sec\n", elapsed);

    std::free(args);
    std::free(input);
    return 0;
#endif
}
