#ifndef RABBITBAM_CRC32_BENCH_H
#define RABBITBAM_CRC32_BENCH_H

#include <stddef.h>
#include <stdint.h>

#define CRC32_BENCH_MAX_CPES 64
#define CRC32_BENCH_VARIANTS 5

enum Crc32BenchVariant {
    CRC32_BENCH_SLICE8 = 0,
    CRC32_BENCH_SLICE8_X2 = 1,
    CRC32_BENCH_SLICE8_LDM = 2,
    CRC32_BENCH_SLICE16_LDM = 3,
    CRC32_BENCH_SIMD_FOLDING = 4
};

enum Crc32BenchStatus {
    CRC32_BENCH_OK = 0,
    CRC32_BENCH_MISMATCH = 1,
    CRC32_BENCH_UNSUPPORTED = 2,
    CRC32_BENCH_NO_LDM = 3
};

typedef struct Crc32BenchArgs {
    uint8_t *input;
    size_t len;
    size_t stride;
    int iterations;
    int warmup_iterations;
    int active_cpes;
    int misalign;
    uint32_t expected_crc;
    uint64_t cycles[CRC32_BENCH_MAX_CPES][CRC32_BENCH_VARIANTS];
    uint64_t bytes[CRC32_BENCH_MAX_CPES][CRC32_BENCH_VARIANTS];
    uint32_t crc[CRC32_BENCH_MAX_CPES][CRC32_BENCH_VARIANTS];
    int status[CRC32_BENCH_MAX_CPES][CRC32_BENCH_VARIANTS];
} Crc32BenchArgs;

#ifdef __cplusplus
extern "C" {
#endif

void slave_crc32_bench(Crc32BenchArgs *args);

#ifdef __cplusplus
}
#endif

#endif
