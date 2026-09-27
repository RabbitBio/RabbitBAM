#include "crc32_bench.h"

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#include <crts.h>
#endif

#include "slave_zlib/lib_common.h"
#include "slave_zlib/crc32_tables.h"

#define CRC32_LDM_SLICE8_BYTES  (8 * 256 * sizeof(u32))
#define CRC32_LDM_SLICE16_BYTES (16 * 256 * sizeof(u32))

static u32 *g_crc32_slice8_ldm[CRC32_BENCH_MAX_CPES];
static u32 *g_crc32_slice16_ldm[CRC32_BENCH_MAX_CPES];
static signed char g_crc32_slice8_ldm_state[CRC32_BENCH_MAX_CPES];
static signed char g_crc32_slice16_ldm_state[CRC32_BENCH_MAX_CPES];

static inline uint64_t crc32_bench_cycle_now(void)
{
#ifdef PLATFORM_SUNWAY
    unsigned long counter = 0;
    asm volatile("rcsr %0, 4" : "=r"(counter));
    return (uint64_t)counter;
#else
    return 0;
#endif
}

static inline u32 crc32_update_byte(const u32 *table, u32 crc, u8 byte)
{
    return (crc >> 8) ^ table[(u8)crc ^ byte];
}

static u32 crc32_slice8_raw(const u32 *table, u32 crc, const u8 *p, size_t len)
{
    const u8 * const end = p + len;
    const u8 *end64;

    for (; ((uintptr_t)p & 7) && p != end; p++)
        crc = crc32_update_byte(table, crc, *p);

    end64 = p + ((end - p) & ~(size_t)7);
    for (; p != end64; p += 8) {
        u32 v1 = le32_bswap(*(const u32 *)(p + 0));
        u32 v2 = le32_bswap(*(const u32 *)(p + 4));

        crc = table[0x700 + (u8)((crc ^ v1) >> 0)] ^
              table[0x600 + (u8)((crc ^ v1) >> 8)] ^
              table[0x500 + (u8)((crc ^ v1) >> 16)] ^
              table[0x400 + (u8)((crc ^ v1) >> 24)] ^
              table[0x300 + (u8)(v2 >> 0)] ^
              table[0x200 + (u8)(v2 >> 8)] ^
              table[0x100 + (u8)(v2 >> 16)] ^
              table[0x000 + (u8)(v2 >> 24)];
    }

    for (; p != end; p++)
        crc = crc32_update_byte(table, crc, *p);

    return crc;
}

static u32 crc32_slice8_x2_raw(const u32 *table, u32 crc, const u8 *p, size_t len)
{
    const u8 * const end = p + len;
    const u8 *end128;

    for (; ((uintptr_t)p & 7) && p != end; p++)
        crc = crc32_update_byte(table, crc, *p);

    end128 = p + ((end - p) & ~(size_t)15);
    for (; p != end128; p += 16) {
        crc = crc32_slice8_raw(table, crc, p, 8);
        crc = crc32_slice8_raw(table, crc, p + 8, 8);
    }

    return crc32_slice8_raw(table, crc, p, (size_t)(end - p));
}

static void crc32_build_slice16_table(u32 *table)
{
    memcpy(table, crc32_slice8_table, CRC32_LDM_SLICE8_BYTES);
    for (int slice = 8; slice < 16; slice++) {
        u32 *dst = table + slice * 256;
        const u32 *prev = table + (slice - 1) * 256;
        for (int i = 0; i < 256; i++) {
            u32 crc = prev[i];
            dst[i] = (crc >> 8) ^ table[(u8)crc];
        }
    }
}

static u32 crc32_slice16_raw(const u32 *table, u32 crc, const u8 *p, size_t len)
{
    const u8 * const end = p + len;
    const u8 *end128;

    for (; ((uintptr_t)p & 7) && p != end; p++)
        crc = crc32_update_byte(table, crc, *p);

    end128 = p + ((end - p) & ~(size_t)15);
    for (; p != end128; p += 16) {
        u32 v0 = le32_bswap(*(const u32 *)(p + 0));
        u32 v1 = le32_bswap(*(const u32 *)(p + 4));
        u32 v2 = le32_bswap(*(const u32 *)(p + 8));
        u32 v3 = le32_bswap(*(const u32 *)(p + 12));

        crc = table[0xf00 + (u8)((crc ^ v0) >> 0)] ^
              table[0xe00 + (u8)((crc ^ v0) >> 8)] ^
              table[0xd00 + (u8)((crc ^ v0) >> 16)] ^
              table[0xc00 + (u8)((crc ^ v0) >> 24)] ^
              table[0xb00 + (u8)(v1 >> 0)] ^
              table[0xa00 + (u8)(v1 >> 8)] ^
              table[0x900 + (u8)(v1 >> 16)] ^
              table[0x800 + (u8)(v1 >> 24)] ^
              table[0x700 + (u8)(v2 >> 0)] ^
              table[0x600 + (u8)(v2 >> 8)] ^
              table[0x500 + (u8)(v2 >> 16)] ^
              table[0x400 + (u8)(v2 >> 24)] ^
              table[0x300 + (u8)(v3 >> 0)] ^
              table[0x200 + (u8)(v3 >> 8)] ^
              table[0x100 + (u8)(v3 >> 16)] ^
              table[0x000 + (u8)(v3 >> 24)];
    }

    return crc32_slice8_raw(table, crc, p, (size_t)(end - p));
}

static int crc32_get_slice8_ldm(int id, const u32 **table)
{
    if (g_crc32_slice8_ldm_state[id] < 0) return CRC32_BENCH_NO_LDM;
    if (g_crc32_slice8_ldm[id] == NULL) {
#ifdef PLATFORM_SUNWAY
        u32 *ldm_table = (u32 *)ldm_malloc(CRC32_LDM_SLICE8_BYTES);
        if (ldm_table == NULL) {
            g_crc32_slice8_ldm_state[id] = -1;
            return CRC32_BENCH_NO_LDM;
        }
        memcpy(ldm_table, crc32_slice8_table, CRC32_LDM_SLICE8_BYTES);
        g_crc32_slice8_ldm[id] = ldm_table;
        g_crc32_slice8_ldm_state[id] = 1;
#else
        return CRC32_BENCH_UNSUPPORTED;
#endif
    }
    *table = g_crc32_slice8_ldm[id];
    return CRC32_BENCH_OK;
}

static int crc32_get_slice16_ldm(int id, const u32 **table)
{
    if (g_crc32_slice16_ldm_state[id] < 0) return CRC32_BENCH_NO_LDM;
    if (g_crc32_slice16_ldm[id] == NULL) {
#ifdef PLATFORM_SUNWAY
        u32 *ldm_table = (u32 *)ldm_malloc(CRC32_LDM_SLICE16_BYTES);
        if (ldm_table == NULL) {
            g_crc32_slice16_ldm_state[id] = -1;
            return CRC32_BENCH_NO_LDM;
        }
        crc32_build_slice16_table(ldm_table);
        g_crc32_slice16_ldm[id] = ldm_table;
        g_crc32_slice16_ldm_state[id] = 1;
#else
        return CRC32_BENCH_UNSUPPORTED;
#endif
    }
    *table = g_crc32_slice16_ldm[id];
    return CRC32_BENCH_OK;
}

static uint32_t crc32_finish(u32 crc)
{
    return ~crc;
}

static u32 crc32_begin(uint32_t crc)
{
    return ~crc;
}

static uint32_t run_variant(int variant, int id, uint32_t crc,
                            const u8 *data, size_t len, int *status)
{
    const u32 *table = crc32_slice8_table;
    u32 raw = crc32_begin(crc);

    *status = CRC32_BENCH_OK;
    switch (variant) {
    case CRC32_BENCH_SLICE8:
        raw = crc32_slice8_raw(crc32_slice8_table, raw, data, len);
        break;
    case CRC32_BENCH_SLICE8_X2:
        raw = crc32_slice8_x2_raw(crc32_slice8_table, raw, data, len);
        break;
    case CRC32_BENCH_SLICE8_LDM:
        *status = crc32_get_slice8_ldm(id, &table);
        if (*status != CRC32_BENCH_OK) return crc;
        raw = crc32_slice8_raw(table, raw, data, len);
        break;
    case CRC32_BENCH_SLICE16_LDM:
        *status = crc32_get_slice16_ldm(id, &table);
        if (*status != CRC32_BENCH_OK) return crc;
        raw = crc32_slice16_raw(table, raw, data, len);
        break;
    case CRC32_BENCH_SIMD_FOLDING:
        *status = CRC32_BENCH_UNSUPPORTED;
        return crc;
    default:
        *status = CRC32_BENCH_UNSUPPORTED;
        return crc;
    }

    return crc32_finish(raw);
}

void slave_crc32_bench(Crc32BenchArgs *args)
{
#ifdef PLATFORM_SUNWAY
    int id = _PEN;
#else
    int id = 0;
#endif
    if (id < 0 || id >= CRC32_BENCH_MAX_CPES) return;
    if (id >= args->active_cpes) return;

    const u8 *data = args->input + (size_t)id * args->stride + args->misalign;
    for (int variant = 0; variant < CRC32_BENCH_VARIANTS; variant++) {
        int status = CRC32_BENCH_OK;
        uint32_t crc = 0;

        for (int i = 0; i < args->warmup_iterations; i++) {
            crc = run_variant(variant, id, crc, data, args->len, &status);
            if (status != CRC32_BENCH_OK) break;
        }

        crc = 0;
        uint64_t t0 = crc32_bench_cycle_now();
        if (status == CRC32_BENCH_OK) {
            for (int i = 0; i < args->iterations; i++) {
                crc = run_variant(variant, id, crc, data, args->len, &status);
                if (status != CRC32_BENCH_OK) break;
            }
        }
        uint64_t t1 = crc32_bench_cycle_now();

        if (status == CRC32_BENCH_OK && crc != args->expected_crc)
            status = CRC32_BENCH_MISMATCH;

        args->cycles[id][variant] = (status == CRC32_BENCH_OK) ? (t1 - t0) : 0;
        args->bytes[id][variant] = (status == CRC32_BENCH_OK) ?
            (uint64_t)args->len * (uint64_t)args->iterations : 0;
        args->crc[id][variant] = crc;
        args->status[id][variant] = status;
    }
}
