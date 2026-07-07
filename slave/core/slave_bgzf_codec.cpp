#include "swbam/cpe_codec.h"
#include "../libdeflate.h"

#include <cstdio>
#include <cstring>

#include <htslib/bgzf.h>
#include <htslib/hts_endian.h>

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

#ifndef RABBITBAM_MPI_LEVEL1_NICE_LEN
#define RABBITBAM_MPI_LEVEL1_NICE_LEN 4
#endif

#ifndef RABBITBAM_MPI_LEVEL1_SINGLE_PROBE
#define RABBITBAM_MPI_LEVEL1_SINGLE_PROBE 1
#endif

#ifndef RABBITBAM_MPI_LEVEL1_STRIDE2_PROBE
#define RABBITBAM_MPI_LEVEL1_STRIDE2_PROBE 1
#endif

#if defined(PLATFORM_SUNWAY) && defined(RABBITBAM_ENABLE_SUNWAY_CRC16_LDM)
extern "C" void libdeflate_crc32_sunway_set_ldm_enabled(int enabled);
#endif

namespace {

struct CompressorCache {
    struct libdeflate_compressor *items[64];
    int levels[64];
};

CompressorCache g_default_compressors = {};
CompressorCache g_mpi_compressors = {};
struct libdeflate_decompressor *g_decompressors[64] = {};

void ConfigureMpiLevel1(struct libdeflate_compressor *compressor,
                        int level) {
    if (level != 1 || !compressor) return;
    libdeflate_set_level1_nice_match_length(
        compressor, RABBITBAM_MPI_LEVEL1_NICE_LEN);
#if RABBITBAM_MPI_LEVEL1_SINGLE_PROBE
    libdeflate_set_level1_single_probe(compressor, 1);
#endif
#if RABBITBAM_MPI_LEVEL1_STRIDE2_PROBE
    libdeflate_set_level1_stride2_probe(compressor, 1);
#endif
}

int InflateBgzf(uint8_t *dst, size_t *dst_len,
                const uint8_t *src, size_t src_len,
                uint32_t expected_crc,
                struct libdeflate_decompressor *decompressor,
                uint64_t *inflate_cycles,
                uint64_t *crc_cycles) {
    if (!decompressor) return -1;

    const unsigned long inflate_t0 = swbam_cpe_cycle_now();
    const int ret = libdeflate_deflate_decompress(
        decompressor, src, src_len, dst, *dst_len, dst_len);
    if (inflate_cycles) {
        *inflate_cycles +=
            (uint64_t)(swbam_cpe_cycle_now() - inflate_t0);
    }
    if (ret != 0) return -1;

    const unsigned long crc_t0 = swbam_cpe_cycle_now();
    const uint32_t crc = libdeflate_crc32(0, dst, *dst_len);
    if (crc_cycles) {
        *crc_cycles += (uint64_t)(swbam_cpe_cycle_now() - crc_t0);
    }
    return crc == expected_crc ? 0 : -2;
}

} // namespace

struct libdeflate_decompressor *swbam_cpe_get_decompressor(
        int cpe_id, uint64_t *alloc_cycles) {
    if (cpe_id >= 0 && cpe_id < 64) {
        struct libdeflate_decompressor *decompressor =
            g_decompressors[cpe_id];
        if (decompressor) return decompressor;
        const unsigned long t0 = swbam_cpe_cycle_now();
        decompressor = libdeflate_alloc_decompressor();
        if (alloc_cycles) {
            *alloc_cycles += (uint64_t)(swbam_cpe_cycle_now() - t0);
        }
        g_decompressors[cpe_id] = decompressor;
        return decompressor;
    }

    const unsigned long t0 = swbam_cpe_cycle_now();
    struct libdeflate_decompressor *decompressor =
        libdeflate_alloc_decompressor();
    if (alloc_cycles) {
        *alloc_cycles += (uint64_t)(swbam_cpe_cycle_now() - t0);
    }
    return decompressor;
}

struct libdeflate_compressor *swbam_cpe_get_compressor(
        int cpe_id, int level, int mpi_tuned, uint64_t *alloc_cycles) {
    CompressorCache *cache = mpi_tuned
        ? &g_mpi_compressors : &g_default_compressors;
    if (cpe_id >= 0 && cpe_id < 64) {
        struct libdeflate_compressor *compressor = cache->items[cpe_id];
        if (compressor && cache->levels[cpe_id] == level) return compressor;
        const unsigned long t0 = swbam_cpe_cycle_now();
        if (compressor) libdeflate_free_compressor(compressor);
        compressor = libdeflate_alloc_compressor(level);
        if (mpi_tuned) ConfigureMpiLevel1(compressor, level);
        if (alloc_cycles) {
            *alloc_cycles += (uint64_t)(swbam_cpe_cycle_now() - t0);
        }
        cache->items[cpe_id] = compressor;
        cache->levels[cpe_id] = compressor ? level : 0;
        return compressor;
    }

    const unsigned long t0 = swbam_cpe_cycle_now();
    struct libdeflate_compressor *compressor =
        libdeflate_alloc_compressor(level);
    if (mpi_tuned) ConfigureMpiLevel1(compressor, level);
    if (alloc_cycles) {
        *alloc_cycles += (uint64_t)(swbam_cpe_cycle_now() - t0);
    }
    return compressor;
}

extern "C" void slave_mpi_common_release_caches(void *unused) {
    (void)unused;
    const int cpe_id = _PEN;
    if (cpe_id < 0 || cpe_id >= 64) return;

    CompressorCache *caches[2] = {
        &g_default_compressors, &g_mpi_compressors
    };
    for (int i = 0; i < 2; ++i) {
        if (caches[i]->items[cpe_id]) {
            libdeflate_free_compressor(caches[i]->items[cpe_id]);
            caches[i]->items[cpe_id] = nullptr;
            caches[i]->levels[cpe_id] = 0;
        }
    }
    if (g_decompressors[cpe_id]) {
        libdeflate_free_decompressor(g_decompressors[cpe_id]);
        g_decompressors[cpe_id] = nullptr;
    }
}

int swbam_cpe_decode_bgzf(bam_block *compressed,
                          bam_block *decoded,
                          struct libdeflate_decompressor *decompressor,
                          uint64_t *inflate_cycles,
                          uint64_t *crc_cycles) {
    if (!compressed || !decoded ||
        compressed->length < BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH) {
        return -1;
    }
    decoded->pos = 0;
    decoded->length = BGZF_MAX_BLOCK_SIZE;
    const uint32_t crc = le_to_u32(
        (uint8_t *)compressed->data + compressed->length - 8);
    size_t decoded_len = BGZF_MAX_BLOCK_SIZE;
    const int ret = InflateBgzf(
        decoded->data, &decoded_len,
        compressed->data + BLOCK_HEADER_LENGTH,
        compressed->length - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH,
        crc, decompressor, inflate_cycles, crc_cycles);
    decoded->length = (unsigned int)decoded_len;
    if (ret != 0) decoded->errcode |= BGZF_ERR_ZLIB;
    return ret;
}

int swbam_cpe_compress_bgzf(void *dst_data, size_t *dst_len,
                            const void *src, size_t src_len,
                            int level,
                            struct libdeflate_compressor *compressor,
                            uint64_t *deflate_cycles,
                            uint64_t *footer_cycles) {
    if (!dst_data || !dst_len || (!src && src_len != 0)) return -1;
    if (src_len == 0) return 0;
    if (level != 0 && !compressor) return -1;
    uint8_t *dst = (uint8_t *)dst_data;

    if (level == 0) {
        if (*dst_len < src_len + 5 + BLOCK_HEADER_LENGTH +
                           BLOCK_FOOTER_LENGTH) {
            return -1;
        }
        dst[BLOCK_HEADER_LENGTH] = 1;
        packInt16(&dst[BLOCK_HEADER_LENGTH + 1], (uint16_t)src_len);
        packInt16(&dst[BLOCK_HEADER_LENGTH + 3],
                  (uint16_t)~(uint16_t)src_len);
        memcpy(dst + BLOCK_HEADER_LENGTH + 5, src, src_len);
        *dst_len = src_len + 5 + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    } else {
        const unsigned long deflate_t0 = swbam_cpe_cycle_now();
        const size_t compressed_len = libdeflate_deflate_compress(
            compressor, src, src_len, dst + BLOCK_HEADER_LENGTH,
            *dst_len - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH);
        if (deflate_cycles) {
            *deflate_cycles +=
                (uint64_t)(swbam_cpe_cycle_now() - deflate_t0);
        }
        if (compressed_len == 0) return -1;
        *dst_len = compressed_len + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    }

    const unsigned long footer_t0 = swbam_cpe_cycle_now();
    memcpy(dst, g_magic, BLOCK_HEADER_LENGTH);
    packInt16(&dst[16], *dst_len - 1);
#if defined(PLATFORM_SUNWAY) && defined(RABBITBAM_ENABLE_SUNWAY_CRC16_LDM)
    libdeflate_crc32_sunway_set_ldm_enabled(level != 0 && compressor != nullptr);
#endif
    const uint32_t crc = libdeflate_crc32(0, src, src_len);
#if defined(PLATFORM_SUNWAY) && defined(RABBITBAM_ENABLE_SUNWAY_CRC16_LDM)
    libdeflate_crc32_sunway_set_ldm_enabled(0);
#endif
    packInt32(&dst[*dst_len - 8], crc);
    packInt32(&dst[*dst_len - 4], src_len);
    if (footer_cycles) {
        *footer_cycles += (uint64_t)(swbam_cpe_cycle_now() - footer_t0);
    }
    return 0;
}
