#ifndef SWBAM_CPE_CODEC_H
#define SWBAM_CPE_CODEC_H

#include "BamTools.h"

#include <cstddef>
#include <cstdint>

struct libdeflate_compressor;
struct libdeflate_decompressor;

struct SwbamCpeDecodePara {
    bam_block *compressed;
    bam_block *decoded;
    int status;
    uint64_t alloc_cycles;
    uint64_t inflate_cycles;
    uint64_t crc_cycles;
    uint64_t total_cycles;
};

struct SwbamCpeCompressPara {
    bam_block *uncompressed;
    bam_block *compressed;
    int level;
    int status;
    int output_size;
    uint64_t alloc_cycles;
    uint64_t deflate_cycles;
    uint64_t footer_cycles;
    uint64_t total_cycles;
};

inline unsigned long swbam_cpe_cycle_now() {
#ifdef PLATFORM_SUNWAY
    unsigned long rpcc = 0;
    asm volatile("rcsr %0, 4" : "=r"(rpcc));
    return rpcc;
#else
    return 0;
#endif
}

struct libdeflate_decompressor *swbam_cpe_get_decompressor(
    int cpe_id, uint64_t *alloc_cycles);

struct libdeflate_compressor *swbam_cpe_get_compressor(
    int cpe_id, int level, int mpi_tuned, uint64_t *alloc_cycles);

int swbam_cpe_decode_bgzf(bam_block *compressed,
                          bam_block *decoded,
                          struct libdeflate_decompressor *decompressor,
                          uint64_t *inflate_cycles,
                          uint64_t *crc_cycles);

int swbam_cpe_compress_bgzf(void *dst, size_t *dst_len,
                            const void *src, size_t src_len,
                            int level,
                            struct libdeflate_compressor *compressor,
                            uint64_t *deflate_cycles,
                            uint64_t *footer_cycles);

extern "C" void slave_mpi_common_release_caches(void *unused);

#endif
