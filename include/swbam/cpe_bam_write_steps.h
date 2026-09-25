#ifndef SWBAM_CPE_BAM_WRITE_STEPS_H
#define SWBAM_CPE_BAM_WRITE_STEPS_H

#include "swbam/cpe_codec.h"

// Shared BGZF compression step. The caller chooses whether it first serializes
// bam1_t records or supplies an already packed uncompressed block.
inline int swbam_cpe_compress_bam_block(
        const bam_block *uncompressed, size_t uncompressed_size,
        bam_block *compressed, int level, int cpe_id, int mpi_tuned,
        uint64_t *alloc_cycles, uint64_t *deflate_cycles,
        uint64_t *footer_cycles) {
    if (!uncompressed || !compressed ||
        (level != 0 && level != 1 && level != 6) ||
        uncompressed_size > BGZF_BLOCK_SIZE) return -1;
    struct libdeflate_compressor *compressor = nullptr;
    if (level != 0) {
        compressor = swbam_cpe_get_compressor(
            cpe_id, level, mpi_tuned, alloc_cycles);
        if (!compressor) return -1;
    }
    size_t output_size = BGZF_MAX_BLOCK_SIZE;
    if (swbam_cpe_compress_bgzf(
            compressed->data, &output_size, uncompressed->data,
            uncompressed_size, level, compressor,
            deflate_cycles, footer_cycles) != 0) return -1;
    compressed->length = static_cast<unsigned int>(output_size);
    return 0;
}

#endif
