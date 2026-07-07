#include "swbam/cpe_codec.h"

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

extern "C" void slave_swbam_compress_bgzf(SwbamCpeCompressPara paras[64]) {
    const int id = _PEN;
    SwbamCpeCompressPara *para = &paras[id];
    para->alloc_cycles = 0;
    para->deflate_cycles = 0;
    para->footer_cycles = 0;
    para->total_cycles = 0;
    para->output_size = 0;
    if (!para->uncompressed) return;

    const unsigned long total_t0 = swbam_cpe_cycle_now();
    struct libdeflate_compressor *compressor = nullptr;
    if (para->level != 0) {
        compressor = swbam_cpe_get_compressor(
            id, para->level, 1, &para->alloc_cycles);
    }
    if (para->level != 0 && !compressor) {
        para->status = -2;
        para->total_cycles =
            (uint64_t)(swbam_cpe_cycle_now() - total_t0);
        return;
    }

    size_t output_size = BGZF_MAX_BLOCK_SIZE;
    const int ret = swbam_cpe_compress_bgzf(
        para->compressed->data, &output_size,
        para->uncompressed->data, para->uncompressed->length,
        para->level, compressor,
        &para->deflate_cycles, &para->footer_cycles);
    if (ret != 0) {
        para->status = -2;
        para->total_cycles =
            (uint64_t)(swbam_cpe_cycle_now() - total_t0);
        return;
    }
    para->compressed->length = (unsigned int)output_size;
    para->compressed->pos = 0;
    para->compressed->block_id = para->uncompressed->block_id;
    para->output_size = (int)output_size;
    para->status = 0;
    para->total_cycles =
        (uint64_t)(swbam_cpe_cycle_now() - total_t0);
}
