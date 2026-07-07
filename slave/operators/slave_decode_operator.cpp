#include "swbam/cpe_codec.h"

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

extern "C" void slave_swbam_decode_bgzf(SwbamCpeDecodePara paras[64]) {
    const int id = _PEN;
    SwbamCpeDecodePara *para = &paras[id];
    para->alloc_cycles = 0;
    para->inflate_cycles = 0;
    para->crc_cycles = 0;
    para->total_cycles = 0;
    if (!para->compressed) return;

    const unsigned long total_t0 = swbam_cpe_cycle_now();
    struct libdeflate_decompressor *decompressor =
        swbam_cpe_get_decompressor(id, &para->alloc_cycles);
    if (!decompressor || !para->decoded) {
        para->status = -2;
        para->total_cycles =
            (uint64_t)(swbam_cpe_cycle_now() - total_t0);
        return;
    }
    if (swbam_cpe_decode_bgzf(
            para->compressed, para->decoded, decompressor,
            &para->inflate_cycles, &para->crc_cycles) != 0) {
        para->status = -2;
        para->total_cycles =
            (uint64_t)(swbam_cpe_cycle_now() - total_t0);
        return;
    }
    para->status = 0;
    para->total_cycles =
        (uint64_t)(swbam_cpe_cycle_now() - total_t0);
}
