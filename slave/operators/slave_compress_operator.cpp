#include "swbam/cpe_bam_write_steps.h"

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
    if (swbam_cpe_compress_bam_block(
            para->uncompressed, para->uncompressed->length,
            para->compressed, para->level, id, 1,
            &para->alloc_cycles, &para->deflate_cycles,
            &para->footer_cycles) != 0) {
        para->status = -2;
        para->total_cycles =
            (uint64_t)(swbam_cpe_cycle_now() - total_t0);
        return;
    }
    para->compressed->pos = 0;
    para->compressed->block_id = para->uncompressed->block_id;
    para->output_size = para->compressed->length;
    para->status = 0;
    para->total_cycles =
        (uint64_t)(swbam_cpe_cycle_now() - total_t0);
}
