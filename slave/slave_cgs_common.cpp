#include "BamTools.h"

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

int writeBam1_to_block(bam_block *&write_block, bam1_t *b, int is_be);
int block_encode_func(bam_block *un_comp, bam_block *comp, int compress_level);

static inline int cgs_tid_common() {
#ifdef PLATFORM_SUNWAY
    return (int)_CGN * CGS_PES_PER_CG + (int)_PEN;
#else
    return 0;
#endif
}

extern "C" void slave_cgs_compressfunc(Comp_Para paras[CGS_NB]) {
    int id = cgs_tid_common();
    if (id < 0 || id >= CGS_NB) return;

    Comp_Para *para = &paras[id];
    if (para->status != 0 || para->input_records == nullptr || para->n_records == 0) return;

    bam_block *uncompressed = para->un_comp_block;
    bam_block *compressed = para->output_block;
    uncompressed->pos = 0;
    uncompressed->length = 0;
    uncompressed->errcode = 0;
    uncompressed->block_id = 0;
    uncompressed->block_address = 0;

    for (int i = 0; i < para->n_records; ++i) {
        writeBam1_to_block(uncompressed, para->input_records[i], 0);
    }

    compressed->length = block_encode_func(uncompressed, compressed, 1);
    para->output_size = compressed->length;
    para->status = compressed->length > 0 ? 0 : -2;
}
