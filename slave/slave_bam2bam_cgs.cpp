#include <cstring>

#include "BamTools.h"

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

int block_decode_func(struct bam_block *comp, struct bam_block *un_comp);
int read_bam(struct bam_block *fq, bam1_t *b, int is_be);

static inline int cgs_tid() {
#ifdef PLATFORM_SUNWAY
    return (int)_CGN * CGS_PES_PER_CG + (int)_PEN;
#else
    return 0;
#endif
}

extern "C" void slave_cgs_decompress_filterfunc(Bam2BamPara paras[CGS_NB]) {
    int id = cgs_tid();
    if (id < 0 || id >= CGS_NB) return;

    Bam2BamPara *para = &paras[id];
    bam_block *comp = para->input_block;
    bam_block *un_comp = para->un_comp_block;
    if (comp == nullptr) return;

    if (block_decode_func(comp, un_comp) != 0) {
        para->status = -2;
        return;
    }

    int total_count = 0;
    int kept_count = 0;
    int ret = -1;
    while (total_count < (int)MAX_RECORDS_PER_BLOCK) {
        bam1_t *b = para->record_base ? para->record_base + total_count
                                      : para->output_records[total_count];
        ret = read_bam(un_comp, b, 0);
        if (ret < 0) break;

        total_count++;
        if (bam_filter_matches(b, para->filter)) {
            para->output_records[kept_count] = b;
            para->bam_lens[kept_count] = (uint32_t)(b->l_data - b->core.l_extranul + 32);
            kept_count++;
        }
    }

    para->n_total_records = total_count;
    para->n_kept_records = kept_count;
    if (ret < -1) {
        para->status = -2;
        return;
    }
    if (total_count == (int)MAX_RECORDS_PER_BLOCK && un_comp->pos < un_comp->length) {
        para->status = -3;
        return;
    }
    para->status = 0;
}

extern "C" void slave_cgs_decompress_bam2bam_passthrough(Bam2BamPara paras[CGS_NB]) {
    int id = cgs_tid();
    if (id < 0 || id >= CGS_NB) return;

    Bam2BamPara *para = &paras[id];
    bam_block *comp = para->input_block;
    bam_block *un_comp = para->un_comp_block;
    if (comp == nullptr) return;

    if (block_decode_func(comp, un_comp) != 0) {
        para->status = -2;
        return;
    }

    int total_count = 0;
    int ret = -1;
    while (total_count < (int)MAX_RECORDS_PER_BLOCK) {
        bam1_t *b = para->output_records[total_count];
        ret = read_bam(un_comp, b, 0);
        if (ret < 0) break;

        para->bam_lens[total_count] = (uint32_t)(b->l_data - b->core.l_extranul + 32);
        total_count++;
    }

    para->n_total_records = total_count;
    para->n_kept_records = total_count;
    if (ret < -1) {
        para->status = -2;
        return;
    }
    if (total_count == (int)MAX_RECORDS_PER_BLOCK && un_comp->pos < un_comp->length) {
        para->status = -3;
        return;
    }
    para->status = 0;
}
