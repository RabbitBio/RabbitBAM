#include <cstring>

#include "BamTools.h"

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

int block_decode_func(struct bam_block *comp, struct bam_block *un_comp);
int read_bam(struct bam_block *fq, bam1_t *b, int is_be);
int sam_format1_append(const bam_hdr_t *h, const bam1_t *b, kstring_t *str);

static inline int cgs_tid_bam2sam() {
#ifdef PLATFORM_SUNWAY
    return (int)_CGN * CGS_PES_PER_CG + (int)_PEN;
#else
    return 0;
#endif
}

extern "C" void slave_cgs_bam_decodefunc(CgsBamDecodePara paras[CGS_NB]) {
    int id = cgs_tid_bam2sam();
    if (id < 0 || id >= CGS_NB) return;

    CgsBamDecodePara *para = &paras[id];
    bam_block *comp = para->input_block;
    bam_block *un_comp = para->un_comp_block;
    if (comp == nullptr) return;

    if (block_decode_func(comp, un_comp) != 0) {
        para->status = -2;
        return;
    }

    int count = 0;
    int ret = -1;
    while (count < (int)MAX_RECORDS_PER_BLOCK) {
        bam1_t *b = para->output_records[count];
        ret = read_bam(un_comp, b, 0);
        if (ret < 0) break;
        count++;
    }

    para->n_records = count;
    if (ret < -1) {
        para->status = -2;
        return;
    }
    if (count == (int)MAX_RECORDS_PER_BLOCK && un_comp->pos < un_comp->length) {
        para->status = -3;
        return;
    }
    para->status = 0;
}

extern "C" void slave_cgs_sam_format_tile(CgsSamFormatBatch *batch) {
    int id = cgs_tid_bam2sam();
    if (id < 0 || id >= CGS_NB) return;

    kstring_t *out = &batch->core_out_lines[id];
    out->l = 0;
    batch->formatted_records[id] = 0;
    batch->status[id] = 0;

    int start = id * CGS_BAM2SAM_FORMAT_RECORDS_PER_CPE;
    if (start >= batch->total_records) return;
    int end = start + CGS_BAM2SAM_FORMAT_RECORDS_PER_CPE;
    if (end > batch->total_records) end = batch->total_records;

    for (int i = start; i < end; ++i) {
        int ret = sam_format1_append(batch->hdr, batch->records[i], out);
        if (ret < 0) {
            batch->status[id] = -2;
            return;
        }
        if (out->l + 1 > out->m) {
            batch->status[id] = -3;
            return;
        }
        out->s[out->l++] = '\n';
        batch->formatted_records[id]++;
    }
}
