#include <cstring>

#include "BamTools.h"

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

static inline int cgs_tid_sam2bam() {
#ifdef PLATFORM_SUNWAY
    return (int)_CGN * CGS_PES_PER_CG + (int)_PEN;
#else
    return 0;
#endif
}

extern "C" void slave_cgs_copy_and_count(CgsSamParseBatch *batch) {
    int id = cgs_tid_sam2bam();
    if (id < 0 || id >= CGS_NB) return;

    CgsSamParseChunk *chunk = &batch->chunks[id];
    if (chunk->src_len == 0) {
        chunk->text_len = 0;
        chunk->count = 0;
        chunk->status = 0;
        return;
    }

    memcpy(chunk->text_buf, chunk->src_ptr, chunk->src_len);
    chunk->text_buf[chunk->src_len] = '\0';
    chunk->text_len = chunk->src_len;

    int count = 0;
    char *p = chunk->text_buf;
    char *end = chunk->text_buf + chunk->text_len;
    while (p < end) {
        char *line_end = p;
        while (line_end < end && *line_end != '\n' && *line_end != '\0') line_end++;

        int line_len = (int)(line_end - p);
        if (line_len > 0 && p[line_len - 1] == '\r') line_len--;
        if (line_len > 0) {
            count++;
            if (count > CGS_SAM2BAM_MAX_BAMS_PER_CHUNK) {
                chunk->count = count;
                chunk->status = -3;
                return;
            }
        }

        if (line_end < end && *line_end == '\n') line_end++;
        p = line_end;
    }
    chunk->count = count;
    chunk->status = 0;
}

extern "C" void slave_cgs_sam_parse_chunk(CgsSamParseBatch *batch) {
    int id = cgs_tid_sam2bam();
    if (id < 0 || id >= CGS_NB) return;

    CgsSamParseChunk *chunk = &batch->chunks[id];
    if (chunk->text_len == 0 || chunk->count == 0) {
        chunk->count = 0;
        chunk->status = 0;
        return;
    }

    int valid_count = 0;
    char *ptr = chunk->text_buf;
    char *end = chunk->text_buf + chunk->text_len;

    while (ptr < end) {
        char *eol = ptr;
        while (eol < end && *eol != '\n' && *eol != '\0') eol++;

        int line_len = (int)(eol - ptr);
        kstring_t ks;
        ks.s = ptr;
        ks.l = line_len;
        ks.m = line_len + 1;

        if (ks.l > 0 && ks.s[ks.l - 1] == '\r') {
            ks.l--;
        }

        if (ks.l > 0) {
            if (valid_count >= CGS_SAM2BAM_MAX_BAMS_PER_CHUNK) {
                chunk->status = -3;
                return;
            }
            char saved_char = ks.s[ks.l];
            ks.s[ks.l] = '\0';

            int ret = sam_parse1(&ks, (sam_hdr_t *)batch->hdr, chunk->bams[valid_count]);

            ks.s[ks.l] = saved_char;

            if (ret >= 0) {
                bam1_t *b = chunk->bams[valid_count];
                chunk->bam_lens[valid_count] = (uint32_t)(b->l_data - b->core.l_extranul + 32);
                valid_count++;
            }
        }

        if (eol < end && *eol == '\n') eol++;
        ptr = eol;
    }
    chunk->count = valid_count;
    chunk->status = 0;
}
