#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include "swbam/bam_types.h"
#include <climits>
#include <cstring>
#include <htslib/hts_endian.h>

#define RB_LIKELY(x) __builtin_expect(!!(x), 1)
#define RB_UNLIKELY(x) __builtin_expect(!!(x), 0)
#include "swbam/cpe_bam_write_steps.h"
#include <slave.h>
#define slave_cycle_now swbam_cpe_cycle_now

int rabbit_bgzf_mul_write(bam_block *&write_block, const void *data, size_t length) {
    const uint8_t *input = (const uint8_t *) data;
    ssize_t remaining = length;
    uint8_t *buffer = (uint8_t *) write_block->data;
    int copy_length = BGZF_BLOCK_SIZE - write_block->pos;
    if (copy_length > remaining) copy_length = remaining;

    memcpy(buffer + write_block->pos, input, copy_length);

    write_block->pos += copy_length;
    input += copy_length;
    remaining -= copy_length;
    return length - remaining;
}

//解析并把一条记录写入未压缩区中
int writeBam1_to_block(bam_block *&write_block, bam1_t *b , int is_be) {
    const bam1_core_t *c = &b->core;
    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
    int i, ok;
    if (c->l_qname - c->l_extranul > 255) {
        hts_log_error("QNAME \"%s\" is longer than 254 characters", bam_get_qname(b));
        errno = EOVERFLOW;
        return -1;
    }
    if (c->n_cigar > 0xffff) block_len += 16; // "16" for "CGBI", 4-byte tag length and 8-byte fake CIGAR


    if (c->pos > INT_MAX ||
        c->mpos > INT_MAX ||
        c->isize < INT_MIN || c->isize > INT_MAX) {
        hts_log_error("Positional data is too large for BAM format");
        return -1;
    }
    x[0] = c->tid;
    x[1] = c->pos;
    x[2] = (uint32_t) c->bin << 16 | c->qual << 8 | (c->l_qname - c->l_extranul);
    if (c->n_cigar > 0xffff) x[3] = (uint32_t) c->flag << 16 | 2;
    else x[3] = (uint32_t) c->flag << 16 | (c->n_cigar & 0xffff);
    x[4] = c->l_qseq;
    x[5] = c->mtid;
    x[6] = c->mpos;
    x[7] = c->isize;

    //这里空间一定是够的，ok要赋值为1！！！！
    ok = 1;
    //ok = (rabbit_bgzf_mul_flush_try(fp, bam_write_compress, write_block, 4 + block_len) >= 0);
    if (is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
        y = block_len;
        if (ok) ok = (rabbit_bgzf_mul_write(write_block, ed_swap_4p(&y), 4) >= 0);
        swap_data(c, b->l_data, b->data, 1);
    } else {
        if (ok) {
            ok = (rabbit_bgzf_mul_write(write_block, &block_len, 4) >= 0);
        }
    }
    if (ok) {
        ok = (rabbit_bgzf_mul_write(write_block, x, 32) >= 0);
    }
    if (ok) ok = (rabbit_bgzf_mul_write(write_block, b->data, c->l_qname - c->l_extranul) >= 0);
    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
        if (ok)
            ok = (rabbit_bgzf_mul_write(write_block, b->data + c->l_qname,
                                        b->l_data - c->l_qname) >= 0);
    } else { // with long CIGAR, insert a fake CIGAR record and move the real CIGAR to the CG:B,I tag
        uint8_t buf[8];
        uint32_t cigar_st, cigar_en, cigar[2];
        hts_pos_t cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
        if (cigreflen >= (1 << 28)) {
            // Length of reference covered is greater than the biggest
            // CIGAR operation currently allowed.
            hts_log_error("Record %s with %d CIGAR ops and ref length %"
            PRIhts_pos
            " cannot be written in BAM.  Try writing SAM or CRAM instead.\n",
                    bam_get_qname(b), c->n_cigar, cigreflen);
            return -1;
        }
        cigar_st = (uint8_t *) bam_get_cigar(b) - b->data;
        cigar_en = cigar_st + c->n_cigar * 4;
        cigar[0] = (uint32_t) c->l_qseq << 4 | BAM_CSOFT_CLIP;
        cigar[1] = (uint32_t) cigreflen << 4 | BAM_CREF_SKIP;
        u32_to_le(cigar[0], buf);
        u32_to_le(cigar[1], buf + 4);
        if (ok) {
            ok = (rabbit_bgzf_mul_write(write_block, buf, 8) >=
                  0); // write cigar: <read_length>S<ref_length>N
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write(write_block, &b->data[cigar_en],
                                        b->l_data - cigar_en) >= 0); // write data after CIGAR
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write(write_block, "CGBI", 4) >= 0); // write CG:B,I
        }
        u32_to_le(c->n_cigar, buf);
        if (ok) {
            ok = (rabbit_bgzf_mul_write(write_block, buf, 4) >=
                  0); // write the true CIGAR length
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write(write_block, &b->data[cigar_st], c->n_cigar * 4) >=
                  0); // write the real CIGAR
        }
    }
    if (is_be) swap_data(c, b->l_data, b->data, 0);
    return ok ? 4 + block_len : -1;
}

static inline int writeBam1_to_block_mpi_fast(bam_block *write_block, bam1_t *b) {
    const bam1_core_t *c = &b->core;
    int raw_l_qname = c->l_qname - c->l_extranul;
    uint32_t block_len;
    uint32_t x[8];
    uint8_t *dst;
    uint32_t rest_len;

    if (RB_UNLIKELY(raw_l_qname <= 0 || raw_l_qname > 255 ||
                    c->n_cigar > 0xffff ||
                    c->pos > INT_MAX ||
                    c->mpos > INT_MAX ||
                    c->isize < INT_MIN || c->isize > INT_MAX)) {
        bam_block *wb = write_block;
        return writeBam1_to_block(wb, b, 0);
    }

    block_len = b->l_data - c->l_extranul + 32;
    if (RB_UNLIKELY((uint64_t)write_block->pos + 4 + block_len > BGZF_BLOCK_SIZE))
        return -1;

    x[0] = c->tid;
    x[1] = c->pos;
    x[2] = (uint32_t)c->bin << 16 | c->qual << 8 | (uint32_t)raw_l_qname;
    x[3] = (uint32_t)c->flag << 16 | (c->n_cigar & 0xffff);
    x[4] = c->l_qseq;
    x[5] = c->mtid;
    x[6] = c->mpos;
    x[7] = c->isize;

    dst = (uint8_t *)write_block->data + write_block->pos;
    memcpy(dst, &block_len, 4);
    memcpy(dst + 4, x, 32);
    memcpy(dst + 36, b->data, (size_t)raw_l_qname);
    rest_len = (uint32_t)(b->l_data - c->l_qname);
    memcpy(dst + 36 + raw_l_qname, b->data + c->l_qname, rest_len);
    write_block->pos += 4 + block_len;
    return 4 + block_len;
}

extern "C" void slave_mpi_compressfunc(Comp_Para paras[64]) {
    int id = _PEN;
    Comp_Para* para = &paras[id];

    if (para->status != 0 || para->input_records == nullptr || para->n_records == 0) return;
    int compress_level = para->compress_level;
    if (compress_level != 0 && compress_level != 1 && compress_level != 6) {
        para->status = -2;
        return;
    }

    bam_block* uncompressed = para->un_comp_block;
    bam_block* compressed = para->output_block;
    uncompressed->pos = 0;
    uncompressed->length = 0;
    uncompressed->errcode = 0;
    uncompressed->block_id = 0;
    uncompressed->block_address = 0;
    para->compress_serialize_cycles = 0;
    para->compress_alloc_cycles = 0;
    para->compress_deflate_cycles = 0;
    para->compress_footer_cycles = 0;
    para->compress_total_cycles = 0;
    para->compress_level = compress_level;

    unsigned long total_t0 = slave_cycle_now();
    unsigned long serialize_t0 = slave_cycle_now();
    for (int i = 0; i < para->n_records; i++) {
        bam1_t* b = para->input_records[i];
        if (RB_UNLIKELY(writeBam1_to_block_mpi_fast(uncompressed, b) < 0)) {
            para->status = -2;
            para->compress_serialize_cycles = (uint64_t)(slave_cycle_now() - serialize_t0);
            para->compress_total_cycles = (uint64_t)(slave_cycle_now() - total_t0);
            return;
        }
    }
    para->compress_serialize_cycles = (uint64_t)(slave_cycle_now() - serialize_t0);

    if (swbam_cpe_compress_bam_block(
            uncompressed, uncompressed->pos, compressed,
            compress_level, id, 1, &para->compress_alloc_cycles,
            &para->compress_deflate_cycles,
            &para->compress_footer_cycles) != 0) {
        compressed->length = -1;
        para->status = -2;
        para->compress_total_cycles = (uint64_t)(slave_cycle_now() - total_t0);
        return;
    }
    para->output_size = compressed->length;
    para->compress_total_cycles = (uint64_t)(slave_cycle_now() - total_t0);
    para->status = 0;
}
