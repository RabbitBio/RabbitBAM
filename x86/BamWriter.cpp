#include "BamWriter.h"

int rabbit_write_deflate_block(BGZF *fp, bam_write_block *write_block) {
    size_t comp_size = BGZF_MAX_BLOCK_SIZE;
    int ret;
    if (!fp->is_gzip)
        ret = rabbit_bgzf_compress(write_block->compressed_data, &comp_size, write_block->uncompressed_data,
                                   write_block->block_offset, fp->compress_level);
    else
        ret = rabbit_bgzf_gzip_compress(fp, write_block->compressed_data, &comp_size, write_block->uncompressed_data,
                                        write_block->block_offset, fp->compress_level);

    if (ret != 0) {
        hts_log_debug("Compression error %d", ret);
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    return comp_size;
}

int rabbit_bgzf_flush(BGZF *fp, bam_write_block *write_block) {
    while (write_block->block_offset > 0) {
        int block_length;
        printf("Write Block Offset : %d\n", write_block->block_offset);
        block_length = rabbit_write_deflate_block(fp, write_block);
        if (block_length < 0) {
            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block_length, NULL));
            return -1;
        }
        if (write_block == nullptr) {
            printf("Write Game Over!!!\n");
        }
        if (hwrite(fp->fp, write_block->compressed_data, block_length) != block_length) {
            printf("Write Failed\n");
            hts_log_error("File write failed (wrong size)");
            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
            return -1;
        }


        write_block->block_offset = 0;
        fp->block_address += block_length;
    }
    write_block->block_offset = 0;
    return 0;
}

int rabbit_bgzf_mul_flush(BGZF *fp, BamWriteCompress *bam_write_compress, bam_write_block *&write_block) {

    bam_write_compress->inputUnCompressData(write_block);
    write_block = bam_write_compress->getEmpty();
    return 0;
}

int rabbit_bgzf_write(BGZF *fp, bam_write_block *&write_block, const void *data, size_t length) {
    const uint8_t *input = (const uint8_t *) data;
    ssize_t remaining = length;
    while (remaining > 0) {
        uint8_t *buffer = (uint8_t *) write_block->uncompressed_data;
        int copy_length = BGZF_BLOCK_SIZE - write_block->block_offset;
        if (copy_length > remaining) copy_length = remaining;
        memcpy(buffer + write_block->block_offset, input, copy_length);
        write_block->block_offset += copy_length;
        input += copy_length;
        remaining -= copy_length;
        if (write_block->block_offset == BGZF_BLOCK_SIZE) {
            if (rabbit_bgzf_flush(fp, write_block) != 0) return -1;
        }
    }
    return length - remaining;
}

int rabbit_bgzf_mul_write_fast(BGZF *fp, BamWriteCompress *bam_write_compress, bam_write_block *&write_block,
                               const void *data, size_t length) {
    const uint8_t *input = (const uint8_t *) data;
    ssize_t remaining = length;
    while (remaining > 0) {
        uint8_t *buffer = (uint8_t *) write_block->uncompressed_data;
        int copy_length = BGZF_BLOCK_SIZE - write_block->block_offset;
        if (copy_length > remaining) copy_length = remaining;
        memcpy(buffer + write_block->block_offset, input, copy_length);
        write_block->block_offset += copy_length;
        input += copy_length;
        remaining -= copy_length;
        //if (write_block->block_offset == BGZF_BLOCK_SIZE) {
        //    if (rabbit_bgzf_mul_flush(fp, bam_write_compress, write_block) != 0) return -1;
        //}
    }
    return length - remaining;
}


int
rabbit_bgzf_mul_write(BGZF *fp, BamWriteCompress *bam_write_compress, bam_write_block *&write_block, const void *data,
                      size_t length) {
    const uint8_t *input = (const uint8_t *) data;
    ssize_t remaining = length;
    while (remaining > 0) {
        uint8_t *buffer = (uint8_t *) write_block->uncompressed_data;
        int copy_length = BGZF_BLOCK_SIZE - write_block->block_offset;
        if (copy_length > remaining) copy_length = remaining;
        memcpy(buffer + write_block->block_offset, input, copy_length);
        write_block->block_offset += copy_length;
        input += copy_length;
        remaining -= copy_length;
        if (write_block->block_offset == BGZF_BLOCK_SIZE) {
            if (rabbit_bgzf_mul_flush(fp, bam_write_compress, write_block) != 0) return -1;
        }
    }
    return length - remaining;
}

int rabbit_bgzf_flush_try(BGZF *fp, bam_write_block *write_block, ssize_t size) {
    if (write_block->block_offset + size > BGZF_BLOCK_SIZE) {
        return rabbit_bgzf_flush(fp, write_block);
    }
    return 0;
}

int
rabbit_bgzf_mul_flush_try(BGZF *fp, BamWriteCompress *bam_write_compress, bam_write_block *&write_block, ssize_t size) {
    if (write_block->block_offset + size > BGZF_BLOCK_SIZE) {
        return rabbit_bgzf_mul_flush(fp, bam_write_compress, write_block);
    }
    return 0;
}

int bam_write_pack(BGZF *fp, BamWriteCompress *bam_write_compress) {
    bam_write_block *block;
    while (1) {

        block = bam_write_compress->getCompressData();
        if (block == nullptr) {
            return 0;
        }
        int block_offset = block->block_offset;
        if (block->block_length < 0) {
            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block->block_length, NULL));
            return -1;
        }
        if (hwrite(fp->fp, block->compressed_data, block->block_length) != block->block_length) {
            hts_log_error("File write failed (wrong size)");
            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
            return -1;
        }
        block->block_offset = 0;
        fp->block_address += block->block_length;
        bam_write_compress->backEmpty(block);
    }

}

void bam_write_compress_pack(BGZF *fp, BamWriteCompress *bam_write_compress) {
    bam_write_block *block;
    while (1) {
        block = bam_write_compress->getUnCompressData();
        if (block == nullptr) {
            break;
        }
        block->block_length = rabbit_write_deflate_block(fp, block);
        bam_write_compress->inputCompressData(block);
    }
    bam_write_compress->CompressThreadComplete();
}

int rabbit_bam_write_test(BGZF *fp, bam_write_block *write_block, bam1_t *b) {
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
    ok = (rabbit_bgzf_flush_try(fp, write_block, 4 + block_len) >= 0);
    if (fp->is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
        y = block_len;
        if (ok) ok = (rabbit_bgzf_write(fp, write_block, ed_swap_4p(&y), 4) >= 0);
        swap_data(c, b->l_data, b->data, 1);
    } else {
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, &block_len, 4) >= 0);
        }
    }
    if (ok) {
        ok = (rabbit_bgzf_write(fp, write_block, x, 32) >= 0);
    }
    if (ok) ok = (rabbit_bgzf_write(fp, write_block, b->data, c->l_qname - c->l_extranul) >= 0);
    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
        if (ok) ok = (rabbit_bgzf_write(fp, write_block, b->data + c->l_qname, b->l_data - c->l_qname) >= 0);
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
            ok = (rabbit_bgzf_write(fp, write_block, buf, 8) >= 0); // write cigar: <read_length>S<ref_length>N
        }
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, &b->data[cigar_en], b->l_data - cigar_en) >=
                  0); // write data after CIGAR
        }
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, "CGBI", 4) >= 0); // write CG:B,I
        }
        u32_to_le(c->n_cigar, buf);
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, buf, 4) >= 0); // write the true CIGAR length
        }
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, &b->data[cigar_st], c->n_cigar * 4) >= 0); // write the real CIGAR
        }
    }
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    return ok ? 4 + block_len : -1;
}

int
rabbit_bam_write_mul_test(BGZF *fp, BamWriteCompress *bam_write_compress, bam_write_block *&write_block, bam1_t *b) {
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
    ok = (rabbit_bgzf_mul_flush_try(fp, bam_write_compress, write_block, 4 + block_len) >= 0);
    if (fp->is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
        y = block_len;
        if (ok) ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, ed_swap_4p(&y), 4) >= 0);
        swap_data(c, b->l_data, b->data, 1);
    } else {
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, &block_len, 4) >= 0);
        }
    }
    if (ok) {
        ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, x, 32) >= 0);
    }
    if (ok) ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, b->data, c->l_qname - c->l_extranul) >= 0);
    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
        if (ok)
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, b->data + c->l_qname,
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
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, buf, 8) >=
                  0); // write cigar: <read_length>S<ref_length>N
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, &b->data[cigar_en],
                                        b->l_data - cigar_en) >= 0); // write data after CIGAR
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, "CGBI", 4) >= 0); // write CG:B,I
        }
        u32_to_le(c->n_cigar, buf);
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, buf, 4) >=
                  0); // write the true CIGAR length
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, &b->data[cigar_st], c->n_cigar * 4) >=
                  0); // write the real CIGAR
        }
    }
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    return ok ? 4 + block_len : -1;
}


void benchmark_write_pack(BamCompleteBlock *completeBlock, samFile *output, sam_hdr_t *hdr, int level) {

    uint8_t *compress_block_test = new uint8_t[BGZF_BLOCK_SIZE];
    uint8_t *uncompress_block_test = new uint8_t[BGZF_BLOCK_SIZE];
    output->fp.bgzf->block_offset = 0;
    output->fp.bgzf->uncompressed_block = uncompress_block_test;
    output->fp.bgzf->compressed_block = compress_block_test;
    output->fp.bgzf->compress_level = level;
    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    if (sam_hdr_write(output, hdr) != 0) {
        printf("HDR Write False!\n");
        return;
    }
    bam_complete_block *un_comp;
    long long ans = 0;
    long long res = 0;
    bam_write_block *write_block = new bam_write_block();
    write_block->block_offset = 0;
    write_block->uncompressed_data = new uint8_t[BGZF_BLOCK_SIZE];
    write_block->compressed_data = new uint8_t[BGZF_BLOCK_SIZE];
    write_block->status = 0;
    while (1) {
        un_comp = completeBlock->getCompleteBlock();
        if (un_comp == nullptr) {
            break;
        }
        int ret;
        while ((ret = (read_bam(un_comp, b, 0))) >= 0) {
            rabbit_bam_write_test(output->fp.bgzf, write_block, b);
            ans++;
        }
        res++;
        completeBlock->backEmpty(un_comp);
    }
    rabbit_bgzf_flush(output->fp.bgzf, write_block);
    printf("Bam1_t Number is %lld\n", ans);
    printf("Block  Number is %lld\n", res);
}

void benchmark_write_mul_pack(BamCompleteBlock *completeBlock, BamWriteCompress *bam_write_compress, samFile *output,
                              sam_hdr_t *hdr, int level) {

    output->fp.bgzf->block_offset = 0;
    output->fp.bgzf->compress_level = level;
    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    if (sam_hdr_write(output, hdr) != 0) {
        printf("HDR Write False!\n");
        return;
    }
    bam_complete_block *un_comp;
    long long ans = 0;
    long long res = 0;
    bam_write_block *write_block = bam_write_compress->getEmpty();
    int bam_num = 1;
    while (1) {
        un_comp = completeBlock->getCompleteBlock();
        if (un_comp == nullptr) {
            break;
        }
        int ret;
        while ((ret = (read_bam(un_comp, b, 0))) >= 0) {
            rabbit_bam_write_mul_test(output->fp.bgzf, bam_write_compress, write_block, b);
            ans++;
        }
        res++;
        completeBlock->backEmpty(un_comp);
    }
    if (write_block->block_offset > 0) {
        bam_write_compress->inputUnCompressData(write_block);
    }
    bam_write_compress->WriteComplete();
}

int
rabbit_bam_write_mul_parallel(BGZF *fp, BamWriteCompress *bam_write_compress, bam_write_block *&write_block, bam1_t *b,
                              std::vector<bam_write_block *> &block_vec) {
    const bam1_core_t *c = &b->core;
    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
    int i, ok = 1;
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
    if (write_block->block_offset + 4 + block_len >= BGZF_BLOCK_SIZE) {
        block_vec.push_back(write_block);
        write_block = bam_write_compress->getEmpty();
    }
    //ok = (rabbit_bgzf_mul_flush_try(fp, bam_write_compress, write_block, 4 + block_len) >= 0);
    if (fp->is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
        y = block_len;
        if (ok) ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, ed_swap_4p(&y), 4) >= 0);
        swap_data(c, b->l_data, b->data, 1);
    } else {
        if (ok) {
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, &block_len, 4) >= 0);
        }
    }
    if (ok) {
        ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, x, 32) >= 0);
    }
    if (ok)
        ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, b->data, c->l_qname - c->l_extranul) >=
              0);
    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
        if (ok)
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, b->data + c->l_qname,
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
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, buf, 8) >=
                  0); // write cigar: <read_length>S<ref_length>N
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, &b->data[cigar_en],
                                             b->l_data - cigar_en) >= 0); // write data after CIGAR
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, "CGBI", 4) >= 0); // write CG:B,I
        }
        u32_to_le(c->n_cigar, buf);
        if (ok) {
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, buf, 4) >=
                  0); // write the true CIGAR length
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, &b->data[cigar_st], c->n_cigar * 4) >=
                  0); // write the real CIGAR
        }
    }
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    return ok ? 4 + block_len : -1;
}

void BamWriter::write_parallel(std::vector<bam1_t *> b_vec) {

#define THREAD_NUM_PW 4
    int vec_size = b_vec.size();
    std::vector < bam_write_block * > todo_push_blocks[THREAD_NUM_PW];
    for (int i = 0; i < THREAD_NUM_PW; i++) {
        blocks[i] = bam_write_compress->getEmpty();
    }

#pragma omp parallel for num_threads(THREAD_NUM_PW) schedule(static)
    for (int i = 0; i < vec_size; i++) {
        int tid = omp_get_thread_num();
        rabbit_bam_write_mul_parallel(output->fp.bgzf, bam_write_compress, blocks[tid], b_vec[i],
                                      todo_push_blocks[tid]);
    }

    for (int i = 0; i < THREAD_NUM_PW; i++) {
        for (auto item: todo_push_blocks[i]) {
            bam_write_compress->inputUnCompressData(item);
        }
        bam_write_compress->inputUnCompressData(blocks[i]);
    }


}


void BamWriter::write(bam1_t *b) {
    rabbit_bam_write_mul_test(output->fp.bgzf, bam_write_compress, write_block, b);

}

void BamWriter::bam_write(bam1_t *b) {
    rabbit_bam_write_mul_test(output->fp.bgzf, bam_write_compress, write_block, b);
}


BamWriter::BamWriter(int threadNumber, int level, int BufferSize) {

    n_thread_write = threadNumber;
    bam_write_compress = new BamWriteCompress(BufferSize, n_thread_write);


}


BamWriter::BamWriter(std::string filename, int threadNumber, int level, int BufferSize, bool is_tgs) {

    if ((output = sam_open(filename.c_str(), "wb")) == NULL) {
        printf("Can`t open this file!\n");
    }

    n_thread_write = threadNumber;
    bam_write_compress = new BamWriteCompress(BufferSize, n_thread_write);


    write_compress_thread = new std::thread *[n_thread_write];
    for (int i = 0; i < n_thread_write; i++)
        write_compress_thread[i] = new std::thread(&bam_write_compress_pack, output->fp.bgzf, bam_write_compress);

    write_output_thread = new std::thread(&bam_write_pack, output->fp.bgzf, bam_write_compress);


    output->fp.bgzf->block_offset = 0;
    output->fp.bgzf->compress_level = level;

    if(is_tgs) {
        write_block=bam_write_compress->getEmpty();
    }

//#ifdef use_parallel_write
//#else
//    write_block=bam_write_compress->getEmpty();
//#endif
    //for(int i = 0; i < THREAD_NUM_PW; i++) {
    //    blocks[i] = bam_write_compress->getEmpty();
    //}

}

BamWriter::BamWriter(std::string filename, sam_hdr_t *hdr, int threadNumber, int level, int BufferSize, bool is_tgs) {

    if ((output = sam_open(filename.c_str(), "wb")) == NULL) {
        printf("Can`t open this file!\n");
    }
    if (sam_hdr_write(output, hdr) != 0) {
        printf("HDR Write False!\n");
    }

    n_thread_write = threadNumber;
    bam_write_compress = new BamWriteCompress(BufferSize, n_thread_write);


    write_compress_thread = new std::thread *[n_thread_write];
    for (int i = 0; i < n_thread_write; i++)
        write_compress_thread[i] = new std::thread(&bam_write_compress_pack, output->fp.bgzf, bam_write_compress);

    write_output_thread = new std::thread(&bam_write_pack, output->fp.bgzf, bam_write_compress);


    output->fp.bgzf->block_offset = 0;
    output->fp.bgzf->compress_level = level;


    if(is_tgs) {
        write_block=bam_write_compress->getEmpty();
    }

//#ifdef use_parallel_write
//#else
//    write_block=bam_write_compress->getEmpty();
//#endif
    //for(int i = 0; i < THREAD_NUM_PW; i++) {
    //    blocks[i] = bam_write_compress->getEmpty();
    //}


}

void BamWriter::hdr_write(sam_hdr_t *hdr) {
    if (sam_hdr_write(output, hdr) != 0) {
        printf("HDR Write False!\n");
    }
}


void BamWriter::set_output(samFile *output, bool is_tgs) {
    this->output = output;
    write_compress_thread = new std::thread *[n_thread_write];
    for (int i = 0; i < n_thread_write; i++)
        write_compress_thread[i] = new std::thread(&bam_write_compress_pack, output->fp.bgzf, bam_write_compress);

    write_output_thread = new std::thread(&bam_write_pack, output->fp.bgzf, bam_write_compress);


    output->fp.bgzf->block_offset = 0;
    output->fp.bgzf->compress_level = 6;

    if(is_tgs) {
        write_block=bam_write_compress->getEmpty();
    }

//#ifdef use_parallel_write
//#else
//    write_block=bam_write_compress->getEmpty();
//#endif

    //for(int i = 0; i < THREAD_NUM_PW; i++) {
    //    blocks[i] = bam_write_compress->getEmpty();
    //}


}

void BamWriter::over_parallel() {

    //for(int i = 0; i < THREAD_NUM_PW; i++) {
    //    if (blocks[i]->block_offset > 0) {
    //        bam_write_compress->inputUnCompressData(blocks[i]);
    //    }
    //}
    bam_write_compress->WriteComplete();
    for (int i = 0; i < n_thread_write; i++) write_compress_thread[i]->join();
    write_output_thread->join();
    int ret = sam_close(output);
    if (ret < 0) {
        fprintf(stderr, "Error closing output.\n");
        //exit_code = EXIT_FAILURE;
    }

}


void BamWriter::over() {

    if (write_block->block_offset > 0) {
        bam_write_compress->inputUnCompressData(write_block);
        write_block = bam_write_compress->getEmpty();
    }
    bam_write_compress->WriteComplete();
    for (int i = 0; i < n_thread_write; i++) write_compress_thread[i]->join();
    write_output_thread->join();
    int ret = sam_close(output);
    if (ret < 0) {
        fprintf(stderr, "Error closing output.\n");
        //exit_code = EXIT_FAILURE;
    }
}

