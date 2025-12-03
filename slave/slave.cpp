
 #include <stdint.h>
 #include <cstdlib>
// #include <mutex>
//#include <libdeflate.h>
// #include <thread>
// #include <atomic>
// #include <condition_variable>
 #include <vector>
// #include <omp.h>
 #include <cmath>
 #include <zlib.h>
 #include <stdio.h>
 #include <stdlib.h>

 #include <sunway/BamTools.h>


 //----------------------------------------------------------------------------------------------
 //bam to sam
// 进行标准的解压缩
int bgzf_uncompress(uint8_t *dst, size_t *dlen,
                    const uint8_t *src, size_t slen,
                    uint32_t expected_crc) {

    struct libdeflate_decompressor *z = libdeflate_alloc_decompressor();
    if (!z) {
        hts_log_error("Call to libdeflate_alloc_decompressor failed");
        return -1;
    }

    int ret = libdeflate_deflate_decompress(z, src, slen, dst, *dlen, dlen);
    libdeflate_free_decompressor(z);

    if (ret != 0) {
        hts_log_error("Inflate operation failed: %d", ret);
        return -1;
    }

    uint32_t crc = libdeflate_crc32(0, (unsigned char *) dst, *dlen);
    if (crc != expected_crc) {
        hts_log_error("CRC32 checksum mismatch");
        return -2;
    }

    return 0;
}

// 进行一个BGZF块的解压缩,处理一个 BGZF 块的头尾
int block_decode_func(struct bam_block *comp, struct bam_block *un_comp) {
    un_comp->pos = 0;
    un_comp->length = BGZF_MAX_BLOCK_SIZE;
    uint32_t crc = le_to_u32((uint8_t *) comp->data + comp->length - 8);                             //le_to_u32  hts_endian.h
    int ret = bgzf_uncompress(un_comp->data, reinterpret_cast<size_t *>(&un_comp->length),
                              comp->data + 18, comp->length - 18, crc);
    if (ret != 0) un_comp->errcode |= BGZF_ERR_ZLIB;
    return ret;
}

// 从解压缩后的块中解析出一个 bam1_t 记录
int read_bam(struct bam_block *fq, bam1_t *b, int is_be) {
    if (fq->pos >= fq->length) return -1;
    bam1_core_t *c = &b->core;
    int32_t block_len, ret, i;
    uint32_t x[8], new_l_data;

    b->l_data = 0;

    if ((ret = Rabbit_bgzf_read(fq, &block_len, 4)) != 4) { //读取四个字节，转换成int
        if (ret == 0) return -1; // normal end-of-file
        else return -2; // truncated
    }

    if (is_be) ed_swap_4p(&block_len);
    if (block_len < 32) return -4;  // block_len includes core data
    if (Rabbit_bgzf_read(fq, x, 32) != 32) return -3; //读取32个字节
    if (is_be) { for (i = 0; i < 8; ++i) ed_swap_4p(x + i); }
    c->tid = x[0];
    c->pos = (int32_t) x[1];
    c->bin = x[2] >> 16;
    c->qual = x[2] >> 8 & 0xff;
    c->l_qname = x[2] & 0xff;
    c->l_extranul = (c->l_qname % 4 != 0) ? (4 - c->l_qname % 4) : 0;
    c->flag = x[3] >> 16;
    c->n_cigar = x[3] & 0xffff;
    c->l_qseq = x[4];
    c->mtid = x[5];
    c->mpos = (int32_t) x[6];
    c->isize = (int32_t) x[7];

    new_l_data = block_len - 32 + c->l_extranul;//block_len + c->l_extranul
    if (new_l_data > INT_MAX || c->l_qseq < 0 || c->l_qname < 1) return -4;
    if (((uint64_t) c->n_cigar << 2) + c->l_qname + c->l_extranul
        + (((uint64_t) c->l_qseq + 1) >> 1) + c->l_qseq > (uint64_t) new_l_data)
        return -4;
    //在从核内存上开辟空间！！！！！！
    if (realloc_bam_data(b, new_l_data) < 0) return -4;
    b->l_data = new_l_data;

    if (Rabbit_bgzf_read(fq, b->data, c->l_qname) != c->l_qname) return -4;
    if (b->data[c->l_qname - 1] != '\0') { // Try to fix missing NUL termination
        if (fixup_missing_qname_nul(b) < 0) return -4;
    }
    for (i = 0; i < c->l_extranul; ++i) b->data[c->l_qname + i] = '\0';
    c->l_qname += c->l_extranul;
    if (b->l_data < c->l_qname ||
        Rabbit_bgzf_read(fq, b->data + c->l_qname, b->l_data - c->l_qname) != b->l_data - c->l_qname)
        return -4;
    if (is_be) swap_data(c, b->l_data, b->data, 0);
    if (bam_tag2cigar(b, 0, 0) < 0) return -4;
    if (c->n_cigar > 0) { // recompute "bin" and check CIGAR-qlen consistency
        hts_pos_t rlen, qlen;
        bam_cigar2rqlens(c->n_cigar, bam_get_cigar(b), &rlen, &qlen);
        if ((b->core.flag & BAM_FUNMAP) || rlen == 0) rlen = 1;
        b->core.bin = hts_reg2bin(b->core.pos, b->core.pos + rlen, 14, 5);
        // Sanity check for broken CIGAR alignments
        if (c->l_qseq > 0 && !(c->flag & BAM_FUNMAP) && qlen != c->l_qseq) {
            hts_log_error("CIGAR and query sequence lengths differ for %s",
                          bam_get_qname(b));
            return -4;
        }
    }

    return 4 + block_len;
}

int Rabbit_bgzf_read(struct bam_block *fq, void *data, unsigned int length) {
    if (length <= 0) return -1;
    if (length > fq->length - fq->pos) printf("One Block Is Small\n");
    length = fq->pos + length > fq->length ? fq->length - fq->pos : length;
    memcpy((uint8_t *) data, fq->data + fq->pos, length);
    fq->pos += length;
    return length;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------
extern "C" void decompressfunc(Para paras[64]) {

    int id = _PEN;             // 从核号（0~63）
    Para* para = &paras[id];

    // 1. 读取压缩块数据
    bam_block* comp = para->input_block;
    bam_block* un_comp = para->un_comp_block;
    if (comp == NULL) {
        return;
    }

    // 2. 解压缩该块
    block_decode_func(comp,un_comp);

    // 3. 解析解压后的SAM/BAM记录
    int count = 0;
    bam1_t* b = NULL;
    b = para->output_records[count]; 
    while(read_bam(un_comp, b, 0)>=0){
        para->l_data_list[count] = b->l_data;
        para->data_list[count] = b->data;
        count++;

        b = para->output_records[count];
    }

    //bam_destroy1(b);

    para->n_records = count;
    para->status = 0;


}

extern "C" void copyfunc(Para paras[64]) {

    int id = _PEN;             // 从核号（0~63）
    Para* para = &paras[id];

    if(para->input_block == NULL) {
        return;
    }

    for (int i = 0; i < para->n_records; i++) {
        // 复制数据到新分配的data区域
        memcpy(para->output_records[i]->data, para->data_list[i], para->l_data_list[i]);
    }

    para->status = 1;  // copy 完成标记

}


//--------------------------------------------------------------------------------------------------------------
//sam to bam

//数据的实际写入函数
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

    //这里空间一定是够的
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

int rabbit_bgzf_compress(void *_dst, size_t *dlen, const void *src, size_t slen, int level) {
    if (slen == 0) {
        // EOF block
        if (*dlen < 28) return -1;
        memcpy(_dst, "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", 28);
        *dlen = 28;
        return 0;
    }

    uint8_t *dst = (uint8_t *) _dst;

    if (level == 0) {
        // Uncompressed data
//        printf("Compress 1\n");
        if (*dlen < slen + 5 + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH) return -1;
//        printf("Compress 2\n");
        dst[BLOCK_HEADER_LENGTH] = 1; // BFINAL=1, BTYPE=00; see RFC1951
//        printf("Compress 3\n");
        u16_to_le(slen, &dst[BLOCK_HEADER_LENGTH + 1]); // length
//        printf("Compress 4\n");
        u16_to_le(~slen, &dst[BLOCK_HEADER_LENGTH + 3]); // ones-complement length
//        printf("Compress 5\n");
        memcpy(dst + BLOCK_HEADER_LENGTH + 5, src, slen);
//        printf("Compress 6\n");
        *dlen = slen + 5 + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
//        printf("Compress 7\n");
    } else {
        level = level > 0 ? level : 6; // libdeflate doesn't honour -1 as default
        // NB levels go up to 12 here.
        struct libdeflate_compressor *z = libdeflate_alloc_compressor(level);
        if (!z) return -1;

        // Raw deflate
        size_t clen =
                libdeflate_deflate_compress(z, src, slen,
                                            dst + BLOCK_HEADER_LENGTH,
                                            *dlen - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH);

        if (clen <= 0) {
            hts_log_error("Call to libdeflate_deflate_compress failed");
            libdeflate_free_compressor(z);
            return -1;
        }

        *dlen = clen + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;

        libdeflate_free_compressor(z);
    }

    // write the header
    memcpy(dst, g_magic, BLOCK_HEADER_LENGTH); // the last two bytes are a place holder for the length of the block
    packInt16(&dst[16], *dlen - 1); // write the compressed length; -1 to fit 2 bytes

    // write the footer
    uint32_t crc = libdeflate_crc32(0, src, slen);
    packInt32((uint8_t * ) & dst[*dlen - 8], crc);
    packInt32((uint8_t * ) & dst[*dlen - 4], slen);
    return 0;
}

int block_encode_func(bam_block *un_comp, bam_block *comp , int compress_level) {
    //int rabbit_write_deflate_block(BGZF *fp, bam_write_block *write_block) 
    size_t comp_size = BGZF_MAX_BLOCK_SIZE;
    int ret;
    ret = rabbit_bgzf_compress(comp->data, &comp_size, un_comp->data, un_comp->pos,
                                     compress_level);
    // if (!fp->is_gzip)
    //     ret = rabbit_bgzf_compress(write_block->compressed_data, &comp_size, write_block->uncompressed_data,
    //                                write_block->block_offset, fp->compress_level);
    // else
    //     ret = rabbit_bgzf_gzip_compress(fp, write_block->compressed_data, &comp_size, write_block->uncompressed_data,
    //                                     write_block->block_offset, fp->compress_level);

    if (ret != 0) {
        hts_log_debug("Compression error %d", ret);
        //fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    return comp_size;
}


extern "C" void slave_compressfunc(Comp_Para paras[64]) {

    int id = _PEN;             // 从核号（0~63）
    int compress_level = 6; // 默认压缩等级
    Comp_Para* para = &paras[id];

    // 忽略空任务
    if (para->status != 0 || para->input_records == nullptr || para->n_records == 0) return;

    bam_block* uncompressed = para->un_comp_block;
    bam_block* compressed = para->output_block;
    uncompressed->pos = 0;     //压缩时看这个实际的数据长度
    uncompressed->length = 0;  //这里好像不用管


    //1. 将 bam1_t 记录 逐个解析并写入 到 uncompressed 中
    for (int i = 0; i < para->n_records; i++) {
        bam1_t* b = para->input_records[i];
        writeBam1_to_block(uncompressed, b , 0);
    }

    //2. 对 uncompressed 进行压缩
    compressed->length = block_encode_func(uncompressed, compressed , compress_level); //TODO para->fp

    para->output_size = compressed->length;
    para->status = 0; // success

}



//-------------------------------------------------------------------------------------------------------------------------------------------
//BamTools-slave



