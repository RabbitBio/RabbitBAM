//
// Created by 赵展 on 2022/1/21.
//

#include "BamTools.h"

int sam_realloc_bam_data(bam1_t *b, size_t desired) {
    uint32_t new_m_data;
    uint8_t *new_data;
    new_m_data = desired;
    kroundup32(new_m_data);
    if (new_m_data < desired) {
        errno = ENOMEM; // Not strictly true but we can't store the size
        return -1;
    }
    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) == 0) {
        new_data = (uint8_t *) (realloc(b->data, new_m_data));
    } else {
        if ((new_data = (uint8_t *) (malloc(new_m_data))) != NULL) {
            if (b->l_data > 0)
                memcpy(new_data, b->data,
                       b->l_data < b->m_data ? b->l_data : b->m_data);
            bam_set_mempolicy(b, bam_get_mempolicy(b) & (~BAM_USER_OWNS_DATA));
        }
    }
    if (!new_data) return -1;
    b->data = new_data;
    b->m_data = new_m_data;
    return 0;
}

inline int realloc_bam_data(bam1_t *b, size_t desired) {
    if (desired <= b->m_data) return 0;
    return sam_realloc_bam_data(b, desired);
}

inline int possibly_expand_bam_data(bam1_t *b, size_t bytes) {
    size_t new_len = (size_t) b->l_data + bytes;

    if (new_len > INT32_MAX || new_len < bytes) { // Too big or overflow
        errno = ENOMEM;
        return -1;
    }
    if (new_len <= b->m_data) return 0;
    return sam_realloc_bam_data(b, new_len);
}

int bam_tag2cigar(bam1_t *b, int recal_bin,
                  int give_warning) // return 0 if CIGAR is untouched; 1 if CIGAR is updated with CG
{
    bam1_core_t *c = &b->core;
    uint32_t cigar_st, n_cigar4, CG_st, CG_en, ori_len = b->l_data, *cigar0, CG_len, fake_bytes;
    uint8_t *CG;

    // test where there is a real CIGAR in the CG tag to move
    if (c->n_cigar == 0 || c->tid < 0 || c->pos < 0) return 0;
    cigar0 = bam_get_cigar(b);
    if (bam_cigar_op(cigar0[0]) != BAM_CSOFT_CLIP || bam_cigar_oplen(cigar0[0]) != c->l_qseq) return 0;
    fake_bytes = c->n_cigar * 4;
    int saved_errno = errno;
    CG = bam_aux_get(b, "CG");
    if (!CG) {
        if (errno != ENOENT) return -1;  // Bad aux data
        errno = saved_errno; // restore errno on expected no-CG-tag case
        return 0;
    }
    if (CG[0] != 'B' || CG[1] != 'I') return 0; // not of type B,I
    CG_len = le_to_u32(CG + 2);
    if (CG_len < c->n_cigar || CG_len >= 1U << 29)
        return 0; // don't move if the real CIGAR length is shorter than the fake cigar length

    // move from the CG tag to the right position
    cigar_st = (uint8_t *) cigar0 - b->data;
    c->n_cigar = CG_len;
    n_cigar4 = c->n_cigar * 4;
    CG_st = CG - b->data - 2;
    CG_en = CG_st + 8 + n_cigar4;
    if (possibly_expand_bam_data(b, n_cigar4 - fake_bytes) < 0) return -1;
    b->l_data =
            b->l_data - fake_bytes + n_cigar4; // we need c->n_cigar-fake_bytes bytes to swap CIGAR to the right place
    memmove(b->data + cigar_st + n_cigar4, b->data + cigar_st + fake_bytes,
            ori_len - (cigar_st + fake_bytes)); // insert c->n_cigar-fake_bytes empty space to make room
    memcpy(b->data + cigar_st, b->data + (n_cigar4 - fake_bytes) + CG_st + 8,
           n_cigar4); // copy the real CIGAR to the right place; -fake_bytes for the fake CIGAR
    if (ori_len > CG_en) // move data after the CG tag
        memmove(b->data + CG_st + n_cigar4 - fake_bytes, b->data + CG_en + n_cigar4 - fake_bytes, ori_len - CG_en);
    b->l_data -= n_cigar4 + 8; // 8: CGBI (4 bytes) and CGBI length (4)
    if (recal_bin)
        b->core.bin = hts_reg2bin(b->core.pos, bam_endpos(b), 14, 5);
    if (give_warning)
        hts_log_error("%s encodes a CIGAR with %d operators at the CG tag", bam_get_qname(b), c->n_cigar);
    return 1;
}

void bam_cigar2rqlens(int n_cigar, const uint32_t *cigar,
                      hts_pos_t *rlen, hts_pos_t *qlen) {
    int k;
    *rlen = *qlen = 0;
    for (k = 0; k < n_cigar; ++k) {
        int type = bam_cigar_type(bam_cigar_op(cigar[k]));
        int len = bam_cigar_oplen(cigar[k]);
        if (type & 1) *qlen += len;
        if (type & 2) *rlen += len;
    }
}

void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host) {
    uint32_t *cigar = (uint32_t *) (data + c->l_qname);
    uint32_t i;
    for (i = 0; i < c->n_cigar; ++i) ed_swap_4p(&cigar[i]);
}

inline int unpackInt16(const uint8_t *buffer) {
    return buffer[0] | buffer[1] << 8;
}

int load_block_from_cache(BGZF *fp, int64_t block_address) {
    khint_t k;
    cache_t *p;

    khash_t(cache) *h = fp->cache->h;
    k = kh_get(cache, h, block_address);
    if (k == kh_end(h)) return 0;
    p = &kh_val(h, k);
    if (fp->block_length != 0) fp->block_offset = 0;
    fp->block_address = block_address;
    fp->block_length = p->size;
    memcpy(fp->uncompressed_block, p->block, p->size);
    if (hseek(fp->fp, p->end_offset, SEEK_SET) < 0) {
        // todo: move the error up
        hts_log_error("Could not hseek to %" PRId64, p->end_offset);
        exit(1);
    }
    return p->size;
}

int check_header(const uint8_t *header) {
    if (header[0] != 31 || header[1] != 139 || header[2] != 8) return -2;
    return ((header[3] & 4) != 0
            && unpackInt16((uint8_t *) &header[10]) == 6
            && header[12] == 'B' && header[13] == 'C'
            && unpackInt16((uint8_t *) &header[14]) == 2) ? 0 : -1;
}

const char *bgzf_zerr(int errnum, z_stream *zs) {
    static char buffer[32];

    /* Return zs->msg if available.
       zlib doesn't set this very reliably.  Looking at the source suggests
       that it may get set to a useful message for deflateInit2, inflateInit2
       and inflate when it returns Z_DATA_ERROR. For inflate with other
       return codes, deflate, deflateEnd and inflateEnd it doesn't appear
       to be useful.  For the likely non-useful cases, the caller should
       pass NULL into zs. */

    if (zs && zs->msg) return zs->msg;

    // gzerror OF((gzFile file, int *errnum)
    switch (errnum) {
        case Z_ERRNO:
            return strerror(errno);
        case Z_STREAM_ERROR:
            return "invalid parameter/compression level, or inconsistent stream state";
        case Z_DATA_ERROR:
            return "invalid or incomplete IO";
        case Z_MEM_ERROR:
            return "out of memory";
        case Z_BUF_ERROR:
            return "progress temporarily not possible, or in() / out() returned an error";
        case Z_VERSION_ERROR:
            return "zlib version mismatch";
        case Z_NEED_DICT:
            return "data was compressed using a dictionary";
        case Z_OK: // 0: maybe gzgets error Z_NULL
        default:
            snprintf(buffer, sizeof(buffer), "[%d] unknown", errnum);
            return buffer;  // FIXME: Not thread-safe.
    }
}

int bgzf_uncompress(uint8_t *dst, size_t *dlen,
                    const uint8_t *src, size_t slen,
                    uint32_t expected_crc) {
//    printf("dst: %d   dlen: %d   src: %d   slen: %d   expected_crc: %d\n", dst[0], *dlen, src[0], slen, expected_crc);

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

int fixup_missing_qname_nul(bam1_t *b) {
    bam1_core_t *c = &b->core;

    // Note this is called before c->l_extranul is added to c->l_qname
    if (c->l_extranul > 0) {
        b->data[c->l_qname++] = '\0';
        c->l_extranul--;
    } else {
        if (b->l_data > INT_MAX - 4) return -1;
        if (realloc_bam_data(b, b->l_data + 4) < 0) return -1;
        b->l_data += 4;
        b->data[c->l_qname++] = '\0';
        c->l_extranul = 3;
    }
    return 0;
}

int block_decode_func(struct bam_block *comp, struct bam_block *un_comp) {
    un_comp->pos = 0;
    un_comp->length = BGZF_MAX_BLOCK_SIZE;
    uint32_t crc = le_to_u32((uint8_t *) comp->data + comp->length - 8);
    int ret = bgzf_uncompress(un_comp->data, reinterpret_cast<size_t *>(&un_comp->length),
                              comp->data + 18, comp->length - 18, crc);
    if (ret != 0) un_comp->errcode |= BGZF_ERR_ZLIB;
    return ret;
}

int read_block(BGZF *fp, struct bam_block *j) {
    uint8_t header[BLOCK_HEADER_LENGTH], *compressed_block;
    int count, size = 0, block_length, remaining;

    // NOTE: Guaranteed to be compressed as we block multi-threading in
    // uncompressed mode.  However it may be gzip compression instead
    // of bgzf.

    // Reading compressed file
    int64_t block_address;
    block_address = htell(fp->fp);

    if (fp->cache_size && load_block_from_cache(fp, block_address)) return 0;
    count = hpeek(fp->fp, header, sizeof(header));
    if (count == 0) // no data read
        return -1;
    int ret;
    if (count != sizeof(header) || (ret = check_header(header)) == -2) {
        j->errcode |= BGZF_ERR_HEADER;
        return -1;
    }
    if (ret == -1) {
        j->errcode |= BGZF_ERR_MT;
        return -1;
    }

    count = hread(fp->fp, header, sizeof(header));
    if (count != sizeof(header)) // no data read
        return -1;

    size = count;
    block_length = unpackInt16((uint8_t *) &header[16]) + 1; // +1 because when writing this number, we used "-1"
    if (block_length < BLOCK_HEADER_LENGTH) {
        j->errcode |= BGZF_ERR_HEADER;
        return -1;
    }

    compressed_block = (uint8_t *) j->data;
    memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
    remaining = block_length - BLOCK_HEADER_LENGTH;
    count = hread(fp->fp, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
    if (count != remaining) {
        j->errcode |= BGZF_ERR_IO;
        return -1;
    }
    size += count;
    j->length = block_length;
    j->block_address = block_address;
    j->errcode = 0;
    return 0;
}

int Rabbit_bgzf_read(struct bam_block *fq, void *data, unsigned int length) {
    if (length <= 0) return -1;
    if (length > fq->length - fq->pos) printf("One Block Is Small\n");
    length = fq->pos + length > fq->length ? fq->length - fq->pos : length;
    memcpy((uint8_t *) data, fq->data + fq->pos, length);
    fq->pos += length;
    return length;
}

void Rabbit_memcpy(void *target, unsigned char *source, unsigned int length) {
    memcpy((uint8_t *) target, source, length);
}

int Rabbit_bgzf_read(struct bam_complete_block *fq, void *data, unsigned int length) {
    if (length <= 0) return -1;
    if (length > fq->length - fq->pos) {
        printf("length is %d\n", length);
        printf("fq length is %d\n", fq->length);
        printf("fq pos is %d\n", fq->pos);
        printf("One Block Is Small\n");
    }
    length = fq->pos + length > fq->length ? fq->length - fq->pos : length;
    memcpy((uint8_t *) data, fq->data + fq->pos, length);
    fq->pos += length;
    return length;
}

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

int read_bam(struct bam_complete_block *fq, bam1_t *b, int is_be) {
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
    /*
     * 未处理这一块
     */
    new_l_data = block_len - 32 + c->l_extranul;//block_len + c->l_extranul
    if (new_l_data > INT_MAX || c->l_qseq < 0 || c->l_qname < 1) {
//        printf("in this not this\n");
        return -4;
    }
    if (((uint64_t) c->n_cigar << 2) + c->l_qname + c->l_extranul
        + (((uint64_t) c->l_qseq + 1) >> 1) + c->l_qseq > (uint64_t) new_l_data) {
//        printf("in this not in this not in this\n");
        return -4;
    }


    if (realloc_bam_data(b, new_l_data) < 0) return -4;
    b->l_data = new_l_data;

    if (Rabbit_bgzf_read(fq, b->data, c->l_qname) != c->l_qname) return -4;
    if (b->data[c->l_qname - 1] != '\0') { // Try to fix missing NUL termination
//        printf("inininininini\n");
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

// Fix bad records where qname is not terminated correctly.
//没处理完全懂？
//static int fixup_missing_qname_nul(bam1_t *b) {
//    bam1_core_t *c = &b->core;
//
//    // Note this is called before c->l_extranul is added to c->l_qname
//    if (c->l_extranul > 0) {
//        b->data[c->l_qname++] = '\0';
//        c->l_extranul--;
//    } else {
//        if (b->l_data > INT_MAX - 4) return -1;
//        if (realloc_bam_data(b, b->l_data + 4) < 0) return -1;
//        b->l_data += 4;
//        b->data[c->l_qname++] = '\0';
//        c->l_extranul = 3;
//    }
//    return 0;
//}
int find_divide_pos(bam_block *block, int last_pos) {
    int divide_pos = last_pos;
    int ret = 0;
    uint32_t x[8], new_l_data;
    while (divide_pos < block->length) {
        Rabbit_memcpy(&ret, block->data + divide_pos, 4);
//        printf("ret is %d\n",ret);
//        printf("divide_pos is %d\n",divide_pos);
//        printf("block length is %d\n",block->length);
        if (ret >= 32) {
            if (divide_pos + 4 + 32 > block->length) {
                break;
            }
            Rabbit_memcpy(x, block->data + divide_pos + 4, 32);
            int pos = (int32_t) x[1];
            int l_qname = x[2] & 0xff;
            int l_extranul = (l_qname % 4 != 0) ? (4 - l_qname % 4) : 0;
            int n_cigar = x[3] & 0xffff;
            int l_qseq = x[4];
            new_l_data = ret - 32 + l_extranul;//block_len + c->l_extranul
            if (new_l_data > INT_MAX || l_qseq < 0 || l_qname < 1) {
//                printf("ai 32 我的老天爷啊\n");
                divide_pos += 4 + 32;
                continue;
            }
            if (((uint64_t) n_cigar << 2) + l_qname + l_extranul
                + (((uint64_t) l_qseq + 1) >> 1) + l_qseq > (uint64_t) new_l_data) {
//                printf("ai 32 我的老天爷啊\n");
                divide_pos += 4 + 32;
                continue;
            }
            if (divide_pos + 4 + 32 + l_qname > block->length) {
//                printf("ai lqname Wrong!!!\n");
                break;
            }
            char fg_char;
            Rabbit_memcpy(&fg_char, block->data + divide_pos + 4 + 32 + l_qname - 1, 1);
            if (fg_char != '\0' && l_extranul <= 0 && new_l_data > INT_MAX - 4) {
//                printf("this is wrong\n");
                divide_pos += 4 + 32 + l_qname;
                continue;
            }

            if (divide_pos + 4 + ret > block->length) {
                break;
            }
            divide_pos += 4 + ret;
        } else {
//            printf("BIG WRONG!!!\n");
            if (divide_pos + 4 > block->length) {
                break;
            }
            divide_pos += 4;
        }
//        printf("One Block Size is %d\n",ret);

    }
//    if (block->length!=divide_pos && block->length - divide_pos < 4) printf("BIG WRONG!!!\n");
    return divide_pos;
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
    packInt32((uint8_t *) &dst[*dlen - 8], crc);
    packInt32((uint8_t *) &dst[*dlen - 4], slen);
    return 0;
}

int rabbit_bgzf_gzip_compress(BGZF *fp, void *_dst, size_t *dlen, const void *src, size_t slen, int level) {
    uint8_t *dst = (uint8_t *) _dst;
    z_stream *zs = fp->gz_stream;
    int flush = slen ? Z_PARTIAL_FLUSH : Z_FINISH;
    zs->next_in = (Bytef *) src;
    zs->avail_in = slen;
    zs->next_out = dst;
    zs->avail_out = *dlen;
    int ret = deflate(zs, flush);
    if (ret == Z_STREAM_ERROR) {
        hts_log_error("Deflate operation failed: %s", bgzf_zerr(ret, NULL));
        return -1;
    }
    if (zs->avail_in != 0) {
        hts_log_error("Deflate block too large for output buffer");
        return -1;
    }
    *dlen = *dlen - zs->avail_out;
    return 0;
}

int find_divide_pos(bam_complete_block *block, int last_pos) {
    int divide_pos = last_pos;
    int ret = 0;
    uint32_t x[8], new_l_data;
    while (divide_pos < block->length) {
        Rabbit_memcpy(&ret, block->data + divide_pos, 4);
//        printf("ret is %d\n",ret);
//        printf("divide_pos is %d\n",divide_pos);
//        printf("block length is %d\n",block->length);
        if (ret >= 32) {
            if (divide_pos + 4 + 32 > block->length) {

                break;
            }
            Rabbit_memcpy(x, block->data + divide_pos + 4, 32);
            int pos = (int32_t) x[1];
            int l_qname = x[2] & 0xff;
            int l_extranul = (l_qname % 4 != 0) ? (4 - l_qname % 4) : 0;
            int n_cigar = x[3] & 0xffff;
            int l_qseq = x[4];
            new_l_data = ret - 32 + l_extranul;//block_len + c->l_extranul
            if (new_l_data > INT_MAX || l_qseq < 0 || l_qname < 1) {
//                printf("ai 32 我的老天爷啊\n");
                divide_pos += 4 + 32;
                continue;
            }
            if (((uint64_t) n_cigar << 2) + l_qname + l_extranul
                + (((uint64_t) l_qseq + 1) >> 1) + l_qseq > (uint64_t) new_l_data) {
//                printf("ai 32 我的老天爷啊\n");
                divide_pos += 4 + 32;
                continue;
            }
            if (divide_pos + 4 + 32 + l_qname > block->length) {
//                printf("ai lqname Wrong!!!\n");
                break;
            }
            char fg_char;
            Rabbit_memcpy(&fg_char, block->data + divide_pos + 4 + 32 + l_qname - 1, 1);
            if (fg_char != '\0') {
//                printf("this is wrong\n");
            }
            if (fg_char != '\0' && l_extranul <= 0 && new_l_data > INT_MAX - 4) {

                if (divide_pos + 4 + 32 + l_qname > block->length) {
//                    printf("in this not happy!\n");
                    break;
                }
                divide_pos += 4 + 32 + l_qname;
                continue;
            }

            if (divide_pos + 4 + ret > block->length) {
                break;
            }
            divide_pos += 4 + ret;
        } else {
//            printf("BIG WRONG!!!\n");
            if (divide_pos + 4 > block->length) {
                break;
            }
            divide_pos += 4;
        }
//        printf("One Block Size is %d\n",ret);

    }
//    if (block->length!=divide_pos && block->length - divide_pos < 4) printf("BIG WRONG!!!\n");
    return divide_pos;
}

std::pair<int, int> find_divide_pos_and_get_read_number(bam_block *block, int last_pos) {
    int divide_pos = last_pos;
    int ans = 0;
    int ret = 0;
    uint32_t x[8], new_l_data;
    while (divide_pos < block->length) {
        Rabbit_memcpy(&ret, block->data + divide_pos, 4);
//        printf("ret is %d\n",ret);
//        printf("divide_pos is %d\n",divide_pos);
//        printf("block length is %d\n",block->length);
        if (ret >= 32) {
            if (divide_pos + 4 + 32 > block->length) {

                break;
            }
            Rabbit_memcpy(x, block->data + divide_pos + 4, 32);
            int pos = (int32_t) x[1];
            int l_qname = x[2] & 0xff;
            int l_extranul = (l_qname % 4 != 0) ? (4 - l_qname % 4) : 0;
            int n_cigar = x[3] & 0xffff;
            int l_qseq = x[4];
            new_l_data = ret - 32 + l_extranul;//block_len + c->l_extranul
            if (new_l_data > INT_MAX || l_qseq < 0 || l_qname < 1) {
//                printf("ai 32 我的老天爷啊\n");
                divide_pos += 4 + 32;
                continue;
            }
            if (((uint64_t) n_cigar << 2) + l_qname + l_extranul
                + (((uint64_t) l_qseq + 1) >> 1) + l_qseq > (uint64_t) new_l_data) {
//                printf("ai 32 我的老天爷啊\n");
                divide_pos += 4 + 32;
                continue;
            }
            if (divide_pos + 4 + 32 + l_qname > block->length) {
//                printf("ai lqname Wrong!!!\n");
                break;
            }
            char fg_char;
            Rabbit_memcpy(&fg_char, block->data + divide_pos + 4 + 32 + l_qname - 1, 1);
            if (fg_char != '\0') {
//                printf("this is wrong\n");
            }
            if (fg_char != '\0' && l_extranul <= 0 && new_l_data > INT_MAX - 4) {

                if (divide_pos + 4 + 32 + l_qname > block->length) {
//                    printf("in this not happy!\n");
                    break;
                }
                divide_pos += 4 + 32 + l_qname;
                continue;
            }

            if (divide_pos + 4 + ret > block->length) {
                break;
            }
            divide_pos += 4 + ret;
            ans++;
        } else {
//            printf("BIG WRONG!!!\n");
            if (divide_pos + 4 > block->length) {
                break;
            }
            divide_pos += 4;
        }
//        printf("One Block Size is %d\n",ret);

    }
//    if (block->length!=divide_pos && block->length - divide_pos < 4) printf("BIG WRONG!!!\n");
    return std::pair<int, int>(divide_pos, ans);
}

std::pair<int, int> find_divide_pos_and_get_read_number(bam_complete_block *block, int last_pos) {
    int divide_pos = last_pos;
    int ans = 0;
    int ret = 0;
    uint32_t x[8], new_l_data;
    while (divide_pos < block->length) {
        Rabbit_memcpy(&ret, block->data + divide_pos, 4);
//        printf("ret is %d\n",ret);
//        printf("divide_pos is %d\n",divide_pos);
//        printf("block length is %d\n",block->length);
        if (ret >= 32) {
            if (divide_pos + 4 + 32 > block->length) {

                break;
            }
            Rabbit_memcpy(x, block->data + divide_pos + 4, 32);
            int pos = (int32_t) x[1];
            int l_qname = x[2] & 0xff;
            int l_extranul = (l_qname % 4 != 0) ? (4 - l_qname % 4) : 0;
            int n_cigar = x[3] & 0xffff;
            int l_qseq = x[4];
            new_l_data = ret - 32 + l_extranul;//block_len + c->l_extranul
            if (new_l_data > INT_MAX || l_qseq < 0 || l_qname < 1) {
//                printf("ai 32 我的老天爷啊\n");
                divide_pos += 4 + 32;
                continue;
            }
            if (((uint64_t) n_cigar << 2) + l_qname + l_extranul
                + (((uint64_t) l_qseq + 1) >> 1) + l_qseq > (uint64_t) new_l_data) {
//                printf("ai 32 我的老天爷啊\n");
                divide_pos += 4 + 32;
                continue;
            }
            if (divide_pos + 4 + 32 + l_qname > block->length) {
//                printf("ai lqname Wrong!!!\n");
                break;
            }
            char fg_char;
            Rabbit_memcpy(&fg_char, block->data + divide_pos + 4 + 32 + l_qname - 1, 1);
            if (fg_char != '\0') {
//                printf("this is wrong\n");
            }
            if (fg_char != '\0' && l_extranul <= 0 && new_l_data > INT_MAX - 4) {

                if (divide_pos + 4 + 32 + l_qname > block->length) {
//                    printf("in this not happy!\n");
                    break;
                }
                divide_pos += 4 + 32 + l_qname;
                continue;
            }

            if (divide_pos + 4 + ret > block->length) {
                break;
            }
            divide_pos += 4 + ret;
            ans++;
        } else {
//            printf("BIG WRONG!!!\n");
            if (divide_pos + 4 > block->length) {
                break;
            }
            divide_pos += 4;
        }
//        printf("One Block Size is %d\n",ret);

    }
//    if (block->length!=divide_pos && block->length - divide_pos < 4) printf("BIG WRONG!!!\n");
    return std::pair<int, int>(divide_pos, ans);
}

int change_data_size(bam_complete_block *block) {
    int new_length;
    if (block->data_size < 8 * BGZF_MAX_BLOCK_SIZE) {
//        printf("In this\n");
        new_length = 2 * block->data_size;
    } else {
//        printf("ai why in this\n");
        new_length = block->data_size + 2 * BGZF_MAX_BLOCK_SIZE;
    }
    unsigned char *data_new = new unsigned char[new_length];
    memcpy(data_new, block->data, block->length * sizeof(unsigned char));
//    delete [] block->data;
    block->data = data_new;
    block->data_size = new_length;
    return new_length;
}

//int rabbit_write_deflate_block(BGZF *fp, bam_write_block* write_block){
//    size_t comp_size = BGZF_MAX_BLOCK_SIZE;
//    int ret;
//    if ( !fp->is_gzip )
//        ret = rabbit_bgzf_compress(write_block->compressed_data, &comp_size, write_block->uncompressed_data, write_block->block_offset, fp->compress_level);
//    else
//        ret = rabbit_bgzf_gzip_compress(fp, write_block->compressed_data, &comp_size, write_block->uncompressed_data, write_block->block_offset, fp->compress_level);
//
//    if ( ret != 0 )
//    {
//        hts_log_debug("Compression error %d", ret);
//        fp->errcode |= BGZF_ERR_ZLIB;
//        return -1;
//    }
//    return comp_size;
//}
//int rabbit_bgzf_flush(BGZF *fp,bam_write_block* write_block)
//{
//    //TODO 此处可能会出现问题
//    while (write_block->block_offset > 0) {
//        int block_length;
//        printf("Write Block Offset : %d\n",write_block->block_offset);
//        block_length = rabbit_write_deflate_block(fp, write_block);
//        if (block_length < 0) {
//            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block_length, NULL));
//            return -1;
//        }
//        if (write_block== nullptr){
//            printf("Write Game Over!!!\n");
//        }
//        if (hwrite(fp->fp, write_block->compressed_data, block_length) != block_length) {
//            printf("Write Failed\n");
//            hts_log_error("File write failed (wrong size)");
//            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
//            return -1;
//        }
//
//
//        write_block->block_offset=0;
//        fp->block_address += block_length;
//    }
//    write_block->block_offset=0;
//    return 0;
//}
//int rabbit_bgzf_mul_flush(BGZF *fp,BamWriteCompress *bam_write_compress,bam_write_block* &write_block)
//{
////    printf("Try to input One Uncompressed data\n");
////    InputBlockNum++;
////    printf("Write Block Offset : %d\n",write_block->block_offset);
//    bam_write_compress->inputUnCompressData(write_block);
//    write_block=bam_write_compress->getEmpty();
////    printf("Get Another Empty Block Block Num : %d\n",write_block->block_num);
//    return 0;
//}
//int rabbit_bgzf_write(BGZF *fp,bam_write_block* &write_block,const void *data, size_t length)
//{
//    const uint8_t *input = (const uint8_t*)data;
//    ssize_t remaining = length;
////    assert(fp->is_write);
//    while (remaining > 0) {
//        uint8_t* buffer = (uint8_t*)write_block->uncompressed_data;
//        int copy_length = BGZF_BLOCK_SIZE - write_block->block_offset;
//        if (copy_length > remaining) copy_length = remaining;
//        memcpy(buffer + write_block->block_offset, input, copy_length);
//        write_block->block_offset += copy_length;
//        input += copy_length;
//        remaining -= copy_length;
//        if (write_block->block_offset == BGZF_BLOCK_SIZE) {
//            //BUG  write block 不是 引用，所以没有改变
//            if (rabbit_bgzf_flush(fp,write_block) != 0) return -1;
//        }
//    }
//    return length - remaining;
//}
//int rabbit_bgzf_mul_write(BGZF *fp, BamWriteCompress *bam_write_compress,bam_write_block* &write_block,const void *data, size_t length)
//{
//    const uint8_t *input = (const uint8_t*)data;
//    ssize_t remaining = length;
////    assert(fp->is_write);
//    while (remaining > 0) {
//        uint8_t* buffer = (uint8_t*)write_block->uncompressed_data;
////        printf("In rabbit bgzf mul write block Block Num : %d\n",write_block->block_num);
//        int copy_length = BGZF_BLOCK_SIZE - write_block->block_offset;
//        if (copy_length > remaining) copy_length = remaining;
////        printf("In rabbit bgzf mul write block Block Num : %d\n",write_block->block_num);
//        memcpy(buffer + write_block->block_offset, input, copy_length);
//        write_block->block_offset += copy_length;
////        printf("In rabbit bgzf mul write block Block Num : %d\n",write_block->block_num);
//        input += copy_length;
//        remaining -= copy_length;
//        if (write_block->block_offset == BGZF_BLOCK_SIZE) {
//            if (rabbit_bgzf_mul_flush(fp,bam_write_compress,write_block) != 0) return -1;
//        }
//    }
//    return length - remaining;
//}
//int rabbit_bgzf_flush_try(BGZF *fp, bam_write_block* write_block,ssize_t size)
//{
//    if (write_block->block_offset + size > BGZF_BLOCK_SIZE) {
//        return rabbit_bgzf_flush(fp,write_block);
//    }
//    return 0;
//}
//int rabbit_bgzf_mul_flush_try(BGZF *fp,BamWriteCompress* bam_write_compress,bam_write_block* &write_block,ssize_t size)
//{
//    if (write_block->block_offset + size > BGZF_BLOCK_SIZE) {
//        return rabbit_bgzf_mul_flush(fp,bam_write_compress,write_block);
//    }
//    return 0;
//}
//int bam_write_pack(BGZF *fp,BamWriteCompress *bam_write_compress){
//    bam_write_block* block;
//    while(1){
////        printf("Try to Get Compress Data\n");
//        block=bam_write_compress->getCompressData();
////        printf("Has Get One Compress Data\n");
//        if (block == nullptr){
//            return 0;
//        }
////        std::this_thread::sleep_for(std::chrono::nanoseconds(5));
//        if (block->block_length < 0) {
//            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block->block_length, NULL));
//            return -1;
//        }
//        if (hwrite(fp->fp, block->compressed_data, block->block_length) != block->block_length) {
////            printf("Write Failed\n");
//            hts_log_error("File write failed (wrong size)");
//            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
//            return -1;
//        }
////        printf("Has write One Block\n");
//        block->block_offset=0;
//        fp->block_address += block->block_length;
//        bam_write_compress->backEmpty(block);
//    }
//}
//void bam_write_compress_pack(BGZF *fp,BamWriteCompress *bam_write_compress){
////    printf("Start Compress\n");
//    bam_write_block * block;
//    while (1){
//        // fg = getRead(comp);
//        //printf("%d is not get One compressed data\n",id);
////        printf("Has Start Try to Get One Uncompress\n");
//        block=bam_write_compress->getUnCompressData();
////        printf("Has get One Uncompress data\n");
//
//        //printf("%d is get One compressed data\n",id);
//        if (block == nullptr) {
//            //printf("%d is Over\n",id);
//            break;
//        }
////        printf("This Uncompress data block num : %d\n",block->block_num);
//        /*
//         * 压缩
//         */
////        int block_num=block->block_num;
//        block->block_length = rabbit_write_deflate_block(fp, block);
////        printf("Has Compress One Block\n");
//        bam_write_compress->inputCompressData(block);
////        printf("Can,t input Compress Data\n");
////        while (!compress->tryinputUnCompressData(un_comp,comp.second)){
////            std::this_thread::sleep_for(std::chrono::milliseconds(1));
////        }
//    }
////    printf("One Compress Thread has been over!\n");
//    bam_write_compress->CompressThreadComplete();
//}
//
//int rabbit_bam_write_test(BGZF *fp,bam_write_block* write_block,bam1_t *b){
//    const bam1_core_t *c = &b->core;
//    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
//    int i, ok;
//    if (c->l_qname - c->l_extranul > 255) {
//        hts_log_error("QNAME \"%s\" is longer than 254 characters", bam_get_qname(b));
//        errno = EOVERFLOW;
//        return -1;
//    }
//    if (c->n_cigar > 0xffff) block_len += 16; // "16" for "CGBI", 4-byte tag length and 8-byte fake CIGAR
//    if (c->pos > INT_MAX ||
//        c->mpos > INT_MAX ||
//        c->isize < INT_MIN || c->isize > INT_MAX) {
//        hts_log_error("Positional data is too large for BAM format");
//        return -1;
//    }
//    x[0] = c->tid;
//    x[1] = c->pos;
//    x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | (c->l_qname - c->l_extranul);
//    if (c->n_cigar > 0xffff) x[3] = (uint32_t)c->flag << 16 | 2;
//    else x[3] = (uint32_t)c->flag << 16 | (c->n_cigar & 0xffff);
//    x[4] = c->l_qseq;
//    x[5] = c->mtid;
//    x[6] = c->mpos;
//    x[7] = c->isize;
//    ok = (rabbit_bgzf_flush_try(fp, write_block, 4 + block_len) >= 0);
//    if (fp->is_be) {
//        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
//        y = block_len;
//        if (ok) ok = (rabbit_bgzf_write(fp, write_block,ed_swap_4p(&y), 4) >= 0);
//        swap_data(c, b->l_data, b->data, 1);
//    } else {
//        if (ok) {
//            ok = (rabbit_bgzf_write(fp, write_block, &block_len, 4) >= 0);
//        }
//    }
//    if (ok) {
//        ok = (rabbit_bgzf_write(fp, write_block, x, 32) >= 0);
//    }
//    if (ok) ok = (rabbit_bgzf_write(fp, write_block, b->data, c->l_qname - c->l_extranul) >= 0);
//    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
//        if (ok) ok = (rabbit_bgzf_write(fp, write_block, b->data + c->l_qname, b->l_data - c->l_qname) >= 0);
//    } else { // with long CIGAR, insert a fake CIGAR record and move the real CIGAR to the CG:B,I tag
//        uint8_t buf[8];
//        uint32_t cigar_st, cigar_en, cigar[2];
//        hts_pos_t cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
//        if (cigreflen >= (1<<28)) {
//            // Length of reference covered is greater than the biggest
//            // CIGAR operation currently allowed.
//            hts_log_error("Record %s with %d CIGAR ops and ref length %" PRIhts_pos
//                                  " cannot be written in BAM.  Try writing SAM or CRAM instead.\n",
//                          bam_get_qname(b), c->n_cigar, cigreflen);
//            return -1;
//        }
//        cigar_st = (uint8_t*)bam_get_cigar(b) - b->data;
//        cigar_en = cigar_st + c->n_cigar * 4;
//        cigar[0] = (uint32_t)c->l_qseq << 4 | BAM_CSOFT_CLIP;
//        cigar[1] = (uint32_t)cigreflen << 4 | BAM_CREF_SKIP;
//        u32_to_le(cigar[0], buf);
//        u32_to_le(cigar[1], buf + 4);
//        if (ok) {
//            ok = (rabbit_bgzf_write(fp, write_block, buf, 8) >= 0); // write cigar: <read_length>S<ref_length>N
//        }
//        if (ok) {
//            ok = (rabbit_bgzf_write(fp, write_block, &b->data[cigar_en], b->l_data - cigar_en) >= 0); // write data after CIGAR
//        }
//        if (ok) {
//            ok = (rabbit_bgzf_write(fp, write_block, "CGBI", 4) >= 0); // write CG:B,I
//        }
//        u32_to_le(c->n_cigar, buf);
//        if (ok) {
//            ok = (rabbit_bgzf_write(fp, write_block, buf, 4) >= 0); // write the true CIGAR length
//        }
//        if (ok) {
//            ok = (rabbit_bgzf_write(fp, write_block, &b->data[cigar_st], c->n_cigar * 4) >= 0); // write the real CIGAR
//        }
//    }
//    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
//    return ok? 4 + block_len : -1;
//}
//int rabbit_bam_write_mul_test(BGZF *fp,BamWriteCompress *bam_write_compress,bam_write_block* &write_block,bam1_t *b){
//    const bam1_core_t *c = &b->core;
//    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
//    int i, ok;
//    if (c->l_qname - c->l_extranul > 255) {
//        hts_log_error("QNAME \"%s\" is longer than 254 characters", bam_get_qname(b));
//        errno = EOVERFLOW;
//        return -1;
//    }
//    if (c->n_cigar > 0xffff) block_len += 16; // "16" for "CGBI", 4-byte tag length and 8-byte fake CIGAR
//    if (c->pos > INT_MAX ||
//        c->mpos > INT_MAX ||
//        c->isize < INT_MIN || c->isize > INT_MAX) {
//        hts_log_error("Positional data is too large for BAM format");
//        return -1;
//    }
//    x[0] = c->tid;
//    x[1] = c->pos;
//    x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | (c->l_qname - c->l_extranul);
//    if (c->n_cigar > 0xffff) x[3] = (uint32_t)c->flag << 16 | 2;
//    else x[3] = (uint32_t)c->flag << 16 | (c->n_cigar & 0xffff);
//    x[4] = c->l_qseq;
//    x[5] = c->mtid;
//    x[6] = c->mpos;
//    x[7] = c->isize;
//    ok = (rabbit_bgzf_mul_flush_try(fp,bam_write_compress, write_block, 4 + block_len) >= 0);
//    if (fp->is_be) {
//        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
//        y = block_len;
//        if (ok) ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block,ed_swap_4p(&y), 4) >= 0);
//        swap_data(c, b->l_data, b->data, 1);
//    } else {
//        if (ok) {
//            ok = (rabbit_bgzf_mul_write(fp,bam_write_compress, write_block, &block_len, 4) >= 0);
//        }
//    }
//    if (ok) {
//        ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, x, 32) >= 0);
//    }
//    if (ok) ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, b->data, c->l_qname - c->l_extranul) >= 0);
//    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
//        if (ok) ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, b->data + c->l_qname, b->l_data - c->l_qname) >= 0);
//    } else { // with long CIGAR, insert a fake CIGAR record and move the real CIGAR to the CG:B,I tag
//        uint8_t buf[8];
//        uint32_t cigar_st, cigar_en, cigar[2];
//        hts_pos_t cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
//        if (cigreflen >= (1<<28)) {
//            // Length of reference covered is greater than the biggest
//            // CIGAR operation currently allowed.
//            hts_log_error("Record %s with %d CIGAR ops and ref length %" PRIhts_pos
//                                  " cannot be written in BAM.  Try writing SAM or CRAM instead.\n",
//                          bam_get_qname(b), c->n_cigar, cigreflen);
//            return -1;
//        }
//        cigar_st = (uint8_t*)bam_get_cigar(b) - b->data;
//        cigar_en = cigar_st + c->n_cigar * 4;
//        cigar[0] = (uint32_t)c->l_qseq << 4 | BAM_CSOFT_CLIP;
//        cigar[1] = (uint32_t)cigreflen << 4 | BAM_CREF_SKIP;
//        u32_to_le(cigar[0], buf);
//        u32_to_le(cigar[1], buf + 4);
//        if (ok) {
//            ok = (rabbit_bgzf_mul_write(fp,bam_write_compress, write_block, buf, 8) >= 0); // write cigar: <read_length>S<ref_length>N
//        }
//        if (ok) {
//            ok = (rabbit_bgzf_mul_write(fp,bam_write_compress,write_block, &b->data[cigar_en], b->l_data - cigar_en) >= 0); // write data after CIGAR
//        }
//        if (ok) {
//            ok = (rabbit_bgzf_mul_write(fp,bam_write_compress,write_block, "CGBI", 4) >= 0); // write CG:B,I
//        }
//        u32_to_le(c->n_cigar, buf);
//        if (ok) {
//            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, buf, 4) >= 0); // write the true CIGAR length
//        }
//        if (ok) {
//            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, &b->data[cigar_st], c->n_cigar * 4) >= 0); // write the real CIGAR
//        }
//    }
//    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
//    return ok? 4 + block_len : -1;
//}


















































