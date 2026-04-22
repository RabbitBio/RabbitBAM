#include "BamTools.h"


void print_bam1(const bam1_t *b)
{
    if (!b) {
        printf("bam1_t pointer is NULL\n");
        return;
    }

    printf("bam1_t at %p {\n", (void*)b);

    // core 信息
    printf("  core:\n");
    printf("    pos=%d, tid=%d, bin=%u, qual=%u, l_extranul=%u\n",
           b->core.pos, b->core.tid, b->core.bin, b->core.qual, b->core.l_extranul);
    printf("    flag=%u, l_qname=%u, n_cigar=%u, l_qseq=%d\n",
           b->core.flag, b->core.l_qname, b->core.n_cigar, b->core.l_qseq);
    printf("    mtid=%d, mpos=%d, isize=%d\n",
           b->core.mtid, b->core.mpos, b->core.isize);

    // id
    printf("  id          : %lu\n", b->id);

    // data 指针
    printf("  data ptr    : %p\n", (void*)b->data);

    // data 内容前 16 字节
    if (b->data && b->l_data > 0) {
        int n = b->l_data < 16 ? b->l_data : 16;
        printf("  data content: ");
        for (int i = 0; i < n; i++) {
            printf("%02X ", b->data[i]);
        }
        if (b->l_data > 16) printf("... (%d bytes total)", b->l_data);
        printf("\n");
    } else {
        printf("  data content: NULL or empty\n");
    }

    // 长度信息
    printf("  l_data      : %d\n", b->l_data);
    printf("  m_data      : %u\n", b->m_data);

    // 内存策略
    printf("  mempolicy   : %u\n", b->mempolicy);

    printf("}\n");
}

void print_bam_block(struct bam_block *blk) {
    if (!blk) {
        printf("bam_block pointer is NULL\n");
        return;
    }

    printf("bam_block at %p {\n", (void*)blk);
    printf("  errcode      : %u\n", blk->errcode);
    printf("  data ptr     : %p\n", (void*)blk->data);
    if (blk->data != NULL && blk->length > 0) {
        // 打印前 16 字节或实际长度
        unsigned int print_len = blk->length < 16 ? blk->length : 16;
        printf("  data content : ");
        for (unsigned int i = 0; i < print_len; i++) {
            printf("%02X ", blk->data[i]);
        }
        if (blk->length > 16) {
            printf("... (%u bytes total)", blk->length);
        }
        printf("\n");
    }
    printf("  length       : %u\n", blk->length);
    printf("  pos          : %u\n", blk->pos);
    printf("  block_address: %lld\n", (long long)blk->block_address);
    printf("  block_id     : %d\n", blk->block_id);
    printf("}\n");
}

void bam_destroy1_sw(bam1_t *b){
    if (b == NULL) return;

    // 1. 释放 data 区域
    if (b->data) {
        // 如果是用户自定义的对齐内存，则用 aligned_free_custom
        if (bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) {
            aligned_free_custom(b->data);
        }
        else {
            // 普通 HTSlib 分配（malloc/realloc）
            free(b->data);
        }

        // 清空，避免悬空指针
        b->data = NULL;
        b->m_data = 0;
        b->l_data = 0;
    }

    // 2. 释放 bam1_t 结构体本身
    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_STRUCT) == 0) {
        free(b);
    }
}

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
        new_data = (uint8_t * )(realloc(b->data, new_m_data));
    } else {
        if ((new_data = (uint8_t * )(malloc(new_m_data))) != NULL) {
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

int realloc_bam_data(bam1_t *b, size_t desired) {
    if (desired <= b->m_data) return 0;
    return sam_realloc_bam_data(b, desired);
}


inline int unpackInt16(const uint8_t *buffer) {
    return buffer[0] | buffer[1] << 8;
}

int load_block_from_cache(BGZF *fp, int64_t block_address) {
    khint_t k;
    cache_t * p;

    khash_t(cache) * h = fp->cache->h;
    k = kh_get(cache, h, block_address);
    if (k == kh_end(h)) return 0;
    p = &kh_val(h, k);
    if (fp->block_length != 0) fp->block_offset = 0;
    fp->block_address = block_address;
    fp->block_length = p->size;
    memcpy(fp->uncompressed_block, p->block, p->size);
    if (hseek(fp->fp, p->end_offset, SEEK_SET) < 0) {
        // todo: move the error up
        hts_log_error("Could not hseek to %"
        PRId64, p->end_offset);
        exit(1);
    }
    return p->size;
}

int check_header(const uint8_t *header) {
    if (header[0] != 31 || header[1] != 139 || header[2] != 8) return -2;
    return ((header[3] & 4) != 0
            && unpackInt16((uint8_t * ) & header[10]) == 6
            && header[12] == 'B' && header[13] == 'C'
            && unpackInt16((uint8_t * ) & header[14]) == 2) ? 0 : -1;
}

//读取一个压缩块，放到j中
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
    block_length = unpackInt16((uint8_t * ) & header[16]) + 1; // +1 because when writing this number, we used "-1"
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

int sam_write1_sw(samFile *fp, const sam_hdr_t *h, const bam1_t *b){

    fp->format.category = sequence_data;
    fp->format.format = sam;

    if (sam_format1(h, b, &fp->line) < 0) return -1;
    kputc('\n', &fp->line);

    if ( hwrite(fp->fp.hfile, fp->line.s, fp->line.l) != fp->line.l ) return -1;     

    return fp->line.l;

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

int sam_read1_sw(samFile *fp, sam_hdr_t *h,  bam1_t *b){
    int ret = sam_read1(fp, h, b);
    return ret;
}

void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host) {
    uint32_t *cigar = (uint32_t * )(data + c->l_qname);
    uint32_t i;
    for (i = 0; i < c->n_cigar; ++i) ed_swap_4p(&cigar[i]);
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



int check_bam_cross_block_ex(const char *bam_path, BamCrossBlockStats *stats, bool verbose) {
    samFile *in = sam_open(bam_path, "r");
    if (!in) {
        if (verbose) fprintf(stderr, "[check_bam_cross_block] ERROR: cannot open '%s'\n", bam_path);
        return -1;
    }
    if (in->format.format != bam) {
        if (verbose) fprintf(stderr, "[check_bam_cross_block] ERROR: '%s' is not a BAM file\n", bam_path);
        sam_close(in);
        return -1;
    }
    sam_hdr_t *hdr = sam_hdr_read(in);
    if (!hdr) {
        if (verbose) fprintf(stderr, "[check_bam_cross_block] ERROR: failed to read BAM header\n");
        sam_close(in);
        return -1;
    }
    sam_hdr_destroy(hdr);

    BGZF *fp = in->fp.bgzf;
    long long total = 0, cross = 0;

    int64_t skip_remaining = 0;  // 上一块末尾跨块 record 还需跳过的字节数
    uint8_t split_hdr[4];        // 4 字节头部被切分时已收集的部分字节
    int     split_hdr_bytes = 0; // 已收集的头部字节数（>0 表示头部被切分）

    while (true) {
        const uint8_t *blk = (const uint8_t *)fp->uncompressed_block;
        int blen = fp->block_length;
        int pos  = fp->block_offset;  

        //1.处理上一块的跨块遗留
        if (split_hdr_bytes > 0) {
            // 4 字节头部被切分，继续从本块收集剩余字节
            while (pos < blen && split_hdr_bytes < 4)
                split_hdr[split_hdr_bytes++] = blk[pos++];
            if (split_hdr_bytes < 4) goto advance;  // 极罕见：头部跨 3 块以上

            int32_t bsize = (int32_t)((uint32_t)split_hdr[0]       |
                                      (uint32_t)split_hdr[1] <<  8  |
                                      (uint32_t)split_hdr[2] << 16  |
                                      (uint32_t)split_hdr[3] << 24);
            split_hdr_bytes = 0;
            if (bsize < 0) break;  // 非法值，退出
            int avail = blen - pos;
            if (bsize <= avail)    pos += bsize;          // body 在本块内
            else                   { skip_remaining = bsize - avail; pos = blen; }
        }

        if (skip_remaining > 0) {
            // 跨块 record 的 body 还需继续跳过
            int avail = blen - pos;
            if ((int64_t)avail >= skip_remaining) {
                pos += (int)skip_remaining;
                skip_remaining = 0;
            } else {
                skip_remaining -= avail;
                pos = blen;
            }
        }

        //2.逐个扫描本块内的 BAM 记录
        while (pos < blen) {
            if (blen - pos < 4) {
                // 4 字节头部跨块（头部的 1~3 字节在本块末尾）
                cross++; total++;
                memcpy(split_hdr, blk + pos, blen - pos);
                split_hdr_bytes = blen - pos;
                pos = blen;
                break;
            }
            int32_t bsize = (int32_t)((uint32_t)blk[pos]       |
                                       (uint32_t)blk[pos+1] <<  8 |
                                       (uint32_t)blk[pos+2] << 16 |
                                       (uint32_t)blk[pos+3] << 24);
            if (bsize <= 0) break;  // 块末尾填充 / 非法值，停止扫描

            if (pos + 4 + bsize > blen) {
                // record body 跨块
                cross++; total++;
                skip_remaining = (int64_t)(pos + 4 + bsize) - blen;
                pos = blen;
                break;
            }
            pos += 4 + bsize;
            total++;
        }

advance:
        //3.读取下一块
        fp->block_offset = blen;  
        if (bgzf_read_block(fp) != 0 || fp->block_length == 0) break;  
    }

    sam_close(in);

    if (stats) {
        stats->total_records = total;
        stats->cross_block_records = cross;
    }

    if (verbose) {
        printf("[check_bam_cross_block] file: %s\n", bam_path);
        printf("  total records  : %lld\n", total);
        printf("  cross-block    : %lld\n", cross);
        if (cross > 0) {
            printf("  result: CROSS-BLOCK RECORDS EXIST — bam2sam fast-path NOT safe\n");
        } else {
            printf("  result: NO cross-block records — bam2sam fast-path is safe\n");
        }
    }

    return cross > 0 ? 1 : 0;
}

// 检测 BAM 文件中是否存在跨 BGZF 块的记录
int check_bam_cross_block(const char *bam_path) {
    return check_bam_cross_block_ex(bam_path, nullptr, true);
}
