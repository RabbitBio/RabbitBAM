
#include <cstring>
#include <climits>
#include <cassert>

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "BamTools.h"
#include "libdeflate.h"
#include <htslib/hts_endian.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

 //----------------------------------------------------------------------------------------------
 //bam to sam 剩下需要的

 static char get_severity_tag(enum htsLogLevel severity)
{
    switch (severity) {
    case HTS_LOG_ERROR:
        return 'E';
    case HTS_LOG_WARNING:
        return 'W';
    case HTS_LOG_INFO:
        return 'I';
    case HTS_LOG_DEBUG:
        return 'D';
    case HTS_LOG_TRACE:
        return 'T';
    default:
        break;
    }

    return '*';
}

/*! Logs an event with severity HTS_LOG_ERROR and default context. Parameters: format, ... */
#define hts_log_error(...) hts_log(HTS_LOG_ERROR, __func__, __VA_ARGS__)

void hts_log(enum htsLogLevel severity, const char *context, const char *format, ...)
{
    int save_errno = errno;
    if (severity <= hts_verbose) {
        va_list argptr;

        fprintf(stderr, "[%c::%s] ", get_severity_tag(severity), context);

        va_start(argptr, format);
        vfprintf(stderr, format, argptr);
        va_end(argptr);

        fprintf(stderr, "\n");
    }
    errno = save_errno;
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

void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host) {
    uint32_t *cigar = (uint32_t * )(data + c->l_qname);
    uint32_t i;
    for (i = 0; i < c->n_cigar; ++i) ed_swap_4p(&cigar[i]);
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

static inline int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

static inline uint8_t *skip_aux(uint8_t *s, uint8_t *end)
{
    int size;
    uint32_t n;
    if (s >= end) return end;
    size = aux_type2size(*s); ++s; // skip type
    switch (size) {
    case 'Z':
    case 'H':
        while (s < end && *s) ++s;
        return s < end ? s + 1 : end;
    case 'B':
        if (end - s < 5) return NULL;
        size = aux_type2size(*s); ++s;
        n = le_to_u32(s);
        s += 4;
        if (size == 0 || end - s < size * n) return NULL;
        return s + size * n;
    case 0:
        return NULL;
    default:
        if (end - s < size) return NULL;
        return s + size;
    }
}

uint8_t *bam_aux_first(const bam1_t *b)
{
    uint8_t *s = bam_get_aux(b);
    uint8_t *end = b->data + b->l_data;
    if (end - s <= 2) { errno = ENOENT; return NULL; }
    return s+2;
}

uint8_t *bam_aux_next(const bam1_t *b, const uint8_t *s)
{
    uint8_t *end = b->data + b->l_data;
    uint8_t *next = s? skip_aux((uint8_t *) s, end) : end;
    if (next == NULL) goto bad_aux;
    if (end - next <= 2) { errno = ENOENT; return NULL; }
    return next+2;

 bad_aux:
    hts_log_error("Corrupted aux data for read %s flag %d",
                  bam_get_qname(b), b->core.flag);
    errno = EINVAL;
    return NULL;
}

uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
{
    uint8_t *s;
    for (s = bam_aux_first(b); s; s = bam_aux_next(b, s))
        if (s[-2] == tag[0] && s[-1] == tag[1]) {
            // Check the tag value is valid and complete
            uint8_t *e = skip_aux(s, b->data + b->l_data);
            if (e == NULL) goto bad_aux;
            if ((*s == 'Z' || *s == 'H') && *(e - 1) != '\0') goto bad_aux;

            return s;
        }

    // errno now as set by bam_aux_first()/bam_aux_next()
    return NULL;

 bad_aux:
    hts_log_error("Corrupted aux data for read %s flag %d",
                  bam_get_qname(b), b->core.flag);
    errno = EINVAL;
    return NULL;
}

hts_pos_t bam_endpos(const bam1_t *b)
{
    hts_pos_t rlen = (b->core.flag & BAM_FUNMAP)? 0 : bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
    if (rlen == 0) rlen = 1;
    return b->core.pos + rlen;
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




  //----------------------------------------------------------------------------------------------
 // sam to bam 剩下需要的

 #define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3) // bit 1: consume query; bit 2: consume reference
 #define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
 #define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
 hts_pos_t bam_cigar2rlen(int n_cigar, const uint32_t *cigar)
{
    int k;
    hts_pos_t l;
    for (k = l = 0; k < n_cigar; ++k)
        if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
            l += bam_cigar_oplen(cigar[k]);
    return l;
}



 //----------------------------------------------------------------------------------------------
 //bam to sam
// 进行标准的解压缩
int bgzf_uncompress(uint8_t *dst, size_t *dlen,
                    const uint8_t *src, size_t slen,
                    uint32_t expected_crc) {

    //确保输入输出缓冲区都是 64 字节对齐的!不然会报错
    //好像注释掉也没报错？？？
    // const uintptr_t ALIGN_MASK = 63;
    // if (((uintptr_t)dst & ALIGN_MASK) != 0) {
    //     printf("The dst data is not aligned\n");
    //     return -1;
    // }
    // if (((uintptr_t)src & ALIGN_MASK) != 0 ) {
    //     printf("The src data is not aligned\n");
    //     return -1;
    // }


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

int Rabbit_bgzf_read(struct bam_block *fq, void *data, unsigned int length) {
    if (length <= 0) return -1;
    if (length > fq->length - fq->pos) printf("One Block Is Small\n");
    length = fq->pos + length > fq->length ? fq->length - fq->pos : length;
    memcpy((uint8_t *) data, fq->data + fq->pos, length);
    fq->pos += length;
    return length;
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
    //这里如果已经分配了最够大小的data空间，就不会再开辟了
    //printf("The new_l_data is %u\n", new_l_data);
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

// 打印 bam1_t
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

// 打印 bam_block 内容的函数
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

//-------------------------------------------------------------------------------------------------------------------------------------------------
extern "C" void decompressfunc(Para paras[64]) {
    //printf("decompressfunc started with _PEN = %d\n", _PEN);

    int id = _PEN;             // 从核号（0~63）
    Para* para = &paras[id];

    // 1. 读取压缩块数据
    bam_block* comp = para->input_block;
    bam_block* un_comp = para->un_comp_block;
    if (comp == NULL) {
        return;
    }

    printf("decompressfunc started with _PEN = %d\n", _PEN);

    //#define TEST_ALIGN
    #ifdef TEST_ALIGN
        unsigned char *p = comp->data;  // 64B aligned
        unsigned char *u = p + 18;      //没对齐

        //test1
        printf("u = %p\n", u);          // 成功
        printf("%d\n", *u);             // 成功
        if (((uintptr_t)u & 63) != 0) {       // 成功
            printf("u is NOT aligned: addr=%p\n", u);
        }
        uint32_t x = *(uint32_t*)u;     // 失败
        printf("misaligned value = %u\n", x);

        //test2
        //成功与失败要看具体地址情况
        // p = 0x500001878040
        // u = 0x500001878052  (misaligned)
        printf("p = %p\n", p);           
        printf("u    = %p  (misaligned)\n", u); 
        // 1) 测试 1 字节 load —— 成功
        uint8_t b1 = *u;
        printf("1-byte load OK: %u\n", b1);
        // 2) 测试 2 字节 load —— 成功
        uint16_t h = *(uint16_t*)u;  
        printf("2-byte load = %u\n", h);
        // 3) 测试 4 字节 load —— 失败 
        uint32_t w = *(uint32_t*)u;  
        printf("4-byte load = %u\n", w);
        // 4) 测试 8 字节 load ——  
        uint64_t g = *(uint64_t*)u;  
        printf("8-byte load = %u\n", g);
    #endif
    
    //#define TEST_MALLOC
    #ifdef TEST_MALLOC
        uint8_t *p = NULL;
        size_t size = 32;

        // 测试 malloc
        p = (uint8_t *)malloc(size);
        if (!p) {
            printf("CPE: malloc failed\n");
            return;
        }
        printf("CPE: malloc(%zu) succeeded, p=%p\n", size, p);

        // 写入测试数据
        for (size_t i = 0; i < size; i++) p[i] = (uint8_t)i;

        // 测试 realloc
        size_t new_size = 64;
        uint8_t *q = (uint8_t *)realloc(p, new_size);
        if (!q) {
            printf("CPE: realloc failed\n");
            free(p); // 原始内存释放
            return;
        }
        printf("CPE: realloc(%zu) succeeded, q=%p\n", new_size, q);

        // 检查数据是否完整
        for (size_t i = 0; i < size; i++) {
            if (q[i] != i) {
                printf("CPE: data corrupted at index %zu\n", i);
                break;
            }
        }

        free(q);
        printf("CPE: free done\n");
    #endif


    // 2. 解压缩该块
    block_decode_func(comp,un_comp);
    print_bam_block(comp);
    print_bam_block(un_comp);
    printf("Complete the decompression!\n");

    // 3. 解析解压后的SAM/BAM记录
    int count = 0;
    bam1_t* b = NULL;
    b = para->output_records[count]; 
    // print_bam1(b);
    // printf("111\n");
    while(read_bam(un_comp, b, 0)>=0){
        // printf("222\n");
        // para->l_data_list[count] = b->l_data;
        // para->data_list[count] = b->data;
        count++;
        //print_bam1(b);

        b = para->output_records[count];
    }

    printf("Complete parsing with %d bam1_t!\n", count);

    para->n_records = count;
    para->status = 0;


}

//-------------------------------------------------------------------------------------------------------------------------------------------------
extern "C" void copyfunc(Para paras[64]) {

    int id = _PEN;             // 从核号（0~63）
    Para* para = &paras[id];

    if(para->input_block == NULL) {
        return;
    }

    printf("copyfunc started with _PEN = %d\n", _PEN);

    for (int i = 0; i < para->n_records; i++) {
        // 复制数据到新分配的data区域
        memcpy(para->output_records[i]->data, para->data_list[i], para->l_data_list[i]);
    }

    para->status = 0;  // copy 完成标记

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

//-------------------------------------------------------------------------------------------------------------------------------------------------
extern "C" void slave_compressfunc(Comp_Para paras[64]) {

    int id = _PEN;             // 从核号（0~63）
    int compress_level = 1; // 默认压缩等级
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





