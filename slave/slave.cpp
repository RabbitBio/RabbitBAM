
#include <cstring>
#include <climits>
#include <cmath>
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
#include <htslib/khash.h>

#include "test.h"

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
    //这里如果记录跨块会输出One Block Is Small----------------
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

//sam_format(void *arg)需要的-----------------------------------------------------------------------------------------------------------

int kputd(double d, kstring_t *s) {
	int len = 0;
	char buf[21], *cp = buf+20, *ep;
	if (d == 0) {
		if (std::signbit(d)) {
			kputsn("-0",2,s);
			return 2;
		} else {
			kputsn("0",1,s);
			return 1;
		}
	}

	if (d < 0) {
		kputc('-',s);
		len = 1;
		d=-d;
	}
	if (!(d >= 0.0001 && d <= 999999)) {
		if (ks_resize(s, s->l + 50) < 0)
			return EOF;
		// We let stdio handle the exponent cases
		int s2 = snprintf(s->s + s->l, s->m - s->l, "%g", d);
		len += s2;
		s->l += s2;
		return len;
	}

	// Correction for rounding - rather ugly
	// Optimised for small numbers.

	uint32_t i;
	if (d<0.001)         i = rint(d*1000000000), cp -= 1;
	else if (d < 0.01)   i = rint(d*100000000),  cp -= 2;
	else if (d < 0.1)    i = rint(d*10000000),   cp -= 3;
	else if (d < 1)      i = rint(d*1000000),    cp -= 4;
	else if (d < 10)     i = rint(d*100000),     cp -= 5;
	else if (d < 100)    i = rint(d*10000),      cp -= 6;
	else if (d < 1000)   i = rint(d*1000),       cp -= 7;
	else if (d < 10000)  i = rint(d*100),        cp -= 8;
	else if (d < 100000) i = rint(d*10),         cp -= 9;
	else                 i = rint(d),            cp -= 10;

	// integer i is always 6 digits, so print it 2 at a time.
	static const char kputuw_dig2r[] =
		"00010203040506070809"
		"10111213141516171819"
		"20212223242526272829"
		"30313233343536373839"
		"40414243444546474849"
		"50515253545556575859"
		"60616263646566676869"
		"70717273747576777879"
		"80818283848586878889"
		"90919293949596979899";

	memcpy(cp-=2, &kputuw_dig2r[2*(i%100)], 2); i /= 100;
	memcpy(cp-=2, &kputuw_dig2r[2*(i%100)], 2); i /= 100;
	memcpy(cp-=2, &kputuw_dig2r[2*(i%100)], 2);

	// Except when it rounds up (d=0.009999999 is i=1000000)
	if (i >= 100)
		*--cp = '0' + (i/100);


	int p = buf+20-cp;
	if (p <= 10) { /* d < 1 */
		// 0.00123 is 123, so add leading zeros and 0.
		ep = cp+5; // 6 precision
		while (p < 10) { // aka d < 1
			*--cp = '0';
			p++;
		}
		*--cp = '.';
		*--cp = '0';
	} else {
		// 123.001 is 123001 with p==13, so move 123 down and add "."
		// Equiv to memmove(cp-1, cp, p-10); cp--;
		char *xp = --cp;
		ep = cp+6;
		while (p > 10) {
			xp[0] = xp[1];
			xp++;
			p--;
		}
		xp[0] = '.';
	}

	// Cull trailing zeros
	while (*ep == '0' && ep > cp)
		ep--;

	// End can be 1 out due to the mostly-6 but occasionally 7 (i==1) case.
	// Also code with "123." which should be "123"
	if (*ep && *ep != '.')
		ep++;
	*ep = 0;

	int sl = ep-cp;
	len += sl;
	kputsn(cp, sl, s);
	return len;
}

int kvsprintf(kstring_t *s, const char *fmt, va_list ap)
{
	va_list args;
	int l;
	va_copy(args, ap);

	if (fmt[0] == '%' && fmt[1] == 'g' && fmt[2] == 0) {
		double d = va_arg(args, double);
		l = kputd(d, s);
		va_end(args);
		return l;
	}

	if (!s->s) {
		const size_t sz = 64;
		s->s =  (char*)malloc(sz);
		if (!s->s)
			return -1;
		s->m = sz;
		s->l = 0;
	}

	l = vsnprintf(s->s + s->l, s->m - s->l, fmt, args); // This line does not work with glibc 2.0. See `man snprintf'.
	va_end(args);
	if (l + 1 > s->m - s->l) {
		if (ks_resize(s, s->l + l + 2) < 0)
			return -1;
		va_copy(args, ap);
		l = vsnprintf(s->s + s->l, s->m - s->l, fmt, args);
		va_end(args);
	}
	s->l += l;
	return l;
}

int ksprintf(kstring_t *s, const char *fmt, ...)
{
	va_list ap;
	int l;
	va_start(ap, fmt);
	l = kvsprintf(s, fmt, ap);
	va_end(ap);
	return l;
}

// With gcc, -O3 or -ftree-loop-vectorize is really key here as otherwise
// this code isn't vectorised and runs far slower than is necessary (even
// with the restrict keyword being used).
static inline void HTS_OPT3
add33(uint8_t *a, const uint8_t * b, int32_t len) {
    uint32_t i;
    for (i = 0; i < len; i++)
        a[i] = b[i]+33;
}

static inline void nibble2base(uint8_t *nib, char *seq, int len) {
    static const char code2base[513] =
        "===A=C=M=G=R=S=V=T=W=Y=H=K=D=B=N"
        "A=AAACAMAGARASAVATAWAYAHAKADABAN"
        "C=CACCCMCGCRCSCVCTCWCYCHCKCDCBCN"
        "M=MAMCMMMGMRMSMVMTMWMYMHMKMDMBMN"
        "G=GAGCGMGGGRGSGVGTGWGYGHGKGDGBGN"
        "R=RARCRMRGRRRSRVRTRWRYRHRKRDRBRN"
        "S=SASCSMSGSRSSSVSTSWSYSHSKSDSBSN"
        "V=VAVCVMVGVRVSVVVTVWVYVHVKVDVBVN"
        "T=TATCTMTGTRTSTVTTTWTYTHTKTDTBTN"
        "W=WAWCWMWGWRWSWVWTWWWYWHWKWDWBWN"
        "Y=YAYCYMYGYRYSYVYTYWYYYHYKYDYBYN"
        "H=HAHCHMHGHRHSHVHTHWHYHHHKHDHBHN"
        "K=KAKCKMKGKRKSKVKTKWKYKHKKKDKBKN"
        "D=DADCDMDGDRDSDVDTDWDYDHDKDDDBDN"
        "B=BABCBMBGBRBSBVBTBWBYBHBKBDBBBN"
        "N=NANCNMNGNRNSNVNTNWNYNHNKNDNBNN";

    int i, len2 = len/2;
    seq[0] = 0;

    for (i = 0; i < len2; i++)
        // Note size_t cast helps gcc optimiser.
        memcpy(&seq[i*2], &code2base[(size_t)nib[i]*2], 2);

    if ((i *= 2) < len)
        seq[i] = seq_nt16_str[bam_seqi(nib, i)];
}

static int sam_format1_append(const bam_hdr_t *h, const bam1_t *b, kstring_t *str)
{
    int i, r = 0;
    uint8_t *s, *end;
    const bam1_core_t *c = &b->core;

    if (c->l_qname == 0)
        return -1;
    r |= kputsn_(bam_get_qname(b), c->l_qname-1-c->l_extranul, str);
    r |= kputc_('\t', str); // query name
    r |= kputw(c->flag, str); r |= kputc_('\t', str); // flag
    if (c->tid >= 0) { // chr
        r |= kputs(h->target_name[c->tid] , str);
        r |= kputc_('\t', str);
    } else r |= kputsn_("*\t", 2, str);
    r |= kputll(c->pos + 1, str); r |= kputc_('\t', str); // pos
    r |= kputw(c->qual, str); r |= kputc_('\t', str); // qual
    if (c->n_cigar) { // cigar
        uint32_t *cigar = bam_get_cigar(b);
        for (i = 0; i < c->n_cigar; ++i) {
            r |= kputw(bam_cigar_oplen(cigar[i]), str);
            r |= kputc_(bam_cigar_opchr(cigar[i]), str);
        }
    } else r |= kputc_('*', str);
    r |= kputc_('\t', str);
    if (c->mtid < 0) r |= kputsn_("*\t", 2, str); // mate chr
    else if (c->mtid == c->tid) r |= kputsn_("=\t", 2, str);
    else {
        r |= kputs(h->target_name[c->mtid], str);
        r |= kputc_('\t', str);
    }
    r |= kputll(c->mpos + 1, str); r |= kputc_('\t', str); // mate pos
    r |= kputll(c->isize, str); r |= kputc_('\t', str); // template len
    if (c->l_qseq) { // seq and qual
        uint8_t *s = bam_get_seq(b);
        if (ks_resize(str, str->l+2+2*c->l_qseq) < 0) goto mem_err;
        char *cp = str->s + str->l;

        // Sequence, 2 bases at a time
        nibble2base(s, cp, c->l_qseq);
        cp[c->l_qseq] = '\t';
        cp += c->l_qseq+1;

        // Quality
        s = bam_get_qual(b);
        i = 0;
        if (s[0] == 0xff) {
            cp[i++] = '*';
        } else {
            add33((uint8_t *)cp, s, c->l_qseq); // cp[i] = s[i]+33;
            i = c->l_qseq;
        }
        cp[i] = 0;
        cp += i;
        str->l = cp - str->s;
    } else r |= kputsn_("*\t*", 3, str);

    s = bam_get_aux(b); // aux
    end = b->data + b->l_data;

    while (end - s >= 4) {
        r |= kputc_('\t', str);
        if ((s = (uint8_t *)sam_format_aux1(s, s[2], s+3, end, str)) == NULL)
            goto bad_aux;
    }
    r |= kputsn("", 0, str); // nul terminate
    if (r < 0) goto mem_err;

    return str->l;

 bad_aux:
    hts_log_error("Corrupted aux data for read %.*s flag %d",
                  b->core.l_qname, bam_get_qname(b), b->core.flag);
    errno = EINVAL;
    return -1;

 mem_err:
    hts_log_error("Out of memory");
    errno = ENOMEM;
    return -1;
}

int sam_format1(const bam_hdr_t *h, const bam1_t *b, kstring_t *str)
{
    str->l = 0;
    return sam_format1_append(h, b, str);
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

    //printf("decompressfunc started with _PEN = %d\n", _PEN);

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
    //print_bam_block(comp);
    //print_bam_block(un_comp);
    //printf("Complete the decompression!\n");

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

    //printf("Complete parsing with %d bam1_t!\n", count);

    para->n_records = count;
    para->status = 0;


}

extern "C" void sam_format(void *arg) {
    SamFormatBatch *batch = (SamFormatBatch *)arg;
    int cid = _PEN;

    int start = cid * BATCH_PER_CORE;
    int end   = start + BATCH_PER_CORE;

    if (start >= batch->count) return;
    if (end > batch->count) end = batch->count;

    for (int i = start; i < end; i++) {
        kstring_t *ks = &batch->sam_lines[i];
        ks->l = 0;

        sam_format1(batch->hdr, batch->bams[i], ks);

        kputc('\n', ks);
    }
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

int rabbit_bgzf_compress(void *_dst, size_t *dlen, const void *src, size_t slen, int level) {
    //写入一个EOF块，但是从核这里用不到，可以在主核写入
    if (slen == 0) {
        // EOF block
        if (*dlen < 28) return -1;
        memcpy(_dst, "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", 28);
        *dlen = 28;
        return 0;
    }

    // const uintptr_t ALIGN_MASK = 63;
    // if (((uintptr_t)src & ALIGN_MASK) != 0 || 
    //     ((uintptr_t)_dst & ALIGN_MASK) != 0) {
    //         printf("NOT ALIGN!\n" );
    //     return -1;
    // }

    uint8_t *dst = (uint8_t *) _dst;

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
    int compress_level = 1;    // 默认压缩等级为1
    Comp_Para* para = &paras[id];

    // 忽略空任务
    if (para->status != 0 || para->input_records == nullptr || para->n_records == 0) return;

    bam_block* uncompressed = para->un_comp_block;
    //print_bam_block(uncompressed);
    bam_block* compressed = para->output_block;
    //print_bam_block(compressed);
    //这里只用到了data和pos两个字段
    uncompressed->pos = 0;     //压缩时看这个作为实际的数据长度
    uncompressed->length = 0;  
    uncompressed->errcode = 0;
    uncompressed->block_id = 0;
    uncompressed->block_address = 0;


    //1. 将 bam1_t 记录 逐个解析并写入 到 uncompressed 中
    for (int i = 0; i < para->n_records; i++) {
        bam1_t* b = para->input_records[i];
        //print_bam1(b);
        writeBam1_to_block(uncompressed, b , 0);
    }
    //print_bam_block(uncompressed);
    //printf("Complete parsing with %d bam1_t!\n", para->n_records);

    //2. 对 uncompressed 进行压缩
    compressed->length = block_encode_func(uncompressed, compressed , compress_level); 
    //print_bam_block(uncompressed);
    //print_bam_block(compressed);
    //printf("Complete the compression!\n");

    para->output_size = compressed->length;
    para->status = 0; // success

}


extern "C" void slave_sam_parse(void *arg) {

    SamFormatBatch *batch = (SamFormatBatch *)arg;
    int tid = _PEN;  // 0~63

    //均匀划分
    for (int i = tid; i < batch->count; i += 64) {
        kstring_t *ks = &batch->sam_lines[i];

        sam_parse1(ks, (sam_hdr_t *)batch->hdr, batch->bams[i]);
        ks->l = 0;
    }
}

