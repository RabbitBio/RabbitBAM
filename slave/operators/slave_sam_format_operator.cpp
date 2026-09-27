#include "BamTools.h"
#include <climits>
#include <cstring>
#include <htslib/hts_endian.h>
#include "sam_parse.h"
#include <cmath>

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

int sam_format1_append(const bam_hdr_t *h, const bam1_t *b, kstring_t *str)
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

extern "C" void sam_format(void *arg) {
    SamFormatBatch *batch = (SamFormatBatch *)arg;
    int cid = _PEN;

    kstring_t *ks_out = &batch->core_out_lines[cid];
    char *start_s = ks_out->s;
    size_t start_m = ks_out->m;
    ks_out->l = 0;
    batch->status[cid] = 0;
    batch->required_capacity[cid] = start_m;

    int total_tasks = batch->count;
    if (total_tasks == 0) return;

    int start, end;
    int base_tasks = total_tasks >> 6;
    int remainder = total_tasks & 63;

    if (cid < remainder) {
        start = cid * (base_tasks + 1);
        end = start + base_tasks + 1;
    } else {
        start = remainder + cid * base_tasks;
        end = start + base_tasks;
    }

    if (start >= total_tasks || start >= end) return;

    for (int i = start; i < end; i++) {
        sam_format1_append(batch->hdr, batch->bams[i], ks_out);
        kputc('\n', ks_out);
        if (ks_out->s != start_s || ks_out->m > start_m) {
            batch->status[cid] = BOUNDS_LIMIT_MAX_SAM_FORMAT_CORE_BUFFER_SIZE;
            batch->required_capacity[cid] = ks_out->m;
            return;
        }
    }
}

