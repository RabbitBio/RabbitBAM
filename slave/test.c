#include "test.h"

//test用------------------------------------------------------------------------------------------------------------
// int test_add(int a, int b)
// {
//     return a + b;
// }

// void test_print(const char *msg)
// {
//     printf("from C: %s\n", msg);
// }


//slave_sam_parse需要的------------------------------------------------------------------------------------------------------------

int sam_realloc_bam_data(bam1_t *b, size_t desired)
{
    uint32_t new_m_data;
    uint8_t *new_data;
    new_m_data = desired;
    kroundup32(new_m_data);
    if (new_m_data < desired) {
        errno = ENOMEM; // Not strictly true but we can't store the size
        return -1;
    }
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (new_m_data > FUZZ_ALLOC_LIMIT) {
        errno = ENOMEM;
        return -1;
    }
#endif
    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) == 0) {
        new_data = realloc(b->data, new_m_data);
    } else {
        if ((new_data = malloc(new_m_data)) != NULL) {
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

static inline int possibly_expand_bam_data(bam1_t *b, size_t bytes) {
    size_t new_len = (size_t) b->l_data + bytes;

    if (new_len > INT32_MAX || new_len < bytes) { // Too big or overflow
        errno = ENOMEM;
        return -1;
    }
    if (new_len <= b->m_data) return 0;
    return sam_realloc_bam_data(b, new_len);
}


static inline int64_t hts_str2int(const char *in, char **end, int bits,
                                    int *failed) {
    uint64_t n = 0, limit = (1ULL << (bits - 1)) - 1;
    uint32_t fast = (bits - 1) * 1000 / 3322 + 1; // log(10)/log(2) ~= 3.322
    const unsigned char *v = (const unsigned char *) in;
    const unsigned int ascii_zero = '0'; // Prevents conversion to signed
    unsigned int d;

    int neg;
    switch(*v) {
    case '-':
        limit++;
        neg=1;
        v++;
        // See "dup" comment below
        while (--fast && *v>='0' && *v<='9')
            n = n*10 + *v++ - ascii_zero;
        break;

    case '+':
        v++;
        // fall through

    default:
        neg = 0;
        // dup of above.  This is somewhat unstable and mainly for code
        // size cheats to prevent instruction cache lines spanning 32-byte
        // blocks in the sam_parse_B_vals calling code.  It's been tested
        // on gcc7, gcc13, clang10 and clang16 with -O2 and -O3.  While
        // not exhaustive, this code duplication gives stable fast results
        // while a single copy does not.
        // (NB: system was "seq4d", so quite old)
        while (--fast && *v>='0' && *v<='9')
            n = n*10 + *v++ - ascii_zero;
        break;
    }

    // NB gcc7 is slow with (unsigned)(*v - ascii_zero) < 10,
    // while gcc13 prefers it.
    if (*v>='0' && !fast) { // rejects ',' and tab
        uint64_t limit_d_10 = limit / 10;
        uint64_t limit_m_10 = limit - 10 * limit_d_10;
        while ((d = *v - ascii_zero) < 10) {
            if (n < limit_d_10 || (n == limit_d_10 && d <= limit_m_10)) {
                n = n*10 + d;
                v++;
            } else {
                do { v++; } while (*v - ascii_zero < 10);
                n = limit;
                *failed = 1;
                break;
            }
        }
    }

    *end = (char *)v;

    return neg ? (int64_t)-n : (int64_t)n;
}

static inline int isprint_c(char c) { return isprint((unsigned char) c); }

/* Utility function for printing possibly malicious text data
 */
const char *
hts_strprint(char *buf, size_t buflen, char quote, const char *s, size_t len)
{
    const char *slim = (len < SIZE_MAX)? &s[len] : NULL;
    char *t = buf, *bufend = buf + buflen;

    size_t qlen = quote? 1 : 0;
    if (quote) *t++ = quote;

    for (; slim? (s < slim) : (*s); s++) {
        char c;
        size_t clen;
        switch (*s) {
        case '\n': c = 'n'; clen = 2; break;
        case '\r': c = 'r'; clen = 2; break;
        case '\t': c = 't'; clen = 2; break;
        case '\0': c = '0'; clen = 2; break;
        case '\\': c = '\\'; clen = 2; break;
        default:
            c = *s;
            if (c == quote) clen = 2;
            else clen = isprint_c(c)? 1 : 4;
            break;
        }

        if (t-buf + clen + qlen >= buflen) {
            while (t-buf + 3 + qlen >= buflen) t--;
            if (quote) *t++ = quote;
            strcpy(t, "...");
            return buf;
        }

        if (clen == 4) {
            snprintf(t, bufend - t, "\\x%02X", (unsigned char) c);
            t += clen;
        }
        else {
            if (clen == 2) *t++ = '\\';
            *t++ = c;
        }
    }

    if (quote) *t++ = quote;
    *t = '\0';
    return buf;
}

static inline uint64_t hts_str2uint(const char *in, char **end, int bits,
                                    int *failed) {
    uint64_t n = 0, limit = (bits < 64 ? (1ULL << bits) : 0) - 1;
    const unsigned char *v = (const unsigned char *) in;
    const unsigned int ascii_zero = '0'; // Prevents conversion to signed
    uint32_t fast = bits * 1000 / 3322 + 1; // log(10)/log(2) ~= 3.322
    unsigned int d;

    if (*v == '+')
        v++;

    while (--fast && *v>='0' && *v<='9')
        n = n*10 + *v++ - ascii_zero;

    if ((unsigned)(*v - ascii_zero) < 10 && !fast) {
        uint64_t limit_d_10 = limit / 10;
        uint64_t limit_m_10 = limit - 10 * limit_d_10;
        while ((d = *v - ascii_zero) < 10) {
            if (n < limit_d_10 || (n == limit_d_10 && d <= limit_m_10)) {
                n = n*10 + d;
                v++;
            } else {
                do { v++; } while (*v - ascii_zero < 10);
                n = limit;
                *failed = 1;
                break;
            }
        }
    }

    *end = (char *)v;
    return n;
}

static inline unsigned int parse_sam_flag(char *v, char **rv, int *overflow) {
    if (*v >= '1' && *v <= '9') {
        return hts_str2uint(v, rv, 16, overflow);
    }
    else if (*v == '0') {
        // handle single-digit "0" directly; otherwise it's hex or octal
        if (v[1] == '\t') { *rv = v+1; return 0; }
        else {
            unsigned long val = strtoul(v, rv, 0);
            if (val > 65535) { *overflow = 1; return 65535; }
            return val;
        }
    }
    else {
        // TODO implement symbolic flag letters
        *rv = v;
        return 0;
    }
}

#define skip_to_comma_(q) do { while (*(q) > '\t' && *(q) != ',') (q)++; } while (0)

static inline int64_t grow_B_array(bam1_t *b, uint32_t *n, size_t size) {
    // Avoid overflow on 32-bit platforms, but it breaks BAM anyway
    if (*n > INT32_MAX*0.666) {
        errno = ENOMEM;
        return -1;
    }

    size_t bytes = (size_t)size * (size_t)(*n>>1);
    if (possibly_expand_bam_data(b, bytes) < 0) {
        hts_log_error("Out of memory");
        return -1;
    }

    (*n)+=*n>>1;
    return 0;
}


HTS_ALIGN32
static char *sam_parse_Bc_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 1) < 0)
                return NULL;
        }
        *(b->data + b->l_data) = hts_str2int(q + 1, &q, 8, overflow);
        b->l_data++;
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_BC_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 1) < 0)
                return NULL;
        }
        if (q[1] != '-') {
            *(b->data + b->l_data) = hts_str2uint(q + 1, &q, 8, overflow);
            b->l_data++;
        } else {
            *overflow = 1;
            q++;
            skip_to_comma_(q);
        }
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_Bs_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 2) < 0)
                return NULL;
        }
        i16_to_le(hts_str2int(q + 1, &q, 16, overflow),
                  b->data + b->l_data);
        b->l_data += 2;
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_BS_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 2) < 0)
                return NULL;
        }
        if (q[1] != '-') {
            u16_to_le(hts_str2uint(q + 1, &q, 16, overflow),
                      b->data + b->l_data);
            b->l_data += 2;
        } else {
            *overflow = 1;
            q++;
            skip_to_comma_(q);
        }
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_Bi_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 4) < 0)
                return NULL;
        }
        i32_to_le(hts_str2int(q + 1, &q, 32, overflow),
                  b->data + b->l_data);
        b->l_data += 4;
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_BI_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 4) < 0)
                return NULL;
        }
        if (q[1] != '-') {
            u32_to_le(hts_str2uint(q + 1, &q, 32, overflow),
                      b->data + b->l_data);
            b->l_data += 4;
        } else {
            *overflow = 1;
            q++;
            skip_to_comma_(q);
        }
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_Bf_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 4) < 0)
                return NULL;
        }
        float_to_le(strtod(q + 1, &q), b->data + b->l_data);
        b->l_data += 4;
    }
    return q;
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

HTS_ALIGN32
static int sam_parse_B_vals_r(char type, uint32_t nalloc, char *in,
                              char **end, bam1_t *b,
                              int *ctr) {
    // Protect against infinite recursion when dealing with invalid input.
    // An example string is "XX:B:C,-".  The lack of a number means min=0,
    // but it overflowed due to "-" and so we repeat ad-infinitum.
    //
    // Loop detection is the safest solution incase there are other
    // strange corner cases with malformed inputs.
    if (++(*ctr) > 2) {
        hts_log_error("Malformed data in B:%c array", type);
        return -1;
    }

    int orig_l = b->l_data;
    char *q = in;
    int32_t size;
    size_t bytes;
    int overflow = 0;

    size = aux_type2size(type);
    if (size <= 0 || size > 4) {
        hts_log_error("Unrecognized type B:%c", type);
        return -1;
    }

    // Ensure space for type + values.
    // The first pass through here we don't know the number of entries and
    // nalloc == 0.  We start with a small working set and then parse the
    // data, growing as needed.
    //
    // If we have a second pass through we do know the number of entries
    // and nalloc is already known.  We have no need to expand the bam data.
    if (!nalloc)
         nalloc=7;

    // Ensure allocated memory is big enough (for current nalloc estimate)
    bytes = (size_t) nalloc * (size_t) size;
    if (bytes / size != nalloc
        || possibly_expand_bam_data(b, bytes + 2 + sizeof(uint32_t))) {
        hts_log_error("Out of memory");
        return -1;
    }

    uint32_t nused = 0;

    b->data[b->l_data++] = 'B';
    b->data[b->l_data++] = type;
    // 32-bit B-array length is inserted later once we know it.
    int b_len_idx = b->l_data;
    b->l_data += sizeof(uint32_t);

    if (type == 'c') {
        if (!(q = sam_parse_Bc_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'C') {
        if (!(q = sam_parse_BC_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 's') {
        if (!(q = sam_parse_Bs_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'S') {
        if (!(q = sam_parse_BS_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'i') {
        if (!(q = sam_parse_Bi_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'I') {
        if (!(q = sam_parse_BI_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'f') {
        if (!(q = sam_parse_Bf_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    }
    if (*q != '\t' && *q != '\0') {
        // Unknown B array type or junk in the numbers
        hts_log_error("Malformed B:%c", type);
        return -1;
    }
    i32_to_le(nused, b->data + b_len_idx);

    if (!overflow) {
        *end = q;
        return 0;
    } else {
        int64_t max = 0, min = 0, val;
        // Given type was incorrect.  Try to rescue the situation.
        char *r = q;
        q = in;
        overflow = 0;
        b->l_data = orig_l;
        // Find out what range of values is present
        while (q < r) {
            val = hts_str2int(q + 1, &q, 64, &overflow);
            if (max < val) max = val;
            if (min > val) min = val;
            skip_to_comma_(q);
        }
        // Retry with appropriate type
        if (!overflow) {
            if (min < 0) {
                if (min >= INT8_MIN && max <= INT8_MAX) {
                    return sam_parse_B_vals_r('c', nalloc, in, end, b, ctr);
                } else if (min >= INT16_MIN && max <= INT16_MAX) {
                    return sam_parse_B_vals_r('s', nalloc, in, end, b, ctr);
                } else if (min >= INT32_MIN && max <= INT32_MAX) {
                    return sam_parse_B_vals_r('i', nalloc, in, end, b, ctr);
                }
            } else {
                if (max < UINT8_MAX) {
                    return sam_parse_B_vals_r('C', nalloc, in, end, b, ctr);
                } else if (max <= UINT16_MAX) {
                    return sam_parse_B_vals_r('S', nalloc, in, end, b, ctr);
                } else if (max <= UINT32_MAX) {
                    return sam_parse_B_vals_r('I', nalloc, in, end, b, ctr);
                }
            }
        }
        // If here then at least one of the values is too big to store
        hts_log_error("Numeric value in B array out of allowed range");
        return -1;
    }
#undef skip_to_comma_
}

HTS_ALIGN32
static int sam_parse_B_vals(char type, char *in, char **end, bam1_t *b)
{
    int ctr = 0;
    uint32_t nalloc = 0;
    return sam_parse_B_vals_r(type, nalloc, in, end, b, &ctr);
}

static inline int isspace_c(char c) { return isspace((unsigned char) c); }

KHASH_SET_INIT_INT(tag)

// Parse tag line and append to bam object b.
// Shared by both SAM and FASTQ parsers.
//
// The difference between the two is how lenient we are to recognising
// non-compliant strings.  The FASTQ parser glosses over arbitrary
// non-SAM looking strings.
static inline int aux_parse(char *start, char *end, bam1_t *b, int lenient,
                            khash_t(tag) *tag_whitelist) {
    int overflow = 0;
    int checkpoint;
    char logbuf[40];
    char *q = start, *p = end;

#define _parse_err(cond, ...)                   \
    do {                                        \
        if (cond) {                             \
            if (lenient) {                      \
                while (q < p && !isspace_c(*q))   \
                    q++;                        \
                while (q < p && isspace_c(*q))    \
                    q++;                        \
                b->l_data = checkpoint;         \
                goto loop;                      \
            } else {                            \
                hts_log_error(__VA_ARGS__);     \
                goto err_ret;                   \
            }                                   \
        }                                       \
    } while (0)

    while (q < p) loop: {
        char type;
        checkpoint = b->l_data;
        if (p - q < 5) {
            if (lenient) {
                break;
            } else {
                hts_log_error("Incomplete aux field");
                goto err_ret;
            }
        }
        _parse_err(q[0] < '!' || q[1] < '!', "invalid aux tag id");

        if (lenient && (q[2] | q[4]) != ':') {
            while (q < p && !isspace_c(*q))
                q++;
            while (q < p && isspace_c(*q))
                q++;
            continue;
        }

        if (tag_whitelist) {
            int tt = q[0]*256 + q[1];
            if (kh_get(tag, tag_whitelist, tt) == kh_end(tag_whitelist)) {
                while (q < p && *q != '\t')
                    q++;
                continue;
            }
        }

        // Copy over id
        if (possibly_expand_bam_data(b, 2) < 0) goto err_ret;
        memcpy(b->data + b->l_data, q, 2); b->l_data += 2;
        q += 3; type = *q++; ++q; // q points to value
        if (type != 'Z' && type != 'H') // the only zero length acceptable fields
            _parse_err(*q <= '\t', "incomplete aux field");

        // Ensure enough space for a double + type allocated.
        if (possibly_expand_bam_data(b, 16) < 0) goto err_ret;

        if (type == 'A' || type == 'a' || type == 'c' || type == 'C') {
            b->data[b->l_data++] = 'A';
            b->data[b->l_data++] = *q++;
        } else if (type == 'i' || type == 'I') {
            if (*q == '-') {
                int32_t x = hts_str2int(q, &q, 32, &overflow);
                if (x >= INT8_MIN) {
                    b->data[b->l_data++] = 'c';
                    b->data[b->l_data++] = x;
                } else if (x >= INT16_MIN) {
                    b->data[b->l_data++] = 's';
                    i16_to_le(x, b->data + b->l_data);
                    b->l_data += 2;
                } else {
                    b->data[b->l_data++] = 'i';
                    i32_to_le(x, b->data + b->l_data);
                    b->l_data += 4;
                }
            } else {
                uint32_t x = hts_str2uint(q, &q, 32, &overflow);
                if (x <= UINT8_MAX) {
                    b->data[b->l_data++] = 'C';
                    b->data[b->l_data++] = x;
                } else if (x <= UINT16_MAX) {
                    b->data[b->l_data++] = 'S';
                    u16_to_le(x, b->data + b->l_data);
                    b->l_data += 2;
                } else {
                    b->data[b->l_data++] = 'I';
                    u32_to_le(x, b->data + b->l_data);
                    b->l_data += 4;
                }
            }
        } else if (type == 'f') {
            b->data[b->l_data++] = 'f';
            float_to_le(strtod(q, &q), b->data + b->l_data);
            b->l_data += sizeof(float);
        } else if (type == 'd') {
            b->data[b->l_data++] = 'd';
            double_to_le(strtod(q, &q), b->data + b->l_data);
            b->l_data += sizeof(double);
        } else if (type == 'Z' || type == 'H') {
            char *end = strchr(q, '\t');
            if (!end) end = q + strlen(q);
            _parse_err(type == 'H' && ((end-q)&1) != 0,
                       "hex field does not have an even number of digits");
            b->data[b->l_data++] = type;
            if (possibly_expand_bam_data(b, end - q + 1) < 0) goto err_ret;
            memcpy(b->data + b->l_data, q, end - q);
            b->l_data += end - q;
            b->data[b->l_data++] = '\0';
            q = end;
        } else if (type == 'B') {
            type = *q++; // q points to the first ',' following the typing byte
            _parse_err(*q && *q != ',' && *q != '\t',
                       "B aux field type not followed by ','");

            if (sam_parse_B_vals(type, q, &q, b) < 0)
                goto err_ret;
        } else _parse_err(1, "unrecognized type %s", hts_strprint(logbuf, sizeof logbuf, '\'', &type, 1));

        while (*q > '\t') { q++; } // Skip any junk to next tab
        q++;
    }

    _parse_err(!lenient && overflow != 0, "numeric value out of allowed range");
#undef _parse_err

    return 0;

err_ret:
    return -2;
}

static inline int isdigit_c(char c) { return isdigit((unsigned char) c); }

static uint32_t read_ncigar(const char *q) {
    uint32_t n_cigar = 0;
    for (; *q && *q != '\t'; ++q)
        if (!isdigit_c(*q)) ++n_cigar;
    if (!n_cigar) {
        hts_log_error("No CIGAR operations");
        return 0;
    }
    if (n_cigar >= 2147483647) {
        hts_log_error("Too many CIGAR operations");
        return 0;
    }

    return n_cigar;
}

static int parse_cigar(const char *in, uint32_t *a_cigar, uint32_t n_cigar) {
    int i, overflow = 0;
    const char *p = in;
    for (i = 0; i < n_cigar; i++) {
        uint32_t len;
        int op;
        char *q;
        len = hts_str2uint(p, &q, 28, &overflow)<<BAM_CIGAR_SHIFT;
        if (q == p) {
            hts_log_error("CIGAR length invalid at position %d (%s)", (int)(i+1), p);
            return 0;
        }
        if (overflow) {
            hts_log_error("CIGAR length too long at position %d (%.*s)", (int)(i+1), (int)(q-p+1), p);
            return 0;
        }
        p = q;
        op = bam_cigar_table[(unsigned char)*p++];
        if (op < 0) {
            hts_log_error("Unrecognized CIGAR operator");
            return 0;
        }
        a_cigar[i] = len;
        a_cigar[i] |= op;
    }

    return p-in;
}

ssize_t bam_parse_cigar(const char *in, char **end, bam1_t *b) {
    size_t n_cigar = 0;
    int diff;

    if (!in || !b) {
        hts_log_error("NULL pointer arguments");
        return -1;
    }
    if (end) *end = (char *)in;

    n_cigar = (*in == '*') ? 0 : read_ncigar(in);
    if (!n_cigar && b->core.n_cigar == 0) {
        if (end) *end = (char *)in+1;
        return 0;
    }

    ssize_t cig_diff = n_cigar - b->core.n_cigar;
    if (cig_diff > 0 &&
        possibly_expand_bam_data(b, cig_diff * sizeof(uint32_t)) < 0) {
        hts_log_error("Memory allocation error");
        return -1;
    }

    uint32_t *cig = bam_get_cigar(b);
    if ((uint8_t *)cig != b->data + b->l_data) {
        // Modifying an BAM existing BAM record
        uint8_t  *seq = bam_get_seq(b);
        memmove(cig + n_cigar, seq, (b->data + b->l_data) - seq);
    }

    if (n_cigar) {
        if (!(diff = parse_cigar(in, cig, n_cigar)))
            return -1;
    } else {
        diff = 1; // handle "*"
    }

    b->l_data += cig_diff * sizeof(uint32_t);
    b->core.n_cigar = n_cigar;
    if (end) *end = (char *)in + diff;

    return n_cigar;
}

hts_pos_t bam_cigar2qlen(int n_cigar, const uint32_t *cigar)
{
    int k;
    hts_pos_t l;
    for (k = l = 0; k < n_cigar; ++k)
        if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
            l += bam_cigar_oplen(cigar[k]);
    return l;
}


typedef struct {
    char *str;
    size_t used;
} string_t;

typedef struct {
    size_t max_length;
    size_t nstrings;
    size_t max_strings;
    string_t *strings;
} string_alloc_t;


typedef struct {
    void   *pool;
    size_t  used;
} pool_t;

typedef struct {
    size_t dsize;
    size_t psize;
    size_t npools;
    pool_t *pools;
    void *free;
} pool_alloc_t;

typedef struct sam_hrec_tag_s {
    struct sam_hrec_tag_s *next;
    const char *str;
    int   len;
} sam_hrec_tag_t;
typedef struct sam_hrec_type_s {
    struct sam_hrec_type_s *next; // circular list of this type
    struct sam_hrec_type_s *prev; // circular list of this type
    struct sam_hrec_type_s *global_next; // circular list of all lines
    struct sam_hrec_type_s *global_prev; // circular list of all lines
    sam_hrec_tag_t *tag;          // first tag
    khint32_t type;               // Two-letter type code as an int
} sam_hrec_type_t;

KHASH_MAP_INIT_INT(sam_hrecs_t, sam_hrec_type_t*)
KHASH_MAP_INIT_STR(m_s2i, int)

// Hash table for removing multiple lines from the header
KHASH_SET_INIT_STR(rm)
// Used for long refs in SAM files

//KHASH_DECLARE(s2i, kh_cstr_t, int64_t)
KHASH_INIT2(s2i,, kh_cstr_t, int64_t, 1, kh_str_hash_func, kh_str_hash_equal)

typedef khash_t(rm) rmhash_t;

/*! Parsed \@SQ lines */
typedef struct {
    const char *name;
    hts_pos_t len;
    sam_hrec_type_t *ty;
} sam_hrec_sq_t;

typedef struct {
    const char *name;
    sam_hrec_type_t *ty;
    int name_len;
    int id;           // numerical ID
} sam_hrec_rg_t;

typedef struct {
    const char *name;
    sam_hrec_type_t *ty;
    int name_len;
    int id;           // numerical ID
    int prev_id;      // -1 if none
} sam_hrec_pg_t;

struct sam_hrecs_t {
    khash_t(sam_hrecs_t) *h;
    sam_hrec_type_t *first_line; //!< First line (usually @HD)
    string_alloc_t *str_pool; //!< Pool of sam_hdr_tag->str strings
    pool_alloc_t   *type_pool;//!< Pool of sam_hdr_type structs
    pool_alloc_t   *tag_pool; //!< Pool of sam_hdr_tag structs

    // @SQ lines / references
    int nref;                  //!< Number of \@SQ lines
    int ref_sz;                //!< Number of entries available in ref[]
    sam_hrec_sq_t *ref;        //!< Array of parsed \@SQ lines
    khash_t(m_s2i) *ref_hash;  //!< Maps SQ SN field to ref[] index

    // @RG lines / read-groups
    int nrg;                   //!< Number of \@RG lines
    int rg_sz;                 //!< number of entries available in rg[]
    sam_hrec_rg_t *rg;         //!< Array of parsed \@RG lines
    khash_t(m_s2i) *rg_hash;   //!< Maps RG ID field to rg[] index

    // @PG lines / programs
    int npg;                   //!< Number of \@PG lines
    int pg_sz;                //!< Number of entries available in pg[]
    int npg_end;               //!< Number of terminating \@PG lines
    int npg_end_alloc;         //!< Size of pg_end field
    sam_hrec_pg_t *pg;         //!< Array of parsed \@PG lines
    khash_t(m_s2i) *pg_hash;   //!< Maps PG ID field to pg[] index
    int *pg_end;               //!< \@PG chain termination IDs

    // @cond internal
    char *ID_buf;             // temporary buffer for sam_hdr_pg_id
    uint32_t ID_buf_sz;
    int ID_cnt;
    // @endcond

    int dirty;                // marks the header as modified, so it can be rebuilt
    int refs_changed;         // Index of first changed ref (-1 if unchanged)
    int pgs_changed;          // New PG line added
    int type_count;
    char (*type_order)[3];
};

#define PSIZE 1024*1024

static int next_power_2(unsigned int v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;

    return v;
}

#define MIN(a, b)		((a) <= (b) ? (a) : (b))

pool_alloc_t *pool_create(size_t dsize) {
    pool_alloc_t *p;

    if (NULL == (p = (pool_alloc_t *)malloc(sizeof(*p))))
        return NULL;

    /* Minimum size is a pointer, for free list */
    dsize = (dsize + sizeof(void *) - 1) & ~(sizeof(void *)-1);
    if (dsize < sizeof(void *))
        dsize = sizeof(void *);
    p->dsize = dsize;
    p->psize = MIN(PSIZE, next_power_2(p->dsize*1024));

    p->npools = 0;
    p->pools = NULL;
    p->free  = NULL;

    return p;
}

#define MIN_STR_SIZE 1024

string_alloc_t *string_pool_create(size_t max_length) {
    string_alloc_t *a_str;

    if (NULL == (a_str = (string_alloc_t *)malloc(sizeof(*a_str)))) {
        return NULL;
    }

    if (max_length < MIN_STR_SIZE) max_length = MIN_STR_SIZE;

    a_str->nstrings    = 0;
    a_str->max_strings = 0;
    a_str->max_length  = max_length;
    a_str->strings     = NULL;

    return a_str;
}

static int sam_hrecs_init_type_order(sam_hrecs_t *hrecs, char *type_list) {
    if (!hrecs)
        return -1;

    if (!type_list) {
        hrecs->type_count = 5;
        hrecs->type_order = calloc(hrecs->type_count, 3);
        if (!hrecs->type_order)
            return -1;
        memcpy(hrecs->type_order[0], "HD", 2);
        memcpy(hrecs->type_order[1], "SQ", 2);
        memcpy(hrecs->type_order[2], "RG", 2);
        memcpy(hrecs->type_order[3], "PG", 2);
        memcpy(hrecs->type_order[4], "CO", 2);
    }

    return 0;
}

static void pool_destroy(pool_alloc_t *p) {
    size_t i;

    for (i = 0; i < p->npools; i++) {
        free(p->pools[i].pool);
    }
    free(p->pools);
    free(p);
}

void string_pool_destroy(string_alloc_t *a_str) {
    size_t i;

    for (i = 0; i < a_str->nstrings; i++) {
        free(a_str->strings[i].str);
    }

    free(a_str->strings);
    free(a_str);
}


static string_t *new_string_pool(string_alloc_t *a_str) {
    string_t *str;

    if (a_str->nstrings == a_str->max_strings) {
        size_t new_max = (a_str->max_strings | (a_str->max_strings >> 2)) + 1;
        str = realloc(a_str->strings, new_max * sizeof(*a_str->strings));

        if (NULL == str) return NULL;

        a_str->strings = str;
        a_str->max_strings = new_max;
    }

    str = &a_str->strings[a_str->nstrings];

    str->str = malloc(a_str->max_length);

    if (NULL == str->str) return NULL;

    str->used = 0;
    a_str->nstrings++;

    return str;
}

char *string_alloc(string_alloc_t *a_str, size_t length) {
    string_t *str;
    char *ret;

    if (length <= 0) return NULL;

    // add to last string pool if we have space
    if (a_str->nstrings) {
        str = &a_str->strings[a_str->nstrings - 1];

        if (str->used + length < a_str->max_length) {
            ret = str->str + str->used;
            str->used += length;
            return ret;
        }
    }

    // increase the max length if needs be
    if (length > a_str->max_length) a_str->max_length = length;

    // need a new string pool
    str = new_string_pool(a_str);

    if (NULL == str) return NULL;

    str->used = length;
    return str->str;
}

char *string_ndup(string_alloc_t *a_str, const char *instr, size_t len) {
    char *str = string_alloc(a_str, len + 1);

    if (NULL == str) return NULL;

    memcpy(str, instr, len);
    str[len] = 0;

    return str;
}

char *string_dup(string_alloc_t *a_str, const char *instr) {
    return string_ndup(a_str, instr, strlen(instr));
}


static int sam_hrecs_refs_from_targets_array(sam_hrecs_t *hrecs,
                                             const sam_hdr_t *bh) {
    int32_t tid = 0;

    if (!hrecs || !bh)
        return -1;

    // This should always be called before parsing the text header
    // so the ref array should start off empty, and we don't have to try
    // to reconcile any existing data.
    if (hrecs->nref > 0) {
        hts_log_error("Called with non-empty ref array");
        return -1;
    }

    if (hrecs->ref_sz < bh->n_targets) {
        sam_hrec_sq_t *new_ref = realloc(hrecs->ref,
                                         bh->n_targets * sizeof(*new_ref));
        if (!new_ref)
            return -1;

        hrecs->ref = new_ref;
        hrecs->ref_sz = bh->n_targets;
    }

    for (tid = 0; tid < bh->n_targets; tid++) {
        khint_t k;
        int r;
        hrecs->ref[tid].name = string_dup(hrecs->str_pool, bh->target_name[tid]);
        if (!hrecs->ref[tid].name) goto fail;
        if (bh->target_len[tid] < UINT32_MAX || !bh->sdict) {
            hrecs->ref[tid].len  = bh->target_len[tid];
        } else {
            khash_t(s2i) *long_refs = (khash_t(s2i) *) bh->sdict;
            k = kh_get(s2i, long_refs, hrecs->ref[tid].name);
            if (k < kh_end(long_refs)) {
                hrecs->ref[tid].len = kh_val(long_refs, k);
            } else {
                hrecs->ref[tid].len = UINT32_MAX;
            }
        }
        hrecs->ref[tid].ty   = NULL;
        k = kh_put(m_s2i, hrecs->ref_hash, hrecs->ref[tid].name, &r);
        if (r < 0) goto fail;
        if (r == 0) {
            hts_log_error("Duplicate entry \"%s\" in target list",
                            hrecs->ref[tid].name);
            return -1;
        } else {
            kh_val(hrecs->ref_hash, k) = tid;
        }
    }
    hrecs->nref = bh->n_targets;
    return 0;

 fail: {
        int32_t i;
        hts_log_error("%s", strerror(errno));
        for (i = 0; i < tid; i++) {
            khint_t k;
            if (!hrecs->ref[i].name) continue;
            k = kh_get(m_s2i, hrecs->ref_hash, hrecs->ref[tid].name);
            if (k < kh_end(hrecs->ref_hash)) kh_del(m_s2i, hrecs->ref_hash, k);
        }
        hrecs->nref = 0;
        return -1;
    }
}

sam_hrecs_t *sam_hrecs_new() {
    //sam_hrecs_t *hrecs = calloc(1, sizeof(*hrecs));
    sam_hrecs_t *hrecs = (sam_hrecs_t *)calloc(1, sizeof(*hrecs));


    if (!hrecs)
        return NULL;

    hrecs->h = kh_init(sam_hrecs_t);
    if (!hrecs->h)
        goto err;

    hrecs->ID_cnt = 1;

    hrecs->nref = 0;
    hrecs->ref_sz = 0;
    hrecs->ref  = NULL;
    if (!(hrecs->ref_hash = kh_init(m_s2i)))
        goto err;
    hrecs->refs_changed = -1;

    hrecs->nrg = 0;
    hrecs->rg_sz = 0;
    hrecs->rg  = NULL;
    if (!(hrecs->rg_hash = kh_init(m_s2i)))
        goto err;

    hrecs->npg = 0;
    hrecs->pg_sz = 0;
    hrecs->pg  = NULL;
    hrecs->npg_end = hrecs->npg_end_alloc = 0;
    hrecs->pg_end = NULL;
    if (!(hrecs->pg_hash = kh_init(m_s2i)))
        goto err;

    if (!(hrecs->tag_pool = pool_create(sizeof(sam_hrec_tag_t))))
        goto err;

    if (!(hrecs->type_pool = pool_create(sizeof(sam_hrec_type_t))))
        goto err;

    if (!(hrecs->str_pool = string_pool_create(65536)))
        goto err;

    if (sam_hrecs_init_type_order(hrecs, NULL))
        goto err;

    return hrecs;

err:
    if (hrecs->h)
        kh_destroy(sam_hrecs_t, hrecs->h);

    if (hrecs->tag_pool)
        pool_destroy(hrecs->tag_pool);

    if (hrecs->type_pool)
        pool_destroy(hrecs->type_pool);

    if (hrecs->str_pool)
        string_pool_destroy(hrecs->str_pool);

    free(hrecs);

    return NULL;
}

void sam_hrecs_free(sam_hrecs_t *hrecs) {
    if (!hrecs)
        return;

    if (hrecs->h)
        kh_destroy(sam_hrecs_t, hrecs->h);

    if (hrecs->ref_hash)
        kh_destroy(m_s2i, hrecs->ref_hash);

    if (hrecs->ref)
        free(hrecs->ref);

    if (hrecs->rg_hash)
        kh_destroy(m_s2i, hrecs->rg_hash);

    if (hrecs->rg)
        free(hrecs->rg);

    if (hrecs->pg_hash)
        kh_destroy(m_s2i, hrecs->pg_hash);

    if (hrecs->pg)
        free(hrecs->pg);

    if (hrecs->pg_end)
        free(hrecs->pg_end);

    if (hrecs->type_pool)
        pool_destroy(hrecs->type_pool);

    if (hrecs->tag_pool)
        pool_destroy(hrecs->tag_pool);

    if (hrecs->str_pool)
        string_pool_destroy(hrecs->str_pool);

    if (hrecs->type_order)
        free(hrecs->type_order);

    if (hrecs->ID_buf)
        free(hrecs->ID_buf);

    free(hrecs);
}

#define MAX_ERROR_QUOTE 320 // Prevent over-long error messages
static void sam_hrecs_error(const char *msg, const char *line, size_t len, size_t lno) {
    int j;

    if (len > MAX_ERROR_QUOTE)
        len = MAX_ERROR_QUOTE;
    for (j = 0; j < len && line[j] != '\n'; j++)
        ;
    hts_log_error("%s at line %zd: \"%.*s\"", msg, lno, j, line);
}

static inline int isalpha_c(char c) { return isalpha((unsigned char) c); }

static inline khint32_t TYPEKEY(const char *type) {
    unsigned int u0 = (unsigned char) type[0];
    unsigned int u1 = (unsigned char) type[1];
    return (u0 << 8) | u1;
}

static pool_t *new_pool(pool_alloc_t *p) {
    size_t n = p->psize / p->dsize;
    pool_t *pool;

    pool = realloc(p->pools, (p->npools + 1) * sizeof(*p->pools));
    if (NULL == pool) return NULL;
    p->pools = pool;
    pool = &p->pools[p->npools];

    pool->pool = malloc(n * p->dsize);
    if (NULL == pool->pool) return NULL;

    pool->used = 0;

    p->npools++;

    return pool;
}

static void sam_hrecs_global_list_add(sam_hrecs_t *hrecs,
                                      sam_hrec_type_t *h_type,
                                      sam_hrec_type_t *after) {
    const khint32_t hd_type = TYPEKEY("HD");
    int update_first_line = 0;

    // First line seen
    if (!hrecs->first_line) {
        hrecs->first_line = h_type->global_next = h_type->global_prev = h_type;
        return;
    }

    // @HD goes at the top (unless there's one already)
    if (h_type->type == hd_type && hrecs->first_line->type != hd_type) {
        after = hrecs->first_line->global_prev;
        update_first_line = 1;
    }

    // If no instructions given, put it at the end
    if (!after)
        after = hrecs->first_line->global_prev;

    h_type->global_prev = after;
    h_type->global_next = after->global_next;
    h_type->global_prev->global_next = h_type;
    h_type->global_next->global_prev = h_type;

    if (update_first_line)
        hrecs->first_line = h_type;
}

void *pool_alloc(pool_alloc_t *p) {
    pool_t *pool;
    void *ret;

    /* Look on free list */
    if (NULL != p->free) {
        ret = p->free;
        p->free = *((void **)p->free);
        return ret;
    }

    /* Look for space in the last pool */
    if (p->npools) {
        pool = &p->pools[p->npools - 1];
        if (pool->used + p->dsize < p->psize) {
            ret = ((char *) pool->pool) + pool->used;
            pool->used += p->dsize;
            return ret;
        }
    }

    /* Need a new pool */
    pool = new_pool(p);
    if (NULL == pool) return NULL;

    pool->used = p->dsize;
    return pool->pool;
}

char *kstrtok(const char *str, const char *sep_in, ks_tokaux_t *aux)
{
	const unsigned char *p, *start, *sep = (unsigned char *) sep_in;
	if (sep) { // set up the table
		if (str == 0 && aux->finished) return 0; // no need to set up if we have finished
		aux->finished = 0;
		if (sep[0] && sep[1]) {
			aux->sep = -1;
			aux->tab[0] = aux->tab[1] = aux->tab[2] = aux->tab[3] = 0;
			for (p = sep; *p; ++p) aux->tab[*p>>6] |= 1ull<<(*p&0x3f);
		} else aux->sep = sep[0];
	}
	if (aux->finished) return 0;
	else if (str) start = (unsigned char *) str, aux->finished = 0;
	else start = (unsigned char *) aux->p + 1;
	if (aux->sep < 0) {
		for (p = start; *p; ++p)
			if (aux->tab[*p>>6]>>(*p&0x3f)&1) break;
	} else {
		// Using strchr is fast for next token, but slower for
		// last token due to extra pass from strlen.  Overall
		// on a VCF parse this func was 146% faster with // strchr.
		// Equiv to:
		// for (p = start; *p; ++p) if (*p == aux->sep) break;

		// NB: We could use strchrnul() here from glibc if detected,
		// which is ~40% faster again, but it's not so portable.
		// i.e.   p = (uint8_t *)strchrnul((char *)start, aux->sep);
		uint8_t *p2 = (uint8_t *)strchr((char *)start, aux->sep);
		p = p2 ? p2 : start + strlen((char *)start);
	}
	aux->p = (const char *) p; // end of token
	if (*p == 0) aux->finished = 1; // no more tokens
	return (char*)start;
}

static int sam_hrecs_add_ref_altnames(sam_hrecs_t *hrecs, int nref, const char *list) {
    const char *token;
    ks_tokaux_t aux;

    if (!list)
        return 0;

    for (token = kstrtok(list, ",", &aux); token; token = kstrtok(NULL, NULL, &aux)) {
        if (aux.p == token)
            continue;

        char *name = string_ndup(hrecs->str_pool, token, aux.p - token);
        if (!name)
            return -1;
        int r;
        khint_t k = kh_put(m_s2i, hrecs->ref_hash, name, &r);
        if (r < 0) return -1;

        if (r > 0)
            kh_val(hrecs->ref_hash, k) = nref;
        else if (kh_val(hrecs->ref_hash, k) != nref)
            hts_log_warning("Duplicate entry AN:\"%s\" in sam header", name);
    }

    return 0;
}

sam_hrec_tag_t *sam_hrecs_find_key(sam_hrec_type_t *type,
                                   const char *key,
                                   sam_hrec_tag_t **prev) {
    sam_hrec_tag_t *tag, *p = NULL;
    if (!type)
        return NULL;

    for (tag = type->tag; tag; p = tag, tag = tag->next) {
        if (tag->str[0] == key[0] && tag->str[1] == key[1]) {
            if (prev)
                *prev = p;
            return tag;
        }
    }

    if (prev)
        *prev = p;

    return NULL;
}

int sam_hrecs_vupdate(sam_hrecs_t *hrecs, sam_hrec_type_t *type, va_list ap) {
    if (!hrecs)
        return -1;

    for (;;) {
        char *k, *v, *str;
        sam_hrec_tag_t *tag, *prev = NULL;

        if (!(k = (char *)va_arg(ap, char *)))
            break;
        if (!(v = va_arg(ap, char *)))
            v = "";

        tag = sam_hrecs_find_key(type, k, &prev);
        if (!tag) {
            if (!(tag = pool_alloc(hrecs->tag_pool)))
                return -1;
            if (prev)
                prev->next = tag;
            else
                type->tag = tag;

            tag->next = NULL;
        }

        tag->len = 3 + strlen(v);
        str = string_alloc(hrecs->str_pool, tag->len+1);
        if (!str)
            return -1;

        if (snprintf(str, tag->len+1, "%2.2s:%s", k, v) < 0)
            return -1;

        tag->str = str;
    }

    hrecs->dirty = 1; //mark text as dirty and force a rebuild

    return 0;
}

static int sam_hrecs_update(sam_hrecs_t *hrecs, sam_hrec_type_t *type, ...) {
    va_list args;
    int res;
    va_start(args, type);
    res = sam_hrecs_vupdate(hrecs, type, args);
    va_end(args);
    return res;
}

static int sam_hrecs_update_hashes(sam_hrecs_t *hrecs,
                                   khint32_t type,
                                   sam_hrec_type_t *h_type) {
    /* Add to reference hash? */
    if (type == TYPEKEY("SQ")) {
        sam_hrec_tag_t *tag = h_type->tag;
        int nref = hrecs->nref;
        const char *name = NULL;
        const char *altnames = NULL;
        hts_pos_t len = -1;
        int r;
        khint_t k;

        while (tag) {
            if (tag->str[0] == 'S' && tag->str[1] == 'N') {
                assert(tag->len >= 3);
                name = tag->str+3;
            } else if (tag->str[0] == 'L' && tag->str[1] == 'N') {
                assert(tag->len >= 3);
                len = strtoll(tag->str+3, NULL, 10);
            } else if (tag->str[0] == 'A' && tag->str[1] == 'N') {
                assert(tag->len >= 3);
                altnames = tag->str+3;
            }
            tag = tag->next;
        }

        if (!name) {
            hts_log_error("Header includes @SQ line with no SN: tag");
            return -1; // SN should be present, according to spec.
        }

        if (len == -1) {
            hts_log_error("Header includes @SQ line \"%s\" with no LN: tag",
                          name);
            return -1; // LN should be present, according to spec.
        }

        // Seen already?
        k = kh_get(m_s2i, hrecs->ref_hash, name);
        if (k < kh_end(hrecs->ref_hash)) {
            nref = kh_val(hrecs->ref_hash, k);
            int ref_changed_flag = 0;

            // Check for hash entry added by sam_hrecs_refs_from_targets_array()
            if (hrecs->ref[nref].ty == NULL) {
                // Attach header line to existing stub entry.
                hrecs->ref[nref].ty = h_type;
                // Check lengths match; correct if not.
                if (len != hrecs->ref[nref].len) {
                    char tmp[32];
                    snprintf(tmp, sizeof(tmp), "%" PRIhts_pos,
                             hrecs->ref[nref].len);
                    if (sam_hrecs_update(hrecs, h_type, "LN", tmp, NULL) < 0)
                        return -1;
                    ref_changed_flag = 1;
                }
                if (sam_hrecs_add_ref_altnames(hrecs, nref, altnames) < 0)
                    return -1;

                if (ref_changed_flag && (hrecs->refs_changed < 0 || hrecs->refs_changed > nref))
                    hrecs->refs_changed = nref;
                return 0;
            }

            // Check to see if an existing entry is being updated
            if (hrecs->ref[nref].ty == h_type) {
                if (hrecs->ref[nref].len != len) {
                    hrecs->ref[nref].len = len;
                    ref_changed_flag = 1;
                }
                if (!hrecs->ref[nref].name || strcmp(hrecs->ref[nref].name, name)) {
                    hrecs->ref[nref].name = name;
                    ref_changed_flag = 1;
                }
                if (sam_hrecs_add_ref_altnames(hrecs, nref, altnames) < 0)
                    return -1;

                if (ref_changed_flag && (hrecs->refs_changed < 0 || hrecs->refs_changed > nref))
                    hrecs->refs_changed = nref;
                return 0;
            }

            // If here, the name is a duplicate.
            // Check to see if it matches the SN: tag from the earlier record.
            if (strcmp(hrecs->ref[nref].name, name) == 0) {
                hts_log_error("Duplicate entry \"%s\" in sam header",
                                name);
                return -1;
            }

            // Clash with an already-seen altname
            // As SN: should be preferred to AN: add this as a new
            // record and update the hash entry to point to it.
            hts_log_warning("Ref name SN:\"%s\" is a duplicate of an existing AN key", name);
            nref = hrecs->nref;
        }

        if (nref == hrecs->ref_sz) {
            size_t new_sz = hrecs->ref_sz >= 4 ? hrecs->ref_sz + (hrecs->ref_sz / 4) : 32;
            sam_hrec_sq_t *new_ref = realloc(hrecs->ref, sizeof(*hrecs->ref) * new_sz);
            if (!new_ref)
                return -1;
            hrecs->ref = new_ref;
            hrecs->ref_sz = new_sz;
        }

        hrecs->ref[nref].name = name;
        hrecs->ref[nref].len  = len;
        hrecs->ref[nref].ty = h_type;

        k = kh_put(m_s2i, hrecs->ref_hash, hrecs->ref[nref].name, &r);
        if (-1 == r) return -1;
        kh_val(hrecs->ref_hash, k) = nref;

        if (sam_hrecs_add_ref_altnames(hrecs, nref, altnames) < 0)
            return -1;

        if (hrecs->refs_changed < 0 || hrecs->refs_changed > hrecs->nref)
            hrecs->refs_changed = hrecs->nref;
        hrecs->nref++;
    }

    /* Add to read-group hash? */
    if (type == TYPEKEY("RG")) {
        sam_hrec_tag_t *tag = sam_hrecs_find_key(h_type, "ID", NULL);
        int nrg = hrecs->nrg, r;
        khint_t k;

        if (!tag) {
            hts_log_error("Header includes @RG line with no ID: tag");
            return -1;  // ID should be present, according to spec.
        }
        assert(tag->str && tag->len >= 3);

        // Seen already?
        k = kh_get(m_s2i, hrecs->rg_hash, tag->str + 3);
        if (k < kh_end(hrecs->rg_hash)) {
            nrg = kh_val(hrecs->rg_hash, k);
            assert(hrecs->rg[nrg].ty != NULL);
            if (hrecs->rg[nrg].ty != h_type) {
                hts_log_warning("Duplicate entry \"%s\" in sam header",
                                tag->str + 3);
            } else {
                hrecs->rg[nrg].name = tag->str + 3;
                hrecs->rg[nrg].name_len = tag->len - 3;
            }
            return 0;
        }

        if (nrg == hrecs->rg_sz) {
            size_t new_sz = hrecs->rg_sz >= 4 ? hrecs->rg_sz + hrecs->rg_sz / 4 : 4;
            sam_hrec_rg_t *new_rg = realloc(hrecs->rg, sizeof(*hrecs->rg) * new_sz);
            if (!new_rg)
                return -1;
            hrecs->rg = new_rg;
            hrecs->rg_sz = new_sz;
        }

        hrecs->rg[nrg].name = tag->str + 3;
        hrecs->rg[nrg].name_len = tag->len - 3;
        hrecs->rg[nrg].ty   = h_type;
        hrecs->rg[nrg].id   = nrg;

        k = kh_put(m_s2i, hrecs->rg_hash, hrecs->rg[nrg].name, &r);
        if (-1 == r) return -1;
        kh_val(hrecs->rg_hash, k) = nrg;

        hrecs->nrg++;
    }

    /* Add to program hash? */
    if (type == TYPEKEY("PG")) {
        sam_hrec_tag_t *tag;
        sam_hrec_pg_t *new_pg;
        int npg = hrecs->npg;

        if (npg == hrecs->pg_sz) {
            size_t new_sz = hrecs->pg_sz >= 4 ? hrecs->pg_sz + hrecs->pg_sz / 4 : 4;
            new_pg = realloc(hrecs->pg, sizeof(*hrecs->pg) * new_sz);
            if (!new_pg)
                return -1;
            hrecs->pg = new_pg;
            hrecs->pg_sz = new_sz;
        }

        tag = h_type->tag;
        hrecs->pg[npg].name = NULL;
        hrecs->pg[npg].name_len = 0;
        hrecs->pg[npg].ty  = h_type;
        hrecs->pg[npg].id   = npg;
        hrecs->pg[npg].prev_id = -1;

        while (tag) {
            if (tag->str[0] == 'I' && tag->str[1] == 'D') {
                /* Avoid duplicate ID tags coming from other applications */
                if (!hrecs->pg[npg].name) {
                    assert(tag->len >= 3);
                    hrecs->pg[npg].name = tag->str + 3;
                    hrecs->pg[npg].name_len = tag->len - 3;
                } else {
                    hts_log_warning("PG line with multiple ID tags. The first encountered was preferred - ID:%s", hrecs->pg[npg].name);
                }
            } else if (tag->str[0] == 'P' && tag->str[1] == 'P') {
                // Resolve later if needed
                khint_t k;
                k = kh_get(m_s2i, hrecs->pg_hash, tag->str+3);

                if (k != kh_end(hrecs->pg_hash)) {
                    int p_id = kh_val(hrecs->pg_hash, k);
                    hrecs->pg[npg].prev_id = hrecs->pg[p_id].id;

                    /* Unmark previous entry as a PG termination */
                    if (hrecs->npg_end > 0 &&
                        hrecs->pg_end[hrecs->npg_end-1] == p_id) {
                        hrecs->npg_end--;
                    } else {
                        int i;
                        for (i = 0; i < hrecs->npg_end; i++) {
                            if (hrecs->pg_end[i] == p_id) {
                                memmove(&hrecs->pg_end[i], &hrecs->pg_end[i+1],
                                        (hrecs->npg_end-i-1)*sizeof(*hrecs->pg_end));
                                hrecs->npg_end--;
                            }
                        }
                    }
                } else {
                    hrecs->pg[npg].prev_id = -1;
                }
            }
            tag = tag->next;
        }

        if (hrecs->pg[npg].name) {
            khint_t k;
            int r;
            k = kh_put(m_s2i, hrecs->pg_hash, hrecs->pg[npg].name, &r);
            if (-1 == r) return -1;
            kh_val(hrecs->pg_hash, k) = npg;
        } else {
            return -1; // ID should be present, according to spec.
        }

        /* Add to npg_end[] array. Remove later if we find a PP line */
        if (hrecs->npg_end >= hrecs->npg_end_alloc) {
            int *new_pg_end;
            int  new_alloc = hrecs->npg_end_alloc ? hrecs->npg_end_alloc*2 : 4;

            new_pg_end = realloc(hrecs->pg_end, new_alloc * sizeof(int));
            if (!new_pg_end)
                return -1;
            hrecs->npg_end_alloc = new_alloc;
            hrecs->pg_end = new_pg_end;
        }
        hrecs->pg_end[hrecs->npg_end++] = npg;

        hrecs->npg++;
    }

    return 0;
}

#ifndef SSIZE_MAX /* SSIZE_MAX is POSIX 1 */
#define SSIZE_MAX LONG_MAX
#endif

static int sam_hrecs_parse_lines(sam_hrecs_t *hrecs, const char *hdr, size_t len) {
    size_t i, lno;

    if (!hrecs || len > SSIZE_MAX)
        return -1;

    if (!len)
        len = strlen(hdr);

    if (len < 3) {
        if (len == 0 || *hdr == '\0') return 0;
        sam_hrecs_error("Header line too short", hdr, len, 1);
        return -1;
    }

    for (i = 0, lno = 1; i < len - 3 && hdr[i] != '\0'; i++, lno++) {
        khint32_t type;
        khint_t k;

        int l_start = i, new;
        sam_hrec_type_t *h_type;
        sam_hrec_tag_t *h_tag, *last;

        if (hdr[i] != '@') {
            sam_hrecs_error("Header line does not start with '@'",
                          &hdr[l_start], len - l_start, lno);
            return -1;
        }

        if (!isalpha_c(hdr[i+1]) || !isalpha_c(hdr[i+2])) {
            sam_hrecs_error("Header line does not have a two character key",
                          &hdr[l_start], len - l_start, lno);
            return -1;
        }
        type = TYPEKEY(&hdr[i+1]);

        i += 3;
        if (i == len || hdr[i] == '\n')
            continue;

        // Add the header line type
        if (!(h_type = pool_alloc(hrecs->type_pool)))
            return -1;
        k = kh_put(sam_hrecs_t, hrecs->h, type, &new);
        if (new < 0)
            return -1;

        h_type->type = type;

        // Add to end of global list
        sam_hrecs_global_list_add(hrecs, h_type, NULL);

        // Form the ring, either with self or other lines of this type
        if (!new) {
            sam_hrec_type_t *t = kh_val(hrecs->h, k), *p;
            p = t->prev;

            assert(p->next == t);
            p->next = h_type;
            h_type->prev = p;

            t->prev = h_type;
            h_type->next = t;
        } else {
            kh_val(hrecs->h, k) = h_type;
            h_type->prev = h_type->next = h_type;
        }

        // Parse the tags on this line
        last = NULL;
        if (type == TYPEKEY("CO")) {
            size_t j;

            if (i == len || hdr[i] != '\t') {
                sam_hrecs_error("Missing tab",
                              &hdr[l_start], len - l_start, lno);
                return -1;
            }

            for (j = ++i; j < len && hdr[j] != '\0' && hdr[j] != '\n'; j++)
                ;

            if (!(h_type->tag = h_tag = pool_alloc(hrecs->tag_pool)))
                return -1;
            h_tag->str = string_ndup(hrecs->str_pool, &hdr[i], j-i);
            h_tag->len = j-i;
            h_tag->next = NULL;
            if (!h_tag->str)
                return -1;

            i = j;

        } else {
            do {
                size_t j;

                if (i == len || hdr[i] != '\t') {
                    sam_hrecs_error("Missing tab",
                                  &hdr[l_start], len - l_start, lno);
                    return -1;
                }

                for (j = ++i; j < len && hdr[j] != '\0' && hdr[j] != '\n' && hdr[j] != '\t'; j++)
                    ;

                if (j - i < 3 || hdr[i + 2] != ':') {
                    sam_hrecs_error("Malformed key:value pair",
                                   &hdr[l_start], len - l_start, lno);
                    return -1;
                }

                if (!(h_tag = pool_alloc(hrecs->tag_pool)))
                    return -1;
                h_tag->str = string_ndup(hrecs->str_pool, &hdr[i], j-i);
                h_tag->len = j-i;
                h_tag->next = NULL;
                if (!h_tag->str)
                    return -1;

                if (last)
                    last->next = h_tag;
                else
                    h_type->tag = h_tag;

                last = h_tag;
                i = j;
            } while (i < len && hdr[i] != '\0' && hdr[i] != '\n');
        }

        /* Update RG/SQ hashes */
        if (-1 == sam_hrecs_update_hashes(hrecs, type, h_type))
            return -1;
    }

    return 0;
}

sam_hrec_type_t *sam_hrecs_find_type_id(sam_hrecs_t *hrecs, const char *type,
                                     const char *ID_key, const char *ID_value) {
    if (!hrecs || !type)
        return NULL;
    sam_hrec_type_t *t1, *t2;
    khint_t k;

    /* Special case for types we have prebuilt hashes on */
    if (ID_key) {
        if (!ID_value)
            return NULL;

        if (type[0]   == 'S' && type[1]   == 'Q' &&
            ID_key[0] == 'S' && ID_key[1] == 'N') {
            k = kh_get(m_s2i, hrecs->ref_hash, ID_value);
            return k != kh_end(hrecs->ref_hash)
                ? hrecs->ref[kh_val(hrecs->ref_hash, k)].ty
                : NULL;
        }

        if (type[0]   == 'R' && type[1]   == 'G' &&
            ID_key[0] == 'I' && ID_key[1] == 'D') {
            k = kh_get(m_s2i, hrecs->rg_hash, ID_value);
            return k != kh_end(hrecs->rg_hash)
                ? hrecs->rg[kh_val(hrecs->rg_hash, k)].ty
                : NULL;
        }

        if (type[0]   == 'P' && type[1]   == 'G' &&
            ID_key[0] == 'I' && ID_key[1] == 'D') {
            k = kh_get(m_s2i, hrecs->pg_hash, ID_value);
            return k != kh_end(hrecs->pg_hash)
                ? hrecs->pg[kh_val(hrecs->pg_hash, k)].ty
                : NULL;
        }
    }

    k = kh_get(sam_hrecs_t, hrecs->h, TYPEKEY(type));
    if (k == kh_end(hrecs->h))
        return NULL;

    if (!ID_key)
        return kh_val(hrecs->h, k);

    t1 = t2 = kh_val(hrecs->h, k);
    do {
        sam_hrec_tag_t *tag;
        for (tag = t1->tag; tag; tag = tag->next) {
            if (tag->str[0] == ID_key[0] && tag->str[1] == ID_key[1]) {
                const char *cp1 = tag->str+3;
                const char *cp2 = ID_value;
                while (*cp1 && *cp1 == *cp2)
                    cp1++, cp2++;
                if (*cp2 || *cp1)
                    continue;
                return t1;
            }
        }
        t1 = t1->next;
    } while (t1 != t2);

    return NULL;
}

static int sam_hrecs_vadd(sam_hrecs_t *hrecs, const char *type, va_list ap, ...) {
    va_list args;
    sam_hrec_type_t *h_type;
    sam_hrec_tag_t *h_tag, *last=NULL;
    int new;
    khint32_t type_i = TYPEKEY(type), k;

    if (!strncmp(type, "HD", 2) && (h_type = sam_hrecs_find_type_id(hrecs, "HD", NULL, NULL)))
        return sam_hrecs_vupdate(hrecs, h_type, ap);

    if (!(h_type = pool_alloc(hrecs->type_pool)))
        return -1;
    k = kh_put(sam_hrecs_t, hrecs->h, type_i, &new);
    if (new < 0)
        return -1;

    h_type->type = type_i;

    // Form the ring, either with self or other lines of this type
    if (!new) {
        sam_hrec_type_t *t = kh_val(hrecs->h, k), *p;
        p = t->prev;

        assert(p->next == t);
        p->next = h_type;
        h_type->prev = p;

        t->prev = h_type;
        h_type->next = t;
    } else {
        kh_val(hrecs->h, k) = h_type;
        h_type->prev = h_type->next = h_type;
    }
    h_type->tag = NULL;

    // Add to global line ordering after any existing line of the same type,
    // or at the end if no line of this type exists yet.
    sam_hrecs_global_list_add(hrecs, h_type, !new ? h_type->prev : NULL);

    // Check linked-list invariants
    assert(h_type->prev->next == h_type);
    assert(h_type->next->prev == h_type);
    assert(h_type->global_prev->global_next == h_type);
    assert(h_type->global_next->global_prev == h_type);

    // Any ... varargs
    va_start(args, ap);
    for (;;) {
        char *key, *val = NULL, *str;

        if (!(key = (char *)va_arg(args, char *)))
            break;
        if (strncmp(type, "CO", 2) && !(val = (char *)va_arg(args, char *)))
            break;
        if (*val == '\0')
            continue;

        if (!(h_tag = pool_alloc(hrecs->tag_pool)))
            return -1;

        if (strncmp(type, "CO", 2)) {
            h_tag->len = 3 + strlen(val);
            str = string_alloc(hrecs->str_pool, h_tag->len+1);
            if (!str || snprintf(str, h_tag->len+1, "%2.2s:%s", key, val) < 0)
                return -1;
            h_tag->str = str;
        } else {
            h_tag->len = strlen(key);
            h_tag->str = string_ndup(hrecs->str_pool, key, h_tag->len);
            if (!h_tag->str)
                return -1;
        }

        h_tag->next = NULL;
        if (last)
            last->next = h_tag;
        else
            h_type->tag = h_tag;

        last = h_tag;
    }
    va_end(args);

    // Plus the specified va_list params
    for (;;) {
        char *key, *val = NULL, *str;

        if (!(key = (char *)va_arg(ap, char *)))
            break;
        if (strncmp(type, "CO", 2) && !(val = (char *)va_arg(ap, char *)))
            break;

        if (!(h_tag = pool_alloc(hrecs->tag_pool)))
            return -1;

        if (strncmp(type, "CO", 2)) {
            h_tag->len = 3 + strlen(val);
            str = string_alloc(hrecs->str_pool, h_tag->len+1);
            if (!str || snprintf(str, h_tag->len+1, "%2.2s:%s", key, val) < 0)
                return -1;
            h_tag->str = str;
        } else {
            h_tag->len = strlen(key);
            h_tag->str = string_ndup(hrecs->str_pool, key, h_tag->len);
            if (!h_tag->str)
                return -1;
        }

        h_tag->next = NULL;
        if (last)
            last->next = h_tag;
        else
            h_type->tag = h_tag;

        last = h_tag;
    }

    if (-1 == sam_hrecs_update_hashes(hrecs, TYPEKEY(type), h_type))
        return -1;

    if (!strncmp(type, "PG", 2))
        hrecs->pgs_changed = 1;

    hrecs->dirty = 1;

    return 0;
}

static int sam_hrecs_add(sam_hrecs_t *hrecs, const char *type, ...) {
    va_list args;
    int res;
    va_start(args, type);
    res = sam_hrecs_vadd(hrecs, type, args, NULL);
    va_end(args);
    return res;
}

static int add_stub_ref_sq_lines(sam_hrecs_t *hrecs) {
    int tid;
    char len[32];

    for (tid = 0; tid < hrecs->nref; tid++) {
        if (hrecs->ref[tid].ty == NULL) {
            snprintf(len, sizeof(len), "%"PRIhts_pos, hrecs->ref[tid].len);
            if (sam_hrecs_add(hrecs, "SQ",
                              "SN", hrecs->ref[tid].name,
                              "LN", len, NULL) != 0)
                return -1;

            // Check that the stub has actually been filled
            if(hrecs->ref[tid].ty == NULL) {
                hts_log_error("Reference stub with tid=%d, name=\"%s\", len=%"PRIhts_pos" could not be filled",
                        tid, hrecs->ref[tid].name, hrecs->ref[tid].len);
                return -1;
            }
        }
    }
    return 0;
}

int sam_hdr_update_target_arrays(sam_hdr_t *bh, const sam_hrecs_t *hrecs,
                                 int refs_changed) {
    if (!bh || !hrecs)
        return -1;

    if (refs_changed < 0)
        return 0;

    // Grow arrays if necessary
    if (bh->n_targets < hrecs->nref) {
        char **new_names = realloc(bh->target_name,
                                   hrecs->nref * sizeof(*new_names));
        if (!new_names)
            return -1;
        bh->target_name = new_names;
        uint32_t *new_lens = realloc(bh->target_len,
                                     hrecs->nref * sizeof(*new_lens));
        if (!new_lens)
            return -1;
        bh->target_len = new_lens;
    }

    // Update names and lengths where changed
    // hrecs->refs_changed is the first ref that has been updated, so ones
    // before that can be skipped.
    int i;
    khint_t k;
    khash_t(s2i) *long_refs = (khash_t(s2i) *) bh->sdict;
    for (i = refs_changed; i < hrecs->nref; i++) {
        if (i >= bh->n_targets
            || strcmp(bh->target_name[i], hrecs->ref[i].name) != 0) {
            if (i < bh->n_targets)
                free(bh->target_name[i]);
            bh->target_name[i] = strdup(hrecs->ref[i].name);
            if (!bh->target_name[i])
                return -1;
        }
        if (hrecs->ref[i].len < UINT32_MAX) {
            bh->target_len[i] = hrecs->ref[i].len;

            if (!long_refs)
                continue;

            // Check if we have an old length, if so remove it.
            k = kh_get(s2i, long_refs, bh->target_name[i]);
            if (k < kh_end(long_refs))
                kh_del(s2i, long_refs, k);
        } else {
            bh->target_len[i] = UINT32_MAX;
            if (bh->hrecs != hrecs) {
                // Called from sam_hdr_dup; need to add sdict entries
                if (!long_refs) {
                    if (!(bh->sdict = long_refs = kh_init(s2i)))
                        return -1;
                }

                // Add / update length
                int absent;
                k = kh_put(s2i, long_refs, bh->target_name[i], &absent);
                if (absent < 0)
                    return -1;
                kh_val(long_refs, k) = hrecs->ref[i].len;
            }
        }
    }

    // Free up any names that have been removed
    for (; i < bh->n_targets; i++) {
        if (long_refs) {
            k = kh_get(s2i, long_refs, bh->target_name[i]);
            if (k < kh_end(long_refs))
                kh_del(s2i, long_refs, k);
        }
        free(bh->target_name[i]);
    }

    bh->n_targets = hrecs->nref;
    return 0;
}

static int rebuild_target_arrays(sam_hdr_t *bh) {
    if (!bh || !bh->hrecs)
        return -1;

    sam_hrecs_t *hrecs = bh->hrecs;
    if (hrecs->refs_changed < 0)
        return 0;

    if (sam_hdr_update_target_arrays(bh, hrecs, hrecs->refs_changed) != 0)
        return -1;

    hrecs->refs_changed = -1;
    return 0;
}

int sam_hdr_fill_hrecs(sam_hdr_t *bh) {
    sam_hrecs_t *hrecs = sam_hrecs_new();

    if (!hrecs)
        return -1;

    if (bh->target_name && bh->target_len && bh->n_targets > 0) {
        if (sam_hrecs_refs_from_targets_array(hrecs, bh) != 0) {
            sam_hrecs_free(hrecs);
            return -1;
        }
    }

    // Parse existing header text
    if (bh->text && bh->l_text > 0) {
        if (sam_hrecs_parse_lines(hrecs, bh->text, bh->l_text) != 0) {
            sam_hrecs_free(hrecs);
            return -1;
        }
    }

    if (add_stub_ref_sq_lines(hrecs) < 0) {
        sam_hrecs_free(hrecs);
        return -1;
    }

    bh->hrecs = hrecs;

    if (hrecs->refs_changed >= 0 && rebuild_target_arrays(bh) != 0)
        return -1;

    return 0;
}

int sam_hdr_name2tid(sam_hdr_t *bh, const char *ref) {
    sam_hrecs_t *hrecs;
    khint_t k;

    if (!bh)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -2;
        hrecs = bh->hrecs;
    }

    if (!hrecs->ref_hash)
        return -1;

    k = kh_get(m_s2i, hrecs->ref_hash, ref);
    return k == kh_end(hrecs->ref_hash) ? -1 : kh_val(hrecs->ref_hash, k);
}

static int bam_tag2cigar(bam1_t *b, int recal_bin, int give_warning) // return 0 if CIGAR is untouched; 1 if CIGAR is updated with CG
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
    if (CG[0] != 'B' || !(CG[1] == 'I' || CG[1] == 'i'))
        return 0; // not of type B,I
    CG_len = le_to_u32(CG + 2);
    if (CG_len < c->n_cigar || CG_len >= 1U<<29) return 0; // don't move if the real CIGAR length is shorter than the fake cigar length

    // move from the CG tag to the right position
    cigar_st = (uint8_t*)cigar0 - b->data;
    c->n_cigar = CG_len;
    n_cigar4 = c->n_cigar * 4;
    CG_st = CG - b->data - 2;
    CG_en = CG_st + 8 + n_cigar4;
    if (possibly_expand_bam_data(b, n_cigar4 - fake_bytes) < 0) return -1;
    b->l_data = b->l_data - fake_bytes + n_cigar4; // we need c->n_cigar-fake_bytes bytes to swap CIGAR to the right place
    memmove(b->data + cigar_st + n_cigar4, b->data + cigar_st + fake_bytes, ori_len - (cigar_st + fake_bytes)); // insert c->n_cigar-fake_bytes empty space to make room
    memcpy(b->data + cigar_st, b->data + (n_cigar4 - fake_bytes) + CG_st + 8, n_cigar4); // copy the real CIGAR to the right place; -fake_bytes for the fake CIGAR
    if (ori_len > CG_en) // move data after the CG tag
        memmove(b->data + CG_st + n_cigar4 - fake_bytes, b->data + CG_en + n_cigar4 - fake_bytes, ori_len - CG_en);
    b->l_data -= n_cigar4 + 8; // 8: CGBI (4 bytes) and CGBI length (4)
    if (recal_bin)
        b->core.bin = hts_reg2bin(b->core.pos, bam_endpos(b), 14, 5);
    if (give_warning)
        hts_log_warning("%s encodes a CIGAR with %d operators at the CG tag", bam_get_qname(b), c->n_cigar);
    return 1;
}

int sam_parse1(kstring_t *s, sam_hdr_t *h, bam1_t *b)
{
#define _read_token(_p) (_p); do { char *tab = strchr((_p), '\t'); if (!tab) goto err_ret; *tab = '\0'; (_p) = tab + 1; } while (0)

#if HTS_ALLOW_UNALIGNED != 0 && ULONG_MAX == 0xffffffffffffffff

// Macro that operates on 64-bits at a time.
#define COPY_MINUS_N(to,from,n,l,failed)                        \
    do {                                                        \
        uint64_u *from8 = (uint64_u *)(from);                   \
        uint64_u *to8 = (uint64_u *)(to);                       \
        uint64_t uflow = 0;                                     \
        size_t l8 = (l)>>3, i;                                  \
        for (i = 0; i < l8; i++) {                              \
            to8[i] = from8[i] - (n)*0x0101010101010101UL;       \
            uflow |= to8[i];                                    \
        }                                                       \
        for (i<<=3; i < (l); ++i) {                             \
            to[i] = from[i] - (n);                              \
            uflow |= to[i];                                     \
        }                                                       \
        failed = (uflow & 0x8080808080808080UL) > 0;            \
    } while (0)

#else

// Basic version which operates a byte at a time
#define COPY_MINUS_N(to,from,n,l,failed) do {                \
        uint8_t uflow = 0;                                   \
        for (i = 0; i < (l); ++i) {                          \
            (to)[i] = (from)[i] - (n);                       \
            uflow |= (uint8_t) (to)[i];                      \
        }                                                    \
        failed = (uflow & 0x80) > 0;                         \
    } while (0)

#endif

#define _get_mem(type_t, x, b, l) if (possibly_expand_bam_data((b), (l)) < 0) goto err_ret; *(x) = (type_t*)((b)->data + (b)->l_data); (b)->l_data += (l)
#define _parse_err(cond, ...) do { if (cond) { hts_log_error(__VA_ARGS__); goto err_ret; } } while (0)
#define _parse_warn(cond, ...) do { if (cond) { hts_log_warning(__VA_ARGS__); } } while (0)

    uint8_t *t;

    char *p = s->s, *q;
    int i, overflow = 0;
    char logbuf[40];
    hts_pos_t cigreflen;
    bam1_core_t *c = &b->core;

    b->l_data = 0;
    memset(c, 0, 32);

    // qname
    q = _read_token(p);

    _parse_warn(p - q <= 1, "empty query name");
    _parse_err(p - q > 255, "query name too long");
    // resize large enough for name + extranul
    if (possibly_expand_bam_data(b, (p - q) + 4) < 0) goto err_ret;
    memcpy(b->data + b->l_data, q, p-q); b->l_data += p-q;

    c->l_extranul = (4 - (b->l_data & 3)) & 3;
    memcpy(b->data + b->l_data, "\0\0\0\0", c->l_extranul);
    b->l_data += c->l_extranul;

    c->l_qname = p - q + c->l_extranul;

    // flag
    c->flag = parse_sam_flag(p, &p, &overflow);
    if (*p++ != '\t') goto err_ret; // malformated flag

    // chr
    q = _read_token(p);
    if (strcmp(q, "*")) {
        _parse_err(h->n_targets == 0, "no SQ lines present in the header");
        c->tid = bam_name2id(h, q);
        _parse_err(c->tid < -1, "failed to parse header");
        _parse_warn(c->tid < 0, "unrecognized reference name %s; treated as unmapped", hts_strprint(logbuf, sizeof logbuf, '"', q, SIZE_MAX));
    } else c->tid = -1;

    // pos
    c->pos = hts_str2uint(p, &p, 63, &overflow) - 1;
    if (*p++ != '\t') goto err_ret;
    if (c->pos < 0 && c->tid >= 0) {
        _parse_warn(1, "mapped query cannot have zero coordinate; treated as unmapped");
        c->tid = -1;
    }
    if (c->tid < 0) c->flag |= BAM_FUNMAP;

    // mapq
    c->qual = hts_str2uint(p, &p, 8, &overflow);
    if (*p++ != '\t') goto err_ret;
    // cigar
    if (*p != '*') {
        uint32_t *cigar = NULL;
        int old_l_data = b->l_data;
        int n_cigar = bam_parse_cigar(p, &p, b);
        if (n_cigar < 1 || *p++ != '\t') goto err_ret;
        cigar = (uint32_t *)(b->data + old_l_data);

        // can't use bam_endpos() directly as some fields not yet set up
        cigreflen = (!(c->flag&BAM_FUNMAP))? bam_cigar2rlen(c->n_cigar, cigar) : 1;
        if (cigreflen == 0) cigreflen = 1;
    } else {
        _parse_warn(!(c->flag&BAM_FUNMAP), "mapped query must have a CIGAR; treated as unmapped");
        c->flag |= BAM_FUNMAP;
        q = _read_token(p);
        cigreflen = 1;
    }
    _parse_err(HTS_POS_MAX - cigreflen <= c->pos,
               "read ends beyond highest supported position");
    c->bin = hts_reg2bin(c->pos, c->pos + cigreflen, 14, 5);
    // mate chr
    q = _read_token(p);
    if (strcmp(q, "=") == 0) {
        c->mtid = c->tid;
    } else if (strcmp(q, "*") == 0) {
        c->mtid = -1;
    } else {
        c->mtid = bam_name2id(h, q);
        _parse_err(c->mtid < -1, "failed to parse header");
        _parse_warn(c->mtid < 0, "unrecognized mate reference name %s; treated as unmapped", hts_strprint(logbuf, sizeof logbuf, '"', q, SIZE_MAX));
    }
    // mpos
    c->mpos = hts_str2uint(p, &p, 63, &overflow) - 1;
    if (*p++ != '\t') goto err_ret;
    if (c->mpos < 0 && c->mtid >= 0) {
        _parse_warn(1, "mapped mate cannot have zero coordinate; treated as unmapped");
        c->mtid = -1;
    }
    // tlen
    c->isize = hts_str2int(p, &p, 64, &overflow);
    if (*p++ != '\t') goto err_ret;
    // seq
    q = _read_token(p);
    if (strcmp(q, "*")) {
        _parse_err(p - q - 1 > INT32_MAX, "read sequence is too long");
        c->l_qseq = p - q - 1;
        hts_pos_t ql = bam_cigar2qlen(c->n_cigar, (uint32_t*)(b->data + c->l_qname));
        _parse_err(c->n_cigar && ql != c->l_qseq, "CIGAR and query sequence are of different length");
        i = (c->l_qseq + 1) >> 1;
        _get_mem(uint8_t, &t, b, i);

        unsigned int lqs2 = c->l_qseq&~1, i;
        for (i = 0; i < lqs2; i+=2)
            t[i>>1] = (seq_nt16_table[(unsigned char)q[i]] << 4) | seq_nt16_table[(unsigned char)q[i+1]];
        for (; i < c->l_qseq; ++i)
            t[i>>1] = seq_nt16_table[(unsigned char)q[i]] << ((~i&1)<<2);
    } else c->l_qseq = 0;
    // qual
    _get_mem(uint8_t, &t, b, c->l_qseq);
    if (p[0] == '*' && (p[1] == '\t' || p[1] == '\0')) {
        memset(t, 0xff, c->l_qseq);
        p += 2;
    } else {
        int failed = 0;
        _parse_err(s->l - (p - s->s) < c->l_qseq
                   || (p[c->l_qseq] != '\t' && p[c->l_qseq] != '\0'),
                   "SEQ and QUAL are of different length");
        COPY_MINUS_N(t, p, 33, c->l_qseq, failed);
        _parse_err(failed, "invalid QUAL character");
        p += c->l_qseq + 1;
    }

    // aux
    if (aux_parse(p, s->s + s->l, b, 0 , NULL) < 0)
        goto err_ret;

    if (bam_tag2cigar(b, 1, 1) < 0)
        return -2;
    return 0;

#undef _parse_warn
#undef _parse_err
#undef _get_mem
#undef _read_token
err_ret:
    return -2;
}

