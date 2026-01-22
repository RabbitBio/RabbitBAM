
#ifndef TEST_C_H
#define TEST_C_H

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>


#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <htslib/hts_endian.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/khash.h>

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

// int test_add(int a, int b);
// void test_print(const char *msg);

int sam_parse1(kstring_t *s, sam_hdr_t *h, bam1_t *b);

#ifdef __cplusplus
}
#endif

#endif // TEST_C_H



