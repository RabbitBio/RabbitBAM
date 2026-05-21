
#ifndef SAM_PARSE_C_H
#define SAM_PARSE_C_H

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

#ifndef MPI_SAM_PARSE_DETAIL
#define MPI_SAM_PARSE_DETAIL 1
#endif

typedef struct {
    const sam_hdr_t *hdr[2];
    int name_len[2];
    int tid[2];
    char name[2][64];
} MpiSamParseFastCache;

typedef struct {
    uint64_t core_cycles;
    uint64_t aux_cycles;
    uint64_t cg_cycles;
} MpiSamParseFastTiming;

int sam_parse1(kstring_t *s, sam_hdr_t *h, bam1_t *b);
int sam_parse1_mpi_fast(kstring_t *s, sam_hdr_t *h, bam1_t *b,
                        MpiSamParseFastCache *cache,
                        MpiSamParseFastTiming *timing);

#ifdef __cplusplus
}
#endif

#endif // SAM_PARSE_C_H
