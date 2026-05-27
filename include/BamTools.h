#ifndef BAMTOOLS_H
#define BAMTOOLS_H

#include <stdint.h>
#include <cstdlib>
#include <mutex>
#include <thread>
#include <atomic>
#include <condition_variable>
#include <vector>
#include <cmath>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>


#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/khash.h>
#include "Globals.h"

#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8


const size_t MAX_RECORDS_PER_BLOCK = 1024;   // legacy fast path per-BGZF block record capacity
const size_t MPI_RECORDS_PER_BLOCK = 2048;   // MPI runtime per-BGZF block record capacity
#define MAX_SAM_LINE_SIZE 8192   // 8 KB     //一条sam文本的最大长度
#define BATCH_PER_CORE 1024
#define BATCH_SIZE (64 * BATCH_PER_CORE)     //一个批次的大小
const size_t MAX_SAM_FORMAT_CORE_BUFFER_SIZE = MAX_SAM_LINE_SIZE * ((BATCH_SIZE / 64) + 1); //每一个核要承担的sam格式化输出缓冲区大小，预留一些空间以防万一

const size_t INIT_DATA_SIZE = 1024;           // short-read fast-path estimate; MPI uses arenas
const int FUSED_SAM2BAM_BAM_POOL_SIZE = 720999;  // SAM2BAM subgroup record capacity
const size_t MPI_BAM_BLOCK_ARENA_SIZE = BGZF_MAX_BLOCK_SIZE * 2;


const int CGS_NUM_CGS = 6;
const int CGS_PES_PER_CG = 64;
const int CGS_NB = CGS_NUM_CGS * CGS_PES_PER_CG;
// const int CGS_BAM2SAM_FORMAT_RECORDS_PER_CPE = 192;
const int CGS_BAM2SAM_FORMAT_RECORDS_PER_CPE = 900;
const size_t CGS_SAM_CHUNK_SIZE = 512 * 1024;
const size_t CGS_SAM_CHUNK_BUFFER_SIZE = CGS_SAM_CHUNK_SIZE + MAX_SAM_LINE_SIZE;
const int CGS_SAM2BAM_MAX_BAMS_PER_CHUNK = 4096;
const int CGS_SAM2BAM_RECORD_POOL_SIZE = FUSED_SAM2BAM_BAM_POOL_SIZE;
const size_t CGS_SAM_FORMAT_CORE_BUFFER_SIZE =
    MAX_SAM_LINE_SIZE * (CGS_BAM2SAM_FORMAT_RECORDS_PER_CPE + 1);

#define SAM_CHUNK_SIZE (4 * 1024 * 1024)         
#define CHUNK_BUFFER_SIZE (5 * 1024 * 1024)      
#define MAX_BAMS_PER_CHUNK 12110              // 4MB最多包含的记录数

typedef struct bam_block bam_block;

typedef struct {
    const sam_hdr_t *hdr;
    bam1_t *bams[BATCH_SIZE];
    kstring_t core_out_lines[64];  
    int status[64];
    size_t required_capacity[64];
    int count;   
} SamFormatBatch;

enum BoundsLimitId {
    BOUNDS_LIMIT_NONE = 0,
    BOUNDS_LIMIT_MAX_RECORDS_PER_BLOCK,
    BOUNDS_LIMIT_MAX_SAM_LINE_SIZE,
    BOUNDS_LIMIT_INIT_DATA_SIZE,
    BOUNDS_LIMIT_MAX_BAMS_PER_CHUNK,
    BOUNDS_LIMIT_CHUNK_BUFFER_SIZE,
    BOUNDS_LIMIT_FUSED_SAM2BAM_BAM_POOL_SIZE,
    BOUNDS_LIMIT_BGZF_RECORD_SIZE,
    BOUNDS_LIMIT_MAX_SAM_FORMAT_CORE_BUFFER_SIZE,
    BOUNDS_LIMIT_GENERIC_RUNTIME_ERROR
};

enum MpiFlagstatCounterIdShared {
    RB_FLAGSTAT_TOTAL = 0,
    RB_FLAGSTAT_PRIMARY,
    RB_FLAGSTAT_SECONDARY,
    RB_FLAGSTAT_SUPPLEMENTARY,
    RB_FLAGSTAT_DUPLICATES,
    RB_FLAGSTAT_PRIMARY_DUPLICATES,
    RB_FLAGSTAT_MAPPED,
    RB_FLAGSTAT_PRIMARY_MAPPED,
    RB_FLAGSTAT_PAIRED,
    RB_FLAGSTAT_READ1,
    RB_FLAGSTAT_READ2,
    RB_FLAGSTAT_PROPERLY_PAIRED,
    RB_FLAGSTAT_PAIR_MAPPED,
    RB_FLAGSTAT_SINGLETONS,
    RB_FLAGSTAT_DIFF_CHR,
    RB_FLAGSTAT_DIFF_CHR_MAPQ5,
    RB_FLAGSTAT_COUNTER_COUNT
};

struct MpiFlagstatCountSlice {
    long long values[RB_FLAGSTAT_COUNTER_COUNT][2];
};

struct MpiFlagstatCountPara {
    int block_id;
    bam_block *input_block;
    bam_block *un_comp_block;
    unsigned char *scratch_data;
    size_t scratch_capacity;
    MpiFlagstatCountSlice *counts;
    int n_total_records;
    int status;
    int record_index;
    long long actual_value;
    long long limit_value;
    int limit_id;
    uint64_t decomp_alloc_cycles;
    uint64_t decomp_inflate_cycles;
    uint64_t decomp_crc_cycles;
    uint64_t decomp_parse_cycles;
    uint64_t decomp_total_cycles;
};

const int RB_MPI_STATS_MAX_INSERT_SIZE = 8000;
const int RB_MPI_STATS_INSERT_BINS = RB_MPI_STATS_MAX_INSERT_SIZE + 1;

enum MpiStatsLongIdShared {
    RB_STATS_NREADS_1ST = 0,
    RB_STATS_NREADS_2ND,
    RB_STATS_NREADS_OTHER,
    RB_STATS_NREADS_FILTERED,
    RB_STATS_NREADS_DUP,
    RB_STATS_NREADS_UNMAPPED,
    RB_STATS_NREADS_SINGLE_MAPPED,
    RB_STATS_NREADS_PAIRED_AND_MAPPED,
    RB_STATS_NREADS_PROPERLY_PAIRED,
    RB_STATS_NREADS_PAIRED_TECH,
    RB_STATS_NREADS_ANOMALOUS,
    RB_STATS_NREADS_MQ0,
    RB_STATS_NREADS_QCFAILED,
    RB_STATS_NREADS_SECONDARY,
    RB_STATS_NREADS_SUPPLEMENTARY,
    RB_STATS_TOTAL_LEN,
    RB_STATS_TOTAL_LEN_1ST,
    RB_STATS_TOTAL_LEN_2ND,
    RB_STATS_TOTAL_LEN_DUP,
    RB_STATS_NBASES_MAPPED,
    RB_STATS_NBASES_MAPPED_CIGAR,
    RB_STATS_NBASES_TRIMMED,
    RB_STATS_NMISMATCHES,
    RB_STATS_MAX_LEN,
    RB_STATS_MAX_LEN_1ST,
    RB_STATS_MAX_LEN_2ND,
    RB_STATS_SUM_QUAL,
    RB_STATS_LONG_COUNT
};

enum MpiStatsOrientDiagIdShared {
    RB_ORIENT_DIAG_REF_IN = 0,
    RB_ORIENT_DIAG_REF_OUT,
    RB_ORIENT_DIAG_REF_OTHER,
    RB_ORIENT_DIAG_DOC_IN,
    RB_ORIENT_DIAG_DOC_OUT,
    RB_ORIENT_DIAG_DOC_OTHER,
    RB_ORIENT_DIAG_REF_IN_LT_READ,
    RB_ORIENT_DIAG_REF_IN_LT_2READ,
    RB_ORIENT_DIAG_REF_IN_GE_2READ,
    RB_ORIENT_DIAG_REF_OUT_LT_READ,
    RB_ORIENT_DIAG_REF_OUT_LT_2READ,
    RB_ORIENT_DIAG_REF_OUT_GE_2READ,
    RB_ORIENT_DIAG_POS_NEG_STRAND_NEG,
    RB_ORIENT_DIAG_POS_NEG_STRAND_POS,
    RB_ORIENT_DIAG_POS_ZERO_STRAND_NEG,
    RB_ORIENT_DIAG_POS_ZERO_STRAND_POS,
    RB_ORIENT_DIAG_POS_POS_STRAND_NEG,
    RB_ORIENT_DIAG_POS_POS_STRAND_POS,
    RB_ORIENT_DIAG_READ1_LEFT,
    RB_ORIENT_DIAG_READ1_RIGHT,
    RB_ORIENT_DIAG_READ1_SAME_POS,
    RB_ORIENT_DIAG_READ2_LEFT,
    RB_ORIENT_DIAG_READ2_RIGHT,
    RB_ORIENT_DIAG_READ2_SAME_POS,
    RB_ORIENT_DIAG_ISIZE_NEG,
    RB_ORIENT_DIAG_ISIZE_ZERO,
    RB_ORIENT_DIAG_ISIZE_POS,
    RB_ORIENT_DIAG_COUNT
};

struct MpiStatsBasicCountSlice {
    long long values[RB_STATS_LONG_COUNT];
    long long isize_inward[RB_MPI_STATS_INSERT_BINS];
    long long isize_outward[RB_MPI_STATS_INSERT_BINS];
    long long isize_other[RB_MPI_STATS_INSERT_BINS];
    long long orient_diag[RB_ORIENT_DIAG_COUNT];
};

struct MpiStatsBlockSortState {
    long long has_coord;
    long long sorted;
    long long first_tid;
    long long first_pos;
    long long last_tid;
    long long last_pos;
};

struct MpiStatsBasicCountPara {
    int block_id;
    bam_block *input_block;
    bam_block *un_comp_block;
    unsigned char *scratch_data;
    size_t scratch_capacity;
    MpiStatsBasicCountSlice *counts;
    MpiStatsBlockSortState sort_state;
    int collect_diag;
    int n_total_records;
    int status;
    int record_index;
    long long actual_value;
    long long limit_value;
    int limit_id;
    uint64_t decomp_alloc_cycles;
    uint64_t decomp_inflate_cycles;
    uint64_t decomp_crc_cycles;
    uint64_t decomp_parse_cycles;
    uint64_t decomp_total_cycles;
};

struct BoundsCheckError {
    const char *pipeline;
    const char *stage;
    const char *limit_name;
    long long limit_value;
    long long actual_value;
    int block_id;
    int chunk_id;
    int record_index;
    int core_id;
};


typedef struct {
    const char *src_ptr;                 
    size_t src_len;                      
    char *text_buf;                      
    size_t text_len;                     
    bam1_t *bams[MAX_BAMS_PER_CHUNK];
    uint32_t *bam_lens;             
    int count;                           
} SamParseChunk;

typedef struct {
    const sam_hdr_t *hdr;
    SamParseChunk chunks[64];            
} SamParseBatch; 

struct MpiSamParseChunk {
    const char *src_ptr;
    size_t src_len;
    char *text_buf;
    size_t text_buf_capacity;
    int owns_text_buf;
    size_t text_len;
    bam1_t *bams;
    unsigned char *bam_data;
    size_t bam_data_capacity;
    size_t bam_data_used;
    uint32_t *bam_offsets;
    uint32_t *bam_lens;
    int count;
    int status;
    int record_index;
    long long actual_value;
    long long limit_value;
    int limit_id;
    size_t max_line_len;
    size_t estimated_bam_data;
    int parse_fast_records;
    int parse_fallback_records;
    uint64_t parse_total_cycles;
    uint64_t parse_core_cycles;
    uint64_t parse_aux_cycles;
    uint64_t parse_cg_cycles;
    uint64_t parse_fallback_cycles;
};

struct MpiSamParseBatch {
    const sam_hdr_t *hdr;
    MpiSamParseChunk chunks[64];
};

typedef struct {
    const sam_hdr_t *hdr;
    bam1_t *bams[BATCH_SIZE];
    kstring_t sam_lines[BATCH_SIZE];  
    int count;   
} SamParseByteBatch;



struct MemReader {
    char *base;
    size_t size;
    size_t pos;
};
struct MemWriter {
    char *data;
    size_t size;
    size_t capacity;
};

struct Para {
    int block_id;       

    bam_block* input_block;   
    bam_block* un_comp_block;

    int decompress_size;      
    std::vector<bam1_t*> output_records;  
    int n_records;             
    int status;               
};

struct BamFilterOptions {
    int min_mapq;
    int max_mapq;
    uint32_t require_flag;
    uint32_t exclude_flag;
    int ref_tid;
    int min_read_len;
    int max_read_len;
};

struct Bam2BamPara {
    int block_id;
    bam_block *input_block;
    bam_block *un_comp_block;
    bam1_t **output_records;
    bam1_t *record_base;
    uint32_t *bam_lens;
    int record_capacity;
    unsigned char *data_arena;
    size_t data_arena_capacity;
    size_t data_arena_used;
    BamFilterOptions filter;
    int n_total_records;
    int n_kept_records;
    uint32_t kept_total_len;
    int status;
    int record_index;
    long long actual_value;
    long long limit_value;
    int limit_id;
    uint64_t decomp_alloc_cycles;
    uint64_t decomp_inflate_cycles;
    uint64_t decomp_crc_cycles;
    uint64_t decomp_parse_cycles;
    uint64_t decomp_total_cycles;
};

struct CgsBamDecodePara {
    int block_id;
    bam_block *input_block;
    bam_block *un_comp_block;
    bam1_t **output_records;
    int n_records;
    int status;
};

struct CgsSamFormatBatch {
    const sam_hdr_t *hdr;
    bam1_t **records;
    int total_records;
    kstring_t core_out_lines[CGS_NB];
    kstring_t core_line_bufs[CGS_NB];
    int status[CGS_NB];
    int formatted_records[CGS_NB];
};

struct CgsSamParseChunk {
    const char *src_ptr;
    size_t src_len;
    char *text_buf;
    size_t text_len;
    bam1_t **bams;
    uint32_t *bam_lens;
    int count;
    int status;
};

struct CgsSamParseBatch {
    const sam_hdr_t *hdr;
    CgsSamParseChunk chunks[CGS_NB];
};

struct CheckedPara {
    int block_id;
    bam_block *input_block;
    bam_block *un_comp_block;
    bam1_t **output_records;
    int n_records;
    int status;
    int record_index;
    long long actual_value;
    long long limit_value;
    int limit_id;
};

struct CheckedBam2BamPara {
    int block_id;
    bam_block *input_block;
    bam_block *un_comp_block;
    bam1_t **output_records;
    uint32_t *bam_lens;
    BamFilterOptions filter;
    int n_total_records;
    int n_kept_records;
    int status;
    int record_index;
    long long actual_value;
    long long limit_value;
    int limit_id;
};

struct Comp_Para {
    int block_id;         
    bam1_t **input_records; 
    int n_records;       

    bam_block *un_comp_block;         
    int un_comp_size;            

    bam_block *output_block; 
    int output_size;     
    int status;          
    int compress_level;
    uint64_t compress_serialize_cycles;
    uint64_t compress_alloc_cycles;
    uint64_t compress_deflate_cycles;
    uint64_t compress_footer_cycles;
    uint64_t compress_total_cycles;
};

struct bam_block {
    unsigned int errcode;
    unsigned char *data;
    unsigned int length;
    unsigned int pos;      //记录记录在块中的偏移量，即当前记录的读取位置，在解析时使用
    int64_t block_address; //该块在整个文件中的偏移位置
    int block_id;          //块的编号
};



// Aligned memory allocation helper function
// Allocates size bytes aligned to alignment (must be power of 2)
// Returns aligned pointer, stores original pointer before aligned address for deallocation
inline unsigned char* aligned_alloc_custom(size_t alignment, size_t size) {
    // Allocate extra space for alignment and storing original pointer
    size_t total_size = size + alignment + sizeof(void*);
    unsigned char* raw_ptr = new unsigned char[total_size];
    if (!raw_ptr) {
        return nullptr;
    }
    
    // Calculate aligned address
    uintptr_t raw_addr = (uintptr_t)raw_ptr;
    uintptr_t aligned_addr = (raw_addr + sizeof(void*) + alignment - 1) & ~(uintptr_t)(alignment - 1);
    unsigned char* aligned_ptr = (unsigned char*)aligned_addr;
    
    // Store original pointer before aligned address
    void** ptr_storage = (void**)(aligned_ptr - sizeof(void*));
    *ptr_storage = raw_ptr;
    
    return aligned_ptr;
}

// Free aligned memory allocated by aligned_alloc_custom
inline void aligned_free_custom(unsigned char* aligned_ptr) {
    if (aligned_ptr) {
        // Retrieve original pointer
        void** ptr_storage = (void**)(aligned_ptr - sizeof(void*));
        unsigned char* raw_ptr = (unsigned char*)(*ptr_storage);
        delete[] raw_ptr;
    }
}

void bam_destroy1_sw(bam1_t *b);

typedef struct {
    int size;
    uint8_t *block;
    int64_t end_offset;
} cache_t;
KHASH_MAP_INIT_INT64(cache, cache_t
)
static const uint8_t g_magic[19] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0";

struct bgzf_cache_t {
    khash_t(cache) *h;
    khint_t last_pos;
};

struct BamCrossBlockStats {
    long long total_records;
    long long cross_block_records;
};


#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#define KS_SEP_MAX   2


int check_bam_cross_block(const char *bam_path);
int check_bam_cross_block_ex(const char *bam_path, BamCrossBlockStats *stats, bool verbose);
void print_bam1(const bam1_t *b);
void print_bam_block(struct bam_block *blk) ;


int read_block(BGZF *fp, struct bam_block *j);

int sam_write1_sw(samFile *fp, const sam_hdr_t *h, const bam1_t *b);

int sam_realloc_bam_data(bam1_t *b, size_t desired);

int realloc_bam_data(bam1_t *b, size_t desired);

const char *bgzf_zerr(int errnum, z_stream *zs);

int sam_read1_sw(samFile *fp, sam_hdr_t *h, bam1_t *b);

void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host);

void bam_cigar2rqlens(int n_cigar, const uint32_t *cigar,
                      hts_pos_t *rlen, hts_pos_t *qlen);

inline int possibly_expand_bam_data(bam1_t *b, size_t bytes);

int bam_tag2cigar(bam1_t *b, int recal_bin,
                  int give_warning); 

int fixup_missing_qname_nul(bam1_t *b);

inline bool bam_filter_matches(const bam1_t *b, const BamFilterOptions &filter) {
    const bam1_core_t &core = b->core;

    if (filter.min_mapq >= 0 && core.qual < filter.min_mapq) return false;
    if (filter.max_mapq >= 0 && core.qual > filter.max_mapq) return false;
    if (filter.require_flag != 0 && (core.flag & filter.require_flag) != filter.require_flag) return false;
    if (filter.exclude_flag != 0 && (core.flag & filter.exclude_flag) != 0) return false;
    if (filter.ref_tid != -2 && core.tid != filter.ref_tid) return false;
    if (filter.min_read_len >= 0 && core.l_qseq < filter.min_read_len) return false;
    if (filter.max_read_len >= 0 && core.l_qseq > filter.max_read_len) return false;

    return true;
}

inline bool bam_filter_is_noop(const BamFilterOptions &filter) {
    return filter.min_mapq < 0 &&
           filter.max_mapq < 0 &&
           filter.require_flag == 0 &&
           filter.exclude_flag == 0 &&
           filter.ref_tid == -2 &&
           filter.min_read_len < 0 &&
           filter.max_read_len < 0;
}

inline void packInt16(uint8_t *buffer, uint16_t value) {
    buffer[0] = value;
    buffer[1] = value >> 8;
}

inline void packInt32(uint8_t *buffer, uint32_t value) {
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}


#endif // BAMTOOLS_H
