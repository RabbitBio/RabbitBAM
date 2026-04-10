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


const size_t MAX_RECORDS_PER_BLOCK = 1024;   //每块最大bam1_t数量，暂定1024块
#define MAX_SAM_LINE_SIZE 8192   // 8 KB     //一条sam文本的最大长度
#define BATCH_PER_CORE 1024
#define BATCH_SIZE (64 * BATCH_PER_CORE)     //一个批次的大小

const size_t INIT_DATA_SIZE = 1024;           //bam1_t 的data最大长度 暂定1KB

#define SAM_CHUNK_SIZE (4 * 1024 * 1024)         
#define CHUNK_BUFFER_SIZE (5 * 1024 * 1024)      
#define MAX_BAMS_PER_CHUNK 12110              // 4MB最多包含的记录数

typedef struct {
    const sam_hdr_t *hdr;
    bam1_t *bams[BATCH_SIZE];
    kstring_t core_out_lines[64];  
    int count;   
} SamFormatBatch;


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

typedef struct bam_block bam_block;


struct Para {
    int block_id;       

    bam_block* input_block;   
    bam_block* un_comp_block;

    int decompress_size;      
    std::vector<bam1_t*> output_records;  
    int n_records;             
    int status;               
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


#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#define KS_SEP_MAX   2


int check_bam_cross_block(const char *bam_path);
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