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
//#define BGZF_MAX_BLOCK_SIZE 0x10000
//#define BGZF_MAX_BLOCK_COMPLETE_SIZE 0x40000
//#define BGZF_MAX_BLOCK_COMPLETE_SIZE 0x10000
//#define THREAD_NUM_P 6


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


typedef struct bam_block bam_block;


#define MAX_SAM_LINE_SIZE 8192   // 8 KB（足够覆盖绝大多数）
#define BATCH_PER_CORE 1024
#define BATCH_SIZE (64 * BATCH_PER_CORE)

#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#define KS_SEP_MAX   2

typedef struct {
    const sam_hdr_t *hdr;
    bam1_t *bams[BATCH_SIZE];
    kstring_t sam_lines[BATCH_SIZE];  // 每条 BAM 对应一条 SAM 字符串
    int count;   // 实际 bam 数
} SamFormatBatch;


struct Para {
    int block_id;              // 当前块号
    bam_block* input_block;    // 压缩块指针
    bam_block* un_comp_block;
    //uint8_t* decompress_buf;   // 解压缓冲区
    int decompress_size;       // 解压后大小
    std::vector<bam1_t*> output_records;  // 解析得到的bam1_t数组
    int n_records;             // 解析得到的条目数
    int status;                // 处理状态标志（0=ok，非0=失败）

    int l_data_list[1024];     // 每条记录的长度（假设一块最多1024条）
    uint8_t *data_list[1024];  // 每条记录的 data 指针
};

struct Comp_Para {
    int block_id;         // 本组中索引（0..63）
    bam1_t **input_records; // 指向该块内的 bam1_t* 数组（MPE 分配并传给 CPE）
    int n_records;        // 该块中记录数

    bam_block *un_comp_block;         // 未压缩数据缓存
    int un_comp_size;            // 未压缩数据大小

    bam_block *output_block; // 输出压缩块（MPE 传给 CPE，用于写入压缩数据）
    int output_size;      // CPE 返回：压缩后长度（字节）
    int status;           // 0 = success, <0 = error, -1 = empty
};

struct bam_block {
    unsigned int errcode;
    //unsigned char data[BGZF_MAX_BLOCK_SIZE];//0x1000
    unsigned char *data;
    unsigned int length;
    unsigned int pos;  //记录记录在块中的偏移量，即当前记录的读取位置，在解析时使用
    int64_t block_address; //该块在整个文件中的偏移位置

    int block_id; //块的编号

    //unsigned int split_pos;
    //unsigned int bam_number;
};


// struct bam_write_block {
//     int status = 0; // 0:uncompress 1:compress
//     int block_num = -1; // -1 : not used correctly
//     int block_offset;
//     int block_length;
//     uint8_t *uncompressed_data;
//     uint8_t *compressed_data;
// };


void print_bam1(const bam1_t *b);
void print_bam_block(struct bam_block *blk) ;

//bam to sam
int read_block(BGZF *fp, struct bam_block *j);

int sam_write1_sw(samFile *fp, const sam_hdr_t *h, const bam1_t *b);

int sam_realloc_bam_data(bam1_t *b, size_t desired);

int realloc_bam_data(bam1_t *b, size_t desired);




//sam to bam
const char *bgzf_zerr(int errnum, z_stream *zs);

int sam_read1_sw(samFile *fp, sam_hdr_t *h, bam1_t *b);


//从核使用的
void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host);

void bam_cigar2rqlens(int n_cigar, const uint32_t *cigar,
                      hts_pos_t *rlen, hts_pos_t *qlen);

inline int possibly_expand_bam_data(bam1_t *b, size_t bytes);

int bam_tag2cigar(bam1_t *b, int recal_bin,
                  int give_warning); // return 0 if CIGAR is untouched; 1 if CIGAR is updated with CG

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