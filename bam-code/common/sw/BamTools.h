#ifndef BAMTOOLS_H
#define BAMTOOLS_H


struct Para {
    int block_id;              // 当前块号
    bam_block* input_block;    // 压缩块指针
    bam_block* un_comp_block;
    //uint8_t* decompress_buf;   // 解压缓冲区
    int decompress_size;       // 解压后大小
    bam1_t* output_records;    // 解析得到的bam1_t数组
    int n_records;             // 解析得到的条目数
    int status;                // 处理状态标志（0=ok，非0=失败）

    int l_data_list[8192];     // 每条记录的长度（假设一块最多8192条）
    uint8_t *data_list[8192];  // 每条记录的 data 指针
};

//--------------------------------------------------------------------------------------------------------------
struct Comp_Para {
    int block_id;         // 本组中索引（0..63）
    bam1_t **input_records; // 指向该块内的 bam1_t* 数组（MPE 分配并传给 CPE）
    int n_records;        // 该块中记录数

    bam_block *un_comp_block;         // 未压缩数据缓存
    int un_comp_size;            // 未压缩数据大小

    bam_block *output_block; // 输出压缩块（MPE 传给 CPE，用于写入压缩数据）
    int output_size;      // CPE 返回：压缩后长度（字节）
    int status;           // 0 = success, <0 = error, -1 = empty
}

struct bam_block {
    unsigned int errcode;
    unsigned char data[BGZF_MAX_BLOCK_SIZE];//0x1000
    unsigned int length;
    unsigned int pos;  //当前记录的读取位置，在解析时使用
    int64_t block_address; //在整个文件中的偏移位置

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


int read_block(BGZF *fp, struct bam_block *j);

int sam_write1_sw(samFile *fp, const sam_hdr_t *h, const bam1_t *b);

int sam_read1_sw(samFile *fp, const sam_hdr_t *h, const bam1_t *b);


#endif // BAMTOOLS_H