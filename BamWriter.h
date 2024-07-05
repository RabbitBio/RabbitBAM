//
// Created by 赵展 on 2022/8/17.
//

#ifndef BAMSTATUS_BAMWRITER_H
#define BAMSTATUS_BAMWRITER_Hi
#include "BamTools.h"
#include <thread>
#include "BamCompress.h"
#include "BamWriteCompress.h"
#include "BamCompleteBlock.h"
/*
 *
 * 这个类为对外的输出接口
 *
 * 根据给出的地址，输出BAM文件，并且设置压缩等级
 *
 * 本类使用异步输出，所以需要在析构函数的时候才会完整真正的全部输出
 *
 * 如果需要清空输出，则需要使用over的类内函数进行清空缓存区
 *
 *
 * 需要先输出HDR部分再进行处理
 *
 */


int rabbit_write_deflate_block(BGZF *fp, bam_write_block* write_block);

int rabbit_bgzf_flush(BGZF *fp,bam_write_block* write_block);

int rabbit_bgzf_mul_flush(BGZF *fp,BamWriteCompress *bam_write_compress,bam_write_block* &write_block);
int rabbit_bgzf_write(BGZF *fp,bam_write_block* &write_block,const void *data, size_t length);
int rabbit_bgzf_mul_write(BGZF *fp, BamWriteCompress *bam_write_compress,bam_write_block* &write_block,const void *data, size_t length);
int rabbit_bgzf_flush_try(BGZF *fp, bam_write_block* write_block,ssize_t size);
int rabbit_bgzf_mul_flush_try(BGZF *fp,BamWriteCompress* bam_write_compress,bam_write_block* &write_block,ssize_t size);
int bam_write_pack(BGZF *fp,BamWriteCompress *bam_write_compress);
void bam_write_compress_pack(BGZF *fp,BamWriteCompress *bam_write_compress);

int rabbit_bam_write_test(BGZF *fp,bam_write_block* write_block,bam1_t *b);
int rabbit_bam_write_mul_test(BGZF *fp,BamWriteCompress *bam_write_compress,bam_write_block* &write_block,bam1_t *b);

void benchmark_write_pack(BamCompleteBlock* completeBlock,samFile *output,sam_hdr_t *hdr,int level);

void benchmark_write_mul_pack(BamCompleteBlock* completeBlock,BamWriteCompress* bam_write_compress,samFile *output,sam_hdr_t *hdr,int level);



class BamWriter {

public:
    BamWriter(int threadNumber=1, int level=6,int BufferSize=200);
    BamWriter(std::string file_name, int threadNumber=1, int level=6,int BufferSize=200);
//    BamWriter(std::string file_name,sam_hdr_t *hdr, int level=6);
//    BamWriter(std::string file_name,sam_hdr_t *hdr, int threadNumber, int level=6);
    BamWriter(std::string file_name,sam_hdr_t *hdr, int threadNumber=1, int level=6,int BufferSize=200);
    ~BamWriter();

    void bam_write(bam1_t* b);

    void set_output(samFile *output);

    void hdr_write(sam_hdr_t* hdr);

    void write(bam1_t* b);

    void over();


private:

    std::thread **write_compress_thread;

    std::thread *write_output_thread;

    BamWriteCompress *bam_write_compress;


    bam_write_block *write_block;

    samFile *output;

    int n_thread_write;

    int write_num=0;

};


#endif //BAMSTATUS_BAMWRITER_H

//export PATH=/home/user_home/zz/API/BamUniversalStatus/build:$PATH
//export LD_LIBRARY_PATH/home/user_home/zz/API/BamUniversalStatus/build:$LD_LIBRARY_PATH