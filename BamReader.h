//
// Created by 赵展 on 2022/8/17.
//

#ifndef BAMSTATUS_BAMREADER_H
#define BAMSTATUS_BAMREADER_H
#include "BamTools.h"
#include "BamRead.h"
#include "BamCompleteBlock.h"
#include "BamCompress.h"
#include "Buffer.h"
#include <thread>

/*
 * 对外开放的API，使用该类可以直接进行读取BAM文件
 *
 * 根据文件名称进行BAM文件的读取
 *
 * 给出BAM1_t
 *
 * 并给出多线程获取无规则BAM1——t ？ 是否给出多线程BAM_BLOCK获取
 *
 *
 *
 */

class BamReader {

public:

    /*
     *  简单初始化BAMREADER 只需要提供文件名称
     */

    BamReader(std::string file_name,int n_thread=8);

    /*
     * 初始化BAMREADER，调整三个缓存池的大小
     *
     * 大小调整策略
     *
     * BamBlock 和 BamCompress 使用的用于传递Block块的缓冲池
     *
     * BamCompress 和 BamCompleteBlock 使用的用于重新划分Block的缓冲池
     *
     * 
     */
     BamReader(std::string file_name,int read_block,int compress_block,int compress_complete_block,int n_thread=8);

    /*
     *  销毁这个类
     */
     ~BamReader(){
        read_thread->join();
        for (int i=0;i<n_thread;i++) compress_thread[i]->join();
        assign_thread->join();
     }

     /*
     *  读取 SAMHDR
     */
    sam_hdr_t* getHeader();

    /*
     *  获取 BAM1t
     *
     *  根据输入进来的bam1_t的指针，装入bam1_t
     *
     *  bool
     *
     *
     */

    bool getBam1_t(bam1_t* b);


    /*
     * 获取Bam——Complete——Clock
     *
     * 针对想要高性能的开发者，返回未被分割的bam1_t，但是保证内部存在完整的bam1_t
     */
    bam_complete_block* getBamCompleteClock();
    void backBamCompleteBlock(bam_complete_block *un_comp);


private:
    BamRead *read;
    BamCompress *compress;
    BamCompleteBlock *completeBlock;


    samFile *sin;
    sam_hdr_t *hdr;
    samFile *output;

    std::thread *read_thread;
    std::thread **compress_thread;
    std::thread *assign_thread;


    bam_complete_block* un_comp;

    int n_thread;

};


#endif //BAMSTATUS_BAMREADER_H
