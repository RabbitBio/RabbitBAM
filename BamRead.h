//
// Created by 赵展 on 2022/1/21.
//

#ifndef BAMSTATUS_BAMREAD_H
#define BAMSTATUS_BAMREAD_H

#include "BamTools.h"
#include <thread>
#include <mutex>
/*
 *  支持输出带序号的bamblock
 *  输出格式pair<bam_block* , int> 前面为读取到的bamblock ，后面为该bamblock的编号
 *  编号从0开始
 *  使用指针进行管理
 */


class BamRead {

public:

    /*
     * 只初始化但是什么也不干；
     */
    BamRead();

    /*
     * 初始化函数
     * 使用BufferSize进行管理
     */
    BamRead(int BufferSize);


    /*
     * 重新调整BufferSize
     */
    void resize(int BufferSize);

    /*
     * 获取一个空的bam_block，进行读取
     * 输入：无
     * 输出：bam_block* ：  指向一个可用的bamblock内存空间
     */
    bam_block* getEmpty();
    /*
     * 放入一个读取到的内存块
     * 输入 ： bam_block*
     * 输出 ： 无
     */
    void inputBlock(bam_block* block);
    /*
     * 获取一个读取到的内存块，及其读取到的编号
     * 输入 ： 无
     * 输出 ： pair<bam_block* , int> 前者为指向内存块的指针，后者为该内存块的编号
     */
    std::pair<bam_block*, int> getReadBlock();
    /*
     * 使用完毕之后，归还读取的内存块
     * 输入 ： bam_block* 使用完毕的内存块
     * 输出 ： 无
     */
    void backBlock(bam_block* block);


    /*
     * Read 完成所有读取
     */
    void ReadComplete();

    bool isComplete(){return read_complete;}
private:

    /*
     * 读取使用到的队列，以及管理的节点
     * 循环队列
     * 使用[]区间([闭（开)
     */
    bam_block** readBlock;
    int read_bg;
    int read_ed;
    int readBlockSize;
    /*
     * 输出使用的队列
     * 使用[]区间
     * bg ： 1
     * ed ： 0
     * consumerBlock 大小需要大于 readBlock
     */
    bam_block** consumerBlock;
    int consumer_bg;
    int consumer_ed;

    int consumerBlockSize;

    std::atomic<int> blockNum;
    int blockTot;

//    int blockNum;
    /*
     * 多线程同时使用consumer 获取block
     */
    std::mutex mtx_consumer;

    /*
     * 多线程同时归还bam_block
     */
    std::mutex mtx_read;

    /*
     * 判定是否完成
     */
    bool read_complete;
};


#endif //BAMSTATUS_BAMREAD_H
