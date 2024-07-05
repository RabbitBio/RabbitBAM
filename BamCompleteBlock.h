//
// Created by 赵展 on 2022/1/22.
//

#ifndef BAMSTATUS_BAMCOMPLETEBLOCK_H
#define BAMSTATUS_BAMCOMPLETEBLOCK_H

#include "BamTools.h"
#include <thread>
class BamCompleteBlock {
public:

   /*
    * 只初始化该类
    */
    BamCompleteBlock();

    /*
     * 重新根据BufferSize进行初始化
     */
    void resize(int BufferSize);




    /*
     * 根据输入的BufferSize进行初始化
     *
     */
    BamCompleteBlock(int BufferSize);

    /*
     * 获取一个空的内存块
     * 输入：无
     * 输出：bam_complete_block* ： 输出一个空白内存块的指针
     *
     */
    bam_complete_block* getEmpty();
    /*
     * 输入一个重新划分之后的内存块，到消费队列中
     * 输入： bam_complete_block* : 重新划分的内存块
     * 输出： 无
     *
     */
    void inputCompleteBlock(bam_complete_block* block);
    /*
     * 获取一个重新划分好的内存块
     * 输入：无
     * 输出：bam_complete_block* : 重新划分好的内存块
     */
    bam_complete_block* getCompleteBlock();

    /*
     * 归还使用完毕的内存块
     * 输入：bam_complete_block* ：使用完毕的内存块
     * 输出：无
     */
    void backEmpty(bam_complete_block* block);
    /*
     * 划分线程已经全部划分完毕
     */
    void is_over(){BlockComplete=true;}
private:
    /*
     * 重新划分使用的内存池
     */
    bam_complete_block** complete_data;
    int complete_bg;
    int complete_ed;
    int complete_size;
    std::mutex mtx_complete;

    /*
     *  消费者输出使用的队列
     */
    bam_complete_block** consumer_data;
    int consumer_bg;
    int consumer_ed;
    int consumer_size;
    std::mutex mtx_consumer;

    bool BlockComplete;

};


#endif //BAMSTATUS_BAMCOMPLETEBLOCK_H
