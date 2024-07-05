//
// Created by 赵展 on 2022/1/21.
//

#ifndef BAMSTATUS_BAMCOMPRESS_H
#define BAMSTATUS_BAMCOMPRESS_H
#include "BamTools.h"
#include <thread>
#include <atomic>

/*
 *  给出解压空间，给出有序的解压内容
 */
/*
 * TODO
 * 目前这种写法存在大问题，但是单一队列很难解决，现在考虑对于block_num取模，然后进行处理 一个方向 可以考虑
 */
class BamCompress {
public:

    /*
     * 只初始化函数
     */
    BamCompress();


    /*
     * 重新分配BamCompress
     */
    void resize(int BufferSize,int threadNumber);


    /*
     * 直接按照BufferSize 进行初始化
     */

    BamCompress(int BufferSize,int threadNumber);

    /*
     * 获取一个空白的内存块
     * 输入：无
     * 输出：bam_block* : 指向该内存块的指针
     */
    bam_block* getEmpty();
    /*
     * 按照顺序插入输出队列
     * 输入：  bam_block* : 解压完成的数据
     *        int ：读入的顺序编号
     * 输出：  无
     */
    void inputUnCompressData(bam_block* data,int block_num);
//    /*
//     * 按照顺序尝试插入输出队列
//     * 输入：  bam_block* : 解压完成的数据
//     *        int ：读入的顺序编号
//     * 输出：  bool ： true 插入成功 false 插入失败
//     */
//    bool tryinputUnCompressData(bam_block* data,int block_num);
    /*
     * 按照顺序获取解压完成的数据
     * 输入：无
     * 输出：bam_block* 解压完成的数据
     */
    bam_block* getUnCompressData();
    /*
     * 返还使用完毕的内存块
     *
     */
    void backEmpty(bam_block* data);

    /*
     * 汇报该解压线程完成解压工作
     * 输入：无
     * 输出：无
     */
    void CompressThreadComplete(){
        mtx_compressThread.lock();
        compressThread--;
        mtx_compressThread.unlock();
    }
public:
    /*
     * 记录等待次数
     */
    int wait_num;
private:
    /*
     * 管理插入序列编号
     */
    std::atomic<int> blockNum;

    /*
     *  空白内存的管理部分
     */
    bam_block** compress_data;
    int compress_bg;
    int compress_ed;
    int compress_size;
    std::mutex mtx_compress;
    std::mutex mtx_input;
    /*
     *  输出队列的管理部分
     */
    bam_block** consumer_data;
    int consumer_bg;
    int consumer_ed;
    int consumer_size;
    bool *is_ok;

    /*
     *  管理是否输出完成
     */
    int compressThread;
    std::mutex mtx_compressThread;

};


#endif //BAMSTATUS_BAMCOMPRESS_H
