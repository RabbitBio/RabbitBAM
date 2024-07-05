//
// Created by 赵展 on 2022/7/6.
//

#ifndef BAMSTATUS_BAMWRITECOMPRESS_H
#define BAMSTATUS_BAMWRITECOMPRESS_H
#include "BamTools.h"
#include <thread>
#include <atomic>
/*
 *  功能上为一个对于处理完成的
 */
class BamWriteCompress {
public:
    BamWriteCompress(int BufferSize,int threadNumber);

    /*
     * 获取一个空白的内存块
     * 输入：无
     * 输出：bam_write_block* : 指向该内存块的指针
     */
    bam_write_block* getEmpty();
    /*
     * 插入一个装载数据的未压缩的块
     */
    void inputUnCompressData(bam_write_block* data);

    bam_write_block* getUnCompressData();

    /*
     * 按照顺序插入输出队列
     * 输入：  bam_block* : 解压完成的数据
     *        int ：读入的顺序编号
     * 输出：  无
     */
    void inputCompressData(bam_write_block* data);
//    /*
//     * 按照顺序尝试插入输出队列
//     * 输入：  bam_block* : 解压完成的数据
//     *        int ：读入的顺序编号
//     * 输出：  bool ： true 插入成功 false 插入失败
//     */
//    bool tryinputUnCompressData(bam_block* data,int block_num);
    /*
     * 按照顺序获取压缩完成的数据
     * 输入：无
     * 输出：bam_block* 解压完成的数据
     */
    bam_write_block* getCompressData();
    /*
     * 返还使用完毕的内存块
     *
     */
    void backEmpty(bam_write_block* data);

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
    void WriteComplete(){
        isWriteComplete=true;
//        printf("Input Uncompress Data Over!\n"
//               "Read Complete need compress bg : %d\n"
//               "Read Complete need compress ed : %d\n",need_compress_bg,need_compress_ed);
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
    bam_write_block** compress_data;
    int compress_bg;
    int compress_ed;
    int compress_size;
    std::mutex mtx_compress;
    std::mutex mtx_input;

    /*
     * 需要压缩处理处理的数据队列
     */
    int blockInputNum; // 根据插入时间增加时间戳
    std::atomic<int> blockInputPos;
    bam_write_block** need_compress_data;
    int need_compress_bg;
    int need_compress_ed;
    int need_compress_size;
    std::mutex mtx_need_compress;

    /*
     *  输出队列的管理部分
     */
    bam_write_block** consumer_data;
    int consumer_bg;
    int consumer_ed;
    int consumer_size;
    bool *is_ok;

    /*
     *  管理是否输出完成
     */
    int compressThread;
    std::mutex mtx_compressThread;

    /*
     * 检测需要输出的文件是否已经完全进入队列
     */
    bool isWriteComplete=false;

};


#endif //BAMSTATUS_BAMWRITECOMPRESS_H
