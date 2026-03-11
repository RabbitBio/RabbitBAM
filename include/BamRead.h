#ifndef BAMREAD_H
#define BAMREADER_H

#include <vector> 
#include <unistd.h>
#include <atomic>

#include "BamTools.h"

class BamRead {
public:
    BamRead(int queue_size );  
    ~BamRead();

    bam_block* getEmpty();  // 获取一个空 block（从池中）
    void backBlock(bam_block* block);  // 回收 block

    void inputBlock64(std::vector<bam_block*>& group);  // 推入一组64个
    std::vector<bam_block*> getBlock64();  // 取出一组（64个）

    bool isComplete() const;
    void markComplete();

private:
    std::vector<bam_block*> *consumer_queue_; // 每个元素是一组64个bam_block*
    int con_bg;
    int con_ed;
    int con_queueSizeLim;

    bam_block **readBlock;
    int read_bg;
    int read_ed;
    int readBlockSize;

    bool read_complete;
};



#endif // BAMREAD_H