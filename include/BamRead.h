#ifndef BAMREAD_H
#define BAMREAD_H

#include <vector> 
#include <unistd.h>
#include <atomic>

#include "BamTools.h"

class BamRead {
public:
    BamRead(int blocks_size, int queue_size );  
    ~BamRead();

    bam_block* getEmpty();  
    void backBlock(bam_block* block);  

    void inputBlock64(std::vector<bam_block*>& group);  
    std::vector<bam_block*> getBlock64();  

    bool isComplete() const;
    void markComplete();

private:
    std::vector<bam_block*> *consumer_queue_; 
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