#ifndef BAMSTATUS_BAMREAD_H
#define BAMSTATUS_BAMREAD_H

#include "BamTools.h"
#include <thread>
#include <mutex>


class BamRead {

public:

    BamRead();

    BamRead(int BufferSize);

    void resize(int BufferSize);

    bam_block *getEmpty();

    void inputBlock(bam_block *block);

    std::pair<bam_block *, int> getReadBlock();

    void backBlock(bam_block *block);

    void ReadComplete();

    bool isComplete() { return read_complete; }

private:

    bam_block **readBlock;
    int read_bg;
    int read_ed;
    int readBlockSize;
    bam_block **consumerBlock;
    int consumer_bg;
    int consumer_ed;
    int consumerBlockSize;
    std::atomic<int> blockNum;
    int blockTot;
    std::mutex mtx_consumer;
    std::mutex mtx_read;
    bool read_complete;
};


#endif //BAMSTATUS_BAMREAD_H
