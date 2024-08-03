#ifndef BAMSTATUS_BAMCOMPLETEBLOCK_H
#define BAMSTATUS_BAMCOMPLETEBLOCK_H

#include "BamTools.h"
#include <thread>

class BamCompleteBlock {
public:


    BamCompleteBlock();

    void resize(int BufferSize);

    BamCompleteBlock(int BufferSize);

    bam_complete_block *getEmpty();

    void inputCompleteBlock(bam_complete_block *block);

    bam_complete_block *getCompleteBlock();

    void backEmpty(bam_complete_block *block);

    void is_over() { BlockComplete = true; }

private:
    bam_complete_block **complete_data;
    int complete_bg;
    int complete_ed;
    int complete_size;
    std::mutex mtx_complete;
    bam_complete_block **consumer_data;
    int consumer_bg;
    int consumer_ed;
    int consumer_size;
    std::mutex mtx_consumer;

    bool BlockComplete;

};


#endif //BAMSTATUS_BAMCOMPLETEBLOCK_H
