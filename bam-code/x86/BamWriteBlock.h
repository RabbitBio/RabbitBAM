#ifndef BAMSTATUS_BAMWRITEBLOCK_H
#define BAMSTATUS_BAMWRITEBLOCK_H

#include "BamTools.h"

class BamWriteBlockConfig {
public:
    BamWriteBlockConfig();

    BamWriteBlockConfig(int Buffer_number);

public:
    int Buffer_number;
    int write_number;
    int complete;
};

class BamWriteBlock {
public:
    BamWriteBlock();

    BamWriteBlock(BamWriteBlockConfig *config);

    std::pair<bam_block *, int> getEmpty();

    void inputblock(int id);

    std::pair<bam_block *, int> getCompressdata();

    void backempty(int id);

    bool isComplete();

    void ReadComplete();

public:
    BamWriteBlockConfig *config;
    std::mutex mtx_read;
    std::mutex mtx_compress;
    bam_block **buffer;
    int *compress;
    int compress_bg;
    int compress_ed;
    int *read;
    int read_bg;
    int read_ed;
};


#endif //BAMSTATUS_BAMWRITEBLOCK_H
