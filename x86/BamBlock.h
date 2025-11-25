#ifndef BAM_BAMBLOCK_H
#define BAM_BAMBLOCK_H

#include "BamTools.h"
#include <thread>

using namespace std;

class BamBlockConfig {
public:
    BamBlockConfig();

    BamBlockConfig(int Buffer_number);

public:
    int Buffer_number;
    int write_number;
    int complete;
};

class BamBlock {
public:
    BamBlock();

    BamBlock(BamBlockConfig *config);

    pair<bam_block *, int> getEmpty();

    void inputblock(int id); // 导入未解压的数据
    pair<bam_block *, int> getCompressdata();

    void backempty(int id);

    bool isComplete();

    void ReadComplete();

public:
    BamBlockConfig *config;
    mutex mtx_read;
    mutex mtx_compress;
    bam_block **buffer;
    int *compress;
    int compress_bg;
    int compress_ed;
    int *read;
    int read_bg;
    int read_ed;
};


#endif //BAM_BAMBLOCK_H
