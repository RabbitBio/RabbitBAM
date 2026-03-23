#ifndef BAMWRITE_H
#define BAMWRITE_H

#include <vector>
#include <htslib/sam.h>  
#include <unistd.h>
#include <cstdio>
#include "BamTools.h"

class BamWrite {
public:
    BamWrite(int queue_size);
    ~BamWrite();

    bam1_t* getEmpty();
    void getEmptyBatch(bam1_t** dst, int n);
    void backBam(bam1_t *b);
    void backBamBatch(bam1_t** const records, int n);

    void inputGroup(std::vector<std::vector<bam1_t*>> group);
    std::vector<std::vector<bam1_t*>> getGroup();

    bool isComplete() const;
    void markComplete();

private:

    bam1_t **bamPool;
    int pool_bg, pool_ed, pool_size;

    std::vector<std::vector<bam1_t*>> *groupQueue;
    int q_bg, q_ed, q_size;

    bool write_complete;

    bam1_t  *bam1_struct_pool_;  // 所有 bam1_t 结构体的连续内存
    uint8_t *data_pool_;         // 所有 data 缓冲区的连续 64 字节对齐内存
};

#endif
