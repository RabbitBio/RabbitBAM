#ifndef BAMWRITE_H
#define BAMWRITE_H

#include <vector>
#include <htslib/sam.h>   //bam_init1();  bam_destroy1();
#include <unistd.h>
#include <cstdio>
#include "BamTools.h"

class BamWrite {
public:
    BamWrite(int queue_size);
    ~BamWrite();

    bam1_t* getEmpty();
    void backBam(bam1_t *b);

    void inputGroup(const std::vector<std::vector<bam1_t*>>& group);
    std::vector<std::vector<bam1_t*>> getGroup();

    bool isComplete() const;
    void markComplete();

private:

    const size_t INIT_DATA_SIZE = 1024;  
    const size_t MAX_RECORDS_PER_BLOCK = 1024; 

    bam1_t **bamPool;
    int pool_bg, pool_ed, pool_size;

    std::vector<std::vector<bam1_t*>> *groupQueue;
    int q_bg, q_ed, q_size;

    bool write_complete;
};

#endif
