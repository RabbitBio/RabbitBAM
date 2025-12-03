#ifndef BAM_WRITE_COMPLETE_H
#define BAM_WRITE_COMPLETE_H

#include <unistd.h>
#include <vector>
#include <cstdio>

#include <sunway/BamTools.h>

class BamWriteComplete {
public:
    BamWriteComplete(int queue_size);
    ~BamWriteComplete();

    bam_block* getBuffer(int idx);

    bam_block* getEmpty();        
    void backBlock(bam_block *b);

    void pushCompressedBlock(bam_block* block); 
    bam_block* getCompressedBlock();       

    bool isComplete() const;
    void markComplete();

private:

    std::vector<bam_block*> buffer_pool_;

    bam_block **blockPool;
    int pool_bg, pool_ed, pool_size;

    bam_block **queue;
    int q_bg, q_ed, q_size;

    bool complete_flag;
};

#endif // BAM_WRITE_COMPLETE_H
