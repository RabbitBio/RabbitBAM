#ifndef BAMCOMPLETE_H
#define BAMCOMPLETE_H

#include <vector>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <cstring>
#include <thread>
#include <unistd.h>

#include "BamTools.h"
#include <htslib/sam.h>  //bam_init1();  bam_destroy1();
#include <htslib/hts.h>

class BamComplete {
public:
    BamComplete(int queue_size);
    ~BamComplete();

    bam_block* getBuffer(int idx);
    std::vector<bam1_t*>& getResultBuf(int idx);

    bam1_t* getEmpty();
    void getEmptyBatch(bam1_t** dst, int n);
    void backBam1_t(bam1_t* bam1);
    void backBam1_tBatch(bam1_t** src, int n);

    void inputBam1_t(bam1_t* &bam1);
    bam1_t* getBam1_t();

    void markComplete();
    bool isComplete() const;

private:
    //队列容量（结果队列中bam1_t条数）
    int queue_size_;

    //解压用缓冲区（批量分配）
    bam_block  *bp_flat_;
    uint8_t    *bp_data_;
    std::vector<bam_block*> buffer_pool_;

    //结果缓冲区（批量分配）
    bam1_t  *rp_flat_;
    uint8_t *rp_data_;
    std::vector<std::vector<bam1_t*>> result_pool_;

    //内存池（批量分配）
    bam1_t  *pq_flat_;
    uint8_t *pq_data_;

    bam1_t** consumer_queue_;
    int con_bg;
    int con_ed;
    int con_queueSizeLim;

    bam1_t** producer_queue_;
    int pro_bg;
    int pro_ed;
    int pro_queueSizeLim;

    bool complete_flag;
};

#endif
