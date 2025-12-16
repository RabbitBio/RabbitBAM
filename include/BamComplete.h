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
    void backBam1_t(bam1_t* bam1);  

    void inputBam1_t(bam1_t* &bam1);  
    bam1_t* getBam1_t();  

    void markComplete();
    bool isComplete() const;

private:


    //队列容量（结果队列中bam1_t条数）
    int queue_size_;                
    //bam1_t 的data最大长度 暂定1KB
    const size_t INIT_DATA_SIZE = 1024;  
    //每块最大bam1_t数量，暂定1024块
    const size_t MAX_RECORDS_PER_BLOCK = 1024; 


    //解压用缓冲区
    std::vector<bam_block*> buffer_pool_;
    std::vector<std::vector<bam1_t*>> result_pool_;


    //内存池实现
    bam1_t** consumer_queue_; 
    int con_bg;
    int con_ed;
    int con_queueSizeLim;

    bam1_t ** producer_queue_;
    int pro_bg;
    int pro_ed;
    int pro_queueSizeLim;

    std::atomic_bool complete_flag;


};

#endif
