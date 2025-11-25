#ifndef BAMCOMPLETE_H
#define BAMCOMPLETE_H

#include "BamTools.h"
#include <vector>
#include <atomic>
#include <mutex>
#include <condition_variable>

class BamComplete {
public:
    BamComplete(int queue_size);
    ~BamComplete();

    // 从核缓冲区（解压用）
    bam_block* getBuffer(int idx);
    // 从核解析bam1_t缓冲区
    bam1_t* getResultBuf(int idx);

    // 从核完成一块，放入结果队列（按顺序）
    void pushBlockResults(const bam1_t* records, int n);

    // 主线程逐条读取bam1_t
    bool popRecord(bam1_t& out);

    bam1_t* BamComplete::backRecord();

    void markComplete();
    bool isComplete() const;

private:
    int max_records_per_block_;     // 每块最大bam1_t数量
    int queue_size_;                // 队列容量（bam1_t条数）

    // 解压用缓冲区
    std::vector<bam_block*> buffer_pool_;
    //内存池
    std::vector<bam1_t*> result_pool_;

    bam1_t** p_out_queue_;  
    std::atomic_int p_queueP2;     // 写指针
    std::atomic_int p_queueC2;     // 读指针
    std::atomic_int p_queueNumNow; // 当前队列中条目数
    std::atomic_int p_queueCanRead;   // 当前可读条目数
    int p_queueSizeLim;            // 队列总容量

    std::atomic_bool complete_flag; // 标记是否全部完成
};

#endif
