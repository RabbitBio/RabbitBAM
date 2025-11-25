#include "BamComplete.h"
#include <cstring>
#include <thread>

BamComplete::BamComplete(int queue_size)
{
    max_records_per_block_ = 8192;
    queue_size_ = queue_size;
    p_queueC2 = 0;
    p_queueP2 = 0;
    p_queueNumNow = 0;
    p_queueCanRead = 0;
    p_queueSizeLim = queue_size;
    complete_flag = false;

    // 初始化缓冲区
    buffer_pool_.resize(64);
    result_pool_.resize(64);

    for (int i = 0; i < 64; ++i) {
        buffer_pool_[i] = new bam_block;
        result_pool_[i] = new bam1_t[max_records_per_block_];
        for(int j = 0; j < max_records_per_block_; j++) {
            result_pool_[i][j] = bam_init1();
        }
    }

    // 初始化结果循环队列
    p_out_queue_ = new bam1_t*[queue_size_];
    for(int i = 0; i < queue_size_; i++) {
        p_out_queue_[i] = bam_init1();
    }
}

BamComplete::~BamComplete() {
    for (auto buf : buffer_pool_) delete buf;
    for (auto res : result_pool_) bam_destroy1(res);

    for (int i = 0; i < queue_size_; ++i) {
        bam_destroy1(p_out_queue_[i]);
    }
    delete[] p_out_queue_;
}

bam_block* BamComplete::getBuffer(int idx) {
    return buffer_pool_[idx];
}

bam1_t* BamComplete::getResultBuf(int idx) {
    return result_pool_[idx];
}

void BamComplete::pushBlockResults(const bam1_t* records, int n) {
    if (!records || n <= 0) return;

    for (int i = 0; i < n; ++i) {
        // 等待队列有空间
        while (p_queueNumNow >= p_queueSizeLim) {
            usleep(10);
        }

        bam1_t* dst = p_out_queue_[p_queueP2];
        bam_copy1(dst, &records[i]);
        p_queueP2 = (p_queueP2 + 1) % p_queueSizeLim;
        p_queueNumNow++;
        p_queueCanRead++;
    }
}

bam1_t* BamComplete::popRecord() {
    while (p_queueCanRead <= 0 ) {
        if(complete_flag) return NULL;  // 已经不可能读到新的了
        usleep(10);
    }
    bam1_t* out;
    out = p_out_queue_[p_queueC2];
    p_queueC2 = (p_queueC2 + 1) % p_queueSizeLim;
    p_queueCanRead--;
    //p_queueNumNow--;
    return out;
}

//这里需要按顺序返回
bool BamComplete::backRecord(bam1_t* record) {
    if (record == nullptr) return false;
    p_queueNumNow--;
    return true;
}

//仅仅是处理完成
void BamComplete::markComplete() {
    complete_flag = true;
}

//是否全部完成
bool BamComplete::isComplete() const {
    return complete_flag;
}
