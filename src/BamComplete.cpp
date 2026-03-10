#include "BamComplete.h"

BamComplete::BamComplete(int queue_size)
{

    // 初始化缓冲区-----------------------------
    buffer_pool_.resize(64);
    result_pool_.resize(64);

    for (int i = 0; i < 64; ++i) {
        buffer_pool_[i] = new bam_block;
         
        // Allocate data buffer with 64-byte alignment for Sunway slave cores
        buffer_pool_[i]->data = aligned_alloc_custom(64, BGZF_MAX_BLOCK_SIZE);
        if (!buffer_pool_[i]->data) {
            printf("Allocation failed - handle error appropriately\n");
            // Allocation failed - handle error appropriately
            delete buffer_pool_[i];
            buffer_pool_[i] = nullptr;
        }


        result_pool_[i].resize(MAX_RECORDS_PER_BLOCK);
        for (int j = 0; j < MAX_RECORDS_PER_BLOCK; ++j) {
            // result_pool_[i][j] = bam_init1();  // 指针赋值没问题
            
            //分配足够大小的空间
            bam1_t *b = bam_init1();
            // 分配 1KB 对齐内存
            b->data = (uint8_t*)aligned_alloc_custom(64, INIT_DATA_SIZE);
            b->m_data = INIT_DATA_SIZE;
            b->l_data = 0;
            b->mempolicy = BAM_USER_OWNS_DATA;
            result_pool_[i][j] = b;

        }
    }

    // 初始化内存池------------------------------------------
    pro_queueSizeLim = queue_size + 1; // 多申请一个，防止边界问题
    producer_queue_ = new bam1_t *[pro_queueSizeLim];
    pro_bg = 0;
    pro_ed = queue_size - 1;

    for (int i = pro_bg; i <= pro_ed; i++) {
        //为data字段分配足够大小的空间
        bam1_t *b = bam_init1();
        // 分配 1KB 对齐内存
        b->data = (uint8_t*)aligned_alloc_custom(64, INIT_DATA_SIZE);
        b->m_data = INIT_DATA_SIZE;
        b->l_data = 0;
        b->mempolicy = BAM_USER_OWNS_DATA;
        producer_queue_[i] = b;
    }

    con_queueSizeLim = queue_size + 5; // 多留一些空间
    con_bg = 1;
    con_ed = 0;
    consumer_queue_ = new bam1_t* [con_queueSizeLim];

    complete_flag = false;
}

BamComplete::~BamComplete() {
    // 释放 buffer_pool_ 中的 bam_block 及其 data
    for (auto buf : buffer_pool_) {
        if (buf) {
            if (buf->data) {
                aligned_free_custom(buf->data);
                buf->data = nullptr;
            }
            delete buf;
        }
    }
    buffer_pool_.clear();

    // 释放 result_pool_ 中的 bam1_t*
    for (auto& row : result_pool_) {
        for (auto ptr : row) {
            if (ptr) {
                bam_destroy1_sw(ptr);
            }
        }
        row.clear();
    }
    result_pool_.clear();

    //释放消费者队列
    if (consumer_queue_) {
        delete[] consumer_queue_;
        consumer_queue_ = nullptr;
    }

    //释放生产者队列
    if (producer_queue_) {
        for (int i = 0; i < pro_queueSizeLim - 1; ++i) {
            if (producer_queue_[i]) {
                bam_destroy1_sw(producer_queue_[i]);
                producer_queue_[i] = nullptr;
            }
        }

        delete[] producer_queue_;
        producer_queue_ = nullptr;
    }
}


bam_block* BamComplete::getBuffer(int idx) {
    return buffer_pool_[idx];
}

std::vector<bam1_t*>& BamComplete::getResultBuf(int idx) {
    return result_pool_[idx];
}


bam1_t* BamComplete::getEmpty() {
    //等待直到有空闲 bam1_t
    while ((pro_ed + 1) % pro_queueSizeLim == pro_bg) {
        usleep(10); 
        //  std::this_thread::sleep_for(std::chrono::nanoseconds(1));
    }
    //先取再+1
    int num = pro_bg;
    pro_bg = (pro_bg + 1) % pro_queueSizeLim;
    return producer_queue_[num];
}

void BamComplete::backBam1_t(bam1_t *bam1) {
    //先+1再放回
    producer_queue_[(pro_ed + 1) % pro_queueSizeLim] = bam1;
    pro_ed = (pro_ed + 1) % pro_queueSizeLim;
}

void BamComplete::inputBam1_t(bam1_t*& bam1) {
    //先+1再放入
    consumer_queue_[(con_ed + 1) % con_queueSizeLim] = bam1;
    con_ed = (con_ed + 1) % con_queueSizeLim;
}

bam1_t* BamComplete::getBam1_t() {
    //等待直到有可用组
    while ((con_ed + 1) % con_queueSizeLim == con_bg) {
        usleep(10); 
        //  std::this_thread::sleep_for(std::chrono::nanoseconds(1));
        if (complete_flag && (con_ed + 1) % con_queueSizeLim == con_bg) return nullptr;
    }
    //先取再+1
    int num = con_bg;
    bam1_t* bam1 = consumer_queue_[con_bg];
    con_bg = (con_bg + 1) % con_queueSizeLim;

    return bam1;
}

void BamComplete::markComplete() {
    complete_flag = true;
}

bool BamComplete::isComplete() const {
    return complete_flag;
}
