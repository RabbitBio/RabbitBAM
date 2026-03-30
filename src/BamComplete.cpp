#include "BamComplete.h"

BamComplete::BamComplete(int queue_size)
{
    //批量分配 buffer_pool_（64块解压缓冲区）
    bp_flat_ = (bam_block*)calloc(64, sizeof(bam_block));
    posix_memalign((void**)&bp_data_, 64, (size_t)64 * BGZF_MAX_BLOCK_SIZE);

    buffer_pool_.resize(64);
    result_pool_.resize(64);
    for (int i = 0; i < 64; ++i) {
        bp_flat_[i].data = bp_data_ + (size_t)i * BGZF_MAX_BLOCK_SIZE;
        buffer_pool_[i]  = &bp_flat_[i];
    }

    //批量分配 result_pool_（64×MAX_RECORDS_PER_BLOCK 个 bam1_t）
    int RP_TOTAL = 64 * (int)MAX_RECORDS_PER_BLOCK;
    rp_flat_ = (bam1_t*)calloc(RP_TOTAL, sizeof(bam1_t));
    posix_memalign((void**)&rp_data_, 64, (size_t)RP_TOTAL * INIT_DATA_SIZE);

    for (int i = 0; i < 64; ++i) {
        result_pool_[i].resize(MAX_RECORDS_PER_BLOCK);
        for (int j = 0; j < (int)MAX_RECORDS_PER_BLOCK; ++j) {
            int fi = i * (int)MAX_RECORDS_PER_BLOCK + j;
            rp_flat_[fi].data      = rp_data_ + (size_t)fi * INIT_DATA_SIZE;
            rp_flat_[fi].m_data    = INIT_DATA_SIZE;
            rp_flat_[fi].l_data    = 0;
            rp_flat_[fi].mempolicy = BAM_USER_OWNS_DATA;
            result_pool_[i][j]     = &rp_flat_[fi];
        }
    }

    //批量分配 producer_queue_（bam1_t 内存池）
    pro_queueSizeLim = queue_size + 1;
    producer_queue_  = new bam1_t*[pro_queueSizeLim];
    pro_bg = 0;
    pro_ed = queue_size - 1;

    pq_flat_ = (bam1_t*)calloc(queue_size, sizeof(bam1_t));
    posix_memalign((void**)&pq_data_, 64, (size_t)queue_size * INIT_DATA_SIZE);

    for (int i = 0; i < queue_size; ++i) {
        pq_flat_[i].data      = pq_data_ + (size_t)i * INIT_DATA_SIZE;
        pq_flat_[i].m_data    = INIT_DATA_SIZE;
        pq_flat_[i].l_data    = 0;
        pq_flat_[i].mempolicy = BAM_USER_OWNS_DATA;
        producer_queue_[i]    = &pq_flat_[i];
    }

    //消费者队列
    con_queueSizeLim = queue_size + 5;
    con_bg = 1;
    con_ed = 0;
    consumer_queue_ = new bam1_t*[con_queueSizeLim];

    complete_flag = false;
}

BamComplete::~BamComplete() {
    // buffer_pool_: bp_flat_ 整块释放
    free(bp_flat_);
    free(bp_data_);
    buffer_pool_.clear();

    // result_pool_: rp_flat_ 整块释放
    free(rp_flat_);
    free(rp_data_);
    result_pool_.clear();

    // producer_queue_: pq_flat_ 整块释放
    free(pq_flat_);
    free(pq_data_);
    if (producer_queue_) {
        delete[] producer_queue_;
        producer_queue_ = nullptr;
    }

    if (consumer_queue_) {
        delete[] consumer_queue_;
        consumer_queue_ = nullptr;
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
    }
    //先取再+1
    int num = pro_bg;
    pro_bg = (pro_bg + 1) % pro_queueSizeLim;
    return producer_queue_[num];
}

void BamComplete::getEmptyBatch(bam1_t** dst, int n) {
    for (int k = 0; k < n; k++) {
        while ((pro_ed + 1) % pro_queueSizeLim == pro_bg) {
            usleep(10);
        }
        dst[k] = producer_queue_[pro_bg];
        pro_bg  = (pro_bg + 1) % pro_queueSizeLim;
    }
}

void BamComplete::backBam1_t(bam1_t* bam1) {
    producer_queue_[(pro_ed + 1) % pro_queueSizeLim] = bam1;
    pro_ed = (pro_ed + 1) % pro_queueSizeLim;
}

void BamComplete::backBam1_tBatch(bam1_t** src, int n) {
    for (int k = 0; k < n; k++) {
        producer_queue_[(pro_ed + 1) % pro_queueSizeLim] = src[k];
        pro_ed = (pro_ed + 1) % pro_queueSizeLim;
    }
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
