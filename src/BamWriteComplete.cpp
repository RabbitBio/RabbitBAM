#include "BamWriteComplete.h"

BamWriteComplete::BamWriteComplete(int queue_size) {

    // 初始化解压缓冲区
    buffer_pool_.resize(64);
    for (int i = 0; i < 64; ++i) {
        buffer_pool_[i] = new bam_block;
        // Allocate data buffer with 64-byte alignment for Sunway slave cores
        buffer_pool_[i]->data = aligned_alloc_custom(64, BGZF_MAX_BLOCK_SIZE);
    }

    // 初始化空闲块池
    pool_size = queue_size + 1;
    blockPool = new bam_block*[pool_size];
    pool_bg = 0;
    pool_ed = queue_size - 1; 
    for (int i = pool_bg; i <= pool_ed; i++) {
        blockPool[i] = new bam_block;
        // Allocate data buffer with 64-byte alignment for Sunway slave cores
        buffer_pool_[i]->data = aligned_alloc_custom(64, BGZF_MAX_BLOCK_SIZE);
    }

    // 初始化压缩块队列
    q_size = queue_size + 5;
    queue = new bam_block*[q_size];
    q_bg = 1;
    q_ed = 0;

    complete_flag = false;
}

BamWriteComplete::~BamWriteComplete() {
    for (int i = 0; i < pool_size-1; i++)
        delete blockPool[i];
    delete[] blockPool;
    delete[] queue;
}

bam_block* BamWriteComplete::getBuffer(int idx) {
    return buffer_pool_[idx];
}

bam_block* BamWriteComplete::getEmpty() {
    while (pool_bg == (pool_ed + 1) % pool_size) {
        usleep(10);
    }
    bam_block* blk = blockPool[pool_bg];
    pool_bg = (pool_bg + 1) % pool_size;
    return blk;
}

void BamWriteComplete::backBlock(bam_block *block) {
    blockPool[(pool_ed + 1) % pool_size] = block;
    pool_ed = (pool_ed + 1) % pool_size;
}

void BamWriteComplete::pushCompressedBlock(bam_block* block) {
    queue[ (q_ed + 1) % q_size ] = block;
    q_ed = (q_ed + 1) % q_size;
}

bam_block* BamWriteComplete::getCompressedBlock() {
    while ((q_ed + 1) % q_size == q_bg) {
        usleep(10);
        if (complete_flag && (q_ed + 1) % q_size == q_bg) return {};
    }

    int num = q_bg;
    auto blk = queue[q_bg];
    q_bg = (q_bg + 1) % q_size;
    return blk;
}

void BamWriteComplete::markComplete() {
    complete_flag = true;
}

bool BamWriteComplete::isComplete() const {
    return complete_flag;
}
