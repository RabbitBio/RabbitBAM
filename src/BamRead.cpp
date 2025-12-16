#include "BamRead.h"

BamRead::BamRead(int queue_size) {

    int total_blocks = queue_size * 64;
    readBlockSize = total_blocks + 1; // 多申请一个，防止边界问题
    readBlock = new bam_block *[readBlockSize];
    read_bg = 0;
    read_ed = total_blocks - 1;


    //for (int i = read_bg; i <= read_ed; i++) readBlock[i] = new bam_block;
    for (int i = read_bg; i <= read_ed; i++) {
        readBlock[i] = new bam_block;
        // Allocate data buffer with 64-byte alignment for Sunway slave cores
        readBlock[i]->data = aligned_alloc_custom(64, BGZF_MAX_BLOCK_SIZE);
        if (!readBlock[i]->data) {
            printf("Allocation failed - handle error appropriately\n");
            // Allocation failed - handle error appropriately
            delete readBlock[i];
            readBlock[i] = nullptr;
        }
    }

    con_queueSizeLim = queue_size + 5; // 多留一些空间
    con_bg = 1;
    con_ed = 0;
    consumer_queue_ = new std::vector<bam_block*>[con_queueSizeLim];

    read_complete = false;
}

BamRead:: ~BamRead() {
        // 释放 readBlock 内存
        if (readBlock) {
            for (int i = 0; i < readBlockSize; i++) {
                if (readBlock[i]) {
                    // 释放对齐的 data
                    if (readBlock[i]->data) {
                        aligned_free_custom(readBlock[i]->data);
                        readBlock[i]->data = nullptr;
                    }
                    delete readBlock[i];
                    readBlock[i] = nullptr;
                }
            }
            delete[] readBlock;
            readBlock = nullptr;
        }

        // 释放 consumer_queue_
        if (consumer_queue_) {
            // 不需要手动清理 vector 内指针，假设外部会处理 bam_block 指针生命周期
            delete[] consumer_queue_;
            consumer_queue_ = nullptr;
        }
}

bam_block* BamRead::getEmpty() {
    //等待直到有空闲 block
    while ((read_ed + 1) % readBlockSize == read_bg) {
        usleep(10); 
    }
    //先取再+1
    int num = read_bg;
    read_bg = (read_bg + 1) % readBlockSize;
    return readBlock[num];
}

void BamRead::backBlock(bam_block *block) {
    //先+1再放回
    readBlock[(read_ed + 1) % readBlockSize] = block;
    read_ed = (read_ed + 1) % readBlockSize;
}

void BamRead::inputBlock64(std::vector<bam_block*>& group) {
    //先+1再放入
    consumer_queue_[(con_ed + 1) % con_queueSizeLim] = group;
    con_ed = (con_ed + 1) % con_queueSizeLim;
}

std::vector<bam_block*> BamRead::getBlock64() {
    //等待直到有可用组
    while ((con_ed + 1) % con_queueSizeLim == con_bg) {
        usleep(10); 
        if (read_complete && (con_ed + 1) % con_queueSizeLim == con_bg) return std::vector<bam_block*>();
    }
    //先取再+1
    int num = con_bg;
    std::vector<bam_block*> group = consumer_queue_[con_bg];
    con_bg = (con_bg + 1) % con_queueSizeLim;

    return group;
}

bool BamRead::isComplete() const {
    return read_complete;
}

void BamRead::markComplete() {
    read_complete = true;
}
