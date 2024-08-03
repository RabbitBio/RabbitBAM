#include "BamRead.h"

BamRead::BamRead() {

}

void BamRead::resize(int BufferSize) {

    readBlockSize = BufferSize + 1;
    readBlock = new bam_block *[readBlockSize];
    read_bg = 0;
    read_ed = BufferSize - 1;
    for (int i = read_bg; i <= read_ed; i++) readBlock[i] = new bam_block;

    consumerBlockSize = 2 * BufferSize + 5;
    consumerBlock = new bam_block *[consumerBlockSize];
    consumer_bg = 1;
    consumer_ed = 0;
    blockNum = 0;
    blockTot = 0;

    read_complete = false;
}


BamRead::BamRead(int BufferSize) {
    readBlockSize = BufferSize + 1;
    readBlock = new bam_block *[readBlockSize];
    read_bg = 0;
    read_ed = BufferSize - 1;
    for (int i = read_bg; i <= read_ed; i++) readBlock[i] = new bam_block;

    consumerBlockSize = 2 * BufferSize + 5;
    consumerBlock = new bam_block *[consumerBlockSize];
    consumer_bg = 1;
    consumer_ed = 0;
    blockNum = 0;
    blockTot = 0;

    read_complete = false;
}

bam_block *BamRead::getEmpty() {
    while ((read_ed + 1) % readBlockSize == read_bg) {
        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
    }
    int num = read_bg;
    read_bg = (read_bg + 1) % readBlockSize;
    return readBlock[num];
}

void BamRead::inputBlock(bam_block *block) {
    consumerBlock[blockTot % consumerBlockSize] = block;
    blockTot += 1;
}

std::pair<bam_block *, int> BamRead::getReadBlock() {
    while (1) {
        int num = blockTot;
        while (blockNum.compare_exchange_strong(num, num, std::memory_order_relaxed)) {
            std::this_thread::sleep_for(std::chrono::nanoseconds(1));
            if (read_complete && blockNum.load(std::memory_order_relaxed) == blockTot)
                return std::pair<bam_block *, int>(NULL, -1);
            num = blockTot;
        }
        if (num < blockTot && blockNum.compare_exchange_strong(num, num + 1, std::memory_order_relaxed)) {
            bam_block *res = consumerBlock[num % consumerBlockSize];
            return std::pair<bam_block *, int>(res, num);
        }
    }

}

void BamRead::backBlock(bam_block *block) {
    mtx_read.lock();
    readBlock[(read_ed + 1) % readBlockSize] = block;
    read_ed = (read_ed + 1) % readBlockSize;
    mtx_read.unlock();
}

void BamRead::ReadComplete() {
    read_complete = true;
}
