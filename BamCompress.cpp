#include "BamCompress.h"


BamCompress::BamCompress() {

};


void BamCompress::resize(int BufferSize, int threadNumber) {

    blockNum = 0;

    compress_bg = 0;
    compress_ed = BufferSize - 1;
    compress_size = BufferSize + 1;
    compress_data = new bam_block *[compress_size];
    for (int i = compress_bg; i <= compress_ed; i++) compress_data[i] = new bam_block;


    consumer_bg = 1;
    consumer_ed = 0;
    consumer_size = BufferSize + 5;
    consumer_data = new bam_block *[consumer_size];
    is_ok = new bool[consumer_size];
    for (int i = 0; i < consumer_size; i++) is_ok[i] = false;


    compressThread = threadNumber;

    wait_num = 0;

}


BamCompress::BamCompress(int BufferSize, int threadNumber) {
    blockNum = 0;

    compress_bg = 0;
    compress_ed = BufferSize - 1;
    compress_size = BufferSize + 1;
    compress_data = new bam_block *[compress_size];
    for (int i = compress_bg; i <= compress_ed; i++) compress_data[i] = new bam_block;


    consumer_bg = 1;
    consumer_ed = 0;
    consumer_size = BufferSize + 5;
    consumer_data = new bam_block *[consumer_size];


    compressThread = threadNumber;

    wait_num = 0;
}


bam_block *BamCompress::getEmpty() {
    mtx_compress.lock();
    while ((compress_ed + 1) % compress_size == compress_bg) {
        mtx_compress.unlock();
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        mtx_compress.lock();
    }
    int num = compress_bg;
    compress_bg = (compress_bg + 1) % compress_size;
    bam_block *res = compress_data[num];
    mtx_compress.unlock();
    return res;
}


void BamCompress::inputUnCompressData(bam_block *data, int block_num) {

    while (block_num != blockNum.load(std::memory_order_acq_rel)) {
        wait_num += 1;
        std::this_thread::sleep_for(std::chrono::nanoseconds(block_num - blockNum) / 8);
    }

    consumer_data[(consumer_ed + 1) % consumer_size] = data;
    consumer_ed = (consumer_ed + 1) % consumer_size;
    blockNum.store(blockNum.load(std::memory_order_acq_rel) + 1, std::memory_order_acq_rel);
}

bam_block *BamCompress::getUnCompressData() {
    while ((consumer_ed + 1) % consumer_size == consumer_bg) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        if (compressThread == 0 && (consumer_ed + 1) % consumer_size == consumer_bg) return nullptr;
    }
    int num = consumer_bg;
    bam_block *res = consumer_data[consumer_bg];
    consumer_bg = (consumer_bg + 1) % consumer_size;

    return res;
}

void BamCompress::backEmpty(bam_block *data) {
    compress_data[(compress_ed + 1) % compress_size] = data;
    compress_ed = (compress_ed + 1) % compress_size;
}


