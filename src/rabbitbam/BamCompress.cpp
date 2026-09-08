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
    int spins = 0;
    mtx_compress.lock();
    while ((compress_ed + 1) % compress_size == compress_bg) {
        mtx_compress.unlock();
        // Was sleep_for(5ms): a fixed 5 ms stall whenever the empty-block ring
        // is momentarily drained throttled the decompression threads. rb_backoff
        // yields first (refill picked up within microseconds, same as the old
        // bare yield), then parks only if persistently starved so a starved
        // decode thread stops thrashing the scheduler under oversubscription.
        // Pure scheduling change — the same blocks flow through in the same order.
        rb_backoff(spins);
        mtx_compress.lock();
    }
    int num = compress_bg;
    compress_bg = (compress_bg + 1) % compress_size;
    bam_block *res = compress_data[num];
    mtx_compress.unlock();
    return res;
}


void BamCompress::inputUnCompressData(bam_block *data, int block_num) {

    int spins = 0;
    while (block_num != blockNum.load(std::memory_order_acq_rel)) {
        wait_num += 1;
        rb_backoff(spins);
    }

    consumer_data[(consumer_ed + 1) % consumer_size] = data;
    consumer_ed = (consumer_ed + 1) % consumer_size;
    blockNum.store(blockNum.load(std::memory_order_acq_rel) + 1, std::memory_order_acq_rel);
}

bam_block *BamCompress::getUnCompressData() {
    int spins = 0;
    while ((consumer_ed + 1) % consumer_size == consumer_bg) {
        // Was sleep_for(1ms): the single assign thread polled decompressed
        // blocks only 1000×/s, capping whole-file throughput. rb_backoff picks
        // up a freshly produced block within microseconds, then parks if the
        // producers fall behind. Same data, same order.
        rb_backoff(spins);
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


