#ifndef BAMSTATUS_BAMCOMPRESS_H
#define BAMSTATUS_BAMCOMPRESS_H

#include "BamTools.h"
#include <thread>
#include <atomic>


class BamCompress {
public:


    BamCompress();

    void resize(int BufferSize, int threadNumber);

    BamCompress(int BufferSize, int threadNumber);

    bam_block *getEmpty();

    void inputUnCompressData(bam_block *data, int block_num);

    bam_block *getUnCompressData();

    void backEmpty(bam_block *data);

    void CompressThreadComplete() {
        mtx_compressThread.lock();
        compressThread--;
        mtx_compressThread.unlock();
    }

public:

    int wait_num;
private:

    std::atomic<int> blockNum;
    bam_block **compress_data;
    int compress_bg;
    int compress_ed;
    int compress_size;
    std::mutex mtx_compress;
    std::mutex mtx_input;
    bam_block **consumer_data;
    int consumer_bg;
    int consumer_ed;
    int consumer_size;
    bool *is_ok;
    int compressThread;
    std::mutex mtx_compressThread;

};


#endif //BAMSTATUS_BAMCOMPRESS_H
