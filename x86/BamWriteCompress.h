#ifndef BAMSTATUS_BAMWRITECOMPRESS_H
#define BAMSTATUS_BAMWRITECOMPRESS_H

#include "BamTools.h"
#include <thread>
#include <atomic>

class BamWriteCompress {
public:
    BamWriteCompress(int BufferSize, int threadNumber);


    bam_write_block *getEmpty();

    void inputUnCompressData(bam_write_block *data);

    bam_write_block *getUnCompressData();

    void inputCompressData(bam_write_block *data);


    bam_write_block *getCompressData();

    void backEmpty(bam_write_block *data);


    void CompressThreadComplete() {
        mtx_compressThread.lock();
        compressThread--;
        mtx_compressThread.unlock();
    }

    void WriteComplete() {
        isWriteComplete = true;
    }

public:

    int wait_num;
private:

    std::atomic<int> blockNum;
    bam_write_block **compress_data;
    int compress_bg;
    int compress_ed;
    int compress_size;
    std::mutex mtx_compress;
    std::mutex mtx_input;
    int blockInputNum;
    std::atomic<int> blockInputPos;
    bam_write_block **need_compress_data;
    int need_compress_bg;
    int need_compress_ed;
    int need_compress_size;
    std::mutex mtx_need_compress;
    bam_write_block **consumer_data;
    int consumer_bg;
    int consumer_ed;
    int consumer_size;
    bool *is_ok;
    int compressThread;
    std::mutex mtx_compressThread;
    bool isWriteComplete = false;

};


#endif //BAMSTATUS_BAMWRITECOMPRESS_H
