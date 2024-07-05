//
// Created by 赵展 on 2022/1/21.
//

#include "BamRead.h"

BamRead::BamRead(){

}

void BamRead::resize(int BufferSize){

    readBlockSize = BufferSize+1;
    readBlock = new bam_block*[readBlockSize];
    read_bg = 0;read_ed = BufferSize-1;
    for (int i=read_bg;i<=read_ed;i++) readBlock[i] = new bam_block;

    consumerBlockSize=2*BufferSize+5;
    consumerBlock = new bam_block*[consumerBlockSize];
    consumer_bg = 1;consumer_ed = 0;
    blockNum = 0;
    blockTot = 0;

    read_complete = false;
}




BamRead::BamRead(int BufferSize){
    readBlockSize = BufferSize+1;
    readBlock = new bam_block*[readBlockSize];
    read_bg = 0;read_ed = BufferSize-1;
    for (int i=read_bg;i<=read_ed;i++) readBlock[i] = new bam_block;

    consumerBlockSize=2*BufferSize+5;
    consumerBlock = new bam_block*[consumerBlockSize];
    consumer_bg = 1;consumer_ed = 0;
    blockNum = 0;
    blockTot = 0;

    read_complete = false;
}

bam_block* BamRead::getEmpty() {
    while ((read_ed+1)%readBlockSize == read_bg){
        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
    }
    int num = read_bg;
    read_bg=(read_bg+1)%readBlockSize;
    return readBlock[num];
}

void BamRead::inputBlock(bam_block* block){
//    consumerBlock[(consumer_ed+1)%consumerBlockSize] = block;
//    consumer_ed = (consumer_ed+1)%consumerBlockSize;
    consumerBlock[blockTot%consumerBlockSize] = block;
    blockTot+=1;
}

std::pair<bam_block*, int> BamRead::getReadBlock(){
//    mtx_consumer.lock();
//    while ((consumer_ed+1)%consumerBlockSize == consumer_bg){
//        mtx_consumer.unlock();
//        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
//        if (read_complete && (consumer_ed+1)%consumerBlockSize == consumer_bg ) return std::pair<bam_block *,int>(NULL,-1);
//        mtx_consumer.lock();
//    }
//    int num_point = consumer_bg;
//    consumer_bg = (consumer_bg+1)%consumerBlockSize;
//    bam_block* res = consumerBlock[num_point];
//    int num_block = blockNum++;
//    mtx_consumer.unlock();
////    return std::pair<bam_block*, int>(consumerBlock[num_point],num_block);
//    return std::pair<bam_block*, int>(res,num_block);
//    while (1){
//        while (blockNum.load(std::memory_order_acq_rel) == blockTot){
//            std::this_thread::sleep_for(std::chrono::nanoseconds(1000));
//            if (read_complete && blockNum.load(std::memory_order_acq_rel) == blockTot )return std::pair<bam_block *,int>(NULL,-1);
//        }
//        int num = blockNum.load(std::memory_order_acq_rel);
//        if (num<blockTot && blockNum.compare_exchange_weak(num,num+1,std::memory_order_acq_rel)){
//            bam_block* res = consumerBlock[num%consumerBlockSize];
//            return std::pair<bam_block*, int>(res,num);
//        }
//    }
    while (1){
        int num = blockTot;
        while (blockNum.compare_exchange_strong(num,num,std::memory_order_relaxed)){
            std::this_thread::sleep_for(std::chrono::nanoseconds(1));
            if (read_complete && blockNum.load(std::memory_order_relaxed) == blockTot )return std::pair<bam_block *,int>(NULL,-1);
            num = blockTot;
        }
//        int num = blockNum.load(std::memory_order_acq_rel);
        if (num<blockTot && blockNum.compare_exchange_strong(num,num+1,std::memory_order_relaxed)){
            bam_block* res = consumerBlock[num%consumerBlockSize];
            return std::pair<bam_block*, int>(res,num);
        }
    }

}

void BamRead::backBlock(bam_block* block){
    mtx_read.lock();
    readBlock[(read_ed+1)%readBlockSize] = block;
    read_ed = (read_ed+1)%readBlockSize;
    mtx_read.unlock();
}

void BamRead::ReadComplete(){
    read_complete=true;
}
