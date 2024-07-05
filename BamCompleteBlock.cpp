//
// Created by 赵展 on 2022/1/22.
//

#include "BamCompleteBlock.h"


BamCompleteBlock::BamCompleteBlock() {
}


void BamCompleteBlock::resize(int BufferSize){

    BlockComplete=false;

    complete_bg=0;
    complete_ed=BufferSize-1;
    complete_size=BufferSize+1;
    complete_data = new bam_complete_block*[complete_size];
    for (int i=complete_bg;i<=complete_ed;i++) {
        complete_data[i] = new bam_complete_block;
        complete_data[i]->data_size= BGZF_MAX_BLOCK_COMPLETE_SIZE;
        //        printf("dadakjbkjansdkjandasdasd\n");
        complete_data[i]->data = new unsigned char[complete_data[i]->data_size]; //????
        complete_data[i]->pos=0;
        complete_data[i]->length=0;
    }



    consumer_bg = 1;
    consumer_ed = 0;
    consumer_size = BufferSize+5;

    consumer_data = new bam_complete_block*[consumer_size];
}

BamCompleteBlock::BamCompleteBlock(int BufferSize){

    BlockComplete=false;

    complete_bg=0;
    complete_ed=BufferSize-1;
    complete_size=BufferSize+1;
    complete_data = new bam_complete_block*[complete_size];
    for (int i=complete_bg;i<=complete_ed;i++) {
        complete_data[i] = new bam_complete_block;
        complete_data[i]->data_size= BGZF_MAX_BLOCK_COMPLETE_SIZE;
        //        printf("dadakjbkjansdkjandasdasd\n");
        complete_data[i]->data = new unsigned char[complete_data[i]->data_size]; //????
        complete_data[i]->pos=0;
        complete_data[i]->length=0;
    }



    consumer_bg = 1;
    consumer_ed = 0;
    consumer_size = BufferSize+5;
    consumer_data = new bam_complete_block*[consumer_size];

}




bam_complete_block* BamCompleteBlock::getEmpty(){
    while ((complete_ed+1)%complete_size == complete_bg){
        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
    }
    int num = complete_bg;
    complete_bg=(complete_bg+1)%complete_size;
    return complete_data[num];
}
/*
 * 输入一个重新划分之后的内存块，到消费队列中
 * 输入： bam_complete_block* : 重新划分的内存块
 * 输出： 无
 *
 */
void BamCompleteBlock::inputCompleteBlock(bam_complete_block* block){
    consumer_data[(consumer_ed+1)%consumer_size] = block;
    consumer_ed = (consumer_ed+1)%consumer_size;
}
/*
 * 获取一个重新划分好的内存块
 * 输入：无
 * 输出：bam_complete_block* : 重新划分好的内存块
 */
bam_complete_block* BamCompleteBlock::getCompleteBlock(){
    mtx_consumer.lock();
    while ((consumer_ed+1)%consumer_size == consumer_bg){
        mtx_consumer.unlock();
        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
        if (BlockComplete && (consumer_ed+1)%consumer_size == consumer_bg ) return nullptr;
        mtx_consumer.lock();
    }
    int num_point = consumer_bg;
    consumer_bg = (consumer_bg+1)%consumer_size;
    mtx_consumer.unlock();
    return consumer_data[num_point];
}

/*
 * 归还使用完毕的内存块
 * 输入：bam_complete_block* ：使用完毕的内存块
 * 输出：无
 */
void BamCompleteBlock::backEmpty(bam_complete_block* block){
    mtx_complete.lock();
    block->pos=0;
    block->length=0;
    complete_data[(complete_ed+1)%complete_size] = block;
    complete_ed = (complete_ed+1)%complete_size;
    mtx_complete.unlock();
}