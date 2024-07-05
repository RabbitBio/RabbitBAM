//
// Created by 赵展 on 2022/7/7.
//

#include "BamWriteBlock.h"
BamWriteBlockConfig::BamBlockConfig(){};
BamWriteBlockConfig::BamBlockConfig(int Buffer_number){
    this->Buffer_number=Buffer_number;
    this->write_number=Buffer_number+10;
    this->complete=0;
}

BamWriteBlock::BamBlock(){};
BamWriteBlock::BamBlock(BamBlockConfig *config){
    this->config= config;
    this->buffer=new bam_block*[this->config->Buffer_number];
    for (int i=0;i<this->config->Buffer_number;++i) this->buffer[i]=new bam_block;
    this->compress=new int[this->config->write_number];
    this->compress_bg=0;
    this->compress_ed=0;
    this->read=new int[this->config->write_number];
    this->read_bg=0;
    this->read_ed=this->config->Buffer_number;
    for (int i=this->read_bg;i<this->read_ed;i++) this->read[i]=i;
}
pair<bam_block *,int> BamWriteBlock::getEmpty(){
    while (read_bg==read_ed){
        this_thread::sleep_for(chrono::milliseconds(5));
        //this_thread::yield();
    }
    int num=read_bg;
    read_bg=(read_bg+1)%config->write_number;
    return pair<bam_block *,int>(this->buffer[read[num]],read[num]);

}
void BamWriteBlock::inputblock(int id){
    compress[compress_ed]=id;
    compress_ed=(compress_ed+1)%config->write_number;
}
pair<bam_block *,int> BamWriteBlock::getCompressdata(){
    mtx_compress.lock();
    while (compress_ed==compress_bg){
        mtx_compress.unlock();
        this_thread::sleep_for(chrono::milliseconds(3));
        //this_thread::yield();
        if (config->complete) return pair<bam_block *,int>(NULL,-1);
        mtx_compress.lock();
    }
    int num=compress_bg;
    compress_bg=(compress_bg+1)%config->write_number;
    mtx_compress.unlock();
    return pair<bam_block *,int>(this->buffer[compress[num]],compress[num]);
}
void BamWriteBlock::backempty(int id){
    mtx_read.lock();
    read[read_ed]=id;
    read_ed=(read_ed+1)%config->write_number;
    mtx_read.unlock();
}
bool BamWriteBlock::isComplete(){
    return config->complete;
};
void BamWriteBlock::ReadComplete() {
    config->complete=1;
}