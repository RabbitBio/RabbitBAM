//
// Created by 赵展 on 2021/3/10.
//
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <zlib.h>
#include <htslib/khash.h>
#include <htslib/thread_pool.h>
#include <stdint.h>
#include <chrono>
#include "config.h"
#include "BamBlock.h"
#include "Buffer.h"
#include "BamStatus.h"
#include "Duplicate.h"
#include "Overrepresent.h"
#include "CLI/CLI.hpp"
#include <sched.h>
#include <unistd.h>
#include "BamRead.h"
#include "BamCompress.h"
#include "BamCompleteBlock.h"
#include "BamTools.h"
#include "BamWriteCompress.h"
#include "BamWriter.h"
#include "BamReader.h"
#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8


#define DEBUG 0


typedef std::chrono::high_resolution_clock Clock;
#ifdef CYCLING
#define TDEF(x_) static unsigned long long int x_##_t0, x_##_t1;
    #define TSTART(x_) x_##_t0 = __rdtsc();
    #define TEND(x_) x_##_t1 = __rdtsc();
    #define TPRINT(x_, str) printf("%-20s \t%.6f\t M cycles\n", str, (double)(x_##_t1 - x_##_t0)/1e6);
#elif defined TIMING
#define TDEF(x_) chrono::high_resolution_clock::time_point x_##_t0, x_##_t1;
#define TSTART(x_) x_##_t0 = Clock::now();
#define TEND(x_) x_##_t1 = Clock::now();
#define TPRINT(x_, str) printf("%-20s \t%.6f\t sec\n", str, chrono::duration_cast<chrono::microseconds>(x_##_t1 - x_##_t0).count()/1e6);
#else
#define TDEF(x_)
#define TSTART(x_)
#define TEND(x_)
#define TPRINT(x_, str)
#endif
//const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
//int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
//uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
//uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};



#include <iostream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>



const size_t PACK_SIZE = 1024;

const size_t QUEUE_SIZE = 1024 * 1024 * 1024;

// 定义 bam_pack 结构体，用于存储打包的 bam1_t 对象
struct bam_pack {
    std::vector<bam1_t*> bams;

    bam_pack() {
        bams.reserve(PACK_SIZE);
    }

    ~bam_pack() {
        for (auto& bam : bams) {
            bam_destroy1(bam);
        }
    }
};

std::queue<bam_pack*> bam_queue;
std::mutex queue_mutex;
std::condition_variable not_full;
std::condition_variable not_empty;

void read_thread_task(BamReader* reader) {
    long num = 0;
    bam_pack* pack = new bam_pack();
    for (size_t i = 0; i < PACK_SIZE; ++i) {
        bam1_t* b = bam_init1();
        if (b == NULL) {
            std::cerr << "[E::" << __func__ << "] Out of memory allocating BAM struct.\n";
            return;
        }
        pack->bams.push_back(b);
    }
    while (true) {
        int cnt = 0;
        for (size_t i = 0; i < PACK_SIZE; ++i) {
            if (!reader->getBam1_t(pack->bams[i])) {
                break;
            }
            num++;
            cnt++;
        }

        if (cnt == 0) {
            //delete pack;
            break;
        }

        std::unique_lock<std::mutex> lock(queue_mutex);
        not_full.wait(lock, [] { return bam_queue.size() < QUEUE_SIZE; });

        bam_queue.push(pack);
        not_empty.notify_one();
    }
    printf("num is %lld\n", num);

    // 通知写线程不再有新数据
    std::unique_lock<std::mutex> lock(queue_mutex);
    bam_queue.push(nullptr);
    not_empty.notify_one();
}

void read_thread_task2(BamReader* reader) {
    long num = 0;
	while (true) {
        bam_pack* pack = new bam_pack();
        for (size_t i = 0; i < PACK_SIZE; ++i) {
            bam1_t* b = bam_init1();
            if (b == NULL) {
                std::cerr << "[E::" << __func__ << "] Out of memory allocating BAM struct.\n";
                delete pack;
                return;
            }

            if (!reader->getBam1_t(b)) {
                bam_destroy1(b);
                delete pack;
                break;
            }
            pack->bams.push_back(b);
        }

        if (pack->bams.empty()) {
            delete pack;
            break;
        }

        std::unique_lock<std::mutex> lock(queue_mutex);
        not_full.wait(lock, [] { return bam_queue.size() < PACK_SIZE; });

        bam_queue.push(pack);
        not_empty.notify_one();
    }

    // 通知写线程不再有新数据
    std::unique_lock<std::mutex> lock(queue_mutex);
    bam_queue.push(nullptr);
    not_empty.notify_one();
	printf("num is %d\n", num);
}

void write_thread_task(BamWriter* writer) {
    printf("queue size %lld\n", bam_queue.size());
    double t0 = GetTime();
    while (true) {
        std::unique_lock<std::mutex> lock(queue_mutex);
        not_empty.wait(lock, [] { return !bam_queue.empty(); });

        bam_pack* pack = bam_queue.front();
        bam_queue.pop();
        //printf("queue size %lld\n", bam_queue.size());

        not_full.notify_one();

        if (pack == nullptr) {
            break;
        }

        for (auto& b : pack->bams) {
            writer->write(b);
        }
        //delete pack;
    }
    printf("main write part cost %lf\n", GetTime() - t0);
    writer->over();
}

void write_thread_task2(BamWriter* writer) {
    double t0 = GetTime();
	while (true) {
        std::unique_lock<std::mutex> lock(queue_mutex);
        not_empty.wait(lock, [] { return !bam_queue.empty(); });

        bam_pack* pack = bam_queue.front();
        bam_queue.pop();

        not_full.notify_one();

        if (pack == nullptr) {
            break;
        }

        for (auto& b : pack->bams) {
            writer->write(b);
        }
        delete pack;
    }
    printf("main write part cost %lf\n", GetTime() - t0);
    writer->over();
}

int InputBlockNum=0;

long long NUM_N[100]={0};
long long NUM_M[100]={0};
long long NUM_TID[100][1000]={0};
void read_pack(BGZF *fp,BamRead *read){
    bam_block * b;
    b=read->getEmpty();
    int count=0;
    while(read_block(fp,b)==0){
        read->inputBlock(b);
//        printf("read block is %d\n",++count);
        b=read->getEmpty();
    }
    read->ReadComplete();
}
void write_pack(Buffer *buffer){
    while(!buffer->is_complete()){
        std::this_thread::sleep_for(chrono::milliseconds(10));
        buffer->output();
    }
}
void compress_pack(BamRead *read,BamCompress *compress){
    pair<bam_block *,int> comp;
    bam_block *un_comp = compress->getEmpty();
    while (1){
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
        comp=read->getReadBlock();
        //printf("%d is get One compressed data\n",id);
        if (comp.second<0) {
            //printf("%d is Over\n",id);
            break;
        }
        block_decode_func(comp.first,un_comp);
        read->backBlock(comp.first);

        std::pair<int,int> tmp_pair= find_divide_pos_and_get_read_number(un_comp);
        un_comp->split_pos=tmp_pair.first,un_comp->bam_number=tmp_pair.second;
        compress->inputUnCompressData(un_comp,comp.second);
//        while (!compress->tryinputUnCompressData(un_comp,comp.second)){
//            std::this_thread::sleep_for(std::chrono::milliseconds(1));
//        }
        un_comp = compress->getEmpty();
    }
    compress->CompressThreadComplete();
}
void assign_pack(BamCompress* compress,BamCompleteBlock* completeBlock){
    bam_block *un_comp = nullptr;
    bam_complete_block *assign_block = completeBlock->getEmpty();
    int need_block_len=0,input_length=0;
    int last_use_block_length=0;
    bool isclean = true;
    int ret = -1;
    while (1){
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
        if ( isclean && un_comp!=nullptr) {
            compress->backEmpty(un_comp);
        }
//        printf("here?\n");
        un_comp = compress->getUnCompressData();
//        printf("here??\n");
        if (un_comp == nullptr) {
            break;
        }
        /*
         *  放满一整个 bam_complete_block
         */
//        printf("here???\n");
//        printf("last use block len is %d\n",last_use_block_length);

        ret = un_comp -> split_pos;
//        ret = find_divide_pos(un_comp);
//        printf("ret number is %d\n",ret);
//        if (ret < 0){
//            printf("unsigned int is wrong\n");
//        }
//            Rabbit_memcpy(&need_block_len,un_comp->data+last_use_block_length,4);
        need_block_len=ret;
//        printf("need block len is %d\n",need_block_len);
//        printf("now_push_length is %d\n",now_push_length);
//        printf("un comp length is %d\n",un_comp->length);
        if (assign_block->length + need_block_len > assign_block->data_size){
            completeBlock->inputCompleteBlock(assign_block);
            assign_block = completeBlock->getEmpty();
        }
        if (ret!=un_comp->length){
//            printf("Input This\n");

//            printf("un comp length is %d\n",un_comp->length);
            memcpy(assign_block->data+assign_block->length, un_comp->data,ret*sizeof(char));
            assign_block->length+=ret;
            completeBlock->inputCompleteBlock(assign_block);
            assign_block = completeBlock->getEmpty();

            memcpy(assign_block->data+assign_block->length, un_comp->data+ret,(un_comp->length - ret)*sizeof(char));
            assign_block->length += (un_comp->length - ret);
            compress->backEmpty(un_comp);
            un_comp = compress->getUnCompressData();
            memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
            assign_block->length += un_comp->length;



            int divide_pos = 0;
//            int ans=0;
            int ret = 0;
            uint32_t x[8], new_l_data;
            while (divide_pos<assign_block->length){
                Rabbit_memcpy(&ret,assign_block->data+divide_pos,4);
//        printf("ret is %d\n",ret);
//        printf("divide_pos is %d\n",divide_pos);
//        printf("block length is %d\n",block->length);
                if (ret>=32){
                    if (divide_pos + 4 + 32 > assign_block->length){
                        compress->backEmpty(un_comp);
                        un_comp = compress->getUnCompressData();
                        if (assign_block->length+un_comp->length > assign_block->data_size){
                            change_data_size(assign_block);
                        }
                        memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
                        assign_block->length += un_comp->length;
//                        break;
                    }
                    Rabbit_memcpy(x,assign_block->data+divide_pos+4,32);
                    int pos = (int32_t)x[1];
                    int l_qname = x[2]&0xff;
                    int l_extranul = (l_qname%4 != 0)? (4 - l_qname%4) : 0;
                    int n_cigar = x[3]&0xffff;
                    int l_qseq = x[4];
                    new_l_data = ret - 32 + l_extranul;//block_len + c->l_extranul
                    if (new_l_data > INT_MAX || l_qseq < 0 || l_qname < 1) {
//                printf("ai 32 我的老天爷啊\n");
                        divide_pos+=4+32;
                        continue;
                    }
                    if (((uint64_t) n_cigar << 2) + l_qname + l_extranul
                        + (((uint64_t) l_qseq + 1) >> 1) + l_qseq > (uint64_t) new_l_data){
//                printf("ai 32 我的老天爷啊\n");
                        divide_pos+=4+32;
                        continue;
                    }
                    while (divide_pos + 4 + 32 + l_qname > assign_block->length){
                        compress->backEmpty(un_comp);
                        un_comp = compress->getUnCompressData();
                        if (assign_block->length+un_comp->length > assign_block->data_size){
                            change_data_size(assign_block);
                        }
                        memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
                        assign_block->length += un_comp->length;
                    }
                    char fg_char;
                    Rabbit_memcpy(&fg_char,assign_block->data+divide_pos+4+32+l_qname-1,1);
                    if (fg_char != '\0') {
//                printf("this is wrong\n");
                    }
                    if (fg_char != '\0' && l_extranul <=0 && new_l_data > INT_MAX -4 ){

                        while (divide_pos + 4 + 32 + l_qname > assign_block->length){
                            compress->backEmpty(un_comp);
                            un_comp = compress->getUnCompressData();
                            if (assign_block->length+un_comp->length > assign_block->data_size){
                                change_data_size(assign_block);
                            }
                            memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
                            assign_block->length += un_comp->length;
                        }
                        divide_pos+=4+32+l_qname;
                        continue;
                    }

                    while (divide_pos + 4 + ret > assign_block->length) {
                        compress->backEmpty(un_comp);
                        un_comp = compress->getUnCompressData();
                        if (assign_block->length+un_comp->length > assign_block->data_size){
                            change_data_size(assign_block);
                        }
                        memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
                        assign_block->length += un_comp->length;
                    }
                    divide_pos+=4+ret;
//                    ans++;
                }else {
//            printf("BIG WRONG!!!\n");
                    if (divide_pos+4 > assign_block->length) {
                        compress->backEmpty(un_comp);
                        un_comp = compress->getUnCompressData();
                        if (assign_block->length+un_comp->length > assign_block->data_size){
                            change_data_size(assign_block);
                        }
                        memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
                        assign_block->length += un_comp->length;
                    }
                    divide_pos+=4;
                }
//        printf("One Block Size is %d\n",ret);

            }
//            while (find_divide_pos(assign_block) != assign_block->length){
////                printf("find divide pos is %d\n",find_divide_pos(assign_block));
////                printf("assign block length is %d\n",assign_block->length);
//                compress->backEmpty(un_comp);
//                un_comp = compress->getUnCompressData();
////                printf("assign block length is %d\n",assign_block->length);
////                printf("assign block data size is %d\n",assign_block->data_size);
////                printf("un comp length is %d\n",un_comp->length);
//                if (assign_block->length+un_comp->length > assign_block->data_size){
//                    change_data_size(assign_block);
//                }
//                memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
//                assign_block->length += un_comp->length;
////                printf("assign_block->length is oooooooo %d\n",assign_block->length);
//            }
//            printf("OK is Over");
            completeBlock->inputCompleteBlock(assign_block);
            assign_block = completeBlock->getEmpty();
        } else {
//            printf("Input Two?\n");
            if (ret != un_comp->length ) {
//                printf("ai nan ding\n");
                break;
            }
            memcpy(assign_block->data+assign_block->length, un_comp->data,ret*sizeof(char));
            assign_block->length += ret ;
            last_use_block_length = 0;
            isclean = true;
        }
    }
//    printf("here??????????????????????\n");
    if (assign_block->length != 0){
        completeBlock->inputCompleteBlock(assign_block);
    }
    completeBlock->is_over();
    //        if (assign_block->length + un_comp->length >BGZF_MAX_BLOCK_COMPLETE_SIZE){
//            completeBlock->inputCompleteBlock(assign_block);
//            assign_block = completeBlock->getEmpty();
//        }

//        int ret = find_divide_pos(un_comp);
//        if (ret != un_comp->length){
//            printf("ret == %d  block length == %d\n",ret,un_comp->length);
//        }
//        memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
//        assign_block->length += un_comp->length;
//        compress->backEmpty(un_comp);

}
void compress_test_pack(BamCompress* compress){
    bam_block *un_comp = nullptr;
    while (1){
        un_comp = compress->getUnCompressData();
        if (un_comp == nullptr) {
            break;
        }
        compress->backEmpty(un_comp);

    }
}
void benchmark_pack(BamCompress* compress,BamCompleteBlock* completeBlock){
    bam_block *un_comp = nullptr;
    bam_complete_block *assign_block = completeBlock->getEmpty();
    int need_block_len=0,input_length=0;
    int last_use_block_length=0;
    bool isclean = true;
    int ret = -1;
    int bam_number = 0;
    while (1){
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
        if ( isclean && un_comp!=nullptr) {
            compress->backEmpty(un_comp);
        }
//        printf("here?\n");
        un_comp = compress->getUnCompressData();
//        printf("here??\n");
        if (un_comp == nullptr) {
            break;
        }
        /*
         *  放满一整个 bam_complete_block
         */
//        printf("here???\n");
//        printf("last use block len is %d\n",last_use_block_length);

        ret = un_comp -> split_pos;
//        ret = find_divide_pos(un_comp);
//        printf("ret number is %d\n",ret);
//        if (ret < 0){
//            printf("unsigned int is wrong\n");
//        }
//            Rabbit_memcpy(&need_block_len,un_comp->data+last_use_block_length,4);
        need_block_len=ret;
//        printf("need block len is %d\n",need_block_len);
//        printf("now_push_length is %d\n",now_push_length);
//        printf("un comp length is %d\n",un_comp->length);
        if (assign_block->length + need_block_len > assign_block->data_size){
            completeBlock->backEmpty(assign_block);
            assign_block = completeBlock->getEmpty();
        }
        //if (ret!=un_comp->length){ // 该分支未经测试
if (0){
	
	
	//            printf("Input This\n");

//            printf("un comp length is %d\n",un_comp->length);
            memcpy(assign_block->data+assign_block->length, un_comp->data,ret*sizeof(char));
            assign_block->length+=ret;
            completeBlock->backEmpty(assign_block);
            bam_number+=find_divide_pos_and_get_read_number(assign_block).second;
            assign_block = completeBlock->getEmpty();

            memcpy(assign_block->data+assign_block->length, un_comp->data+ret,(un_comp->length - ret)*sizeof(char));
            assign_block->length += (un_comp->length - ret);
            compress->backEmpty(un_comp);
            un_comp = compress->getUnCompressData();
            memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
            assign_block->length += un_comp->length;
            while (find_divide_pos(assign_block) != assign_block->length){
//                printf("find divide pos is %d\n",find_divide_pos(assign_block));
//                printf("assign block length is %d\n",assign_block->length);
                compress->backEmpty(un_comp);
                un_comp = compress->getUnCompressData();
//                printf("assign block length is %d\n",assign_block->length);
//                printf("assign block data size is %d\n",assign_block->data_size);
//                printf("un comp length is %d\n",un_comp->length);
                if (assign_block->length+un_comp->length > assign_block->data_size){
                    change_data_size(assign_block);
                }
                memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
                assign_block->length += un_comp->length;
//                printf("assign_block->length is oooooooo %d\n",assign_block->length);
            }

//            printf("OK is Over");
            bam_number+=find_divide_pos_and_get_read_number(assign_block).second;
            completeBlock->backEmpty(assign_block);
            assign_block = completeBlock->getEmpty();
        } else {
//            printf("Input Two?\n");
            if (ret != un_comp->length ) {
                printf("ai nan ding\n");
                break;
            }
//            printf("1\n");
//            printf("assign block length is %d\n",assign_block->length);
//            printf("assign block data size is %d\n",assign_block->data_size);
//            printf("un comp length is %d\n",un_comp->length);
            //memcpy(assign_block->data+assign_block->length, un_comp->data,ret*sizeof(unsigned char));

//            printf("2\n");
            assign_block->length += ret ;
            last_use_block_length = 0;
            isclean = true;
            bam_number+=un_comp->bam_number;
        }
    }
//    printf("here??????????????????????\n");
    if (assign_block->length != 0){
        completeBlock->backEmpty(assign_block);
    }
    completeBlock->is_over();
    printf("Bam number is %d\n",bam_number);
    //        if (assign_block->length + un_comp->length >BGZF_MAX_BLOCK_COMPLETE_SIZE){
//            completeBlock->inputCompleteBlock(assign_block);
//            assign_block = completeBlock->getEmpty();
//        }

//        int ret = find_divide_pos(un_comp);
//        if (ret != un_comp->length){
//            printf("ret == %d  block length == %d\n",ret,un_comp->length);
//        }
//        memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
//        assign_block->length += un_comp->length;
//        compress->backEmpty(un_comp);

}
void benchmark_bam_pack(BamCompleteBlock* completeBlock,samFile *output,sam_hdr_t *hdr){

    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    hts_set_threads(output, 20);
    bam_complete_block* un_comp;
    long long ans = 0;
    long long res = 0;
    while (1){
        un_comp = completeBlock->getCompleteBlock();
        if (un_comp == nullptr){
            break;
        }
//        printf("assign over block length is %d\n",un_comp->length);
        int ret;
        while ((ret=(read_bam(un_comp,b,0)))>=0) {
//            printf("One Bam1_t Size is %d\n",ret);
//            printf("This Bam1_t Char Number is %d\n",b->core.l_qseq);


//            sam_write1(output,hdr,b);
            ans++;
        }
        res++;

        completeBlock->backEmpty(un_comp);
    }
    printf("Bam1_t Number is %lld\n",ans);
    printf("Block  Number is %lld\n",res);
}




//int rabbit_write_deflate_block(BGZF *fp, bam_write_block* write_block){
//    size_t comp_size = BGZF_MAX_BLOCK_SIZE;
//    int ret;
//   if ( !fp->is_gzip )
//        ret = rabbit_bgzf_compress(write_block->compressed_data, &comp_size, write_block->uncompressed_data, write_block->block_offset, fp->compress_level);
//    else
//        ret = rabbit_bgzf_gzip_compress(fp, write_block->compressed_data, &comp_size, write_block->uncompressed_data, write_block->block_offset, fp->compress_level);
//
//    if ( ret != 0 )
//    {
//        hts_log_debug("Compression error %d", ret);
//        fp->errcode |= BGZF_ERR_ZLIB;
//        return -1;
//    }
//    return comp_size;
//}
//int rabbit_bgzf_flush(BGZF *fp,bam_write_block* write_block)
//{
//    //TODO 此处可能会出现问题
//    while (write_block->block_offset > 0) {
//        int block_length;
//        printf("Write Block Offset : %d\n",write_block->block_offset);
//        block_length = rabbit_write_deflate_block(fp, write_block);
//        if (block_length < 0) {
//            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block_length, NULL));
//            return -1;
//        }
//        if (write_block== nullptr){
//            printf("Write Game Over!!!\n");
//        }
//        if (hwrite(fp->fp, write_block->compressed_data, block_length) != block_length) {
//            printf("Write Failed\n");
//            hts_log_error("File write failed (wrong size)");
//            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
//            return -1;
//        }
//
//
//        write_block->block_offset=0;
//        fp->block_address += block_length;
//    }
//    write_block->block_offset=0;
//    return 0;
//}
//int rabbit_bgzf_mul_flush(BGZF *fp,BamWriteCompress *bam_write_compress,bam_write_block* &write_block)
//{
////    printf("Try to input One Uncompressed data\n");
////    InputBlockNum++;
////    printf("Write Block Offset : %d\n",write_block->block_offset);
//    bam_write_compress->inputUnCompressData(write_block);
//    write_block=bam_write_compress->getEmpty();
////    printf("Get Another Empty Block Block Num : %d\n",write_block->block_num);
//    return 0;
//}
//int rabbit_bgzf_write(BGZF *fp,bam_write_block* &write_block,const void *data, size_t length)
//{
//    const uint8_t *input = (const uint8_t*)data;
//    ssize_t remaining = length;
////    assert(fp->is_write);
//    while (remaining > 0) {
//        uint8_t* buffer = (uint8_t*)write_block->uncompressed_data;
//        int copy_length = BGZF_BLOCK_SIZE - write_block->block_offset;
//        if (copy_length > remaining) copy_length = remaining;
//        memcpy(buffer + write_block->block_offset, input, copy_length);
//        write_block->block_offset += copy_length;
//        input += copy_length;
//        remaining -= copy_length;
//        if (write_block->block_offset == BGZF_BLOCK_SIZE) {
//            //BUG  write block 不是 引用，所以没有改变
//            if (rabbit_bgzf_flush(fp,write_block) != 0) return -1;
//        }
//    }
//    return length - remaining;
//}
//int rabbit_bgzf_mul_write(BGZF *fp, BamWriteCompress *bam_write_compress,bam_write_block* &write_block,const void *data, size_t length)
//{
//    const uint8_t *input = (const uint8_t*)data;
//    ssize_t remaining = length;
////    assert(fp->is_write);
//    while (remaining > 0) {
//        uint8_t* buffer = (uint8_t*)write_block->uncompressed_data;
////        printf("In rabbit bgzf mul write block Block Num : %d\n",write_block->block_num);
//        int copy_length = BGZF_BLOCK_SIZE - write_block->block_offset;
//        if (copy_length > remaining) copy_length = remaining;
////        printf("In rabbit bgzf mul write block Block Num : %d\n",write_block->block_num);
//        memcpy(buffer + write_block->block_offset, input, copy_length);
//        write_block->block_offset += copy_length;
////        printf("In rabbit bgzf mul write block Block Num : %d\n",write_block->block_num);
//        input += copy_length;
//        remaining -= copy_length;
//        if (write_block->block_offset == BGZF_BLOCK_SIZE) {
//            if (rabbit_bgzf_mul_flush(fp,bam_write_compress,write_block) != 0) return -1;
//        }
//    }
//    return length - remaining;
//}
//int rabbit_bgzf_flush_try(BGZF *fp, bam_write_block* write_block,ssize_t size)
//{
//    if (write_block->block_offset + size > BGZF_BLOCK_SIZE) {
//        return rabbit_bgzf_flush(fp,write_block);
//    }
//    return 0;
//}
//int rabbit_bgzf_mul_flush_try(BGZF *fp,BamWriteCompress* bam_write_compress,bam_write_block* &write_block,ssize_t size)
//{
//    if (write_block->block_offset + size > BGZF_BLOCK_SIZE) {
//        return rabbit_bgzf_mul_flush(fp,bam_write_compress,write_block);
//    }
//    return 0;
//}
//int bam_write_pack(BGZF *fp,BamWriteCompress *bam_write_compress){
//    bam_write_block* block;
//    while(1){
////        printf("Try to Get Compress Data\n");
//        block=bam_write_compress->getCompressData();
////        printf("Has Get One Compress Data\n");
//        if (block == nullptr){
//            return 0;
//        }
////        std::this_thread::sleep_for(std::chrono::nanoseconds(5));
//        if (block->block_length < 0) {
//            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block->block_length, NULL));
//            return -1;
//        }
//        if (hwrite(fp->fp, block->compressed_data, block->block_length) != block->block_length) {
////            printf("Write Failed\n");
//            hts_log_error("File write failed (wrong size)");
//            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
//            return -1;
//        }
////        printf("Has write One Block\n");
//        block->block_offset=0;
//        fp->block_address += block->block_length;
//        bam_write_compress->backEmpty(block);
//    }
//}
//void bam_write_compress_pack(BGZF *fp,BamWriteCompress *bam_write_compress){
////    printf("Start Compress\n");
//    bam_write_block * block;
//    while (1){
//        // fg = getRead(comp);
//        //printf("%d is not get One compressed data\n",id);
////        printf("Has Start Try to Get One Uncompress\n");
//        block=bam_write_compress->getUnCompressData();
////        printf("Has get One Uncompress data\n");
//
//        //printf("%d is get One compressed data\n",id);
//        if (block == nullptr) {
//            //printf("%d is Over\n",id);
//            break;
//        }
////        printf("This Uncompress data block num : %d\n",block->block_num);
//        /*
//         * 压缩
//         */
////        int block_num=block->block_num;
//        block->block_length = rabbit_write_deflate_block(fp, block);
////        printf("Has Compress One Block\n");
//        bam_write_compress->inputCompressData(block);
////        printf("Can,t input Compress Data\n");
////        while (!compress->tryinputUnCompressData(un_comp,comp.second)){
////            std::this_thread::sleep_for(std::chrono::milliseconds(1));
////        }
//    }
////    printf("One Compress Thread has been over!\n");
//    bam_write_compress->CompressThreadComplete();
//}
//
//int rabbit_bam_write_test(BGZF *fp,bam_write_block* write_block,bam1_t *b){
//    const bam1_core_t *c = &b->core;
//    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
//    int i, ok;
//    if (c->l_qname - c->l_extranul > 255) {
//        hts_log_error("QNAME \"%s\" is longer than 254 characters", bam_get_qname(b));
//        errno = EOVERFLOW;
//        return -1;
//    }
//    if (c->n_cigar > 0xffff) block_len += 16; // "16" for "CGBI", 4-byte tag length and 8-byte fake CIGAR
//    if (c->pos > INT_MAX ||
//        c->mpos > INT_MAX ||
//        c->isize < INT_MIN || c->isize > INT_MAX) {
//        hts_log_error("Positional data is too large for BAM format");
//        return -1;
//    }
//    x[0] = c->tid;
//    x[1] = c->pos;
//    x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | (c->l_qname - c->l_extranul);
//    if (c->n_cigar > 0xffff) x[3] = (uint32_t)c->flag << 16 | 2;
//    else x[3] = (uint32_t)c->flag << 16 | (c->n_cigar & 0xffff);
//    x[4] = c->l_qseq;
//    x[5] = c->mtid;
//    x[6] = c->mpos;
//    x[7] = c->isize;
//    ok = (rabbit_bgzf_flush_try(fp, write_block, 4 + block_len) >= 0);
//    if (fp->is_be) {
//        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
//        y = block_len;
//        if (ok) ok = (rabbit_bgzf_write(fp, write_block,ed_swap_4p(&y), 4) >= 0);
//        swap_data(c, b->l_data, b->data, 1);
//    } else {
//        if (ok) {
//            ok = (rabbit_bgzf_write(fp, write_block, &block_len, 4) >= 0);
//        }
//    }
//    if (ok) {
//        ok = (rabbit_bgzf_write(fp, write_block, x, 32) >= 0);
//    }
//    if (ok) ok = (rabbit_bgzf_write(fp, write_block, b->data, c->l_qname - c->l_extranul) >= 0);
//    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
//        if (ok) ok = (rabbit_bgzf_write(fp, write_block, b->data + c->l_qname, b->l_data - c->l_qname) >= 0);
//    } else { // with long CIGAR, insert a fake CIGAR record and move the real CIGAR to the CG:B,I tag
//        uint8_t buf[8];
//        uint32_t cigar_st, cigar_en, cigar[2];
//        hts_pos_t cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
//        if (cigreflen >= (1<<28)) {
//            // Length of reference covered is greater than the biggest
//            // CIGAR operation currently allowed.
//            hts_log_error("Record %s with %d CIGAR ops and ref length %" PRIhts_pos
//                                  " cannot be written in BAM.  Try writing SAM or CRAM instead.\n",
//                          bam_get_qname(b), c->n_cigar, cigreflen);
//            return -1;
//        }
//        cigar_st = (uint8_t*)bam_get_cigar(b) - b->data;
//        cigar_en = cigar_st + c->n_cigar * 4;
//        cigar[0] = (uint32_t)c->l_qseq << 4 | BAM_CSOFT_CLIP;
//        cigar[1] = (uint32_t)cigreflen << 4 | BAM_CREF_SKIP;
//        u32_to_le(cigar[0], buf);
//        u32_to_le(cigar[1], buf + 4);
//        if (ok) {
//            ok = (rabbit_bgzf_write(fp, write_block, buf, 8) >= 0); // write cigar: <read_length>S<ref_length>N
//        }
//        if (ok) {
//            ok = (rabbit_bgzf_write(fp, write_block, &b->data[cigar_en], b->l_data - cigar_en) >= 0); // write data after CIGAR
//        }
//        if (ok) {
//            ok = (rabbit_bgzf_write(fp, write_block, "CGBI", 4) >= 0); // write CG:B,I
//        }
//        u32_to_le(c->n_cigar, buf);
//        if (ok) {
//            ok = (rabbit_bgzf_write(fp, write_block, buf, 4) >= 0); // write the true CIGAR length
//        }
//        if (ok) {
//            ok = (rabbit_bgzf_write(fp, write_block, &b->data[cigar_st], c->n_cigar * 4) >= 0); // write the real CIGAR
//        }
//    }
//    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
//    return ok? 4 + block_len : -1;
//}
//int rabbit_bam_write_mul_test(BGZF *fp,BamWriteCompress *bam_write_compress,bam_write_block* &write_block,bam1_t *b){
//    const bam1_core_t *c = &b->core;
//    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
//    int i, ok;
//    if (c->l_qname - c->l_extranul > 255) {
//        hts_log_error("QNAME \"%s\" is longer than 254 characters", bam_get_qname(b));
//        errno = EOVERFLOW;
//        return -1;
//    }
//    if (c->n_cigar > 0xffff) block_len += 16; // "16" for "CGBI", 4-byte tag length and 8-byte fake CIGAR
//    if (c->pos > INT_MAX ||
//        c->mpos > INT_MAX ||
//        c->isize < INT_MIN || c->isize > INT_MAX) {
//        hts_log_error("Positional data is too large for BAM format");
//        return -1;
//    }
//    x[0] = c->tid;
//    x[1] = c->pos;
//    x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | (c->l_qname - c->l_extranul);
//    if (c->n_cigar > 0xffff) x[3] = (uint32_t)c->flag << 16 | 2;
//    else x[3] = (uint32_t)c->flag << 16 | (c->n_cigar & 0xffff);
//    x[4] = c->l_qseq;
//    x[5] = c->mtid;
//    x[6] = c->mpos;
//    x[7] = c->isize;
//    ok = (rabbit_bgzf_mul_flush_try(fp,bam_write_compress, write_block, 4 + block_len) >= 0);
//    if (fp->is_be) {
//        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
//        y = block_len;
//        if (ok) ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block,ed_swap_4p(&y), 4) >= 0);
//        swap_data(c, b->l_data, b->data, 1);
//    } else {
//        if (ok) {
//            ok = (rabbit_bgzf_mul_write(fp,bam_write_compress, write_block, &block_len, 4) >= 0);
//        }
//    }
//    if (ok) {
//        ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, x, 32) >= 0);
//    }
//    if (ok) ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, b->data, c->l_qname - c->l_extranul) >= 0);
//    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
//        if (ok) ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, b->data + c->l_qname, b->l_data - c->l_qname) >= 0);
//    } else { // with long CIGAR, insert a fake CIGAR record and move the real CIGAR to the CG:B,I tag
//        uint8_t buf[8];
//        uint32_t cigar_st, cigar_en, cigar[2];
//        hts_pos_t cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
//        if (cigreflen >= (1<<28)) {
//            // Length of reference covered is greater than the biggest
//            // CIGAR operation currently allowed.
//            hts_log_error("Record %s with %d CIGAR ops and ref length %" PRIhts_pos
//                                  " cannot be written in BAM.  Try writing SAM or CRAM instead.\n",
//                          bam_get_qname(b), c->n_cigar, cigreflen);
//            return -1;
//        }
//        cigar_st = (uint8_t*)bam_get_cigar(b) - b->data;
//        cigar_en = cigar_st + c->n_cigar * 4;
//        cigar[0] = (uint32_t)c->l_qseq << 4 | BAM_CSOFT_CLIP;
//        cigar[1] = (uint32_t)cigreflen << 4 | BAM_CREF_SKIP;
//        u32_to_le(cigar[0], buf);
//        u32_to_le(cigar[1], buf + 4);
//        if (ok) {
//            ok = (rabbit_bgzf_mul_write(fp,bam_write_compress, write_block, buf, 8) >= 0); // write cigar: <read_length>S<ref_length>N
//        }
//        if (ok) {
//            ok = (rabbit_bgzf_mul_write(fp,bam_write_compress,write_block, &b->data[cigar_en], b->l_data - cigar_en) >= 0); // write data after CIGAR
//        }
//        if (ok) {
//            ok = (rabbit_bgzf_mul_write(fp,bam_write_compress,write_block, "CGBI", 4) >= 0); // write CG:B,I
//        }
//        u32_to_le(c->n_cigar, buf);
//        if (ok) {
//            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, buf, 4) >= 0); // write the true CIGAR length
//        }
//        if (ok) {
//            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, &b->data[cigar_st], c->n_cigar * 4) >= 0); // write the real CIGAR
//        }
//    }
//    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
//    return ok? 4 + block_len : -1;
//}

//void benchmark_write_pack(BamCompleteBlock* completeBlock,samFile *output,sam_hdr_t *hdr,int level){
//
//    uint8_t* compress_block_test = new uint8_t[BGZF_BLOCK_SIZE];
//    uint8_t* uncompress_block_test = new uint8_t[BGZF_BLOCK_SIZE];
//    output->fp.bgzf->block_offset=0;
//    output->fp.bgzf->uncompressed_block=uncompress_block_test;
//    output->fp.bgzf->compressed_block=compress_block_test;
//    output->fp.bgzf->compress_level=level;
//    bam1_t *b;
//    if ((b = bam_init1()) == NULL) {
//        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
//    }
//    if (sam_hdr_write(output, hdr) != 0){
//        printf("HDR Write False!\n");
//        return ;
//    }
//    bam_complete_block* un_comp;
//    long long ans = 0;
//    long long res = 0;
//    bam_write_block *write_block=new bam_write_block();
//    write_block->block_offset=0;
//    write_block->uncompressed_data=new uint8_t[BGZF_BLOCK_SIZE];
//    write_block->compressed_data=new uint8_t[BGZF_BLOCK_SIZE];
//    write_block->status=0;
//    while (1){
//        un_comp = completeBlock->getCompleteBlock();
//        if (un_comp == nullptr){
//            break;
//        }
////        printf("assign over block length is %d\n",un_comp->length);
//        int ret;
//        while ((ret=(read_bam(un_comp,b,0)))>=0) {
////            printf("One Bam1_t Size is %d\n",ret);
////            printf("This Bam1_t Char Number is %d\n",b->core.l_qseq);
///*
// *  尝试单线程输出
// */
////            sam_write1(output,hdr,b);
//            rabbit_bam_write_test(output->fp.bgzf,write_block,b);
//            ans++;
//        }
//        res++;
//        completeBlock->backEmpty(un_comp);
//    }
//    printf("Bam1_t Number is %lld\n",ans);
//    printf("Block  Number is %lld\n",res);
//}
//
//void benchmark_write_mul_pack(BamCompleteBlock* completeBlock,BamWriteCompress* bam_write_compress,samFile *output,sam_hdr_t *hdr,int level){
//
//    output->fp.bgzf->block_offset=0;
//    output->fp.bgzf->compress_level=level;
//    bam1_t *b;
//    if ((b = bam_init1()) == NULL) {
//        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
//    }
//    if (sam_hdr_write(output, hdr) != 0){
//        printf("HDR Write False!\n");
//        return ;
//    }
//    bam_complete_block* un_comp;
//    long long ans = 0;
//    long long res = 0;
//    bam_write_block *write_block=bam_write_compress->getEmpty();
////    printf("Write Pack Num : %d\n",write_block->block_num);
////    printf("Write Pack Block Offset : %d\n",write_block->block_offset);
////    printf("Has get One Empty!\n");
//    int bam_num=1;
//    while (1){
//        un_comp = completeBlock->getCompleteBlock();
//        if (un_comp == nullptr){
//            break;
//        }
////        printf("assign over block length is %d\n",un_comp->length);
//        int ret;
//        while ((ret=(read_bam(un_comp,b,0)))>=0) {
////            printf("One Bam1_t Size is %d\n",ret);
////            printf("This Bam1_t Char Number is %d\n",b->core.l_qseq);
///*
// *  尝试单线程输出
// */
////            sam_write1(output,hdr,b);
////            printf("Try to output one Bam1_t : %d\n",bam_num++);
//            rabbit_bam_write_mul_test(output->fp.bgzf,bam_write_compress,write_block,b);
//            ans++;
//        }
//        res++;
//        completeBlock->backEmpty(un_comp);
//    }
//    if (write_block->block_offset>0) {
//        InputBlockNum++;
////        rabbit_bgzf_mul_flush(output->fp.bgzf,bam_write_compress,write_block);
//        bam_write_compress->inputUnCompressData(write_block);
//    }
////    printf("The Input Block Num : %d\n",InputBlockNum);
//    bam_write_compress->WriteComplete();
////    printf("Bam1_t Number is %lld\n",ans);
////    printf("Block  Number is %lld\n",res);
//}

void basic_status_pack(BamCompleteBlock* completeBlock,BamStatus *status){
    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        printf("bam1_t is not ready\n");
    }
    bam_complete_block* un_comp;
    while (1){
        un_comp = completeBlock->getCompleteBlock();
        if (un_comp == nullptr){
            break;
        }
        int ret;
        while ((ret=(read_bam(un_comp,b,0)))>=0) {

            status->statusbam(b);
        }
        completeBlock->backEmpty(un_comp);
    }

}

int main(int argc,char* argv[]){





    CLI::App app("RabbitBAM");

    string inputfile;

    string outputfile("./BAMStatus.html");
    int n_thread=1;
    int n_thread_write=1;
    int level = 6;

    CLI::App *bam2fq = app.add_subcommand("bam2fq", "BAM format turn to FastQ format");
    bam2fq->add_option("-i", inputfile, "input File name")->required()->check(CLI::ExistingFile);
    bam2fq->add_option("-o", outputfile, "output File name");
    bam2fq->add_option("-w,-@,-n,--threads",n_thread,"thread number");


    CLI::App *bamstatus = app.add_subcommand("bamstatus", "Analyze BAM files");
    bamstatus->add_option("-i", inputfile, "input File name")->required()->check(CLI::ExistingFile);
    bamstatus->add_option("-o", outputfile, "output File name");
    bamstatus->add_option("-w,-@,-n,--threads",n_thread,"thread number");


    CLI::App *benchmark = app.add_subcommand("benchmark", "Performance Testing");
    benchmark->add_option("-i", inputfile, "input File name")->required()->check(CLI::ExistingFile);
    benchmark->add_option("-o", outputfile, "output File name");
    benchmark->add_option("-w,-@,-n,--threads",n_thread,"thread number");


    CLI::App *htslib_test = app.add_subcommand("htslib_test", "Htslib sam_read API Performance Testing");
    htslib_test->add_option("-i", inputfile, "input File name")->required()->check(CLI::ExistingFile);
    htslib_test->add_option("-o", outputfile, "output File name");
    htslib_test->add_option("-w,-@,-n,--threads",n_thread,"thread number");


    CLI::App *benchmark_count = app.add_subcommand("benchmark_count", "Banchmark Count Performance Testing");
    benchmark_count->add_option("-i", inputfile, "input File name")->required()->check(CLI::ExistingFile);
    benchmark_count->add_option("-o", outputfile, "output File name");
    benchmark_count->add_option("-w,-@,-n,--threads",n_thread,"thread number");



    CLI::App *compress_test = app.add_subcommand("compress_test", "Compress Performance Testing");
    compress_test->add_option("-i", inputfile, "input File name")->required()->check(CLI::ExistingFile);
    compress_test->add_option("-o", outputfile, "output File name");
    compress_test->add_option("-w,-@,-n,--threads",n_thread,"thread number");



    CLI::App *write_test = app.add_subcommand("write_test", "Write Testing One thread");
    write_test->add_option("-i", inputfile, "input File name")->required()->check(CLI::ExistingFile);
    write_test->add_option("-o", outputfile, "output File name");
    write_test->add_option("-w,-@,-n,--threads",n_thread,"thread number");
    write_test->add_option("-l,--level",level,"zip level");


    CLI::App *write_mul_test = app.add_subcommand("write_mul_test", "Write Testing with multi thread");
    write_mul_test->add_option("-i", inputfile, "input File name")->required()->check(CLI::ExistingFile);
    write_mul_test->add_option("-o", outputfile, "output File name");
    write_mul_test->add_option("--nr",n_thread,"Read thread number");
    write_mul_test->add_option("--nw",n_thread_write,"Write thread number");
    write_mul_test->add_option("-l,--level",level,"zip level");


    CLI::App *api_test = app.add_subcommand("api_test","use api to read and write");
    api_test->add_option("-i",inputfile,"input File name")->required()->check(CLI::ExistingFile);
    api_test->add_option("-o",outputfile,"output File name");
    api_test->add_option("--nr",n_thread,"Read thread number");
    api_test->add_option("--nw",n_thread_write,"Write thread number");
    api_test->add_option("-l,--level",level,"zip level");



    CLI11_PARSE(app, argc, argv);
    if (app.get_subcommands().size()>1){
        printf("you should input one command!!!\n");
        return 0;
    }
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "bamstatus")==0){
        TDEF(fq)
        TSTART(fq)
        printf("Starting Running Bam Analyze\n");
        samFile *sin;
        sam_hdr_t *hdr;
        ofstream fout;
        fout.open(outputfile);
        if ((sin=sam_open(inputfile.c_str(),"r"))==NULL){
            printf("Can`t open this file!\n");
            return 0;
        }
        if ((hdr = sam_hdr_read(sin)) == NULL) {
            return  0;
        }
        /*
         *  读取和处理准备
         */
        BamRead read(2000);
        BamCompress compress(8000,n_thread);
        BamCompleteBlock completeBlock(10);

        /*
         * 分析准备
         */
        BamStatus **status=new BamStatus*[n_thread];


        thread **Bam = new thread *[n_thread];



        thread *read_thread = new thread(&read_pack,sin->fp.bgzf,&read);
        thread **compress_thread = new thread *[n_thread];
        for (int i=0;i<n_thread;i++){
            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(i,&cpuset);
            compress_thread[i]=new thread(&compress_pack,&read,&compress);
            int rc =pthread_setaffinity_np(compress_thread[i]->native_handle(),sizeof(cpu_set_t), &cpuset);

        }
        thread *assign_thread = new thread(&assign_pack,&compress,&completeBlock);
//        thread *consumer_thread = new thread(&benchmark_pack,&completeBlock);
        for (int i=0;i<n_thread;i++){
            status[i]=new BamStatus(inputfile);
            Bam[i]=new thread(&basic_status_pack,&completeBlock,status[i]);
        }
        read_thread->join();
        for (int i=0;i<n_thread;i++) compress_thread[i]->join();
        assign_thread->join();

        for (int i=0;i<n_thread;i++) Bam[i]->join();
        cout << "In here ?" << endl;
        for (int i=1;i<n_thread;i++) {
            status[0]->add(status[i]);
        }
        status[0]->statusAll();
        status[0]->reportHTML(&fout);
        sam_close(sin);
        TEND(fq)
        TPRINT(fq,"time is : ");
    }
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "benchmark")==0){
        TDEF(fq)
        TSTART(fq)
        printf("Starting Running Benchmark\n");
        printf("BGZF_MAX_BLOCK_COMPLETE_SIZE is %d\n",BGZF_MAX_BLOCK_COMPLETE_SIZE);
        //if (strcmp(outputfile.substr(outputfile.size()-4).c_str(),"html")==0) outputfile=("./output.fastq");
        samFile *sin;
        sam_hdr_t *hdr;
        samFile *output;
        if ((sin=sam_open(inputfile.c_str(),"r"))==NULL){
            printf("Can`t open this file!\n");
            return 0;
        }
        if ((output=sam_open(outputfile.c_str(),"w"))==NULL){
            printf("Can`t open this file!\n");
            return 0;
        }
        if ((hdr = sam_hdr_read(sin)) == NULL) {
            return  0;
        }
        /*
        *  读取和处理准备
        */
        BamRead read(8000);
        BamCompress compress(4000,n_thread);
        BamCompleteBlock completeBlock(200);

        printf("Malloc Memory is Over\n");
        /*
        * 分析准备
        */
        thread *read_thread = new thread(&read_pack,sin->fp.bgzf,&read);
        thread **compress_thread = new thread *[n_thread];



        for (int i=0;i<n_thread;i++){
            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(i,&cpuset);
            compress_thread[i]=new thread(&compress_pack,&read,&compress);
            int rc =pthread_setaffinity_np(compress_thread[i]->native_handle(),sizeof(cpu_set_t), &cpuset);
        }
        thread *assign_thread = new thread(&assign_pack,&compress,&completeBlock);
        //thread *consumer_thread = new thread(&benchmark_pack,&completeBlock);
        int  consumer_thread_number = 1;
        thread **consumer_thread = new thread*[consumer_thread_number];
        for (int i=0;i<consumer_thread_number;i++) consumer_thread[i] = new thread(&benchmark_bam_pack,&completeBlock,output,hdr);
        read_thread->join();
        for (int i=0;i<n_thread;i++) compress_thread[i]->join();
        assign_thread->join();
        //consumer_thread->join();
        for (int i=0;i<consumer_thread_number;i++) consumer_thread[i]->join();
        sam_close(sin);
        sam_close(output);
        printf("Wait num is %d\n",compress.wait_num);
        TEND(fq)
        TPRINT(fq,"time is : ");
    }
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "htslib_test")==0){
        TDEF(fq)
        TSTART(fq)
        printf("Starting Htslib sam_read API Running Benchmark\n");
        samFile *sin;
        sam_hdr_t *hdr;
        samFile *output;
        if ((sin=sam_open(inputfile.c_str(),"r"))==NULL){
            printf("Can`t open this file!\n");
            return 0;
        }
        if ((output=sam_open(outputfile.c_str(),"w"))==NULL){
            printf("Can`t open this file!\n");
            return 0;
        }
        output->format.compression_level=6;
        if ((hdr = sam_hdr_read(sin)) == NULL) {
            return  0;
        }

        sam_hdr_write(output,hdr);
        bam1_t *b;
        if ((b = bam_init1()) == NULL) {
            fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
        }
        htsThreadPool p = {NULL, 0};
        p.pool = hts_tpool_init(n_thread);
        hts_set_opt(sin,  HTS_OPT_THREAD_POOL, &p);
        hts_set_threads(output, n_thread);
        int num = 0;
//        output->format.format=bam;

        while(sam_read1(sin, hdr, b)>=0){
            num++;
            sam_write1(output,hdr,b);
        }
//        sam_itr_querys()
//        while (sam_itr_next(sin,,b)>=0){
//            num++;
//        }
        printf("Bam Number is %d\n",num);
        sam_close(sin);
        sam_close(output);
        TEND(fq)
        TPRINT(fq,"time is : ");
    }
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "benchmark_count")==0){
        TDEF(fq)
        TSTART(fq)
        printf("Starting Running Benchmark Count\n");
        printf("BGZF_MAX_BLOCK_COMPLETE_SIZE is %d\n",BGZF_MAX_BLOCK_COMPLETE_SIZE);
        if (strcmp(outputfile.substr(outputfile.size()-4).c_str(),"html")==0) outputfile=("./output.fastq");
        samFile *sin;
        sam_hdr_t *hdr;
        ofstream fout;
        fout.open(outputfile);
        if ((sin=sam_open(inputfile.c_str(),"r"))==NULL){
            printf("Can`t open this file!\n");
            return 0;
        }
        if ((hdr = sam_hdr_read(sin)) == NULL) {
            return  0;
        }
        /*
         *  读取和处理准备
         */
        BamRead read(2000);
        BamCompress compress(8000,n_thread);
        BamCompleteBlock completeBlock(10);

        printf("Malloc Memory is Over\n");
        /*
         * 分析准备
         */
        thread *read_thread = new thread(&read_pack,sin->fp.bgzf,&read);
        thread **compress_thread = new thread *[n_thread];

        for (int i=0;i<n_thread;i++){
            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(i,&cpuset);
            compress_thread[i]=new thread(&compress_pack,&read,&compress);
            int rc =pthread_setaffinity_np(compress_thread[i]->native_handle(),sizeof(cpu_set_t), &cpuset);

        }
        thread *assign_thread = new thread(&benchmark_pack,&compress,&completeBlock);
//        thread *consumer_thread = new thread(&benchmark_pack,&completeBlock);
        read_thread->join();
        for (int i=0;i<n_thread;i++) compress_thread[i]->join();
        assign_thread->join();
//        consumer_thread->join();
        sam_close(sin);
        printf("Wait num is %d\n",compress.wait_num);
        TEND(fq)
        TPRINT(fq,"time is : ");
    }
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "compress_test")==0){
        TDEF(fq)
        TSTART(fq)
        printf("Starting Running Compress Test Benchmark\n");
        printf("BGZF_MAX_BLOCK_COMPLETE_SIZE is %d\n",BGZF_MAX_BLOCK_COMPLETE_SIZE);
        if (strcmp(outputfile.substr(outputfile.size()-4).c_str(),"html")==0) outputfile=("./output.fastq");
        samFile *sin;
        sam_hdr_t *hdr;
        ofstream fout;
        fout.open(outputfile);
        if ((sin=sam_open(inputfile.c_str(),"r"))==NULL){
            printf("Can`t open this file!\n");
            return 0;
        }
        if ((hdr = sam_hdr_read(sin)) == NULL) {
            return  0;
        }
        /*
        *  读取和处理准备
        */
        BamRead read(8000);
        BamCompress compress(4000,n_thread);
        BamCompleteBlock completeBlock(200);

        printf("Malloc Memory is Over\n");
        /*
        * 分析准备
        */
        thread *read_thread = new thread(&read_pack,sin->fp.bgzf,&read);
        thread **compress_thread = new thread *[n_thread];



        for (int i=0;i<n_thread;i++){
            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(i,&cpuset);
            compress_thread[i]=new thread(&compress_pack,&read,&compress);
            int rc =pthread_setaffinity_np(compress_thread[i]->native_handle(),sizeof(cpu_set_t), &cpuset);
        }
        thread *assign_thread = new thread(&compress_test_pack,&compress);
        //thread *consumer_thread = new thread(&benchmark_pack,&completeBlock);
//        int  consumer_thread_number = 2;
//        thread **consumer_thread = new thread*[consumer_thread_number];
//        for (int i=0;i<consumer_thread_number;i++) consumer_thread[i] = new thread(&benchmark_bam_pack,&completeBlock);
        read_thread->join();
        for (int i=0;i<n_thread;i++) compress_thread[i]->join();
        assign_thread->join();
        //consumer_thread->join();
//        for (int i=0;i<consumer_thread_number;i++) consumer_thread[i]->join();
        sam_close(sin);
        printf("Wait num is %d\n",compress.wait_num);
        TEND(fq)
        TPRINT(fq,"time is : ");
    }
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "write_test")==0){
        TDEF(fq)
        TSTART(fq)
        printf("Starting Running Write Test\n");
        printf("BGZF_MAX_BLOCK_COMPLETE_SIZE is %d\n",BGZF_MAX_BLOCK_COMPLETE_SIZE);
        printf("output File Name is %s\n",outputfile.c_str());
        //if (strcmp(outputfile.substr(outputfile.size()-4).c_str(),"html")==0) outputfile=("./output.fastq");
        samFile *sin;
        sam_hdr_t *hdr;
        samFile *output;
        if ((sin=sam_open(inputfile.c_str(),"r"))==NULL){
            printf("Can`t open this file!\n");
            return 0;
        }

        if ((output=sam_open(outputfile.c_str(),"wb"))==NULL){
            printf("Can`t open this file!\n");
            return 0;
        }
        if ((hdr = sam_hdr_read(sin)) == NULL) {
            return  0;
        }
        /*
        *  读取和处理准备
        */
        BamRead read(8000);
        BamCompress compress(4000,n_thread);
        BamCompleteBlock completeBlock(200);

        printf("Malloc Memory is Over\n");
        /*
        * 分析准备
        */
        thread *read_thread = new thread(&read_pack,sin->fp.bgzf,&read);
        thread **compress_thread = new thread *[n_thread];


        for (int i=0;i<n_thread;i++){
            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(i,&cpuset);
            compress_thread[i]=new thread(&compress_pack,&read,&compress);
            int rc =pthread_setaffinity_np(compress_thread[i]->native_handle(),sizeof(cpu_set_t), &cpuset);
        }
        thread *assign_thread = new thread(&assign_pack,&compress,&completeBlock);
        //thread *consumer_thread = new thread(&benchmark_pack,&completeBlock);
        int  consumer_thread_number = 1;
        thread **consumer_thread = new thread*[consumer_thread_number];
        for (int i=0;i<consumer_thread_number;i++) consumer_thread[i] = new thread(&benchmark_write_pack,&completeBlock,output,hdr,level);

        read_thread->join();
        for (int i=0;i<n_thread;i++) compress_thread[i]->join();
        assign_thread->join();
        //consumer_thread->join();
        for (int i=0;i<consumer_thread_number;i++) consumer_thread[i]->join();
        sam_close(sin);
        sam_close(output);
        printf("Wait num is %d\n",compress.wait_num);
        TEND(fq)
        TPRINT(fq,"time is : ");
    }
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "write_mul_test")==0){
        TDEF(fq)
        TSTART(fq)
        printf("Starting Running Write Test\n");
        printf("BGZF_MAX_BLOCK_COMPLETE_SIZE is %d\n",BGZF_MAX_BLOCK_COMPLETE_SIZE);
        printf("output File Name is %s\n",outputfile.c_str());
        //if (strcmp(outputfile.substr(outputfile.size()-4).c_str(),"html")==0) outputfile=("./output.fastq");
        samFile *sin;
        sam_hdr_t *hdr;
        samFile *output;
        if ((sin=sam_open(inputfile.c_str(),"r"))==NULL){
            printf("Can`t open this file!\n");
            return 0;
        }

        if ((output=sam_open(outputfile.c_str(),"wb"))==NULL){
            printf("Can`t open this file!\n");
            return 0;
        }
        if ((hdr = sam_hdr_read(sin)) == NULL) {
            return  0;
        }
        /*
        *  读取和处理准备
        */
        BamRead read(8000);
        BamCompress compress(4000,n_thread);
        BamCompleteBlock completeBlock(200);

        printf("Malloc Memory is Over\n");
        /*
        * 分析准备
        */
        thread *read_thread = new thread(&read_pack,sin->fp.bgzf,&read);
        thread **compress_thread = new thread *[n_thread];


        for (int i=0;i<n_thread;i++){
            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(i,&cpuset);
            compress_thread[i]=new thread(&compress_pack,&read,&compress);
            int rc =pthread_setaffinity_np(compress_thread[i]->native_handle(),sizeof(cpu_set_t), &cpuset);
        }
        thread *assign_thread = new thread(&assign_pack,&compress,&completeBlock);
        //thread *consumer_thread = new thread(&benchmark_pack,&completeBlock);
        /*
         * 压缩准备
         */
        BamWriteCompress *bam_write_compress = new BamWriteCompress(4000,n_thread_write);


        int  consumer_thread_number = 1;
        thread **consumer_thread = new thread*[consumer_thread_number];
        for (int i=0;i<consumer_thread_number;i++) consumer_thread[i] = new thread(&benchmark_write_mul_pack,&completeBlock,bam_write_compress,output,hdr,level);

        thread **write_compress_thread = new thread*[n_thread_write];
        for (int i=0;i<n_thread_write;i++) write_compress_thread[i] = new thread(&bam_write_compress_pack,output->fp.bgzf,bam_write_compress);

        thread *write_output_thread = new thread(&bam_write_pack,output->fp.bgzf,bam_write_compress);

        read_thread->join();
        for (int i=0;i<n_thread;i++) compress_thread[i]->join();
        assign_thread->join();

        printf("Read Thread Has Been Over!!!\n");
        for (int i=0;i<consumer_thread_number;i++) consumer_thread[i]->join();
        printf("Consumer Thread Over!!!\n");
        for (int i=0;i<n_thread_write;i++) write_compress_thread[i]->join();
        printf("Write Compress Thread Over!!!\n");
        write_output_thread->join();
        printf("Write Output Thread Over!!!\n");

        sam_close(sin);
        sam_close(output);
        printf("Wait num is %d\n",compress.wait_num);
        TEND(fq)
        TPRINT(fq,"time is : ");
    }
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "api_test")==0){
        TDEF(fq)
        TSTART(fq)

        printf("Starting Running API Test\n");
        printf("BGZF_MAX_BLOCK_COMPLETE_SIZE is %d\n",BGZF_MAX_BLOCK_COMPLETE_SIZE);
        printf("output File Name is %s\n",outputfile.c_str());
        //if (strcmp(outputfile.substr(outputfile.size()-4).c_str(),"html")==0) outputfile=("./output.fastq");
//        samFile *sin;
//        sam_hdr_t *hdr;
//        samFile *output;
//        if ((sin=sam_open(inputfile.c_str(),"r"))==NULL){
//            printf("Can`t open this file!\n");
//            return 0;
//        }
//
//        if ((output=sam_open(outputfile.c_str(),"wb"))==NULL){
//            printf("Can`t open this file!\n");
//            return 0;
//        }
//        if ((hdr = sam_hdr_read(sin)) == NULL) {
//            return  0;
//        }


        /*
         * 开始创建BamRead 和 BamWriter
         */
#define init_test 
//#define debug_test 
//#define queue_test
//#define parallel_test
//#define merge_parallel_test
//tagg

#ifdef queue_test

        long long num=0;
        int big_size = 256;
        //int big_size = 100 << 10;
        //int big_size = 800 << 10;
        double t0 = GetTime();
		BamReader *reader = new BamReader(inputfile, n_thread);
        printf("new reader cost %lf\n", GetTime() - t0);
        
        t0 = GetTime();
		std::thread reader_t(read_thread_task, reader);
		reader_t.join();
        printf("read cost %lf\n", GetTime() - t0);

        printf("============\n");

        t0 = GetTime();
        BamWriter *writer = new BamWriter(outputfile, reader->getHeader(), n_thread_write, level, big_size);
        printf("new writer cost %lf\n", GetTime() - t0);


        
        t0 = GetTime();
		std::thread writer_t(write_thread_task, writer);
		writer_t.join();
        printf("write cost %lf\n", GetTime() - t0);

#endif

#ifdef debug_test

        double t0 = GetTime();
        BamReader *reader = new BamReader(inputfile,n_thread);
        printf("new reader cost %lf\n", GetTime() - t0);

        const int NN = 1000;

        t0 = GetTime();
        bam1_t **b = new bam1_t*[NN];
        for(int i = 0; i < NN; i++) {
            if ((b[i] = bam_init1()) == NULL) {
                fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
            }
        }
        printf("init items cost %lf\n", GetTime() - t0);

        t0 = GetTime();
        long long num=0;
        while (reader->getBam1_t(b[num % NN])){
            num++;
        }
        printf("read items cost %lf\n", GetTime() - t0);


        t0 = GetTime();
        //int big_size = 4 << 20;
        int big_size = 8 << 10;
        BamWriter *writer = new BamWriter(outputfile,reader->getHeader(),n_thread_write,level, big_size);
        printf("new writer cost %lf\n", GetTime() - t0);

        //sleep(5);
        double t1 = GetTime();
        t0 = GetTime();
        for(long long i = 0; i < num; i++) {
            writer->write(b[i % NN]);
        }
        printf("writer cost %lf\n", GetTime() - t0);
        t0 = GetTime();
        writer->over();
        printf("laste consumer cost %lf\n", GetTime() - t0);
        printf("consumer cost %lf\n", GetTime() - t1);

#endif


#ifdef debug_read

        double t0 = GetTime();
        BamReader *reader = new BamReader(inputfile,n_thread);
        printf("new reader cost %lf\n", GetTime() - t0);

        const int NN = 1000;

        t0 = GetTime();
        bam1_t **b = new bam1_t*[NN];
        for(int i = 0; i < NN; i++) {
            if ((b[i] = bam_init1()) == NULL) {
                fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
            }
        }
        printf("init items cost %lf\n", GetTime() - t0);

        t0 = GetTime();
        //int big_size = 4 << 20;
        //int big_size = 256;
        //BamWriter *writer = new BamWriter(outputfile,reader->getHeader(),n_thread_write,level, big_size);
        //printf("new writer cost %lf\n", GetTime() - t0);

        double t1 = GetTime();
        t0 = GetTime();
        long long num=0;
        while (reader->getBam1_t(b[num % NN])){
            num++;
            //writer->write(b[num % NN]);
        }
        printf("new producer cost %lf\n", GetTime() - t0);
        t0 = GetTime();
        //writer->over();
        printf("new laste consumer cost %lf\n", GetTime() - t0);
        printf("new consumer cost %lf\n", GetTime() - t1);


#endif


#ifdef debug_test_all

        double t0 = GetTime();
        BamReader *reader = new BamReader(inputfile,n_thread);
        printf("new reader cost %lf\n", GetTime() - t0);

        const int NN = 1000;

        t0 = GetTime();
        bam1_t **b = new bam1_t*[NN];
        for(int i = 0; i < NN; i++) {
            if ((b[i] = bam_init1()) == NULL) {
                fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
            }
        }
        printf("init items cost %lf\n", GetTime() - t0);

        t0 = GetTime();
        //int big_size = 4 << 20;
        int big_size = 256;
        BamWriter *writer = new BamWriter(outputfile,reader->getHeader(),n_thread_write,level, big_size);
        printf("new writer cost %lf\n", GetTime() - t0);

        double t1 = GetTime();
        t0 = GetTime();
        long long num=0;
        double t_read_bam = 0;
        double t_write_bam = 0;
        double t2;
        while (true) {
            t2 = GetTime();
            int res = reader->getBam1_t(b[num % NN]);
            t_read_bam += GetTime() - t2;

            if(res == 0) break;
            num++;
            t2 = GetTime();
            writer->write(b[num % NN]);
            t_write_bam += GetTime() - t2;
        }
        printf("read: %lf; write: %lf\n", t_read_bam, t_write_bam);
        printf("total producer cost %lf\n", GetTime() - t0);
        t0 = GetTime();
        writer->over();
        printf("total laste consumer cost %lf\n", GetTime() - t0);
        printf("total consumer cost %lf\n", GetTime() - t1);


#endif

#ifdef init_test

        //int big_size = 800 << 10;
        int big_size = 256;
        //int big_size = 32 << 10;
        double t0 = GetTime();
		BamReader *reader = new BamReader(inputfile, n_thread);
        printf("new reader cost %lf\n", GetTime() - t0);

        t0 = GetTime();
        BamWriter *writer = new BamWriter(outputfile, reader->getHeader(), n_thread_write, level, big_size);
        printf("new writer cost %lf\n", GetTime() - t0);

        bam1_t *b;
        if ((b = bam_init1()) == NULL) {
            fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
        }
        long long num = 0;
        t0 = GetTime();
        double t1 = GetTime();
        while (reader->getBam1_t(b)) {
            num++;
            writer->write(b);
//            if (num%1000 == 0) printf("Bam1_t Num is %d\n",num);
        }
        printf("new process 1+2 cost %lf\n", GetTime() - t0);

        writer->over();
        printf("new total process cost %lf\n", GetTime() - t1);
#endif


#ifdef merge_parallel_test 

        //int big_size = 800 << 10;
        int big_size = 256;
        //int big_size = 32 << 10;
        double write_time = 0;

        double t0, t1;

        t0 = GetTime();
		BamReader *reader = new BamReader(inputfile, n_thread);
        printf("new reader cost %lf\n", GetTime() - t0);

        t0 = GetTime();
        BamWriter *writer = new BamWriter(outputfile, reader->getHeader(), n_thread_write, level, big_size);
        printf("new writer cost %lf\n", GetTime() - t0);

        const int vec_N = 32 << 10;

        t0 = GetTime();
        std::vector<bam1_t *> b_vec(vec_N);
        for(int i = 0; i < b_vec.size(); i++) {
            if ((b_vec[i] = bam_init1()) == NULL) {
                fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
            }
        }
        printf("init bams cost %lf\n", GetTime() - t0);

        long long num = 0;

        t0 = GetTime();
        while (reader->getBam1_t(b_vec[num % vec_N])) {
            num++;
        }
        printf("tot read cost %lf\n", GetTime() - t0);

        t0 = GetTime();
        int vec_size = 0;
        for(int i = 0; i < num; i++) {
            vec_size++;
            if(vec_size == vec_N) {
                t1 = GetTime();
                writer->write_parallel(b_vec);
                vec_size = 0;
                write_time += GetTime() - t1;
            }
        }

        t1 = GetTime();
        if(vec_size) {
            b_vec.resize(vec_size);
            writer->write_parallel(b_vec);
        }
        write_time += GetTime() - t1;

        printf("new process cost %lf\n", GetTime() - t0);

        t0 = GetTime();
        writer->over();
        printf("new last process cost %lf\n", GetTime() - t0);
        printf("tot write cost %lf\n", write_time);
#endif


#ifdef parallel_test 

        //int big_size = 800 << 10;
        int big_size = 4 << 10;
        //int big_size = 32 << 10;
        double write_time = 0;
        double read_time = 0;

        double t0, t1;

        t0 = GetTime();
		//BamReader *reader = new BamReader(inputfile, 800 << 10, 800 << 10, 800 << 10, n_thread);
		BamReader *reader = new BamReader(inputfile, n_thread);
        printf("new reader cost %lf\n", GetTime() - t0);

        t0 = GetTime();
        BamWriter *writer = new BamWriter(outputfile, reader->getHeader(), n_thread_write, level, big_size);
        printf("new writer cost %lf\n", GetTime() - t0);

        const int vec_N = 40 << 10;

        t0 = GetTime();
        std::vector<bam1_t*> b_vec[THREAD_NUM_P];
//#pragma omp parallel for num_threads(THREAD_NUM_P) schedule(static)
        for(int i = 0; i < THREAD_NUM_P; i++) {
            for(int j = 0; j < vec_N; j++) {
                bam1_t* item = bam_init1();
                if(item == NULL) {
                    fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
                }
                b_vec[i].push_back(item);
            }
        }
        printf("init bams cost %lf\n", GetTime() - t0);

        long long num = 0;

        t0 = GetTime();

        //bam1_t* AAA = bam_init1();
        
        bool done = 0;
        while (done == 0) {
            t1 = GetTime();
            auto res_vec = reader->getBam1_t_parallel(b_vec);
            int res_vec_size = res_vec.size();
            //printf("reader get vec size %d\n", res_vec_size);
            num += res_vec_size;
            //int r_res = reader->getBam1_t(AAA);
            read_time += GetTime() - t1;

            if(res_vec_size == 0) break;
            //if(r_res == 0) break;

            t1 = GetTime();
#ifdef use_parallel_write
            writer->write_parallel(res_vec);
#else
            for(int i = 0; i < res_vec_size; i++) {
                writer->write(res_vec[i]);
            }
#endif
            write_time += GetTime() - t1;
            
        }

        //for(auto item : mp) {
        //    printf(" == %d %lld\n", item.first, item.second);
        //}

        printf("tot cost %lf\n", GetTime() - t0);
        printf("tot read cost %lf\n", read_time);
        printf("tot write cost %lf\n", write_time);

        t0 = GetTime();
#ifdef use_parallel_write
        writer->over_parallel();
#else
        writer->over();
#endif
        printf("new last process cost %lf\n", GetTime() - t0);
#endif

        cout << "Bam1_t Num : "<< num << endl;


        TEND(fq)
        TPRINT(fq,"time is : ");
        /*
         *  二代数据线程读写比例为 1 ：4
         *  三代数据线程读写比例为 1 ：4
         *
         *  测试下来读写速率为700mb/s 在64线程下
         *
         *  注意需要有额外的四个线程
         */
    }
}



