//
// Created by 赵展 on 2022/8/17.
//

#include "BamReader.h"



/*
 *  From block_mul.cpp
 *
 */
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
void compress_test_pack(BamCompress* compress) {
    bam_block *un_comp = nullptr;
    while (1) {
        un_comp = compress->getUnCompressData();
        if (un_comp == nullptr) {
            break;
        }
        compress->backEmpty(un_comp);

    }
}


BamReader::BamReader(std::string file_name,int n_thread){

   /*
    * 准备Sin等
    */
    if ((sin=sam_open(file_name.c_str(),"r"))==NULL){
        printf("Can`t open this file!\n");

    }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
    }






    //TODO 可以考虑获取用户的内存大小，根据此来获得Reader等的大小
    read=new BamRead(1024);
    //TODO 默认8线程
    this->n_thread = n_thread;
    compress=new BamCompress(1024,n_thread);
    //TODO
    completeBlock=new BamCompleteBlock(256);


    /*
     *
     * 开始创建线程，并且开始读取
     *
     */
    read_thread = new thread(&read_pack,sin->fp.bgzf,read);
    compress_thread = new thread *[n_thread];

    for (int i=0;i<n_thread;i++){
//        cpu_set_t cpuset;
//        CPU_ZERO(&cpuset);
//        CPU_SET(i,&cpuset);
        compress_thread[i]=new thread(&compress_pack,read,compress);
//        int rc =pthread_setaffinity_np(compress_thread[i]->native_handle(),sizeof(cpu_set_t), &cpuset);

    }
    thread *assign_thread = new thread(&assign_pack,compress,completeBlock);
//        thread *consumer_thread = new thread(&benchmark_pack,&completeBlock);


    un_comp = completeBlock->getCompleteBlock();


}


BamReader::BamReader(std::string file_name,int read_block,int compress_block,int compress_complete_block,int n_thread){
    /*
     * 准备Sin等
     */
    if ((sin=sam_open(file_name.c_str(),"r"))==NULL){
        printf("Can`t open this file!\n");

    }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
    }






    //TODO 可以考虑获取用户的内存大小，根据此来获得Reader等的大小
    read=new BamRead(read_block);
    //TODO 默认8线程
    this->n_thread = n_thread;
    compress=new BamCompress(compress_block,n_thread);
    //TODO
    completeBlock=new BamCompleteBlock(compress_block);


    /*
     *
     * 开始创建线程，并且开始读取
     *
     */
    read_thread = new thread(&read_pack,sin->fp.bgzf,read);
    compress_thread = new thread *[n_thread];

    for (int i=0;i<n_thread;i++){
//        cpu_set_t cpuset;
//        CPU_ZERO(&cpuset);
//        CPU_SET(i,&cpuset);
        compress_thread[i]=new thread(&compress_pack,read,compress);
//        int rc =pthread_setaffinity_np(compress_thread[i]->native_handle(),sizeof(cpu_set_t), &cpuset);

    }
    thread *assign_thread = new thread(&assign_pack,compress,completeBlock);
//        thread *consumer_thread = new thread(&benchmark_pack,&completeBlock);


    un_comp = completeBlock->getCompleteBlock();
}



/*
*  读取 SAMHDR
*/
sam_hdr_t* BamReader::getHeader(){
    return hdr;
}

/*
 *  获取 BAM1t
 *
 *  根据输入进来的bam1_t的指针，装入bam1_t
 *
 *  无返回值
 *
 *
 */

bool BamReader::getBam1_t(bam1_t* b){
        int ret;
        while (un_comp!=nullptr){
            if ((ret=(read_bam(un_comp,b,0)))>=0) {
                return true;
            }else{
                completeBlock->backEmpty(un_comp);
                un_comp = completeBlock->getCompleteBlock();
            }
        }
        return false;
}


/*
 * 获取Bam——Complete——Clock
 *
 * 针对想要高性能的开发者，返回未被分割的bam1_t，但是保证内部存在完整的bam1_t
 */
bam_complete_block* BamReader::getBamCompleteClock(){
    return completeBlock->getCompleteBlock();
}

void BamReader::backBamCompleteBlock(bam_complete_block *un_comp){
    completeBlock->backEmpty(un_comp);
}
