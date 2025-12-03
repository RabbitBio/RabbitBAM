
#include "swbam.h"

using namespace std;

extern "C" {
#include <athread.h>
#include <pthread.h>
    void slave_decompressfunc();
    void copyfunc();
    void slave_compressfunc();
}

SwBam::SwBam(CmdInfo *cmd_info1) {
    printf("SwBam::SwBam(CmdInfo *cmd_info1)\n");
    cmd_info_ = cmd_info1;

}

SwBam::~SwBam() {
    printf("SwBam::~SwBam()\n");
   

}

void SwBam::ProducerSwBamTask(BGZF *fp, BamRead *read) {
    printf("void SwBam::ProducerSwBamTask()\n");

    std::vector<bam_block*> tmp_chunks;
    bam_block* block = read->getEmpty();
    int global_block_id = 0;
    int block_num = 0;
    int group_num = 0;
    int ret=-1;

    while (true) {
        ret = read_block(fp, block);
        if (ret < 0) break; //读取错误或结束

        //读到了一块
        block->block_id = global_block_id++;
        block->pos = 0;
        block_num++;
        tmp_chunks.push_back(block);
        block = read->getEmpty();

        if (tmp_chunks.size() == 64) {
            read->inputBlock64(tmp_chunks);
            group_num++;
            tmp_chunks.clear();
        }
    }

    //文件读取结束
    if(tmp_chunks.size()) {
        read->inputBlock64(tmp_chunks);
        group_num++;
        tmp_chunks.clear();
    }

    //标记读取处理结束
    read->markComplete();  

}

void SwBam::ConsumerSwBamTask(BamRead *read, BamComplete *complete) {
    printf("void SwBam::ConsumerSwBamTask()\n");

    athread_init();
 
    while (true) {

        auto tmp_chunks = read->getBlock64();

        //任务处理结束
        if (tmp_chunks.empty()) break; 

        //不足64块的补齐
        for(int i = tmp_chunks.size(); i < 64; i++) {
            tmp_chunks.push_back(NULL);
        }

        //处理从核所需要的参数
        Para degz_paras[64];

        for (int i = 0; i < 64; i++) {
            if(tmp_chunks[i] == NULL) {
                degz_paras[i].input_block = NULL ;
                degz_paras[i].decompress_size = 0;
                degz_paras[i].n_records = 0;
                degz_paras[i].status = -1; //标记为空块
                continue;
            }
            degz_paras[i].block_id = i;
            degz_paras[i].input_block = tmp_chunks[i];
            degz_paras[i].un_comp_block = complete->getBuffer(i);
            degz_paras[i].output_records = complete->getResultBuf(i); 
            degz_paras[i].status = 0;
        }

        //调用从核处理一块：解压缩，解析
        {
            __real_athread_spawn((void *) slave_decompressfunc, degz_paras, 1);
            athread_join();
        }

        //为每条记录分配data区域
        for(int i = 0; i < 64; i++) {
            if (degz_paras[i].status == 0 && degz_paras[i].input_block != NULL) {
                for (int j = 0; j < degz_paras[i].n_records; j++) {
                   int new_l_data = degz_paras[i].l_data_list[j];
                   degz_paras[i].output_records[j]->l_data = 0;
                   realloc_bam_data(degz_paras[i].output_records[j], new_l_data);
                   degz_paras[i].output_records[j]->l_data = new_l_data;
                }
            }
        }

        //调用从核复制数据到新分配的data区域
        {
            __real_athread_spawn((void *) slave_copyfunc, degz_paras, 1);
            athread_join();
        }

        //处理完的结果放入到complete中
        for (int i = 0; i < 64; i++) {
            if (degz_paras[i].status == 0 && degz_paras[i].input_block != NULL) {
                // 依次收集每个块的bam1_t结果
                for (int j = 0; j < degz_paras[i].n_records; j++) {
                   complete->pushBlockResults(degz_paras[i].output_records,degz_paras[i].n_records);
                }
            }
            read->backBlock(tmp_chunks[i]); // 回收块
        }
    }

    complete->markComplete(); // 标记全部处理完成

}


int writeBam1_tToSam(samFile *fp, const sam_hdr_t *h, const bam1_t *b) {
    if (sam_write1_sw(fp, h, b) < 0) {
                fprintf(stderr, "Error writing SAM record\n");
                return -1;
    }
    return 0;
}


//---------------------------------------------------------------------------------------------------------------------

void SwBam:: ProducerSwBamTask2(samFile *fp, BamWrite *write){
    printf("void SwBam::ProducerSwBamTask2()\n");

    // 当前块和当前组的缓存
    std::vector<bam1_t*> cur_block;
    std::vector<std::vector<bam1_t*>> cur_group;

    bam1_t *b = write->getEmpty();
    int ret;
    bam1_core_t *c;
    uint32_t bam_len , total_len = 0;

    while (true) {
        // 读取一条记录，放入块中
        ret = sam_read1_sw(fp, NULL, b);
        if (ret < 0) break;  // 文件结束

        c = &b->core;
        bam_len = b->l_data - c->l_extranul + 32;
        if( bam_len + 4 + total_len <= BGZF_BLOCK_SIZE ) { //压缩前的最大大小
            // 放入当前记录到块中
            cur_block.push_back(b);
            total_len = total_len + bam_len + 4;
        }else{
            // 当前块已满，打包入组，开启新块
            cur_group.push_back(cur_block);
            cur_block.clear();
            cur_block.push_back(b);
            total_len = bam_len + 4;
        }

        b = write->getEmpty();

        // 达到64个块后打包压入队列
        if ((int)cur_group.size() >= 64) {
            write->inputGroup(cur_group);
            cur_group.clear();
        }
    }

    // 收尾：处理未满的块或组
    if (!cur_block.empty()) {
        cur_group.push_back(cur_block);
        cur_block.clear();
    }
    if (!cur_group.empty()) {
        write->inputGroup(cur_group);
        cur_group.clear();
    }

    write->markComplete();

}


void SwBam:: ConsumerSwBamTask2 (BamWrite *write, BamWriteComplete *complete){
    printf("void SwBam::ConsumerSwBamTask2()\n");

    athread_init();
 
    while (true) {

        auto group = write->getGroup();

        //任务处理结束
        if (group.empty()) break; 

        //不足64块的补齐
        for(int i = group.size(); i < 64; i++) {
            group.push_back(std::vector<bam1_t*>{});
        }

        //处理从核所需要的参数
        Comp_Para comp_paras[64];

        for (int i = 0; i < 64; i++) {
            if(group[i].empty()) {
                comp_paras[i].block_id = i;
                comp_paras[i].input_records = nullptr;
                comp_paras[i].n_records = 0;
                comp_paras[i].output_block = nullptr;
                comp_paras[i].un_comp_block = nullptr;
                comp_paras[i].output_size = 0;
                comp_paras[i].status = -1; // empty
                continue;
            }
            comp_paras[i].block_id = i;
            comp_paras[i].input_records = group[i].data();   
            comp_paras[i].n_records = group[i].size();      
            comp_paras[i].output_block = complete->getEmpty();
            comp_paras[i].un_comp_block = complete->getBuffer(i);
            comp_paras[i].output_size = 0;
            comp_paras[i].status = 0;
        }

        //调用从核处理一块：解析，压缩成bgzf块
        {
            __real_athread_spawn((void*)slave_compressfunc, comp_paras, 1);
            athread_join();
        }

        for (int i = 0; i < 64; i++) {
            //处理完的结果放入到complete中
            if (comp_paras[i].status == 0 && comp_paras[i].output_block != nullptr) {
                complete->pushCompressedBlock(comp_paras[i].output_block);
            }

            // 回收输入的bam1_t
            if (comp_paras[i].input_records && comp_paras[i].n_records > 0) {
                for (int r = 0; r < comp_paras[i].n_records; ++r) {
                    bam1_t *rec = comp_paras[i].input_records[r];
                    write->backBam(rec);
                }
            }

        }
    }
    complete->markComplete(); // 标记全部处理完成
}

//一块bgzf写入文件中
int SwBam:: writeBlockTobam(BGZF *fp, bam_block *block) {
    int block_offset = block->pos;
    if (block->length < 0) {
        hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block->length, NULL));
        return -1;
    }
    if (hwrite(fp->fp, block->data, block->length) != block->length) {
        hts_log_error("File write failed (wrong size)");
        fp->errcode |= BGZF_ERR_IO; // possibly truncated file
        return -1;
    }
    block->pos = 0;
    fp->block_address += block->length;    
    return 0;
}


void SwBam::ProcessSwBam() {
    printf("void SwBam::ProcessSwBam()\n");

    //in文件打开
    sin = sam_open(cmd_info_->in_file_name_.c_str(), "r");
    if (sin == NULL) {
        fprintf(stderr, "Error opening input file %s\n", cmd_info_->in_file_name_.c_str());
    }

    //out文件打开
    sout = sam_open(cmd_info_->out_file_name_.c_str(), "wb");
    if (sout == NULL) {
        fprintf(stderr, "Error opening output file %s\n", cmd_info_->out_file_name_.c_str());
    }

    //头部读取
    hdr = sam_hdr_read(sin);
    if (hdr == NULL) {
        fprintf(stderr, "Error reading header from input file %s\n", cmd_info_->in_file_name_.c_str());
    }

    //头部写入
    if (sam_hdr_write(sout, hdr) != 0) {
        fprintf(stderr, "Error writing header to output file %s\n", cmd_info_->out_file_name_.c_str());
    }

    switch (sin->format.format) {
        case bam:{
            //开辟多线程处理 
            //进行队列等的初始化
            read = new BamRead(50);
            complete = new BamComplete(50);
            //生产者线程---
            thread producer(bind(&SwBam::ProducerSwBamTask, this, sin->fp.bgzf, read));

            //消费者线程----
            thread consumer(bind(&SwBam::ConsumerSwBamTask, this, read , complete));
            
            //写入线程
            //thread writer(bind(&SwBam::WriteSwBamTask, this, ));

            //剩下是主线程
            long long num = 0;
            bam1_t *b;
            while (true) {
                b = complete->popRecord();
                if(b == NULL) break;
                num++;
                int ret = writeBam1_tToSam(sout,hdr,b); 
                complete->backRecord(b);
            }

            printf("num %lld\n", num);

            break;
        }
            

        case sam:{
            write = new BamWrite(50);
            writeComplete = new BamWriteComplete(50);
            //生产者线程---
            thread producer2(bind(&SwBam::ProducerSwBamTask2, this, sin, write));
            //消费者线程----
            thread consumer2(bind(&SwBam::ConsumerSwBamTask2, this, write , writeComplete));

            //剩下是主线程
            long long num2 = 0;
            bam_block* comp_block;
            while (true) {
                comp_block = writeComplete->getCompressedBlock();
                if(comp_block == nullptr) break;
                num2++;
                int write_ret = writeBlockTobam(sout->fp.bgzf, comp_block);
                writeComplete->backBlock(comp_block);
            }

            printf("num2 %lld\n", num2);

            break;
        }
            

        default:
            fprintf(stderr, "Unknown file format\n");
            break;
    }


    //释放内存
    sam_hdr_destroy(hdr);
    
}


