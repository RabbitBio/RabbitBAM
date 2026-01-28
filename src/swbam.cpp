#include "swbam.h"

using namespace std;

extern "C" {
#include <athread.h>
#include <pthread.h>
    void slave_decompressfunc();
    void slave_copyfunc();
    void slave_compressfunc();
    void slave_sam_format();
    void slave_sam_parse();
}

// Global mutex to coordinate reader and writer spawn
std::mutex g_athread_spawn_mutex;

double temp1 = 0 , temp2 = 0 ,temp3 = 0 , temp4 = 0;
double t_sam2bam_write = 0 , t_sam2bam_read = 0 , t_sam2bam_slave = 0 , t_sam_parse = 0;
double t_bam2sam_write = 0 , t_bam2sam_read = 0 , t_bam2sam_slave = 0;

SwBam::SwBam(CmdInfo *cmd_info1) {
    cmd_info_ = cmd_info1;

}

SwBam::~SwBam() {

}

void SwBam::ProducerSwBamTask(BGZF *fp, BamRead *read) {
    printf("ProducerSwBamTask started\n");
    double t0 = GetTime();

    std::vector<bam_block*> tmp_chunks;
    bam_block* block = read->getEmpty();
    int global_block_id = 0;
    int block_num = 0;
    int group_num = 0;
    int ret=-1;

    while (true) {

        // temp1 = GetTime();
        ret = read_block(fp, block);
        // t_bam2sam_read += GetTime() - temp1;
        
        if (ret < 0) break; //读取错误或结束

        if(block->length == 28){
            break;
        }

        //读到了一块非EOF块
        block->block_id = ++global_block_id;
        block->pos = 0;
        block_num++;
        //printf("block_id=%d, block_size=%d\n",global_block_id,block->length);
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

    printf("ProducerSwBamTask finished with block_num = %d , group_num = %d , cost %lf\n", block_num , group_num , GetTime() - t0 );

    //标记读取处理结束
    read->markComplete();  

}

void SwBam::ConsumerSwBamTask(BamRead *read, BamComplete *complete) {
    printf("ConsumerSwBamTask started\n");
    double t0 = GetTime();
 
    while (true) {
        //printf("INTO while!\n");

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
        // temp2 = GetTime();
        {
            std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
            __real_athread_spawn((void *) slave_decompressfunc, degz_paras, 1);
            athread_join();
        }
        // t_bam2sam_slave += GetTime() - temp2;


        //为每条记录分配data区域
        // for(int i = 0; i < 64; i++) {
        //     if (degz_paras[i].status == 0 && degz_paras[i].input_block != NULL) {
        //         for (int j = 0; j < degz_paras[i].n_records; j++) {
        //            int new_l_data = degz_paras[i].l_data_list[j];
        //            degz_paras[i].output_records[j]->l_data = 0;
        //            realloc_bam_data(degz_paras[i].output_records[j], new_l_data);
        //            degz_paras[i].output_records[j]->l_data = new_l_data;
        //         }
        //     }
        // }

        //调用从核复制数据到新分配的data区域
        // {
        //     __real_athread_spawn((void *) slave_copyfunc, degz_paras, 1);
        //     athread_join();
        // }

        //处理完的结果放入到complete中,这一部分是可以在从核中实现的，能优化一点性能！！！！！！！！！！！！
        //在这里打印一共有几块，每一块有几条！！
        for (int i = 0; i < 64; i++) {
            if (degz_paras[i].status == 0 && degz_paras[i].input_block != NULL) {
                // 依次收集每个块的bam1_t结果
                //printf("degz_paras[%d].n_records = %d\n",i , degz_paras[i].n_records);

                for (int j = 0; j < degz_paras[i].n_records; j++) {
                    bam1_t* bam1 = complete->getEmpty();
                    //这里需要复制一份，因为缓冲区还要使用！！！
                    (void)bam_copy1(bam1, degz_paras[i].output_records[j]);
                    complete->inputBam1_t(bam1);
                }

                read->backBlock(degz_paras[i].input_block); // 回收块
            }
        }
    }

    complete->markComplete(); // 标记全部处理完成

    printf("ConsumerSwBamTask finished , cost %lf!\n", GetTime() - t0); 
}

int SwBam::writeBam1_tToSam(samFile *fp, const sam_hdr_t *h, const bam1_t *b) {
    if (sam_write1_sw(fp, h, b) < 0) {
        fprintf(stderr, "Error writing SAM record\n");
        return -1;
    }
    return 0;
}

void init_batch_buffers(SamFormatBatch *batch)
{
    batch->count = 0;

    for (int i = 0; i < BATCH_SIZE; ++i) {
        kstring_t *ks = &batch->sam_lines[i];

        ks->l = 0;
        ks->m = MAX_SAM_LINE_SIZE;
        ks->s = (char *)aligned_alloc_custom(64, MAX_SAM_LINE_SIZE);

        if (!ks->s) {
            fprintf(stderr, "Failed to alloc sam buffer %d\n", i);
            abort();
        }
    }
}

void destroy_batch_buffers(SamFormatBatch *batch)
{
    for (int i = 0; i < BATCH_SIZE; ++i) {
        kstring_t *ks = &batch->sam_lines[i];
        if (ks->s) {
            aligned_free_custom((unsigned char*)ks->s);
            ks->s = NULL;
        }
        ks->l = ks->m = 0;
    }
}



//---------------------------------------------------------------------------------------------------------------------

struct MemReader {
    char *base;
    size_t size;
    size_t pos;
};

inline bool mem_getline(MemReader &r, kstring_t *ks)
{
    if (r.pos >= r.size)
        return false;

    size_t start = r.pos;
    while (r.pos < r.size && r.base[r.pos] != '\n')
        r.pos++;

    size_t len = r.pos - start;

    if (len + 2 > ks->m) {
        fprintf(stderr,
                "SAM line too long: %zu > %zu\n",
                len, ks->m);
        abort();
    }

    memcpy(ks->s, r.base + start, len);
    ks->s[len] = '\n';
    ks->s[len + 1] = 0;
    ks->l = len + 1;

    if (r.pos < r.size && r.base[r.pos] == '\n')
        r.pos++;

    return true;
}

void SwBam:: ProducerSwBamTask2_parallel_memory( BamWrite *write, sam_hdr_t *h ,size_t sam_size, char *sam_mem){
    printf("ProducerSwBamTask2_parallel_memory started\n");
    double t0 = GetTime();

    // 当前块和当前组的缓存
    std::vector<bam1_t*> cur_block;
    std::vector<std::vector<bam1_t*>> cur_group;

    bam1_t *b = write->getEmpty();
    int ret;
    bam1_core_t *c;
    uint32_t bam_len , total_len = 0;
    int block_nums=0 , group_nums=0 , bam1_t_nums=0;

    MemReader reader;
    reader.base = sam_mem;
    reader.size = sam_size;
    reader.pos  = 0;

    // SamFormatBatch batch;
    // init_batch_buffers(&batch);
    // batch.hdr = h;
    kstring_t ks_struct;
    kstring_t *ks = &ks_struct;
    ks->l = 0;
    ks->m = MAX_SAM_LINE_SIZE;
    ks->s = (char *)malloc(MAX_SAM_LINE_SIZE);

     while (true) {

        // 读取一条记录，放入块中
        temp1 = GetTime();
        bool ok = mem_getline(reader, ks);
        t_sam2bam_read += GetTime() - temp1;

       if (!ok) break;
        bam1_t_nums++;

        temp4 = GetTime();
        sam_format1(h, b, ks);
        kputc('\n', ks);
        t_sam_parse += GetTime() - temp4;


        c = &b->core;
        bam_len = b->l_data - c->l_extranul + 32;
        if( bam_len + 4 + total_len <= BGZF_BLOCK_SIZE ) { //压缩前的最大大小
            // 放入当前记录到块中
            cur_block.push_back(b);
            total_len = total_len + bam_len + 4;
        }else{
            // 当前块已满，打包入组，开启新块
            block_nums++;
            cur_group.push_back(cur_block);
            cur_block.clear();
            cur_block.push_back(b);
            total_len = bam_len + 4;
        }

        b = write->getEmpty();

        // 达到64个块后打包压入队列
        if ((int)cur_group.size() >= 64) {
            group_nums++;
            write->inputGroup(cur_group);
            cur_group.clear();
        }
    }

    // while (true) {
    //     batch.count = 0;

    //     //1. 主核读取 SAM 文本 ----------
    //     while (batch.count < BATCH_SIZE) {

    //         temp1 = GetTime();
    //         bool ok = mem_getline(reader, &batch.sam_lines[batch.count]);
    //         t_sam2bam_read += GetTime() - temp1;

    //         if (!ok) break;
    //         batch.bams[batch.count] = write->getEmpty();
    //         batch.count++;
    //     }
    //     if (batch.count == 0) break;

    //     //2. 从核并行 sam_parse1 ----------
    //     temp4 = GetTime();
    //     // {
    //     //     std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
    //     //     __real_athread_spawn((void *)slave_sam_parse, &batch, 1);
    //     //     athread_join();
    //     // }
    //     for (int i = 0; i < batch.count; i++) {
    //         kstring_t *ks = &batch.sam_lines[i];
    //         ks->l = 0;
    //         sam_format1(batch.hdr, batch.bams[i], ks);
    //         kputc('\n', ks);
    //     }
    //     t_sam_parse += GetTime() - temp4;

    //     //3. 主核对一个批次进行打包 ---------
    //     for (int i = 0; i < batch.count; i++) {
    //         bam1_t *b = batch.bams[i];
    //         bam1_t_nums++;

    //         c = &b->core;
    //         bam_len = b->l_data - c->l_extranul + 32;
    //         if( bam_len + 4 + total_len <= BGZF_BLOCK_SIZE ) { //压缩前的最大大小
    //             // 放入当前记录到块中
    //             cur_block.push_back(b);
    //             total_len = total_len + bam_len + 4;
    //         }else{
    //             // 当前块已满，打包入组，开启新块
    //             block_nums++;
    //             cur_group.push_back(cur_block);
    //             cur_block.clear();
    //             cur_block.push_back(b);
    //             total_len = bam_len + 4;
    //         }

    //         // 达到64个块后打包压入队列
    //         if ((int)cur_group.size() >= 64) {
    //             group_nums++;
    //             write->inputGroup(cur_group);
    //             cur_group.clear();
    //         }
    //     }
        
    // }

    // 收尾：处理未满的块或组
    if (!cur_block.empty()) {
        block_nums++;
        cur_group.push_back(cur_block);
        cur_block.clear();
    }
    if (!cur_group.empty()) {
        group_nums++;
        write->inputGroup(cur_group);
        cur_group.clear();
    }

    // destroy_batch_buffers(&batch);

    printf("ProducerSwBamTask2 finished with bam1_t_nums = %d , block_num = %d , group_num = %d , cost %lf\n", bam1_t_nums, block_nums , group_nums , GetTime() - t0 );

    write->markComplete();

}

void SwBam:: ProducerSwBamTask2_parallel(samFile *fp, BamWrite *write, sam_hdr_t *h){
    printf("ProducerSwBamTask2_parallel started\n");
    double t0 = GetTime();

    // 当前块和当前组的缓存
    std::vector<bam1_t*> cur_block;
    std::vector<std::vector<bam1_t*>> cur_group;

    bam1_t *b = write->getEmpty();
    int ret;
    bam1_core_t *c;
    uint32_t bam_len , total_len = 0;
    int block_nums=0 , group_nums=0 , bam1_t_nums=0;

    SamFormatBatch batch;
    init_batch_buffers(&batch);
    batch.hdr = h;

    while (true) {
        batch.count = 0;

        // 处理 header 预读行
        if (fp->line.l != 0) {
            printf("The header has a line!!!");
            batch.sam_lines[0] = fp->line;
            batch.bams[0] = write->getEmpty();
            batch.count = 1;
            fp->line.l = 0;
        }

        //1. 主核读取 SAM 文本 ----------
        while (batch.count < BATCH_SIZE) {

            // temp1 = GetTime();
            int ret = hts_getline(fp, KS_SEP_LINE,&batch.sam_lines[batch.count]);
            // t_sam2bam_read += GetTime() - temp1;

            if (ret < 0) break;
            batch.bams[batch.count] = write->getEmpty();
            batch.count++;
        }
        if (batch.count == 0) break;

        //2. 从核并行 sam_parse1 ----------
        // temp4 = GetTime();
        {
            std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
            __real_athread_spawn((void *)slave_sam_parse, &batch, 1);
            athread_join();
        }
        // t_sam_parse += GetTime() - temp4;

        //3. 主核对一个批次进行打包 ---------
        for (int i = 0; i < batch.count; i++) {
            bam1_t *b = batch.bams[i];
            bam1_t_nums++;

            c = &b->core;
            bam_len = b->l_data - c->l_extranul + 32;
            if( bam_len + 4 + total_len <= BGZF_BLOCK_SIZE ) { //压缩前的最大大小
                // 放入当前记录到块中
                cur_block.push_back(b);
                total_len = total_len + bam_len + 4;
            }else{
                // 当前块已满，打包入组，开启新块
                block_nums++;
                cur_group.push_back(cur_block);
                cur_block.clear();
                cur_block.push_back(b);
                total_len = bam_len + 4;
            }

            // 达到64个块后打包压入队列
            if ((int)cur_group.size() >= 64) {
                group_nums++;
                write->inputGroup(cur_group);
                cur_group.clear();
            }
        }
        // t_sam_parse += GetTime() - temp4;
        
    }

    // 收尾：处理未满的块或组
    if (!cur_block.empty()) {
        block_nums++;
        cur_group.push_back(cur_block);
        cur_block.clear();
    }
    if (!cur_group.empty()) {
        group_nums++;
        write->inputGroup(cur_group);
        cur_group.clear();
    }

    destroy_batch_buffers(&batch);

    printf("ProducerSwBamTask2 finished with bam1_t_nums = %d , block_num = %d , group_num = %d , cost %lf\n", bam1_t_nums, block_nums , group_nums , GetTime() - t0 );

    write->markComplete();

}

void SwBam:: ProducerSwBamTask2(samFile *fp, BamWrite *write, sam_hdr_t *h){
    printf("ProducerSwBamTask2 started\n");
    double t0 = GetTime();

    // 当前块和当前组的缓存
    std::vector<bam1_t*> cur_block;
    std::vector<std::vector<bam1_t*>> cur_group;

    bam1_t *b = write->getEmpty();
    int ret;
    bam1_core_t *c;
    uint32_t bam_len , total_len = 0;
    int block_nums=0 , group_nums=0 , bam1_t_nums=0;

    while (true) {

        // 读取一条记录，放入块中
        temp1 = GetTime();
        ret = sam_read1_sw(fp, h, b);
        t_sam2bam_read += GetTime() - temp1;

        if (ret < 0) break;  // 文件结束
        bam1_t_nums++;

        c = &b->core;
        bam_len = b->l_data - c->l_extranul + 32;
        if( bam_len + 4 + total_len <= BGZF_BLOCK_SIZE ) { //压缩前的最大大小
            // 放入当前记录到块中
            cur_block.push_back(b);
            total_len = total_len + bam_len + 4;
        }else{
            // 当前块已满，打包入组，开启新块
            block_nums++;
            cur_group.push_back(cur_block);
            cur_block.clear();
            cur_block.push_back(b);
            total_len = bam_len + 4;
        }

        b = write->getEmpty();

        // 达到64个块后打包压入队列
        if ((int)cur_group.size() >= 64) {
            group_nums++;
            write->inputGroup(cur_group);
            cur_group.clear();
        }
    }

    // 收尾：处理未满的块或组
    if (!cur_block.empty()) {
        block_nums++;
        cur_group.push_back(cur_block);
        cur_block.clear();
    }
    if (!cur_group.empty()) {
        group_nums++;
        write->inputGroup(cur_group);
        cur_group.clear();
    }

    printf("ProducerSwBamTask2 finished with bam1_t_nums = %d , block_num = %d , group_num = %d , cost %lf\n", bam1_t_nums, block_nums , group_nums , GetTime() - t0 );

    write->markComplete();

}

void SwBam:: ConsumerSwBamTask2 (BamWrite *write, BamWriteComplete *complete){
    printf("ConsumerSwBamTask2 started\n");
    double t0 = GetTime();

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
        // temp2 = GetTime();
        {
            std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
            __real_athread_spawn((void*)slave_compressfunc, comp_paras, 1);
            athread_join();
        }
        // t_sam2bam_slave += GetTime() - temp2;

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
    printf("ConsumerSwBamTask2 finished , cost %lf!\n", GetTime() - t0); 
}

//一块bgzf写入文件中
int SwBam:: writeBlockTobam(BGZF *fp, bam_block *block) {
    int block_offset = block->pos;
    if (block->length < 0) {
        hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block->length, NULL));
        return -1;
    }

    // temp3 = GetTime();
    if (hwrite(fp->fp, block->data, block->length) != block->length) {
        hts_log_error("File write failed (wrong size)");
        fp->errcode |= BGZF_ERR_IO; // possibly truncated file
        return -1;
    }
    // t_sam2bam_write += GetTime() - temp3;

    block->pos = 0;
    fp->block_address += block->length;    
    return 0;
}

void SwBam::ProcessSwBam() {

    //从内存中读取数据
    //#define USE_MEMORY
    #ifdef USE_MEMORY

    //in文件打开
    sin = sam_open(cmd_info_->in_file_name_.c_str(), "r");
    if (sin == NULL) {
        fprintf(stderr, "Error opening input file %s\n", cmd_info_->in_file_name_.c_str());
    }
    //out文件打开
    const char *out_mode = get_sam_open_mode(cmd_info_->out_file_name_);
    sout = sam_open(cmd_info_->out_file_name_.c_str(), out_mode);
    if (sout == NULL) {
        fprintf(stderr, "Error opening output file %s\n", cmd_info_->out_file_name_.c_str());
    }

    double t0 = GetTime();
    double tbody;
    size_t sam_size;
    char *sam_mem;

    switch (sin->format.format) {
        case bam:{
        
            break;
        }

        case sam:{

            FILE *f = fopen(cmd_info_->in_file_name_.c_str(), "r");
            if (!f) { perror("fopen"); exit(1); }

            // 获取文件大小
            fseek(f, 0, SEEK_END);
            sam_size = ftell(f);
            fseek(f, 0, SEEK_SET);

            // 分配内存
            sam_mem = (char*)malloc(sam_size);
            if (!sam_mem) { fprintf(stderr, "Failed to allocate memory for SAM file\n"); exit(1); }

            // 读取文件
            size_t read_bytes = fread(sam_mem, 1, sam_size, f);
            if (read_bytes != sam_size) {
                fprintf(stderr, "Failed to read SAM file\n");
                exit(1);
            }
            fclose(f);

            printf("Loaded SAM into memory: %.2f MB\n", sam_size / 1024.0 / 1024.0);
            break;
        }

        default:
            fprintf(stderr, "Unknown file format\n");
            break;
    }

    printf("Complete the memory cost %lf\n", GetTime() - t0);

    #else

    double t0 = GetTime();
    double tbody;

    //in文件打开
    sin = sam_open(cmd_info_->in_file_name_.c_str(), "r");
    if (sin == NULL) {
        fprintf(stderr, "Error opening input file %s\n", cmd_info_->in_file_name_.c_str());
    }

    //out文件打开
    const char *out_mode = get_sam_open_mode(cmd_info_->out_file_name_);
    sout = sam_open(cmd_info_->out_file_name_.c_str(), out_mode);
    if (sout == NULL) {
        fprintf(stderr, "Error opening output file %s\n", cmd_info_->out_file_name_.c_str());
    }
    printf("open the files cost %lf\n", GetTime() - t0);


    t0 = GetTime();
    //头部读取
    hdr = sam_hdr_read(sin);
    if (hdr == NULL) {
        fprintf(stderr, "Error reading header from input file %s\n", cmd_info_->in_file_name_.c_str());
    }

    //头部写入
    if (sam_hdr_write(sout, hdr) != 0) {
        fprintf(stderr, "Error writing header to output file %s\n", cmd_info_->out_file_name_.c_str());
    }
    printf("Complete the head cost %lf\n", GetTime() - t0);

    #endif

    tbody = GetTime();
    switch (sin->format.format) {
        case bam:{
            //这里初始化的性能非常非常差需要后续进行优化！！
            t0 = GetTime();
            //BamRead(x) x*64是bam_block的内存池大小
            //BamComplete(x) x是bam1_t的内存池大小

            //read = new BamRead(50);
            // read = new BamRead(260);
            // complete = new BamComplete(2932000);
            read = new BamRead(5);
            complete = new BamComplete(132000);

            printf("Complete the queue initialization cost %lf\n", GetTime() - t0);

            //生产者线程---
            thread producer(bind(&SwBam::ProducerSwBamTask, this, sin->fp.bgzf, read));
            // producer.join();

            //消费者线程----
            thread consumer(bind(&SwBam::ConsumerSwBamTask, this, read , complete));
            // consumer.join();

            //剩下是主线程
            t0 = GetTime();

            #define SAM_FORMAT_PARALLEL
            #ifdef SAM_FORMAT_PARALLEL
            printf("Enable the SAM_FORMAT_PARALLEL!!!\n");

            SamFormatBatch batch;
            init_batch_buffers(&batch);
            batch.hdr = hdr;

            long long num = 0;
            bam1_t *b;
            while (true) {
                batch.count = 0;

                // temp4 = GetTime();
                // 1. 收集 bam1_t
                while (batch.count < BATCH_SIZE) {
                    b = complete->getBam1_t();
                    if (!b) break;
                    num++;

                    batch.bams[batch.count] = b;
                    kstring_t *ks = &batch.sam_lines[batch.count];
                    ks->l = 0;

                    batch.count++;
                }
                if (batch.count == 0) break;

                // 2. 从核并行 sam_format1
                sout->format.category = sequence_data;
                sout->format.format = sam;

                // temp4 = GetTime();
                {
                    std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
                    __real_athread_spawn((void*)slave_sam_format, &batch, 1);
                    athread_join();
                }
                // t_sam_parse += GetTime() - temp4;

                // 3. 主核按顺序写入文件中
                for (int i = 0; i < batch.count; i++) {
                    kstring_t *ks = &batch.sam_lines[i];

                    // temp3 = GetTime();
                    if (hwrite(sout->fp.hfile, ks->s, ks->l) != ks->l) {
                        fprintf(stderr, "write failed\n");
                    }
                    // t_bam2sam_write += GetTime() - temp3;

                    complete->backBam1_t(batch.bams[i]);
                }
            }

            destroy_batch_buffers(&batch);
            #else
            long long num = 0;
            bam1_t *b;
            while (true) {
                b = complete->getBam1_t();
                if(b == NULL) break;
                
                //print_bam1(b);
                num++;

                temp3 = GetTime();
                int ret = writeBam1_tToSam(sout,hdr,b); 
                t_bam2sam_write += GetTime() - temp3;

                complete->backBam1_t(b);
            }
            #endif

            producer.join();
            consumer.join();
            // printf("The actual time of reading bam cost %lf\n", t_bam2sam_read);
            // printf("The actual time of slave cost %lf\n", t_bam2sam_slave);
            // printf("The actual time of sam parsing cost %lf\n", t_sam_parse);
            // printf("The actual time of writing to sam cost %lf\n", t_bam2sam_write);
            printf("Complete writing to sam cost %lf\n", GetTime() - t0);
            printf("The total bam1_t nums is %lld\n", num);
            break;
        }
            

        case sam:{
            t0 = GetTime();
            // const size_t INIT_DATA_SIZE = 1024;  
            // const size_t MAX_RECORDS_PER_BLOCK = 1024; 
            //BamWrite(x) x*MAX_RECORDS_PER_BLOCK是bam1_t的内存池大小
            //BamWriteComplete(x) x是bam_block的内存池大小

            write = new BamWrite(500);
            writeComplete = new BamWriteComplete(500);
            // write = new BamWrite(2870);
            // writeComplete = new BamWriteComplete(16400);
            // write = new BamWrite(128);
            // writeComplete = new BamWriteComplete(128);
            printf("Complete the queue initialization cost %lf\n", GetTime() - t0);

            //生产者线程---
            #define SAM_PARSE_PARALLEL
            #ifdef SAM_PARSE_PARALLEL
            #ifdef USE_MEMORY
            printf("Enable the SAM_PARSE_PARALLEL!!!\n");
            printf("Enable the USE_MEMORY!!!\n");
            thread producer2(bind(&SwBam::ProducerSwBamTask2_parallel_memory, this, write, hdr, sam_size  , sam_mem));
            #else
            printf("Enable the SAM_PARSE_PARALLEL!!!\n");
            thread producer2(bind(&SwBam::ProducerSwBamTask2_parallel, this, sin, write, hdr));
            #endif
            #else
            thread producer2(bind(&SwBam::ProducerSwBamTask2, this, sin, write, hdr));
            #endif
            // producer2.join();

            //消费者线程----
            thread consumer2(bind(&SwBam::ConsumerSwBamTask2, this, write , writeComplete));
            // consumer2.join();

            //剩下是主线程
            t0 = GetTime();
            long long num2 = 0;
            bam_block* comp_block;
            while (true) {
                comp_block = writeComplete->getCompressedBlock();
                if(comp_block == nullptr) break;
                num2++;
                //print_bam_block(comp_block);

                int write_ret = writeBlockTobam(sout->fp.bgzf, comp_block);
                writeComplete->backBlock(comp_block);
            }

            producer2.join();
            consumer2.join();
            // printf("The actual time of reading sam cost %lf\n", t_sam2bam_read);
            // printf("The actual time of sam parsing cost %lf\n", t_sam_parse);
            // printf("The actual time of slave cost %lf\n", t_sam2bam_slave);
            // printf("The actual time of writing to bam cost %lf\n", t_sam2bam_write);
            printf("Complete writing to bam cost %lf\n", GetTime() - t0);
            printf("The total BGZF nums is %lld\n", num2);
            break;
        }
            

        default:
            fprintf(stderr, "Unknown file format\n");
            break;
    }

    //释放内存
    sam_hdr_destroy(hdr);
    int ret;
    ret = hts_close(sout);
    if (ret < 0) {
        fprintf(stderr, "Error closing output.\n");
    }
    ret = hts_close(sin);
    if (ret < 0) {
        fprintf(stderr, "Error closing input.\n");
    }

    printf("Complete the body cost %lf\n", GetTime() - tbody);

    
}


