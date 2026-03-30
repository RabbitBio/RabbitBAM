#include "swbam.h"

using namespace std;

extern "C" {
#include <athread.h>
#include <pthread.h>
    void slave_decompressfunc();
    void slave_compressfunc();
    void slave_sam_format();
    void slave_sam_parse();
    void slave_sam_parse_chunk();
    void slave_copy_and_count(); 
}

// Global mutex to coordinate reader and writer spawn
std::mutex g_athread_spawn_mutex;

double temp1 = 0 , temp2 = 0 ,temp3 = 0 , temp4 = 0;
double t_sam2bam_write = 0 , t_sam2bam_read = 0 , t_sam2bam_slave = 0 , t_sam_parse = 0 , t_sam2bam_for = 0;
double t_bam2sam_write = 0 , t_bam2sam_read = 0 , t_bam2sam_slave = 0 , t_bam2sam_for = 0;

SwBam::SwBam(CmdInfo *cmd_info1) {
    cmd_info_ = cmd_info1;

}

SwBam::~SwBam() { }

void init_sam_format_batch_buffers(SamFormatBatch *batch)
{
    batch->count = 0;
    
    int batch_per_core = (BATCH_SIZE / 64) + 1; 
    size_t core_buf_size = MAX_SAM_LINE_SIZE * batch_per_core;

    for (int i = 0; i < 64; ++i) {
        kstring_t *ks = &batch->core_out_lines[i];

        ks->l = 0;
        ks->m = core_buf_size;
        ks->s = (char *)aligned_alloc_custom(64, core_buf_size);

        if (!ks->s) {
            fprintf(stderr, "Failed to alloc sam buffer for core %d\n", i);
            abort();
        }
    }
}

void destroy_sam_format_batch_buffers(SamFormatBatch *batch)
{
    for (int i = 0; i < 64; ++i) {
        kstring_t *ks = &batch->core_out_lines[i];
        if (ks->s) {
            aligned_free_custom((unsigned char*)ks->s);
            ks->s = NULL;
        }
        ks->l = ks->m = 0;
    }
}

void init_sam_parse_batch(SamParseBatch *batch) {
    for (int i = 0; i < 64; ++i) {
        batch->chunks[i].text_buf = (char *)aligned_alloc_custom(64, CHUNK_BUFFER_SIZE);
        if (!batch->chunks[i].text_buf) {
            fprintf(stderr, "Failed to allocate 5MB chunk for core %d\n", i);
            abort();
        }
        // bam_lens：64 字节对齐，供从核写入预计算的 bam_len，pack 阶段顺序读取
        batch->chunks[i].bam_lens = (uint32_t *)aligned_alloc_custom(64, MAX_BAMS_PER_CHUNK * sizeof(uint32_t));
        if (!batch->chunks[i].bam_lens) {
            fprintf(stderr, "Failed to allocate bam_lens for core %d\n", i);
            abort();
        }
        batch->chunks[i].text_len = 0;
        batch->chunks[i].count = 0;
    }
}

void destroy_sam_parse_batch(SamParseBatch *batch) {
    for (int i = 0; i < 64; ++i) {
        if (batch->chunks[i].text_buf) {
            aligned_free_custom((unsigned char*)batch->chunks[i].text_buf);
            batch->chunks[i].text_buf = nullptr;
        }
        if (batch->chunks[i].bam_lens) {
            aligned_free_custom((unsigned char*)batch->chunks[i].bam_lens);
            batch->chunks[i].bam_lens = nullptr;
        }
    }
}

void init_batch_buffers(SamParseByteBatch *batch)
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

void destroy_batch_buffers(SamParseByteBatch *batch)
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


//使用内存读写-----------------------------------------------------------------------------------------------------------------------
inline bool mem_getline(MemReader &r, kstring_t *ks)
{
    if (r.pos >= r.size)
        return false;

    size_t start = r.pos;

    while (r.pos < r.size && r.base[r.pos] != '\n')
        r.pos++;

    size_t len = r.pos - start;

    if (len > 0 && r.base[r.pos - 1] == '\r')
        len--;

    //保证空间
    if (len + 1 > ks->m) {
        ks->m = len + 1;
        ks->s = (char*)realloc(ks->s, ks->m);
    }

    memcpy(ks->s, r.base + start, len);

    ks->s[len] = '\0';
    ks->l = len;

    if (r.pos < r.size && r.base[r.pos] == '\n')
        r.pos++;

    return true;
}

void init_mem_writer(MemWriter &w, size_t cap = 64 * 1024 * 1024)
{
    w.data = (char*)malloc(cap);
    w.capacity = cap;
    w.size = 0;
}

void ensure_capacity(MemWriter &w, size_t add)
{
    if (w.size + add > w.capacity) {
        while (w.size + add > w.capacity)
            w.capacity *= 2;
        w.data = (char*)realloc(w.data, w.capacity);
    }
}        

inline void write_block_to_mem(MemWriter &w, bam_block *block)
{
    ensure_capacity(w, block->length);

    memcpy(
        w.data + w.size,
        block->data,
        block->length
    );

    w.size += block->length;
}

inline int mem_read_block(char *base, size_t size, size_t &pos, bam_block *block)
{
    if (pos >= size)
        return -1;

    uint16_t bsize = *(uint16_t*)(base + pos + 16);
    bsize += 1;

    if (pos + bsize > size)
        return -1;

    memcpy(block->data, base + pos, bsize);

    block->length = bsize;
    pos += bsize;

    return bsize;
}

inline void write_sam_to_mem(MemWriter &w, const char *data, size_t len)
{
    ensure_capacity(w, len);

    memcpy(w.data + w.size, data, len);

    w.size += len;
}


//---------------------------------------------------------------------------------------------------------------------
void SwBam::ProducerSwBamTask_memory(BGZF *fp, BamRead *read , char *bam_mem, size_t bam_size) {
    printf("ProducerSwBamTask_memory started\n");
    double t0 = GetTime();

    std::vector<bam_block*> tmp_chunks;
    bam_block* block = read->getEmpty();
    int global_block_id = 0;
    int block_num = 0;
    int group_num = 0;
    int ret=-1;

    size_t pos = fp->block_address;
    tmp_chunks.reserve(64);

    while (true) {

        // temp1 = GetTime();
        ret = mem_read_block(
                bam_mem,
                bam_size,
                pos,
                block);
        // t_bam2sam_read += GetTime() - temp1;
        
        if (ret < 0) break; //读取错误或结束

        if(block->length == 28){
            // printf("Read a EOF block!\n");
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

    printf("ProducerSwBamTask_memory finished with block_num = %d , group_num = %d , cost %lf\n", block_num , group_num , GetTime() - t0 );

    //标记读取处理结束
    read->markComplete();  

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
            // printf("Read a EOF block!\n");
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
            degz_paras[i].status = 0;
            degz_paras[i].output_records = complete->getResultBuf(i); 
        }

        //调用从核进行解压缩，解析
        // temp2 = GetTime();
        {
            std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
            __real_athread_spawn((void *) slave_decompressfunc, degz_paras, 1);
            athread_join();
        }
        // t_bam2sam_slave += GetTime() - temp2;

        // temp2 = GetTime();
        for (int i = 0; i < 64; i++) {
            if (degz_paras[i].status == 0 && degz_paras[i].input_block != NULL) {

                int n_records = degz_paras[i].n_records;
                auto &records = degz_paras[i].output_records; 

                for (int j = 0; j < n_records; j++) {
                    bam1_t* queue_bam = complete->getEmpty(); 
                    
                    bam1_t temp = *queue_bam;
                    *queue_bam = *(records[j]);
                    *(records[j]) = temp;
                
                    complete->inputBam1_t(queue_bam);
                }
                read->backBlock(degz_paras[i].input_block);
            }
        }
        // t_bam2sam_for += GetTime() - temp2;
    }

    complete->markComplete(); // 标记全部处理完成

    printf("ConsumerSwBamTask finished , cost %lf!\n", GetTime() - t0); 
}

// FusedBamToSam: 单线程融合 Producer + Consumer + Writer
void SwBam::FusedBamToSam(BamRead *read, BamComplete *complete,
                           sam_hdr_t *h, MemReader &reader, MemWriter &mem_writer) {
    double t1 = GetTime();
    double t_decomp = 0, t_format = 0 , t_read = 0 , t_collect = 0;
    double ts;
    const int NB = 64;

    // int max_bams_per_block = 0 , max_data_size = 0 , max_line_size = 0;
    Para degz_para[NB];
    for (int b = 0; b < NB; b++)
        degz_para[b].output_records = complete->getResultBuf(b);

    SamFormatBatch *fmt_cur  = new SamFormatBatch;
    SamFormatBatch *fmt_prev = new SamFormatBatch;
    init_sam_format_batch_buffers(fmt_cur);  fmt_cur->hdr  = h;
    init_sam_format_batch_buffers(fmt_prev); fmt_prev->hdr = h;

    bool has_pending_write = false;

    auto flush_pending_write = [&]() {
        if (!has_pending_write) return;
        for (int i = 0; i < NB; i++) {
            kstring_t *ks = &fmt_prev->core_out_lines[i];
            // if(ks->l> max_line_size) { max_line_size = ks->l; }   //边界检测！
            if (ks->l > 0) write_sam_to_mem(mem_writer, ks->s, ks->l);
        }
        has_pending_write = false;
    };

    auto do_read_group = [&](std::vector<bam_block*> &blocks) {
        blocks.clear();
        for (int b = 0; b < NB; b++) {
            bam_block *blk = read->getEmpty();
            int ret = mem_read_block(reader.base, reader.size, reader.pos, blk);
            if (ret < 0 || blk->length == 28) {
                read->backBlock(blk);
                break;
            }
            blk->block_id = b;
            blk->pos      = 0;
            blocks.push_back(blk);
        }
    };

    long long total_records = 0, group_count = 0, block_count = 0;

    std::vector<bam_block*> prev_blocks;
    prev_blocks.reserve(NB);

    // 预读第一组
    std::vector<bam_block*> next_blocks;
    next_blocks.reserve(NB);
    do_read_group(next_blocks);

    while (!next_blocks.empty()) {
        std::vector<bam_block*> cur_blocks = std::move(next_blocks);
        next_blocks.reserve(NB);

        int n_blocks = (int)cur_blocks.size();
        block_count += n_blocks;

        // 设置 Para
        for (int b = 0; b < NB; b++) {
            if (b < n_blocks) {
                degz_para[b].block_id      = b;
                degz_para[b].input_block   = cur_blocks[b];
                degz_para[b].un_comp_block = complete->getBuffer(b);
                degz_para[b].status        = 0;
                degz_para[b].n_records     = 0;
            } else {
                degz_para[b].input_block = nullptr;
                degz_para[b].status      = -1;
                degz_para[b].n_records   = 0;
            }
        }

        ts = GetTime();
        __real_athread_spawn((void*)slave_decompressfunc, degz_para, 1);
        flush_pending_write();
        for (auto blk : prev_blocks) read->backBlock(blk);
        prev_blocks.clear();
        athread_join();
        t_decomp += GetTime() - ts;

        // 收集 bam1_t 到 SamFormatBatch 中
        ts = GetTime();
        fmt_cur->count = 0;
        for (int b = 0; b < n_blocks; b++) {
            int nr = degz_para[b].n_records;
            // if(nr> max_bams_per_block) { max_bams_per_block = nr; }   //边界检测！
            for (int r = 0; r < nr; r++){
                fmt_cur->bams[fmt_cur->count++] = degz_para[b].output_records[r];
                // if(fmt_cur->bams[fmt_cur->count-1]->l_data > max_data_size) { max_data_size = fmt_cur->bams[fmt_cur->count-1]->l_data; }   //边界检测！
            }
            total_records += nr;
        }
        t_collect += GetTime() - ts;

        for (int i = 0; i < NB; i++) fmt_cur->core_out_lines[i].l = 0;
        ts = GetTime();
        __real_athread_spawn((void*)slave_sam_format, fmt_cur, 1);
        do_read_group(next_blocks);
        athread_join();
        t_format += GetTime() - ts;

        std::swap(fmt_cur, fmt_prev);
        prev_blocks = std::move(cur_blocks);
        has_pending_write = (fmt_prev->count > 0);
        group_count++;
    }

    flush_pending_write();
    for (auto blk : prev_blocks) read->backBlock(blk);

    destroy_sam_format_batch_buffers(fmt_cur);
    destroy_sam_format_batch_buffers(fmt_prev);
    delete fmt_cur;
    delete fmt_prev;

    // printf("  max_bams_per_block=%d  max_data_size=%d  max_line_size=%d\n", max_bams_per_block, max_data_size, max_line_size);
    printf("FusedBamToSam finished. blocks=%lld groups=%lld records=%lld cost=%.3f s\n",
           block_count, group_count, total_records, GetTime() - t1);
    printf("  decomp_slave=%.3f  format_slave=%.3f  collect=%.3f \n",
           t_decomp, t_format, t_collect);
}

int SwBam::writeBam1_tToSam(samFile *fp, const sam_hdr_t *h, const bam1_t *b) {
    if (sam_write1_sw(fp, h, b) < 0) {
        fprintf(stderr, "Error writing SAM record\n");
        return -1;
    }
    return 0;
}


//---------------------------------------------------------------------------------------------------------------------

void SwBam::ProducerSwBamTask2_parallel_memory_OP(BamWrite *write, sam_hdr_t *h, MemReader reader) {
    printf("ProducerSwBamTask2_parallel_memory (Read-Forward) started\n");
    double t0 = GetTime();

    std::vector<bam1_t*> cur_block;
    std::vector<std::vector<bam1_t*>> cur_group;
    //预留容量，避免 push_back 时频繁扩容
    cur_block.reserve(MAX_RECORDS_PER_BLOCK); 
    cur_group.reserve(64);

    bam1_core_t *c __attribute__((unused));
    uint32_t bam_len, total_len = 0;
    int block_nums = 0, group_nums = 0, bam1_t_nums = 0;

    SamParseBatch batch;
    init_sam_parse_batch(&batch);
    batch.hdr = h;

    while (reader.pos < reader.size) {
        int active_chunks = 0;

        // 阶段 1：主核只做边界切分，记录 src_ptr/src_len，不再 memcpy
        temp1 = GetTime();
        for (int i = 0; i < 64; i++) {
            SamParseChunk *chunk = &batch.chunks[i];
            if (reader.pos >= reader.size) {
                // 清零未使用 chunk，防止 slave_copy_and_count 处理旧数据
                chunk->src_len = 0;
                chunk->count   = 0;
                continue;
            }

            size_t start_pos = reader.pos;
            size_t end_pos   = start_pos + SAM_CHUNK_SIZE;

            if (end_pos < reader.size) {
                while (end_pos < reader.size && reader.base[end_pos] != '\n') {
                    end_pos++;
                }
                if (end_pos < reader.size && reader.base[end_pos] == '\n') {
                    end_pos++;
                }
            } else {
                end_pos = reader.size;
            }

            size_t read_len = end_pos - start_pos;
            if (read_len >= CHUNK_BUFFER_SIZE) {
                fprintf(stderr, "FATAL: chunk exceeds buffer! read_len=%zu\n", read_len);
                abort();
            }

            // 只记录原始指针和长度，memcpy 留给从核并行完成
            chunk->src_ptr = reader.base + start_pos;
            chunk->src_len = read_len;
            chunk->count   = 0;

            reader.pos = end_pos;
            active_chunks++;
        }
        t_sam2bam_read += GetTime() - temp1;

        if (active_chunks == 0) break;

        // 阶段 2：从核并行 copy + count
        {
            std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
            __real_athread_spawn((void *)slave_copy_and_count, &batch, 1);
            athread_join();
        }

        // 阶段 3：主核分配 bam1_t*（getEmptyBatch 一次 memcpy 取出 count 个，
        for (int i = 0; i < active_chunks; i++) {
            SamParseChunk *chunk = &batch.chunks[i];
            write->getEmptyBatch(chunk->bams, chunk->count);
        }

        // 阶段 4：从核并行解析
        temp4 = GetTime();
        {
            std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
            __real_athread_spawn((void *)slave_sam_parse_chunk, &batch, 1);
            athread_join();
        }
        t_sam_parse += GetTime() - temp4;

        // 阶段 5：主核打包写入队列
        temp4 = GetTime();
        for (int i = 0; i < active_chunks; i++) {
            SamParseChunk *chunk = &batch.chunks[i];
            for (int j = 0; j < chunk->count; j++) {
                bam1_t *b = chunk->bams[j];
                bam1_t_nums++;
                bam_len = (int)chunk->bam_lens[j]; // 从核已在 parse 时预算好，无需再访问 bam1_t

                if (bam_len + 4 + total_len <= BGZF_BLOCK_SIZE) {
                    cur_block.push_back(b);
                    total_len += bam_len + 4;
                } else {
                    block_nums++;
                    // move 而非 copy：避免 ~3200 次/轮次的堆拷贝，move 后 cur_block 变空
                    cur_group.push_back(std::move(cur_block));
                    cur_block.clear();
                    cur_block.reserve(MAX_RECORDS_PER_BLOCK); // 立即重新预留，避免下次 push_back 重分配
                    cur_block.push_back(b);
                    total_len = bam_len + 4;
                }

                if ((int)cur_group.size() >= 64) {
                    group_nums++;
                    // move 传入：inputGroup 内部再 move 进队列，整条路径无深拷贝
                    write->inputGroup(std::move(cur_group));
                    cur_group.clear();
                    cur_group.reserve(64);
                }
            }
        }
        t_sam2bam_for += GetTime() - temp4;
    }

    if (!cur_block.empty()) {
        block_nums++;
        cur_group.push_back(std::move(cur_block));
    }
    if (!cur_group.empty()) {
        group_nums++;
        write->inputGroup(std::move(cur_group));
    }

    write->markComplete();
    printf("ProducerSwBamTask2_parallel_memory finished. bam1_t_nums=%d, blocks=%d, groups=%d, cost %lf\n",
           bam1_t_nums, block_nums, group_nums, GetTime() - t0);
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

    SamParseByteBatch batch;
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

    printf("ProducerSwBamTask2_parallel finished with bam1_t_nums = %d , block_num = %d , group_num = %d , cost %lf\n", bam1_t_nums, block_nums , group_nums , GetTime() - t0 );

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
        temp2 = GetTime();
        {
            std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
            __real_athread_spawn((void*)slave_compressfunc, comp_paras, 1);
            athread_join();
        }
        t_sam2bam_slave += GetTime() - temp2;

        for (int i = 0; i < 64; i++) {
            //处理完的结果放入到complete中
            if (comp_paras[i].status == 0 && comp_paras[i].output_block != nullptr) {
                complete->pushCompressedBlock(comp_paras[i].output_block);
            }

            if (comp_paras[i].input_records && comp_paras[i].n_records > 0) {
                write->backBamBatch(comp_paras[i].input_records, comp_paras[i].n_records);
            }
        }
    }

    complete->markComplete(); 
    printf("ConsumerSwBamTask2 finished , cost %lf!\n", GetTime() - t0); 
}

// FusedSamToBam: 单线程融合 Producer + Consumer + Writer
void SwBam::FusedSamToBam(BamWrite *write, BamWriteComplete *complete,
                          sam_hdr_t *h, MemReader reader, MemWriter &mem_writer) {
    double t0 , t1 = GetTime();
    double t_copy_count = 0, t_parse = 0, t_compress = 0, t_pack = 0;
    double ts;

    SamParseBatch batch;
    init_sam_parse_batch(&batch);
    batch.hdr = h;

    std::vector<bam1_t*> cur_block;
    std::vector<std::vector<bam1_t*>> cur_group;
    cur_block.reserve(MAX_RECORDS_PER_BLOCK);
    cur_group.reserve(64);

    uint32_t bam_len, total_len = 0;
    int block_nums = 0, group_nums = 0, bam1_t_nums = 0;
    long long bgzf_nums = 0;

    // int max_bams_per_chunk = 0 , max_data_size = 0;

    // 双缓冲 Comp_Para：A/B 交替使用，一个给当前 spawn，一个保存上次结果
    Comp_Para comp_buf_A[64], comp_buf_B[64];
    Comp_Para *comp_active  = comp_buf_A;
    Comp_Para *comp_pending = comp_buf_B;
    std::vector<std::vector<bam1_t*>> pending_group;
    bool has_pending = false;

    // 处理上一次 compress 的结果（写盘 + 回收 bam1_t）
    auto flush_pending = [&]() {
        if (!has_pending) return;
        for (int k = 0; k < 64; k++) {
            if (comp_pending[k].status == 0 && comp_pending[k].output_block) {
                write_block_to_mem(mem_writer, comp_pending[k].output_block);
                complete->backBlock(comp_pending[k].output_block);
                bgzf_nums++;
            }
            if (comp_pending[k].input_records && comp_pending[k].n_records > 0) {
                write->backBamBatch(comp_pending[k].input_records, comp_pending[k].n_records);
            }
        }
        pending_group.clear();
        has_pending = false;
    };

    // 压缩一个 group（64 blocks）并将上一次结果在从核工作期间写盘
    auto do_compress = [&]() {
        for (int k = (int)cur_group.size(); k < 64; k++)
            cur_group.push_back({});

        for (int k = 0; k < 64; k++) {
            if (cur_group[k].empty()) {
                comp_active[k].block_id = k;
                comp_active[k].input_records = nullptr;
                comp_active[k].n_records = 0;
                comp_active[k].output_block = nullptr;
                comp_active[k].un_comp_block = nullptr;
                comp_active[k].output_size = 0;
                comp_active[k].status = -1;
                continue;
            }
            comp_active[k].block_id = k;
            comp_active[k].input_records = cur_group[k].data();
            comp_active[k].n_records = (int)cur_group[k].size();
            comp_active[k].output_block = complete->getEmpty();
            comp_active[k].un_comp_block = complete->getBuffer(k);
            comp_active[k].output_size = 0;
            comp_active[k].status = 0;
        }

        ts = GetTime();
        __real_athread_spawn((void*)slave_compressfunc, comp_active, 1);
        flush_pending();   // 主核：在从核压缩期间处理上一轮结果
        athread_join();
        t_compress += GetTime() - ts;

        // 当前结果变为 pending，下一次 spawn 时再处理
        pending_group = std::move(cur_group);
        std::swap(comp_active, comp_pending);
        has_pending = true;

        cur_group.clear();
        cur_group.reserve(64);
    };

    while (reader.pos < reader.size) {
        int active_chunks = 0;

        // 阶段 1：主核边界切分，打包一组数据（64个4M大小的chunk）放入batch中，几乎不花时间
        for (int i = 0; i < 64; i++) {
            SamParseChunk *chunk = &batch.chunks[i];
            if (reader.pos >= reader.size) {
                chunk->src_len = 0;
                chunk->count   = 0;
                continue;
            }
            size_t start_pos = reader.pos;
            size_t end_pos   = start_pos + SAM_CHUNK_SIZE;
            if (end_pos < reader.size) {
                while (end_pos < reader.size && reader.base[end_pos] != '\n') end_pos++;
                if (end_pos < reader.size && reader.base[end_pos] == '\n') end_pos++;
            } else {
                end_pos = reader.size;
            }
            size_t read_len = end_pos - start_pos;
            // if (read_len >= CHUNK_BUFFER_SIZE) {
            //     fprintf(stderr, "FATAL: chunk exceeds buffer! read_len=%zu\n", read_len);
            //     abort();
            // }
            chunk->src_ptr = reader.base + start_pos;
            chunk->src_len = read_len;
            chunk->count   = 0;
            reader.pos = end_pos;
            active_chunks++;
        }
        if (active_chunks == 0) break;

        // 阶段 2：从核 copy+count
        ts = GetTime();
        __real_athread_spawn((void *)slave_copy_and_count, &batch, 1);
        flush_pending();   // 主核：在从核 copy 期间写上一轮压缩结果
        athread_join();
        t_copy_count += GetTime() - ts;

        // 阶段 3：主核 getEmptyBatch，几乎不花时间
        int pre_alloc_count[64];
        for (int i = 0; i < active_chunks; i++) {
            SamParseChunk *chunk = &batch.chunks[i];
            pre_alloc_count[i] = chunk->count;
            // if(chunk->count> max_bams_per_chunk) { max_bams_per_chunk = chunk->count; }   //边界检测！
            write->getEmptyBatch(chunk->bams, chunk->count);
        }

        // 阶段 4：从核并行解析
        ts = GetTime();
        __real_athread_spawn((void *)slave_sam_parse_chunk, &batch, 1);
        athread_join();
        t_parse += GetTime() - ts;

        // 阶段 4.5：回收解析失败的 bam1_t，几乎不花时间
        for (int i = 0; i < active_chunks; i++) {
            SamParseChunk *chunk = &batch.chunks[i];
            for (int j = chunk->count; j < pre_alloc_count[i]; j++)
                write->backBam(chunk->bams[j]);
        }

        // 阶段 5：打包（将全部bam1_t连续打包，64个block为一组） + 压缩（交织执行）
        t0 = GetTime();
        for (int i = 0; i < active_chunks; i++) {
            SamParseChunk *chunk = &batch.chunks[i];
            for (int j = 0; j < chunk->count; j++) {
                bam1_t *b = chunk->bams[j];
                bam1_t_nums++;
                bam_len = (int)chunk->bam_lens[j];
                // if(b->l_data > max_data_size) { max_data_size = b->l_data; }    //边界检测！

                if (bam_len + 4 + total_len <= BGZF_BLOCK_SIZE) {
                    cur_block.push_back(b);
                    total_len += bam_len + 4;
                } else {
                    block_nums++;
                    cur_group.push_back(std::move(cur_block));
                    cur_block.clear();
                    cur_block.reserve(MAX_RECORDS_PER_BLOCK);
                    cur_block.push_back(b);
                    total_len = bam_len + 4;
                }

                if ((int)cur_group.size() >= 64) {
                    group_nums++;
                    do_compress();
                }
            }
        }
        t_pack += GetTime() - t0;
    }

    // 收尾
    if (!cur_block.empty()) {
        block_nums++;
        cur_group.push_back(std::move(cur_block));
    }
    if (!cur_group.empty()) {
        group_nums++;
        do_compress();
    }
    flush_pending();

    destroy_sam_parse_batch(&batch);

    // printf("max_bams_per_chunk=%d, max_data_size=%d\n",max_bams_per_chunk , max_data_size);
    printf("FusedSamToBam finished. bam1_t=%d, blocks=%d, groups=%d, bgzf=%lld, cost %lf\n",
           bam1_t_nums, block_nums, group_nums, bgzf_nums, GetTime() - t1);
    printf("  copy_count_slave=%lf  parse_slave=%lf  pack(+compress)=%lf  compress_slave=%lf\n",
           t_copy_count, t_parse, t_pack, t_compress);
}

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
    double t0 = GetTime();

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


    #define USE_MEMORY
    #ifdef USE_MEMORY

    t0 = GetTime();
    MemReader reader;
    MemWriter mem_writer;
    size_t sam_size;
    char *sam_mem = nullptr;
    size_t bam_size;
    char *bam_mem = nullptr;

    switch (sin->format.format) {
        case bam:{
            printf("MEMORY BAM\n");

            //内存读--
            FILE *f = fopen(cmd_info_->in_file_name_.c_str(), "rb");
            if (!f) { perror("fopen"); exit(1); }
            fseek(f, 0, SEEK_END);
            bam_size = ftell(f);
            fseek(f, 0, SEEK_SET);
            bam_mem = (char*)malloc(bam_size);
            fread(bam_mem, 1, bam_size, f);
            fclose(f);
            printf("Loaded BAM into memory: %.2f MB\n",bam_size / 1024.0 / 1024.0);

            reader.base = bam_mem;
            reader.size = bam_size;
            reader.pos = sin->fp.bgzf->block_address;

            //内存写--
            size_t sam_estimate = bam_size * 5;
            init_mem_writer(mem_writer, sam_estimate);
            printf("Initialized SAM memory writer: %.2f MB\n",sam_estimate / 1024.0 / 1024.0);
        
            break;
        }

        case sam:{
            printf("MEMORY SAM\n");

            //内存读--
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

            // 初始化 reader
            reader.base = sam_mem;
            reader.size = sam_size;
            reader.pos  = 0;

            // 跳过 header
            kstring_t tmp;
            tmp.l = 0;
            tmp.m = MAX_SAM_LINE_SIZE;
            tmp.s = (char*)malloc(MAX_SAM_LINE_SIZE);

            while (true) {
                size_t old_pos = reader.pos;
                if (!mem_getline(reader, &tmp))
                    break;
                if (tmp.l == 0)
                    continue;
                if (tmp.s[0] != '@') {
                    reader.pos = old_pos;   // 回退到该行开始
                    break;
                }
            }
            free(tmp.s);

            //内存写--
            size_t bam_estimate = sam_size / 2;
            init_mem_writer(mem_writer, bam_estimate);
            printf("Initialized BAM memory writer: %.2f MB\n",bam_estimate / 1024.0 / 1024.0);

            break;
        }

        default:
            fprintf(stderr, "Unknown file format\n");
            break;
    }
    printf("Complete the memory cost %lf\n", GetTime() - t0);
    #endif

    double tbody;
    tbody = GetTime();
    switch (sin->format.format) {
        case bam:{
            t0 = GetTime();

            //单线程：融合版本，内存模式，从内存读写数据
            #define USE_FUSED_BAM2SAM
            #if defined(USE_MEMORY) && defined(USE_FUSED_BAM2SAM)
            {
                printf("Enable FUSED BAM2SAM single-thread mode (USE_MEMORY + USE_FUSED_BAM2SAM)!!!\n");
                // BamRead(128,1): 128块压缩block池（两轮各64），queue不用
                // BamComplete(1): 仅需 buffer_pool_（解压缓冲）和 result_pool_（bam1_t）
                t0 = GetTime();
                read     = new BamRead(128, 1);
                complete = new BamComplete(1);
                printf("Complete the queue initialization cost %lf\n", GetTime() - t0);
                FusedBamToSam(read, complete, hdr, reader, mem_writer);
            }

            //流水线：内存模式，从内存读写数据
            #elif defined(USE_MEMORY)
            {
                printf("Enable 3-thread BAM2SAM pipeline (USE_MEMORY)!!!\n");
                //BamRead(x , y) x是blocks的内存池大小，y是打包之后队列的大小
                //BamComplete(x) x是bam1_t的内存池大小

                //串行测试用
                // read = new BamRead(16640 , 260);
                // complete = new BamComplete(2932000);

                read     = new BamRead(320, 5);
                complete = new BamComplete(132000);
                printf("Complete the queue initialization cost %lf\n", GetTime() - t0);

                thread producer(bind(&SwBam::ProducerSwBamTask_memory, this, sin->fp.bgzf, read, bam_mem, bam_size));
                thread consumer(bind(&SwBam::ConsumerSwBamTask, this, read, complete));

                t0 = GetTime();
                SamFormatBatch batch;
                init_sam_format_batch_buffers(&batch);
                batch.hdr = hdr;

                long long num = 0;
                bam1_t *b;
                while (true) {
                    batch.count = 0;
                    while (batch.count < BATCH_SIZE) {
                        b = complete->getBam1_t();
                        if (!b) break;
                        num++;
                        batch.bams[batch.count++] = b;
                    }
                    if (batch.count == 0) break;

                    {
                        std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
                        __real_athread_spawn((void*)slave_sam_format, &batch, 1);
                        athread_join();
                    }

                    for (int i = 0; i < 64; i++) {
                        kstring_t *ks = &batch.core_out_lines[i];
                        write_sam_to_mem(mem_writer, ks->s, ks->l);
                    }
                    complete->backBam1_tBatch(batch.bams, batch.count);
                }

                destroy_sam_format_batch_buffers(&batch);
                producer.join();
                consumer.join();
                printf("Complete main thread writing to sam cost %lf\n", GetTime() - t0);
                printf("The total bam1_t nums is %lld\n", num);
            }

            //流水线：非内存模式，从磁盘读写数据
            #else
            {
                printf("Enable 3-thread BAM2SAM pipeline (disk mode)!!!\n");
                read     = new BamRead(320, 5);
                complete = new BamComplete(132000);
                printf("Complete the queue initialization cost %lf\n", GetTime() - t0);

                thread producer(bind(&SwBam::ProducerSwBamTask, this, sin->fp.bgzf, read));
                thread consumer(bind(&SwBam::ConsumerSwBamTask, this, read, complete));

                t0 = GetTime();
                SamFormatBatch batch;
                init_sam_format_batch_buffers(&batch);
                batch.hdr = hdr;

                long long num = 0;
                bam1_t *b;
                while (true) {
                    batch.count = 0;
                    while (batch.count < BATCH_SIZE) {
                        b = complete->getBam1_t();
                        if (!b) break;
                        num++;
                        batch.bams[batch.count++] = b;
                    }
                    if (batch.count == 0) break;

                    {
                        std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
                        __real_athread_spawn((void*)slave_sam_format, &batch, 1);
                        athread_join();
                    }

                    for (int i = 0; i < 64; i++) {
                        kstring_t *ks = &batch.core_out_lines[i];
                        if (hwrite(sout->fp.hfile, ks->s, ks->l) != (ssize_t)ks->l)
                            fprintf(stderr, "write failed\n");
                    }
                    complete->backBam1_tBatch(batch.bams, batch.count);
                }

                destroy_sam_format_batch_buffers(&batch);
                producer.join();
                consumer.join();
                printf("Complete main thread writing to sam cost %lf\n", GetTime() - t0);
                printf("The total bam1_t nums is %lld\n", num);
            }
            #endif

            break;
        }
            

        case sam:{
            t0 = GetTime();

            //单线程：融合版本，内存模式，从内存读写数据
            #define USE_FUSED_SAM2BAM
            #if defined(USE_MEMORY) && defined(USE_FUSED_SAM2BAM)
            {
                printf("Enable FUSED single-thread mode (USE_MEMORY + USE_FUSED_SAM2BAM)!!!\n");
                //BamWrite(x , y) x是bam1_t的内存池大小 , y是打包的block的size
                //BamWriteComplete(x) x是bam_block的内存池大小
                write = new BamWrite(655360 ,640);
                writeComplete = new BamWriteComplete(130);
                printf("Complete the queue initialization cost %lf\n", GetTime() - t0);

                FusedSamToBam(write, writeComplete, hdr, reader, mem_writer);
            }


            //流水线：内存模式，从内存读写数据
            #elif defined(USE_MEMORY)
            {
                printf("Enable the SAM_PARSE_PARALLEL (3-thread, read and write from memory)!!!\n");
                //BamWrite(x , y) x是bam1_t的内存池大小 , y是打包的block的size
                //BamWriteComplete(x) x是bam_block的内存池大小

                //串行测试用
                // write = new BamWrite(2938880 , 2870);
                // writeComplete = new BamWriteComplete(16400);

                write = new BamWrite(655360 ,640);
                writeComplete = new BamWriteComplete(130);
                printf("Complete the queue initialization cost %lf\n", GetTime() - t0);
                
                thread producer2(bind(&SwBam::ProducerSwBamTask2_parallel_memory_OP, this, write, hdr, reader));
                // producer2.join();
                thread consumer2(bind(&SwBam::ConsumerSwBamTask2, this, write , writeComplete));
                // consumer2.join();

                t0 = GetTime();
                long long num2 = 0;
                bam_block* comp_block;
                while (true) {
                    comp_block = writeComplete->getCompressedBlock();
                    if (comp_block == nullptr) break;
                    num2++;
                    write_block_to_mem(mem_writer, comp_block);
                    writeComplete->backBlock(comp_block);
                }

                producer2.join();
                consumer2.join();
                printf("The actual time of reading sam cost %lf\n", t_sam2bam_read);
                printf("The actual time of sam parsing cost %lf\n", t_sam_parse);
                printf("The actual time of for loops in producer cost %lf\n", t_sam2bam_for);
                printf("The actual time of slave cost %lf\n", t_sam2bam_slave);
                printf("The actual time of writing to bam cost %lf\n", t_sam2bam_write);
                printf("Complete main thread writing to bam cost %lf\n", GetTime() - t0);
                printf("The total BGZF nums is %lld\n", num2);
            }


            //流水线：非内存模式，从磁盘读写数据
            #else
            {
                printf("Enable the SAM_PARSE_PARALLEL (3-thread, read and write from disk)!!!\n");
                //BamWrite(x , y) x是bam1_t的内存池大小 , y是打包的block的size
                //BamWriteComplete(x) x是bam_block的内存池大小

                //串行测试用
                // write = new BamWrite(2938880 , 2870);
                // writeComplete = new BamWriteComplete(16400);

                write = new BamWrite(655360 ,640);
                writeComplete = new BamWriteComplete(130);
                printf("Complete the queue initialization cost %lf\n", GetTime() - t0);
                
                thread producer2(bind(&SwBam::ProducerSwBamTask2_parallel, this, sin, write, hdr));
                thread consumer2(bind(&SwBam::ConsumerSwBamTask2, this, write, writeComplete));
                t0 = GetTime();
                long long num2 = 0;
                bam_block* comp_block;

                while (true) {
                    comp_block = writeComplete->getCompressedBlock();
                    if (comp_block == nullptr) break;
                    num2++;
                    (void)writeBlockTobam(sout->fp.bgzf, comp_block);
                    writeComplete->backBlock(comp_block);
                }

                producer2.join();
                consumer2.join();
                printf("Complete main thread writing to bam cost %lf\n", GetTime() - t0);
                printf("The total BGZF nums is %lld\n", num2);
            }
            #endif
            break;
        }
            
        default:
            fprintf(stderr, "Unknown file format\n");
            break;
    }


    printf("Complete the body1 cost %lf\n", GetTime() - tbody);

    #ifdef USE_MEMORY
    #define DUMP_MEM
    #ifdef DUMP_MEM
    //将内存中的主体数据落盘到sout中
    double dump_t0 = GetTime();
    if (sout->format.format == bam) {
        bgzf_flush(sout->fp.bgzf);
        if (hwrite(sout->fp.bgzf->fp, mem_writer.data, mem_writer.size) != mem_writer.size) {
            fprintf(stderr, "Error dumping memory to BAM output file\n");
        }
    } else if (sout->format.format == sam) {
        hflush(sout->fp.hfile); 
        if (hwrite(sout->fp.hfile, mem_writer.data, mem_writer.size) != mem_writer.size) {
            fprintf(stderr, "Error dumping memory to SAM output file\n");
        }
    }
    // 释放内存池
    if (mem_writer.data) {
        free(mem_writer.data);
        mem_writer.data = nullptr;
    }
    if (sam_mem){
        free(sam_mem);
    }
    if (bam_mem){
        free(bam_mem);
    }
    printf("Dump memory to output file cost %lf\n", GetTime() - dump_t0);
    #endif
    #endif

    //释放内存
    t0 = GetTime();
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
    printf("close the files cost %lf\n", GetTime() - t0);
    printf("Complete the body cost %lf\n", GetTime() - tbody);

}


