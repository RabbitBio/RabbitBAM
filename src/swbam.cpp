#include "swbam.h"

using namespace std;

extern "C" {
#include <athread.h>
#include <pthread.h>
    void slave_decompressfunc();
    void slave_decompress_checked();
    void slave_decompress_filterfunc();
    void slave_decompress_filter_checked();
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

namespace {

bool HasBamFilterRequest(const CmdInfo *cmd_info) {
    return cmd_info->min_mapq_ >= 0 ||
           cmd_info->max_mapq_ >= 0 ||
           cmd_info->require_flag_ != 0 ||
           cmd_info->exclude_flag_ != 0 ||
           !cmd_info->ref_name_.empty() ||
           cmd_info->min_read_len_ >= 0 ||
           cmd_info->max_read_len_ >= 0;
}

BamFilterOptions BuildBamFilterOptions(const CmdInfo *cmd_info) {
    BamFilterOptions filter;
    filter.min_mapq = cmd_info->min_mapq_;
    filter.max_mapq = cmd_info->max_mapq_;
    filter.require_flag = cmd_info->require_flag_;
    filter.exclude_flag = cmd_info->exclude_flag_;
    filter.ref_tid = -2;
    filter.min_read_len = cmd_info->min_read_len_;
    filter.max_read_len = cmd_info->max_read_len_;
    return filter;
}

bool IsSamLikeFormat(int format) {
    return format == sam || format == text_format;
}

bool IsBamLikeFormat(int format) {
    return format == bam || format == binary_format;
}

const char *BoundsLimitName(int limit_id) {
    switch (limit_id) {
        case BOUNDS_LIMIT_MAX_RECORDS_PER_BLOCK: return "MAX_RECORDS_PER_BLOCK";
        case BOUNDS_LIMIT_MAX_SAM_LINE_SIZE: return "MAX_SAM_LINE_SIZE";
        case BOUNDS_LIMIT_INIT_DATA_SIZE: return "INIT_DATA_SIZE";
        case BOUNDS_LIMIT_MAX_BAMS_PER_CHUNK: return "MAX_BAMS_PER_CHUNK";
        case BOUNDS_LIMIT_CHUNK_BUFFER_SIZE: return "CHUNK_BUFFER_SIZE";
        case BOUNDS_LIMIT_FUSED_SAM2BAM_BAM_POOL_SIZE: return "FUSED_SAM2BAM_BAM_POOL_SIZE";
        case BOUNDS_LIMIT_BGZF_RECORD_SIZE: return "BGZF_BLOCK_SIZE";
        case BOUNDS_LIMIT_MAX_SAM_FORMAT_CORE_BUFFER_SIZE: return "MAX_SAM_FORMAT_CORE_BUFFER_SIZE";
        case BOUNDS_LIMIT_GENERIC_RUNTIME_ERROR: return "RUNTIME_STATUS";
        default: return "UNKNOWN_LIMIT";
    }
}

void ClearBoundsError(BoundsCheckError *err) {
    if (!err) return;
    err->pipeline = nullptr;
    err->stage = nullptr;
    err->limit_name = nullptr;
    err->limit_value = 0;
    err->actual_value = 0;
    err->block_id = -1;
    err->chunk_id = -1;
    err->record_index = -1;
    err->core_id = -1;
}

void SetBoundsError(BoundsCheckError *err,
                    const char *pipeline,
                    const char *stage,
                    int limit_id,
                    long long limit_value,
                    long long actual_value,
                    int block_id,
                    int chunk_id,
                    int record_index,
                    int core_id) {
    if (!err) return;
    err->pipeline = pipeline;
    err->stage = stage;
    err->limit_name = BoundsLimitName(limit_id);
    err->limit_value = limit_value;
    err->actual_value = actual_value;
    err->block_id = block_id;
    err->chunk_id = chunk_id;
    err->record_index = record_index;
    err->core_id = core_id;
}

void SetGenericBoundsError(BoundsCheckError *err,
                           const char *pipeline,
                           const char *stage,
                           long long actual_value,
                           int block_id,
                           int chunk_id,
                           int record_index,
                           int core_id) {
    SetBoundsError(err, pipeline, stage, BOUNDS_LIMIT_GENERIC_RUNTIME_ERROR, 0,
                   actual_value, block_id, chunk_id, record_index, core_id);
}

void PrintBoundsError(const BoundsCheckError &err) {
    fprintf(stderr,
            "BOUNDS CHECK FAILED: pipeline=%s stage=%s limit=%s limit_value=%lld actual_value=%lld block_id=%d chunk_id=%d record_index=%d core_id=%d\n",
            err.pipeline ? err.pipeline : "unknown",
            err.stage ? err.stage : "unknown",
            err.limit_name ? err.limit_name : "unknown",
            err.limit_value,
            err.actual_value,
            err.block_id,
            err.chunk_id,
            err.record_index,
            err.core_id);
    fprintf(stderr, "ERROR: boundary validation failed while --validate-bounds was enabled.\n");
}

void SetBoundsErrorFromCheckedPara(BoundsCheckError *err,
                                   const char *pipeline,
                                   const char *stage,
                                   const CheckedPara &para,
                                   int core_id) {
    if (para.limit_id == BOUNDS_LIMIT_NONE) {
        SetGenericBoundsError(err, pipeline, stage, para.status, para.block_id, -1,
                              para.record_index, core_id);
        return;
    }
    SetBoundsError(err, pipeline, stage, para.limit_id, para.limit_value, para.actual_value,
                   para.block_id, -1, para.record_index, core_id);
}

void SetBoundsErrorFromCheckedBam2BamPara(BoundsCheckError *err,
                                          const char *pipeline,
                                          const char *stage,
                                          const CheckedBam2BamPara &para,
                                          int core_id) {
    if (para.limit_id == BOUNDS_LIMIT_NONE) {
        SetGenericBoundsError(err, pipeline, stage, para.status, para.block_id, -1,
                              para.record_index, core_id);
        return;
    }
    SetBoundsError(err, pipeline, stage, para.limit_id, para.limit_value, para.actual_value,
                   para.block_id, -1, para.record_index, core_id);
}

} // namespace

SwBam::SwBam(CmdInfo *cmd_info1) {
    cmd_info_ = cmd_info1;
    hdr = nullptr;
    sin = nullptr;
    sout = nullptr;
    read = nullptr;
    complete = nullptr;
    write = nullptr;
    writeComplete = nullptr;
}

SwBam::~SwBam() {
    delete read;
    delete complete;
    delete write;
    delete writeComplete;
}

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

    // if (w.size + len > w.capacity) {
    //     w.size = 0; 
    // }

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
    double t_decomp = 0, t_format = 0 , t_read = 0 , t_write = 0, t_collect = 0;
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
        // double t0 = GetTime();
        if (!has_pending_write) return;
        for (int i = 0; i < NB; i++) {
            kstring_t *ks = &fmt_prev->core_out_lines[i];
            // if(ks->l> max_line_size) { max_line_size = ks->l; }   //边界检测！
            if (ks->l > 0) write_sam_to_mem(mem_writer, ks->s, ks->l);
        }
        has_pending_write = false;
        // t_write += GetTime() - t0;
    };

    auto do_read_group = [&](std::vector<bam_block*> &blocks) {
        // double t0 = GetTime();
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
        // t_read += GetTime() - t0;
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

        // ts = GetTime();
        // flush_pending_write();
        // for (auto blk : prev_blocks) read->backBlock(blk);
        // prev_blocks.clear();
        // t_write += GetTime() - ts;
        

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

        // ts = GetTime();
        // do_read_group(next_blocks);
        // t_read += GetTime() - ts;

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
    printf("  decomp_slave=%.3f  format_slave=%.3f  collect=%.3f  read=%.3f  write=%.3f \n",
           t_decomp, t_format, t_collect, t_read, t_write);
}

namespace {

int FormatSamRecordsChecked(const sam_hdr_t *h,
                            bam1_t **records,
                            int total_records,
                            MemWriter &mem_writer,
                            BoundsCheckError *bounds_error) {
    kstring_t line = {0, MAX_SAM_LINE_SIZE + 1, (char *)malloc(MAX_SAM_LINE_SIZE + 1)};
    if (!line.s) {
        SetGenericBoundsError(bounds_error, "bam2sam", "sam_format_alloc", errno, -1, -1, -1, -1);
        return -1;
    }

    int base_tasks = total_records >> 6;
    int remainder = total_records & 63;
    for (int cid = 0; cid < 64; ++cid) {
        size_t core_total_len = 0;
        int start = 0;
        int end = 0;
        if (cid < remainder) {
            start = cid * (base_tasks + 1);
            end = start + base_tasks + 1;
        } else {
            start = remainder + cid * base_tasks;
            end = start + base_tasks;
        }

        if (start >= total_records || start >= end) continue;

        for (int idx = start; idx < end; ++idx) {
            int line_len = sam_format1(h, records[idx], &line);
            if (line_len < 0) {
                SetGenericBoundsError(bounds_error, "bam2sam", "sam_format", line_len, -1, -1, idx, cid);
                free(line.s);
                return -1;
            }
            if ((size_t)line_len > MAX_SAM_LINE_SIZE) {
                SetBoundsError(bounds_error, "bam2sam", "sam_format",
                               BOUNDS_LIMIT_MAX_SAM_LINE_SIZE, MAX_SAM_LINE_SIZE, line_len,
                               -1, -1, idx, cid);
                free(line.s);
                return -1;
            }
            if (core_total_len + (size_t)line_len + 1 > MAX_SAM_FORMAT_CORE_BUFFER_SIZE) {
                SetBoundsError(bounds_error, "bam2sam", "sam_format_batch",
                               BOUNDS_LIMIT_MAX_SAM_FORMAT_CORE_BUFFER_SIZE,
                               MAX_SAM_FORMAT_CORE_BUFFER_SIZE, core_total_len + line_len + 1,
                               -1, -1, idx, cid);
                free(line.s);
                return -1;
            }
            if (line.l > 0) write_sam_to_mem(mem_writer, line.s, line.l);
            write_sam_to_mem(mem_writer, "\n", 1);
            core_total_len += (size_t)line_len + 1;
        }
    }

    free(line.s);
    return 0;
}

} // namespace

int SwBam::FusedBamToSamChecked(BamRead *read, BamComplete *complete,
                                sam_hdr_t *h, MemReader &reader, MemWriter &mem_writer,
                                BoundsCheckError *bounds_error) {
    double t0 = GetTime();
    const int NB = 64;
    CheckedPara degz_para[NB];
    for (int b = 0; b < NB; ++b) {
        degz_para[b].output_records = complete->getResultBuf(b).data();
        degz_para[b].input_block = nullptr;
        degz_para[b].un_comp_block = nullptr;
        degz_para[b].n_records = 0;
        degz_para[b].status = -1;
        degz_para[b].record_index = -1;
        degz_para[b].actual_value = 0;
        degz_para[b].limit_value = 0;
        degz_para[b].limit_id = BOUNDS_LIMIT_NONE;
    }

    auto do_read_group = [&](std::vector<bam_block*> &blocks) {
        blocks.clear();
        blocks.reserve(NB);
        for (int b = 0; b < NB; ++b) {
            bam_block *blk = read->getEmpty();
            int ret = mem_read_block(reader.base, reader.size, reader.pos, blk);
            if (ret < 0 || blk->length == 28) {
                read->backBlock(blk);
                break;
            }
            blk->block_id = b;
            blk->pos = 0;
            blocks.push_back(blk);
        }
    };

    std::vector<bam_block*> next_blocks;
    next_blocks.reserve(NB);
    do_read_group(next_blocks);

    long long block_count = 0;
    long long group_count = 0;
    long long total_records = 0;

    while (!next_blocks.empty()) {
        std::vector<bam_block*> cur_blocks = std::move(next_blocks);
        next_blocks.clear();
        next_blocks.reserve(NB);

        const int n_blocks = (int)cur_blocks.size();
        block_count += n_blocks;

        for (int b = 0; b < NB; ++b) {
            degz_para[b].block_id = b;
            degz_para[b].input_block = (b < n_blocks) ? cur_blocks[b] : nullptr;
            degz_para[b].un_comp_block = (b < n_blocks) ? complete->getBuffer(b) : nullptr;
            degz_para[b].n_records = 0;
            degz_para[b].status = (b < n_blocks) ? 0 : -1;
            degz_para[b].record_index = -1;
            degz_para[b].actual_value = 0;
            degz_para[b].limit_value = 0;
            degz_para[b].limit_id = BOUNDS_LIMIT_NONE;
        }

        {
            std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
            __real_athread_spawn((void*)slave_decompress_checked, degz_para, 1);
            athread_join();
        }

        for (int b = 0; b < n_blocks; ++b) {
            if (degz_para[b].status != 0) {
                SetBoundsErrorFromCheckedPara(bounds_error, "bam2sam", "decompress", degz_para[b], b);
                for (auto blk : cur_blocks) read->backBlock(blk);
                return -1;
            }
        }

        std::vector<bam1_t*> records;
        records.reserve((size_t)n_blocks * MAX_RECORDS_PER_BLOCK);
        for (int b = 0; b < n_blocks; ++b) {
            for (int r = 0; r < degz_para[b].n_records; ++r)
                records.push_back(degz_para[b].output_records[r]);
            total_records += degz_para[b].n_records;
        }

        if (!records.empty() &&
            FormatSamRecordsChecked(h, records.data(), (int)records.size(), mem_writer, bounds_error) != 0) {
            for (auto blk : cur_blocks) read->backBlock(blk);
            return -1;
        }

        for (auto blk : cur_blocks) read->backBlock(blk);
        do_read_group(next_blocks);
        group_count++;
    }

    printf("FusedBamToSamChecked finished. blocks=%lld groups=%lld records=%lld cost=%.3f s\n",
           block_count, group_count, total_records, GetTime() - t0);
    return 0;
}

// FusedBamToBam: 单线程融合 Producer + Consumer + Writer
void SwBam::FusedBamToBam(BamRead *read, BamComplete *complete, BamWriteComplete *write_complete,
                          MemReader &reader, MemWriter &mem_writer, const BamFilterOptions &filter) {
    double t0 = GetTime();
    double t_decomp_filter = 0, t_pack = 0, t_compress = 0, t_read = 0, t_write = 0;
    const int NB = 64;

    Bam2BamPara paras[NB];
    for (int b = 0; b < NB; ++b) {
        paras[b].output_records = complete->getResultBuf(b).data();
        paras[b].bam_lens = (uint32_t *)aligned_alloc_custom(64, MAX_RECORDS_PER_BLOCK * sizeof(uint32_t));
        if (!paras[b].bam_lens) {
            fprintf(stderr, "Failed to allocate bam2bam bam_lens for block %d\n", b);
            abort();
        }
        paras[b].input_block = nullptr;
        paras[b].un_comp_block = nullptr;
        paras[b].n_total_records = 0;
        paras[b].n_kept_records = 0;
        paras[b].status = -1;
        paras[b].filter = filter;
    }

    Comp_Para comp_buf_A[NB], comp_buf_B[NB];
    Comp_Para *comp_active = comp_buf_A;
    Comp_Para *comp_pending = comp_buf_B;
    std::vector<std::vector<bam1_t*>> pending_group;
    bool has_pending = false;
    long long bgzf_blocks = 0;

    long long total_records = 0;
    long long kept_records = 0;
    long long dropped_records = 0;
    long long input_blocks = 0;
    long long group_count = 0;

    auto flush_pending = [&]() {
        double flush_t0 = GetTime();

        if (!has_pending) return;
        for (int k = 0; k < NB; ++k) {
            if (comp_pending[k].status == 0 && comp_pending[k].output_block) {
                write_block_to_mem(mem_writer, comp_pending[k].output_block);
                write_complete->backBlock(comp_pending[k].output_block);
                bgzf_blocks++;
            }
        }
        pending_group.clear();
        has_pending = false;

        t_write += GetTime() - flush_t0;
    };

    auto do_read_group = [&](std::vector<bam_block*> &blocks) {
        double read_t0 = GetTime();

        blocks.clear();
        blocks.reserve(NB);
        for (int b = 0; b < NB; ++b) {
            bam_block *blk = read->getEmpty();
            int ret = mem_read_block(reader.base, reader.size, reader.pos, blk);
            if (ret < 0 || blk->length == 28) {
                read->backBlock(blk);
                break;
            }
            blk->block_id = b;
            blk->pos = 0;
            blocks.push_back(blk);
        }

        t_read += GetTime() - read_t0;
    };

    auto do_compress = [&](std::vector<std::vector<bam1_t*>> &group,
                           std::vector<bam_block*> &next_blocks) {
        for (int k = (int)group.size(); k < NB; ++k)
            group.push_back({});

        for (int k = 0; k < NB; ++k) {
            if (group[k].empty()) {
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
            comp_active[k].input_records = group[k].data();
            comp_active[k].n_records = (int)group[k].size();
            comp_active[k].output_block = write_complete->getEmpty();
            comp_active[k].un_comp_block = write_complete->getBuffer(k);
            comp_active[k].output_size = 0;
            comp_active[k].status = 0;
        }

        double compress_t0 = GetTime();
        __real_athread_spawn((void*)slave_compressfunc, comp_active, 1);
        flush_pending();
        do_read_group(next_blocks);
        athread_join();
        t_compress += GetTime() - compress_t0;

        pending_group = std::move(group);
        std::swap(comp_active, comp_pending);
        has_pending = true;

        group.clear();
        group.reserve(NB);
    };

    std::vector<bam_block*> next_blocks;
    next_blocks.reserve(NB);
    do_read_group(next_blocks);

    while (!next_blocks.empty()) {
        std::vector<bam_block*> cur_blocks = std::move(next_blocks);
        next_blocks.clear();
        next_blocks.reserve(NB);

        const int n_blocks = (int)cur_blocks.size();
        input_blocks += n_blocks;

        for (int b = 0; b < NB; ++b) {
            paras[b].filter = filter;
            if (b < n_blocks) {
                paras[b].block_id = b;
                paras[b].input_block = cur_blocks[b];
                paras[b].un_comp_block = complete->getBuffer(b);
                paras[b].n_total_records = 0;
                paras[b].n_kept_records = 0;
                paras[b].status = 0;
            } else {
                paras[b].block_id = b;
                paras[b].input_block = nullptr;
                paras[b].un_comp_block = nullptr;
                paras[b].n_total_records = 0;
                paras[b].n_kept_records = 0;
                paras[b].status = -1;
            }
        }

        double decomp_t0 = GetTime();
        __real_athread_spawn((void*)slave_decompress_filterfunc, paras, 1);
        athread_join();
        t_decomp_filter += GetTime() - decomp_t0;

        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].status != 0) {
                fprintf(stderr, "FATAL: slave_decompress_filterfunc failed on block %d with status %d\n",
                        b, paras[b].status);
                abort();
            }
        }

        std::vector<std::vector<bam1_t*>> cur_group;
        cur_group.reserve(NB);
        std::vector<bam1_t*> cur_block;
        cur_block.reserve(MAX_RECORDS_PER_BLOCK);
        uint32_t total_len = 0;

        double pack_t0 = GetTime();
        for (int b = 0; b < n_blocks; ++b) {
            total_records += paras[b].n_total_records;
            kept_records += paras[b].n_kept_records;
            dropped_records += paras[b].n_total_records - paras[b].n_kept_records;

            for (int r = 0; r < paras[b].n_kept_records; ++r) {
                bam1_t *record = paras[b].output_records[r];
                uint32_t bam_len = paras[b].bam_lens[r];

                if (bam_len + 4 + total_len <= BGZF_BLOCK_SIZE) {
                    cur_block.push_back(record);
                    total_len += bam_len + 4;
                } else {
                    if (cur_block.empty()) {
                        fprintf(stderr, "FATAL: bam2bam pack encountered an oversized BAM record\n");
                        abort();
                    }
                    cur_group.push_back(std::move(cur_block));
                    cur_block.clear();
                    cur_block.reserve(MAX_RECORDS_PER_BLOCK);
                    cur_block.push_back(record);
                    total_len = bam_len + 4;
                }
            }
        }

        if (!cur_block.empty()) {
            cur_group.push_back(std::move(cur_block));
        }

        if (cur_group.size() > (size_t)NB) {
            fprintf(stderr, "FATAL: bam2bam produced %zu output blocks from one input group\n", cur_group.size());
            abort();
        }

        for (auto blk : cur_blocks) read->backBlock(blk);
        t_pack += GetTime() - pack_t0;

        group_count++;
        if (cur_group.empty()) {
            flush_pending();
            do_read_group(next_blocks);
            continue;
        }

        do_compress(cur_group, next_blocks);
    }

    flush_pending();

    for (int b = 0; b < NB; ++b) {
        if (paras[b].bam_lens) {
            aligned_free_custom((unsigned char*)paras[b].bam_lens);
            paras[b].bam_lens = nullptr;
        }
    }

    printf("FusedBamToBam finished. in_blocks=%lld groups=%lld total_records=%lld kept_records=%lld dropped_records=%lld bgzf_blocks=%lld cost=%.3f s\n",
           input_blocks, group_count, total_records, kept_records, dropped_records, bgzf_blocks, GetTime() - t0);
    printf("  decomp_filter_slave=%.3f  pack=%.3f  compress_slave=%.3f  read=%.3f  write=%.3f\n",
           t_decomp_filter, t_pack, t_compress, t_read, t_write);
}

int SwBam::FusedBamToBamChecked(BamRead *read, BamComplete *complete, BamWriteComplete *write_complete,
                                MemReader &reader, MemWriter &mem_writer, const BamFilterOptions &filter,
                                BoundsCheckError *bounds_error) {
    double t0 = GetTime();
    const int NB = 64;

    CheckedBam2BamPara paras[NB];
    for (int b = 0; b < NB; ++b) {
        paras[b].output_records = complete->getResultBuf(b).data();
        paras[b].bam_lens = (uint32_t *)aligned_alloc_custom(64, MAX_RECORDS_PER_BLOCK * sizeof(uint32_t));
        if (!paras[b].bam_lens) {
            SetGenericBoundsError(bounds_error, "bam2bam", "alloc_bam_lens", errno, -1, -1, -1, b);
            for (int j = 0; j < b; ++j) {
                if (paras[j].bam_lens) aligned_free_custom((unsigned char*)paras[j].bam_lens);
            }
            return -1;
        }
        paras[b].status = -1;
        paras[b].record_index = -1;
        paras[b].actual_value = 0;
        paras[b].limit_value = 0;
        paras[b].limit_id = BOUNDS_LIMIT_NONE;
        paras[b].filter = filter;
    }

    Comp_Para comp_buf_A[NB], comp_buf_B[NB];
    Comp_Para *comp_active = comp_buf_A;
    Comp_Para *comp_pending = comp_buf_B;
    std::vector<std::vector<bam1_t*>> pending_group;
    bool has_pending = false;
    long long bgzf_blocks = 0;

    auto release_pending_outputs = [&](Comp_Para *buf) {
        for (int k = 0; k < NB; ++k) {
            if (buf[k].output_block) {
                write_complete->backBlock(buf[k].output_block);
                buf[k].output_block = nullptr;
            }
        }
    };

    auto do_read_group = [&](std::vector<bam_block*> &blocks) {
        blocks.clear();
        blocks.reserve(NB);
        for (int b = 0; b < NB; ++b) {
            bam_block *blk = read->getEmpty();
            int ret = mem_read_block(reader.base, reader.size, reader.pos, blk);
            if (ret < 0 || blk->length == 28) {
                read->backBlock(blk);
                break;
            }
            blk->block_id = b;
            blk->pos = 0;
            blocks.push_back(blk);
        }
    };

    auto flush_pending = [&]() {
        if (!has_pending) return;
        for (int k = 0; k < NB; ++k) {
            if (comp_pending[k].status == 0 && comp_pending[k].output_block) {
                write_block_to_mem(mem_writer, comp_pending[k].output_block);
                write_complete->backBlock(comp_pending[k].output_block);
                comp_pending[k].output_block = nullptr;
                bgzf_blocks++;
            }
        }
        pending_group.clear();
        has_pending = false;
    };

    auto do_compress = [&](std::vector<std::vector<bam1_t*>> &group,
                           std::vector<bam_block*> &next_blocks) -> int {
        for (int k = (int)group.size(); k < NB; ++k)
            group.push_back({});

        for (int k = 0; k < NB; ++k) {
            if (group[k].empty()) {
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
            comp_active[k].input_records = group[k].data();
            comp_active[k].n_records = (int)group[k].size();
            comp_active[k].output_block = write_complete->getEmpty();
            comp_active[k].un_comp_block = write_complete->getBuffer(k);
            comp_active[k].output_size = 0;
            comp_active[k].status = 0;
        }

        __real_athread_spawn((void*)slave_compressfunc, comp_active, 1);
        flush_pending();
        do_read_group(next_blocks);
        athread_join();

        for (int k = 0; k < NB; ++k) {
            if (!group[k].empty() && comp_active[k].status != 0) {
                SetGenericBoundsError(bounds_error, "bam2bam", "compress", comp_active[k].status,
                                      k, -1, -1, k);
                release_pending_outputs(comp_active);
                return -1;
            }
        }

        pending_group = std::move(group);
        std::swap(comp_active, comp_pending);
        has_pending = true;
        group.clear();
        group.reserve(NB);
        return 0;
    };

    std::vector<bam_block*> next_blocks;
    next_blocks.reserve(NB);
    do_read_group(next_blocks);

    long long total_records = 0;
    long long kept_records = 0;
    long long dropped_records = 0;
    long long input_blocks = 0;
    long long group_count = 0;

    while (!next_blocks.empty()) {
        std::vector<bam_block*> cur_blocks = std::move(next_blocks);
        next_blocks.clear();
        next_blocks.reserve(NB);

        const int n_blocks = (int)cur_blocks.size();
        input_blocks += n_blocks;

        for (int b = 0; b < NB; ++b) {
            paras[b].filter = filter;
            paras[b].block_id = b;
            paras[b].input_block = (b < n_blocks) ? cur_blocks[b] : nullptr;
            paras[b].un_comp_block = (b < n_blocks) ? complete->getBuffer(b) : nullptr;
            paras[b].n_total_records = 0;
            paras[b].n_kept_records = 0;
            paras[b].status = (b < n_blocks) ? 0 : -1;
            paras[b].record_index = -1;
            paras[b].actual_value = 0;
            paras[b].limit_value = 0;
            paras[b].limit_id = BOUNDS_LIMIT_NONE;
        }

        __real_athread_spawn((void*)slave_decompress_filter_checked, paras, 1);
        athread_join();

        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].status != 0) {
                SetBoundsErrorFromCheckedBam2BamPara(bounds_error, "bam2bam", "decompress_filter",
                                                     paras[b], b);
                for (auto blk : cur_blocks) read->backBlock(blk);
                release_pending_outputs(comp_active);
                release_pending_outputs(comp_pending);
                for (int j = 0; j < NB; ++j) {
                    aligned_free_custom((unsigned char*)paras[j].bam_lens);
                    paras[j].bam_lens = nullptr;
                }
                return -1;
            }
        }

        std::vector<std::vector<bam1_t*>> cur_group;
        cur_group.reserve(NB);
        std::vector<bam1_t*> cur_block;
        cur_block.reserve(MAX_RECORDS_PER_BLOCK);
        uint32_t total_len = 0;

        for (int b = 0; b < n_blocks; ++b) {
            total_records += paras[b].n_total_records;
            kept_records += paras[b].n_kept_records;
            dropped_records += paras[b].n_total_records - paras[b].n_kept_records;

            for (int r = 0; r < paras[b].n_kept_records; ++r) {
                bam1_t *record = paras[b].output_records[r];
                uint32_t bam_len = paras[b].bam_lens[r];
                if ((uint32_t)(bam_len + 4) > BGZF_BLOCK_SIZE) {
                    SetBoundsError(bounds_error, "bam2bam", "pack",
                                   BOUNDS_LIMIT_BGZF_RECORD_SIZE, BGZF_BLOCK_SIZE, bam_len + 4,
                                   paras[b].block_id, -1, r, b);
                    for (auto blk : cur_blocks) read->backBlock(blk);
                    release_pending_outputs(comp_active);
                    release_pending_outputs(comp_pending);
                    for (int j = 0; j < NB; ++j) {
                        aligned_free_custom((unsigned char*)paras[j].bam_lens);
                        paras[j].bam_lens = nullptr;
                    }
                    return -1;
                }

                if (bam_len + 4 + total_len <= BGZF_BLOCK_SIZE) {
                    cur_block.push_back(record);
                    total_len += bam_len + 4;
                } else {
                    cur_group.push_back(std::move(cur_block));
                    cur_block.clear();
                    cur_block.reserve(MAX_RECORDS_PER_BLOCK);
                    cur_block.push_back(record);
                    total_len = bam_len + 4;
                }
            }
        }

        if (!cur_block.empty()) cur_group.push_back(std::move(cur_block));
        for (auto blk : cur_blocks) read->backBlock(blk);

        group_count++;
        if (cur_group.empty()) {
            flush_pending();
            do_read_group(next_blocks);
            continue;
        }
        if (do_compress(cur_group, next_blocks) != 0) {
            for (int j = 0; j < NB; ++j) {
                aligned_free_custom((unsigned char*)paras[j].bam_lens);
                paras[j].bam_lens = nullptr;
            }
            return -1;
        }
    }

    flush_pending();

    for (int b = 0; b < NB; ++b) {
        if (paras[b].bam_lens) {
            aligned_free_custom((unsigned char*)paras[b].bam_lens);
            paras[b].bam_lens = nullptr;
        }
    }

    printf("FusedBamToBamChecked finished. in_blocks=%lld groups=%lld total_records=%lld kept_records=%lld dropped_records=%lld bgzf_blocks=%lld cost=%.3f s\n",
           input_blocks, group_count, total_records, kept_records, dropped_records, bgzf_blocks, GetTime() - t0);
    return 0;
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
    double t_copy_count = 0, t_parse = 0, t_compress = 0, t_pack = 0 , t_write = 0;
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
        double t0 = GetTime();

        if (!has_pending) return;
        for (int k = 0; k < 64; k++) {
            if (comp_pending[k].status == 0 && comp_pending[k].output_block) {

                // double t0 = GetTime();
                write_block_to_mem(mem_writer, comp_pending[k].output_block);
                // t_write += GetTime() - t0;

                complete->backBlock(comp_pending[k].output_block);
                bgzf_nums++;
            }
            if (comp_pending[k].input_records && comp_pending[k].n_records > 0) {
                write->backBamBatch(comp_pending[k].input_records, comp_pending[k].n_records);
            }
        }
        pending_group.clear();
        has_pending = false;

        t_write += GetTime() - t0;
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
        flush_pending();  
        athread_join();
        t_compress += GetTime() - ts;
        // flush_pending();  

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
        flush_pending();   
        athread_join();
        t_copy_count += GetTime() - ts;
        // flush_pending();   

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
    printf("  copy_count_slave=%lf  parse_slave=%lf  pack(+compress)=%lf  compress_slave=%lf  write=%lf\n",
           t_copy_count, t_parse, t_pack, t_compress, t_write);
}

int SwBam::FusedSamToBamChecked(BamWriteComplete *complete, sam_hdr_t *h,
                                MemReader reader, MemWriter &mem_writer, BoundsCheckError *bounds_error) {
    double t1 = GetTime();

    std::vector<bam1_t*> cur_block;
    std::vector<std::vector<bam1_t*>> cur_group;
    cur_block.reserve(MAX_RECORDS_PER_BLOCK);
    cur_group.reserve(64);

    uint32_t bam_len = 0;
    uint32_t total_len = 0;
    long long block_nums = 0;
    long long group_nums = 0;
    long long bam1_t_nums = 0;
    long long bgzf_nums = 0;
    long long cur_group_records = 0;
    long long pending_records = 0;

    Comp_Para comp_buf_A[64], comp_buf_B[64];
    Comp_Para *comp_active = comp_buf_A;
    Comp_Para *comp_pending = comp_buf_B;
    std::vector<std::vector<bam1_t*>> pending_group;
    bool has_pending = false;

    auto destroy_group_records = [](std::vector<std::vector<bam1_t*>> &group) {
        for (auto &block : group) {
            for (bam1_t *b : block) bam_destroy1(b);
        }
        group.clear();
    };
    auto destroy_block_records = [](std::vector<bam1_t*> &block) {
        for (bam1_t *b : block) bam_destroy1(b);
        block.clear();
    };
    auto release_outputs = [&](Comp_Para *buf) {
        for (int i = 0; i < 64; ++i) {
            if (buf[i].output_block) {
                complete->backBlock(buf[i].output_block);
                buf[i].output_block = nullptr;
            }
        }
    };
    auto cleanup_state = [&]() {
        destroy_block_records(cur_block);
        destroy_group_records(cur_group);
        destroy_group_records(pending_group);
        release_outputs(comp_active);
        release_outputs(comp_pending);
    };

    auto flush_pending = [&]() {
        if (!has_pending) return;
        for (int k = 0; k < 64; ++k) {
            if (comp_pending[k].status == 0 && comp_pending[k].output_block) {
                write_block_to_mem(mem_writer, comp_pending[k].output_block);
                complete->backBlock(comp_pending[k].output_block);
                comp_pending[k].output_block = nullptr;
                bgzf_nums++;
            }
        }
        destroy_group_records(pending_group);
        pending_records = 0;
        has_pending = false;
    };

    auto do_compress = [&]() -> int {
        for (int k = (int)cur_group.size(); k < 64; ++k)
            cur_group.push_back({});

        for (int k = 0; k < 64; ++k) {
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

        __real_athread_spawn((void*)slave_compressfunc, comp_active, 1);
        flush_pending();
        athread_join();

        for (int k = 0; k < 64; ++k) {
            if (!cur_group[k].empty() && comp_active[k].status != 0) {
                SetGenericBoundsError(bounds_error, "sam2bam", "compress", comp_active[k].status,
                                      k, -1, -1, k);
                release_outputs(comp_active);
                return -1;
            }
        }

        pending_group = std::move(cur_group);
        pending_records = cur_group_records;
        std::swap(comp_active, comp_pending);
        has_pending = true;
        cur_group.clear();
        cur_group.reserve(64);
        cur_group_records = 0;
        return 0;
    };

    std::vector<char> line_buf(MAX_SAM_LINE_SIZE + 1);

    while (reader.pos < reader.size) {
        struct ChunkPlan {
            size_t start;
            size_t end;
            int count;
        };
        std::vector<ChunkPlan> plans;
        plans.reserve(64);
        long long batch_record_count = 0;

        for (int chunk_id = 0; chunk_id < 64 && reader.pos < reader.size; ++chunk_id) {
            size_t start_pos = reader.pos;
            size_t end_pos = start_pos + SAM_CHUNK_SIZE;
            if (end_pos < reader.size) {
                while (end_pos < reader.size && reader.base[end_pos] != '\n') end_pos++;
                if (end_pos < reader.size && reader.base[end_pos] == '\n') end_pos++;
            } else {
                end_pos = reader.size;
            }

            size_t read_len = end_pos - start_pos;
            if (read_len + 1 > CHUNK_BUFFER_SIZE) {
                SetBoundsError(bounds_error, "sam2bam", "copy_count",
                               BOUNDS_LIMIT_CHUNK_BUFFER_SIZE, CHUNK_BUFFER_SIZE, read_len + 1,
                               -1, chunk_id, -1, -1);
                cleanup_state();
                return -1;
            }

            int count = 0;
            size_t pos = start_pos;
            while (pos < end_pos) {
                size_t line_start = pos;
                while (pos < end_pos && reader.base[pos] != '\n') pos++;
                size_t line_len = pos - line_start;
                if (line_len > 0 && reader.base[line_start + line_len - 1] == '\r') line_len--;
                if (pos < end_pos && reader.base[pos] == '\n') pos++;
                if (line_len == 0) continue;
                if (line_len > MAX_SAM_LINE_SIZE) {
                    SetBoundsError(bounds_error, "sam2bam", "copy_count",
                                   BOUNDS_LIMIT_MAX_SAM_LINE_SIZE, MAX_SAM_LINE_SIZE, line_len,
                                   -1, chunk_id, count, -1);
                    cleanup_state();
                    return -1;
                }
                count++;
                if (count > MAX_BAMS_PER_CHUNK) {
                    SetBoundsError(bounds_error, "sam2bam", "copy_count",
                                   BOUNDS_LIMIT_MAX_BAMS_PER_CHUNK, MAX_BAMS_PER_CHUNK, count,
                                   -1, chunk_id, count - 1, -1);
                    cleanup_state();
                    return -1;
                }
            }

            plans.push_back({start_pos, end_pos, count});
            batch_record_count += count;
            reader.pos = end_pos;
        }

        if (plans.empty()) break;

        long long estimated_in_flight = pending_records + cur_group_records + (long long)cur_block.size() + batch_record_count;
        if (estimated_in_flight > FUSED_SAM2BAM_BAM_POOL_SIZE) {
            SetBoundsError(bounds_error, "sam2bam", "pool_estimate",
                           BOUNDS_LIMIT_FUSED_SAM2BAM_BAM_POOL_SIZE, FUSED_SAM2BAM_BAM_POOL_SIZE,
                           estimated_in_flight, -1, -1, -1, -1);
            cleanup_state();
            return -1;
        }

        for (int chunk_id = 0; chunk_id < (int)plans.size(); ++chunk_id) {
            size_t pos = plans[chunk_id].start;
            int record_index = 0;
            while (pos < plans[chunk_id].end) {
                size_t line_start = pos;
                while (pos < plans[chunk_id].end && reader.base[pos] != '\n') pos++;
                size_t line_len = pos - line_start;
                if (line_len > 0 && reader.base[line_start + line_len - 1] == '\r') line_len--;
                if (pos < plans[chunk_id].end && reader.base[pos] == '\n') pos++;
                if (line_len == 0) continue;

                memcpy(line_buf.data(), reader.base + line_start, line_len);
                line_buf[line_len] = '\0';
                kstring_t ks = {line_len, line_len + 1, line_buf.data()};
                bam1_t *b = bam_init1();
                if (!b) {
                    SetGenericBoundsError(bounds_error, "sam2bam", "sam_parse_alloc", errno,
                                          -1, chunk_id, record_index, -1);
                    cleanup_state();
                    return -1;
                }
                int ret = sam_parse1(&ks, h, b);
                if (ret < 0) {
                    bam_destroy1(b);
                    SetGenericBoundsError(bounds_error, "sam2bam", "sam_parse", ret,
                                          -1, chunk_id, record_index, -1);
                    cleanup_state();
                    return -1;
                }
                if (b->m_data > INIT_DATA_SIZE || b->l_data > (int)INIT_DATA_SIZE) {
                    long long actual = std::max<long long>(b->m_data, b->l_data);
                    bam_destroy1(b);
                    SetBoundsError(bounds_error, "sam2bam", "sam_parse",
                                   BOUNDS_LIMIT_INIT_DATA_SIZE, INIT_DATA_SIZE, actual,
                                   -1, chunk_id, record_index, -1);
                    cleanup_state();
                    return -1;
                }

                bam_len = (uint32_t)(b->l_data - b->core.l_extranul + 32);
                if ((uint32_t)(bam_len + 4) > BGZF_BLOCK_SIZE) {
                    bam_destroy1(b);
                    SetBoundsError(bounds_error, "sam2bam", "pack",
                                   BOUNDS_LIMIT_BGZF_RECORD_SIZE, BGZF_BLOCK_SIZE, bam_len + 4,
                                   -1, chunk_id, record_index, -1);
                    cleanup_state();
                    return -1;
                }

                bam1_t_nums++;
                if (bam_len + 4 + total_len <= BGZF_BLOCK_SIZE) {
                    cur_block.push_back(b);
                    total_len += bam_len + 4;
                } else {
                    block_nums++;
                    cur_group_records += cur_block.size();
                    cur_group.push_back(std::move(cur_block));
                    cur_block.clear();
                    cur_block.reserve(MAX_RECORDS_PER_BLOCK);
                    cur_block.push_back(b);
                    total_len = bam_len + 4;
                }

                if ((int)cur_group.size() >= 64) {
                    group_nums++;
                    if (do_compress() != 0) {
                        cleanup_state();
                        return -1;
                    }
                }
                record_index++;
            }
        }
    }

    if (!cur_block.empty()) {
        block_nums++;
        cur_group_records += cur_block.size();
        cur_group.push_back(std::move(cur_block));
    }
    if (!cur_group.empty()) {
        group_nums++;
        if (do_compress() != 0) {
            cleanup_state();
            return -1;
        }
    }
    flush_pending();

    printf("FusedSamToBamChecked finished. bam1_t=%lld, blocks=%lld, groups=%lld, bgzf=%lld, cost=%lf\n",
           bam1_t_nums, block_nums, group_nums, bgzf_nums, GetTime() - t1);
    return 0;
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

int SwBam::ProcessSwBam() {
    double t0 = GetTime();
    double t_header = 0, t_total = 0;
    MemReader reader = {};
    MemWriter mem_writer = {};
    BoundsCheckError bounds_error = {};
    size_t sam_size = 0;
    size_t bam_size = 0;
    char *sam_mem = nullptr;
    char *bam_mem = nullptr;
    bool ran_body = false;
    int exit_code = 1;
    int input_format = -1;
    int output_format = -1;
    bool bam_to_bam = false;
    bool bam_to_sam = false;
    bool sam_to_bam = false;
    bool filter_requested = HasBamFilterRequest(cmd_info_);
    BamFilterOptions bam_filter = BuildBamFilterOptions(cmd_info_);
    ClearBoundsError(&bounds_error);

    // in文件打开
    sin = sam_open(cmd_info_->in_file_name_.c_str(), "r");
    if (sin == NULL) {
        fprintf(stderr, "Error opening input file %s\n", cmd_info_->in_file_name_.c_str());
        return exit_code;
    }

    // out文件打开
    const char *out_mode = get_sam_open_mode(cmd_info_->out_file_name_);
    sout = sam_open(cmd_info_->out_file_name_.c_str(), out_mode);
    if (sout == NULL) {
        fprintf(stderr, "Error opening output file %s\n", cmd_info_->out_file_name_.c_str());
        goto cleanup;
    }
    printf("open the files cost %lf---\n", GetTime() - t0);

    t0 = GetTime();
    hdr = sam_hdr_read(sin);
    if (hdr == NULL) {
        fprintf(stderr, "Error reading header from input file %s\n", cmd_info_->in_file_name_.c_str());
        goto cleanup;
    }

    input_format = sin->format.format;
    output_format = sout->format.format;
    if (IsBamLikeFormat(input_format)) input_format = bam;
    else if (IsSamLikeFormat(input_format)) input_format = sam;
    if (IsBamLikeFormat(output_format)) output_format = bam;
    else if (IsSamLikeFormat(output_format)) output_format = sam;
    bam_to_bam = IsBamLikeFormat(input_format) && IsBamLikeFormat(output_format);
    bam_to_sam = IsBamLikeFormat(input_format) && IsSamLikeFormat(output_format);
    sam_to_bam = IsSamLikeFormat(input_format) && IsBamLikeFormat(output_format);

    if (!bam_to_bam && !bam_to_sam && !sam_to_bam) {
        fprintf(stderr, "Unsupported conversion: input format %d -> output format %d\n", input_format, output_format);
        goto cleanup;
    }

    if (filter_requested && !bam_to_bam) {
        fprintf(stderr, "ERROR: BAM filtering options are only supported for BAM -> BAM.\n");
        goto cleanup;
    }

    if (cmd_info_->min_mapq_ > cmd_info_->max_mapq_ && cmd_info_->max_mapq_ >= 0) {
        fprintf(stderr, "ERROR: --min-mapq cannot be greater than --max-mapq.\n");
        goto cleanup;
    }

    if (cmd_info_->min_read_len_ > cmd_info_->max_read_len_ && cmd_info_->max_read_len_ >= 0) {
        fprintf(stderr, "ERROR: --min-read-len cannot be greater than --max-read-len.\n");
        goto cleanup;
    }

    if (!cmd_info_->ref_name_.empty()) {
        int ref_tid = sam_hdr_name2tid(hdr, cmd_info_->ref_name_.c_str());
        if (ref_tid < 0) {
            fprintf(stderr, "ERROR: reference name '%s' does not exist in the BAM header.\n",
                    cmd_info_->ref_name_.c_str());
            goto cleanup;
        }
        bam_filter.ref_tid = ref_tid;
    }

    //这里不做验证，直接放行
    // if (bam_to_bam) {
    //     BamCrossBlockStats cross_stats = {0, 0};
    //     int cross_ret = check_bam_cross_block_ex(cmd_info_->in_file_name_.c_str(), &cross_stats, false);
    //     if (cross_ret < 0) {
    //         fprintf(stderr, "ERROR: failed to pre-check cross-block records for BAM -> BAM.\n");
    //         goto cleanup;
    //     }
    //     if (cross_ret > 0) {
    //         fprintf(stderr,
    //                 "ERROR: input BAM contains %lld cross-block records; BAM -> BAM fused fast-path is disabled for this file.\n",
    //                 cross_stats.cross_block_records);
    //         goto cleanup;
    //     }
    // }

    if (sam_hdr_write(sout, hdr) != 0) {
        fprintf(stderr, "Error writing header to output file %s\n", cmd_info_->out_file_name_.c_str());
        goto cleanup;
    }

    t_header = GetTime() - t0;
    t_total += t_header;
    printf("Complete the head cost %lf\n", t_header);

    #define USE_MEMORY
    #ifdef USE_MEMORY
    t0 = GetTime();
    switch (input_format) {
        case bam: {
            printf("MEMORY BAM\n");

            FILE *f = fopen(cmd_info_->in_file_name_.c_str(), "rb");
            if (!f) {
                perror("fopen");
                goto cleanup;
            }
            fseek(f, 0, SEEK_END);
            bam_size = ftell(f);
            fseek(f, 0, SEEK_SET);
            bam_mem = (char*)malloc(bam_size);
            if (!bam_mem) {
                fclose(f);
                fprintf(stderr, "Failed to allocate memory for BAM file\n");
                goto cleanup;
            }
            if (fread(bam_mem, 1, bam_size, f) != bam_size) {
                fclose(f);
                fprintf(stderr, "Failed to read BAM file\n");
                goto cleanup;
            }
            fclose(f);
            printf("Loaded BAM into memory: %.2f MB\n", bam_size / 1024.0 / 1024.0);

            reader.base = bam_mem;
            reader.size = bam_size;
            reader.pos = sin->fp.bgzf->block_address;

            // size_t GB = 1024LL * 1024LL * 1024LL;
            // size_t sam_estimate = 1 * GB;
            size_t writer_estimate = bam_to_sam ? bam_size * 5 : bam_size;
            if (writer_estimate == 0) writer_estimate = 64 * 1024 * 1024;
            init_mem_writer(mem_writer, writer_estimate);
            printf("Initialized %s memory writer: %.2f MB\n",
                   output_format == bam ? "BAM" : "SAM",
                   writer_estimate / 1024.0 / 1024.0);
            break;
        }

        case sam: {
            printf("MEMORY SAM\n");

            FILE *f = fopen(cmd_info_->in_file_name_.c_str(), "r");
            if (!f) {
                perror("fopen");
                goto cleanup;
            }
            fseek(f, 0, SEEK_END);
            sam_size = ftell(f);
            fseek(f, 0, SEEK_SET);
            sam_mem = (char*)malloc(sam_size);
            if (!sam_mem) {
                fclose(f);
                fprintf(stderr, "Failed to allocate memory for SAM file\n");
                goto cleanup;
            }
            if (fread(sam_mem, 1, sam_size, f) != sam_size) {
                fclose(f);
                fprintf(stderr, "Failed to read SAM file\n");
                goto cleanup;
            }
            fclose(f);
            printf("Loaded SAM into memory: %.2f MB\n", sam_size / 1024.0 / 1024.0);

            reader.base = sam_mem;
            reader.size = sam_size;
            reader.pos = 0;

            kstring_t tmp;
            tmp.l = 0;
            tmp.m = MAX_SAM_LINE_SIZE;
            tmp.s = (char*)malloc(MAX_SAM_LINE_SIZE);
            if (!tmp.s) {
                fprintf(stderr, "Failed to allocate temporary SAM header buffer\n");
                goto cleanup;
            }

            while (true) {
                size_t old_pos = reader.pos;
                if (!mem_getline(reader, &tmp))
                    break;
                if (tmp.l == 0)
                    continue;
                if (tmp.s[0] != '@') {
                    reader.pos = old_pos;
                    break;
                }
            }
            free(tmp.s);

            size_t bam_estimate = sam_size / 2;
            if (bam_estimate == 0) bam_estimate = 64 * 1024 * 1024;
            init_mem_writer(mem_writer, bam_estimate);
            printf("Initialized BAM memory writer: %.2f MB\n", bam_estimate / 1024.0 / 1024.0);
            break;
        }

        default:
            fprintf(stderr, "Unknown input file format\n");
            goto cleanup;
    }
    printf("Complete the memory cost %lf\n", GetTime() - t0);
    #endif

    {
        double tbody = GetTime();

        switch (input_format) {
            case bam: {
                t0 = GetTime();

                if (output_format == sam) {
                    #define USE_FUSED_BAM2SAM
                    #if defined(USE_MEMORY) && defined(USE_FUSED_BAM2SAM)
                    {
                        printf("Enable FUSED BAM2SAM single-thread mode (USE_MEMORY + USE_FUSED_BAM2SAM)!!!\n");
                        //BamRead(x , y) x是blocks的内存池大小，y是打包之后队列的大小
                        //BamComplete(x) x是bam1_t的内存池大小
                        // BamRead(128,1): 128块压缩block池（两轮各64），queue不用
                        // BamComplete(1): 仅需 buffer_pool_（解压缓冲）和 result_pool_（bam1_t）
                        read = new BamRead(128, 1);
                        complete = new BamComplete(1);
                        printf("Complete the queue initialization cost %lf\n", GetTime() - t0);
                        if (cmd_info_->validate_bounds_) {
                            printf("Enable CHECKED BAM2SAM mode (--validate-bounds)!!!\n");
                            if (FusedBamToSamChecked(read, complete, hdr, reader, mem_writer, &bounds_error) != 0) {
                                PrintBoundsError(bounds_error);
                                goto cleanup;
                            }
                        } else {
                            FusedBamToSam(read, complete, hdr, reader, mem_writer);
                        }
                    }
                    #elif defined(USE_MEMORY)
                    {
                        printf("Enable 3-thread BAM2SAM pipeline (USE_MEMORY)!!!\n");

                        read = new BamRead(320, 5);
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
                    #else
                    {
                        printf("Enable 3-thread BAM2SAM pipeline (disk mode)!!!\n");
                        read = new BamRead(320, 5);
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
                } else if (output_format == bam) {
                    #define USE_FUSED_BAM2BAM
                    #if defined(USE_MEMORY) && defined(USE_FUSED_BAM2BAM)
                    {
                        printf("Enable FUSED BAM2BAM single-thread mode (USE_MEMORY + USE_FUSED_BAM2BAM)!!!\n");
                        read = new BamRead(128, 1);
                        complete = new BamComplete(1);
                        writeComplete = new BamWriteComplete(130);
                        printf("Complete the queue initialization cost %lf\n", GetTime() - t0);
                        if (cmd_info_->validate_bounds_) {
                            printf("Enable CHECKED BAM2BAM mode (--validate-bounds)!!!\n");
                            if (FusedBamToBamChecked(read, complete, writeComplete, reader, mem_writer,
                                                     bam_filter, &bounds_error) != 0) {
                                PrintBoundsError(bounds_error);
                                goto cleanup;
                            }
                        } else {
                            FusedBamToBam(read, complete, writeComplete, reader, mem_writer, bam_filter);
                        }
                    }
                    #else
                    {
                        fprintf(stderr, "BAM -> BAM is only supported in fused memory mode in this version.\n");
                        goto cleanup;
                    }
                    #endif
                }

                break;
            }

            case sam: {
                t0 = GetTime();

                #define USE_FUSED_SAM2BAM
                #if defined(USE_MEMORY) && defined(USE_FUSED_SAM2BAM)
                {
                    printf("Enable FUSED single-thread mode (USE_MEMORY + USE_FUSED_SAM2BAM)!!!\n");
                    //BamWrite(x , y) x是bam1_t的内存池大小 , y是打包的block的size
                    //BamWriteComplete(x) x是bam_block的内存池大小
                    writeComplete = new BamWriteComplete(130);
                    printf("Complete the queue initialization cost %lf\n", GetTime() - t0);
                    if (cmd_info_->validate_bounds_) {
                        printf("Enable CHECKED SAM2BAM mode (--validate-bounds)!!!\n");
                        if (FusedSamToBamChecked(writeComplete, hdr, reader, mem_writer, &bounds_error) != 0) {
                            PrintBoundsError(bounds_error);
                            goto cleanup;
                        }
                    } else {
                        write = new BamWrite(FUSED_SAM2BAM_BAM_POOL_SIZE, 1);
                        FusedSamToBam(write, writeComplete, hdr, reader, mem_writer);
                    }
                }
                #elif defined(USE_MEMORY)
                {
                    printf("Enable the SAM_PARSE_PARALLEL (3-thread, read and write from memory)!!!\n");

                    write = new BamWrite(655360 ,640);
                    writeComplete = new BamWriteComplete(130);
                    printf("Complete the queue initialization cost %lf\n", GetTime() - t0);
                    
                    thread producer2(bind(&SwBam::ProducerSwBamTask2_parallel_memory_OP, this, write, hdr, reader));
                    thread consumer2(bind(&SwBam::ConsumerSwBamTask2, this, write , writeComplete));

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
                #else
                {
                    printf("Enable the SAM_PARSE_PARALLEL (3-thread, read and write from disk)!!!\n");

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
                goto cleanup;
        }

        ran_body = true;
        exit_code = 0;
        t_total += GetTime() - tbody;
        printf("Complete the body cost %lf\n", GetTime() - tbody);
        printf("Complete the total cost %lf---\n", t_total);
    }

    #ifdef USE_MEMORY
    #define DUMP_MEM
    #ifdef DUMP_MEM
    if (ran_body && exit_code == 0) {
        double dump_t0 = GetTime();
        if (sout->format.format == bam) {
            if (bgzf_flush(sout->fp.bgzf) != 0) {
                fprintf(stderr, "Error flushing BAM output header/body stream\n");
            }
            if (hwrite(sout->fp.bgzf->fp, mem_writer.data, mem_writer.size) != mem_writer.size) {
                fprintf(stderr, "Error dumping memory to BAM output file\n");
            }
        } else if (sout->format.format == sam) {
            if (hflush(sout->fp.hfile) != 0) {
                fprintf(stderr, "Error flushing SAM output header/body stream\n");
            }
            if (hwrite(sout->fp.hfile, mem_writer.data, mem_writer.size) != mem_writer.size) {
                fprintf(stderr, "Error dumping memory to SAM output file\n");
            }
        }
        printf("Dump memory to output file cost %lf\n", GetTime() - dump_t0);
    }
    #endif
    #endif

cleanup:
    if (mem_writer.data) {
        free(mem_writer.data);
        mem_writer.data = nullptr;
    }
    if (sam_mem) {
        free(sam_mem);
        sam_mem = nullptr;
    }
    if (bam_mem) {
        free(bam_mem);
        bam_mem = nullptr;
    }

    t0 = GetTime();
    if (hdr) {
        sam_hdr_destroy(hdr);
        hdr = nullptr;
    }
    if (sout) {
        int ret = hts_close(sout);
        if (ret < 0) {
            fprintf(stderr, "Error closing output.\n");
        }
        sout = nullptr;
    }
    if (sin) {
        int ret = hts_close(sin);
        if (ret < 0) {
            fprintf(stderr, "Error closing input.\n");
        }
        sin = nullptr;
    }
    printf("close the files cost %lf---\n", GetTime() - t0);
    return exit_code;
}
