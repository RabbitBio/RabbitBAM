#include "swbam_mpi.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

extern "C" {
    void slave_mpi_decompress_bam2bam_passthrough();
    void slave_sam_format();
}

namespace {

struct MpiBamToSamBlockSet {
    bam_block *blocks;
    unsigned char *data;
    int n;
};

struct MpiBamToSamRecordSet {
    bam1_t *records;
    unsigned char *data;
    bam1_t **ptrs;
    uint32_t *bam_lens;
    int capacity;
};

int MpiAllocateBamToSamBlockSet(MpiBamToSamBlockSet *set, int n) {
    set->n = n;
    set->blocks = (bam_block *)aligned_alloc_custom(64, (size_t)n * sizeof(bam_block));
    set->data = aligned_alloc_custom(64, (size_t)n * BGZF_MAX_BLOCK_SIZE);
    if (!set->blocks || !set->data) return -1;
    memset(set->blocks, 0, (size_t)n * sizeof(bam_block));
    for (int i = 0; i < n; ++i) {
        set->blocks[i].data = set->data + (size_t)i * BGZF_MAX_BLOCK_SIZE;
        set->blocks[i].block_id = i;
    }
    return 0;
}

void MpiFreeBamToSamBlockSet(MpiBamToSamBlockSet *set) {
    if (set->blocks) aligned_free_custom((unsigned char *)set->blocks);
    if (set->data) aligned_free_custom(set->data);
    set->blocks = nullptr;
    set->data = nullptr;
    set->n = 0;
}

int MpiAllocateBamToSamRecordSet(MpiBamToSamRecordSet *set, int total_records) {
    set->capacity = total_records;
    set->records = (bam1_t *)aligned_alloc_custom(64, (size_t)total_records * sizeof(bam1_t));
    set->data = aligned_alloc_custom(64, (size_t)total_records * INIT_DATA_SIZE);
    set->ptrs = (bam1_t **)aligned_alloc_custom(64, (size_t)total_records * sizeof(bam1_t *));
    set->bam_lens = (uint32_t *)aligned_alloc_custom(64, (size_t)total_records * sizeof(uint32_t));
    if (!set->records || !set->data || !set->ptrs || !set->bam_lens) return -1;
    memset(set->records, 0, (size_t)total_records * sizeof(bam1_t));
    memset(set->bam_lens, 0, (size_t)total_records * sizeof(uint32_t));
    for (int i = 0; i < total_records; ++i) {
        set->records[i].data = set->data + (size_t)i * INIT_DATA_SIZE;
        set->records[i].m_data = INIT_DATA_SIZE;
        set->records[i].l_data = 0;
        set->records[i].mempolicy = BAM_USER_OWNS_DATA;
        set->ptrs[i] = &set->records[i];
    }
    return 0;
}

void MpiFreeBamToSamRecordSet(MpiBamToSamRecordSet *set) {
    if (set->records) aligned_free_custom((unsigned char *)set->records);
    if (set->data) aligned_free_custom(set->data);
    if (set->ptrs) aligned_free_custom((unsigned char *)set->ptrs);
    if (set->bam_lens) aligned_free_custom((unsigned char *)set->bam_lens);
    set->records = nullptr;
    set->data = nullptr;
    set->ptrs = nullptr;
    set->bam_lens = nullptr;
    set->capacity = 0;
}

int MpiMemReadBamToSamBlock(char *base, size_t size, size_t &pos, bam_block *block) {
    if (pos >= size) return -1;
    if (pos + BLOCK_HEADER_LENGTH > size) return -1;
    uint16_t bsize = *(uint16_t *)(base + pos + 16);
    bsize += 1;
    if (pos + bsize > size) return -1;
    memcpy(block->data, base + pos, bsize);
    block->length = bsize;
    block->pos = 0;
    block->errcode = 0;
    block->block_address = (int64_t)pos;
    pos += bsize;
    return bsize;
}

int MpiInitSamFormatBatchMPI(SamFormatBatch *batch, sam_hdr_t *hdr) {
    memset(batch, 0, sizeof(SamFormatBatch));
    batch->hdr = hdr;
    batch->count = 0;
    for (int i = 0; i < 64; ++i) {
        kstring_t *ks = &batch->core_out_lines[i];
        ks->l = 0;
        ks->m = MAX_SAM_FORMAT_CORE_BUFFER_SIZE;
        ks->s = (char *)aligned_alloc_custom(64, MAX_SAM_FORMAT_CORE_BUFFER_SIZE);
        if (!ks->s) {
            for (int j = 0; j < i; ++j) {
                aligned_free_custom((unsigned char *)batch->core_out_lines[j].s);
                batch->core_out_lines[j].s = nullptr;
                batch->core_out_lines[j].l = batch->core_out_lines[j].m = 0;
            }
            return -1;
        }
    }
    return 0;
}

void MpiResetSamFormatBatchMPI(SamFormatBatch *batch, int count) {
    batch->count = count;
    for (int i = 0; i < 64; ++i) {
        batch->core_out_lines[i].l = 0;
    }
}

void MpiDestroySamFormatBatchMPI(SamFormatBatch *batch) {
    if (!batch) return;
    for (int i = 0; i < 64; ++i) {
        if (batch->core_out_lines[i].s) {
            aligned_free_custom((unsigned char *)batch->core_out_lines[i].s);
            batch->core_out_lines[i].s = nullptr;
        }
        batch->core_out_lines[i].l = batch->core_out_lines[i].m = 0;
    }
}

void MpiInitEmptyBamToSamPara(Bam2BamPara *para, int block_id) {
    para->block_id = block_id;
    para->input_block = nullptr;
    para->un_comp_block = nullptr;
    para->output_records = nullptr;
    para->record_base = nullptr;
    para->bam_lens = nullptr;
    para->filter.min_mapq = -1;
    para->filter.max_mapq = -1;
    para->filter.require_flag = 0;
    para->filter.exclude_flag = 0;
    para->filter.ref_tid = -2;
    para->filter.min_read_len = -1;
    para->filter.max_read_len = -1;
    para->n_total_records = 0;
    para->n_kept_records = 0;
    para->kept_total_len = 0;
    para->status = -1;
    para->decomp_alloc_cycles = 0;
    para->decomp_inflate_cycles = 0;
    para->decomp_crc_cycles = 0;
    para->decomp_parse_cycles = 0;
    para->decomp_total_cycles = 0;
}

void MpiAccumulateDecompDetail(const Bam2BamPara *paras,
                               int active_blocks,
                               double decomp_wall,
                               MpiBamToSamStats *stats) {
    if (!stats || active_blocks <= 0 || decomp_wall <= 0.0) return;

    const Bam2BamPara *critical = nullptr;
    uint64_t critical_total = 0;
    for (int b = 0; b < active_blocks; ++b) {
        if (paras[b].decomp_total_cycles >= critical_total) {
            critical_total = paras[b].decomp_total_cycles;
            critical = &paras[b];
        }
    }
    if (!critical || critical_total == 0) {
        stats->t_decomp_other += decomp_wall;
        return;
    }

    const double scale = decomp_wall / (double)critical_total;
    const double alloc_time = scale * (double)critical->decomp_alloc_cycles;
    const double inflate_time = scale * (double)critical->decomp_inflate_cycles;
    const double crc_time = scale * (double)critical->decomp_crc_cycles;
    const double parse_time = scale * (double)critical->decomp_parse_cycles;
    double other_time = decomp_wall - alloc_time - inflate_time - crc_time - parse_time;
    if (other_time < 0.0) other_time = 0.0;

    stats->t_decomp_alloc += alloc_time;
    stats->t_decomp_inflate += inflate_time;
    stats->t_decomp_crc += crc_time;
    stats->t_decomp_parse += parse_time;
    stats->t_decomp_other += other_time;
}

} // namespace

int FusedBamToSamMPI(MemReader &reader, MemWriter &mem_writer,
                     sam_hdr_t *hdr,
                     MpiBamToSamStats *stats) {
    const int NB = 64;
    MpiBamToSamStats local_stats = {};
    if (!stats) stats = &local_stats;
    double fused_t0 = GetTime();

    double alloc_t0 = GetTime();
    Bam2BamPara paras[NB];
    MpiBamToSamBlockSet input_blocks = {};
    MpiBamToSamBlockSet un_blocks = {};
    MpiBamToSamRecordSet record_set = {};
    SamFormatBatch *fmt_cur = nullptr;
    SamFormatBatch *fmt_prev = nullptr;
    bool has_pending_write = false;
    int ret_code = -1;

    for (int i = 0; i < NB; ++i) MpiInitEmptyBamToSamPara(&paras[i], i);

    fmt_cur = (SamFormatBatch *)aligned_alloc_custom(64, sizeof(SamFormatBatch));
    fmt_prev = (SamFormatBatch *)aligned_alloc_custom(64, sizeof(SamFormatBatch));
    if (fmt_cur) memset(fmt_cur, 0, sizeof(SamFormatBatch));
    if (fmt_prev) memset(fmt_prev, 0, sizeof(SamFormatBatch));
    if (!fmt_cur || !fmt_prev ||
        MpiAllocateBamToSamBlockSet(&input_blocks, NB) != 0 ||
        MpiAllocateBamToSamBlockSet(&un_blocks, NB) != 0 ||
        MpiAllocateBamToSamRecordSet(&record_set, NB * (int)MAX_RECORDS_PER_BLOCK) != 0 ||
        MpiInitSamFormatBatchMPI(fmt_cur, hdr) != 0 ||
        MpiInitSamFormatBatchMPI(fmt_prev, hdr) != 0) {
        fprintf(stderr, "ERROR: failed to allocate MPI bam2sam 1CG workspace.\n");
        MpiDestroySamFormatBatchMPI(fmt_cur);
        MpiDestroySamFormatBatchMPI(fmt_prev);
        if (fmt_cur) aligned_free_custom((unsigned char *)fmt_cur);
        if (fmt_prev) aligned_free_custom((unsigned char *)fmt_prev);
        MpiFreeBamToSamBlockSet(&input_blocks);
        MpiFreeBamToSamBlockSet(&un_blocks);
        MpiFreeBamToSamRecordSet(&record_set);
        stats->t_alloc_init += GetTime() - alloc_t0;
        stats->t_fused_total += GetTime() - fused_t0;
        return -1;
    }
    stats->t_alloc_init += GetTime() - alloc_t0;

    auto flush_pending_write = [&]() -> int {
        if (!has_pending_write) return 0;
        double write_t0 = GetTime();
        for (int i = 0; i < NB; ++i) {
            kstring_t *ks = &fmt_prev->core_out_lines[i];
            if (ks->l > 0 && MpiWriteBytesToMem(mem_writer, ks->s, ks->l) != 0) {
                fprintf(stderr, "ERROR: failed to append MPI SAM text to memory writer.\n");
                return -1;
            }
        }
        has_pending_write = false;
        stats->t_write += GetTime() - write_t0;
        return 0;
    };

    auto do_read_group = [&](int *n_blocks) -> int {
        double read_t0 = GetTime();
        int count = 0;
        for (int b = 0; b < NB; ++b) {
            bam_block *blk = &input_blocks.blocks[b];
            int ret = MpiMemReadBamToSamBlock(reader.base, reader.size, reader.pos, blk);
            if (ret < 0 || blk->length == 28) break;
            blk->block_id = b;
            blk->pos = 0;
            count++;
        }
        *n_blocks = count;
        stats->t_read += GetTime() - read_t0;
        return 0;
    };

    int next_n_blocks = 0;
    if (do_read_group(&next_n_blocks) != 0) goto cleanup;

    while (next_n_blocks > 0) {
        int n_blocks = next_n_blocks;
        next_n_blocks = 0;
        stats->input_blocks += n_blocks;

        for (int b = 0; b < NB; ++b) {
            MpiInitEmptyBamToSamPara(&paras[b], b);
            paras[b].output_records = record_set.ptrs + (size_t)b * MAX_RECORDS_PER_BLOCK;
            paras[b].record_base = record_set.records + (size_t)b * MAX_RECORDS_PER_BLOCK;
            paras[b].bam_lens = record_set.bam_lens + (size_t)b * MAX_RECORDS_PER_BLOCK;
            if (b < n_blocks) {
                paras[b].input_block = &input_blocks.blocks[b];
                paras[b].un_comp_block = &un_blocks.blocks[b];
                paras[b].status = 0;
            }
        }

        double decomp_t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_decompress_bam2bam_passthrough, paras, 1);
        #ifdef ENABLE_MASKING
        if (flush_pending_write() != 0) goto cleanup;
        #endif
        athread_join();
        double decomp_wall = GetTime() - decomp_t0;
        stats->t_decomp += decomp_wall;
        MpiAccumulateDecompDetail(paras, n_blocks, decomp_wall, stats);
        
        #ifndef ENABLE_MASKING
        if (flush_pending_write() != 0) goto cleanup;
        #endif

        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].status != 0) {
                if (paras[b].status == -3) {
                    fprintf(stderr, "ERROR: MPI bam2sam input block %d exceeds MAX_RECORDS_PER_BLOCK or contains an unsupported cross-block record.\n",
                            b);
                } else {
                    fprintf(stderr, "ERROR: MPI bam2sam decompress/parse failed on input block %d with status %d.\n",
                            b, paras[b].status);
                }
                goto cleanup;
            }
        }

        double collect_t0 = GetTime();
        int format_count = 0;
        for (int b = 0; b < n_blocks; ++b) {
            stats->total_records += paras[b].n_total_records;
            for (int r = 0; r < paras[b].n_total_records; ++r) {
                if (format_count >= BATCH_SIZE) {
                    fprintf(stderr, "ERROR: MPI bam2sam format batch capacity exceeded.\n");
                    goto cleanup;
                }
                fmt_cur->bams[format_count++] = paras[b].output_records[r];
            }
        }
        stats->t_collect += GetTime() - collect_t0;

        if (format_count > 0) {
            MpiResetSamFormatBatchMPI(fmt_cur, format_count);
            double format_t0 = GetTime();
            __real_athread_spawn((void *)slave_sam_format, fmt_cur, 1);
            #ifdef ENABLE_MASKING
            if (do_read_group(&next_n_blocks) != 0) goto cleanup;
            #endif
            athread_join();
            stats->t_format += GetTime() - format_t0;

            #ifndef ENABLE_MASKING
            if (do_read_group(&next_n_blocks) != 0) goto cleanup;
            #endif

            std::swap(fmt_cur, fmt_prev);
            has_pending_write = true;
            stats->format_tiles++;
        } else {
            if (do_read_group(&next_n_blocks) != 0) goto cleanup;
        }

        stats->group_count++;
    }

    if (flush_pending_write() != 0) goto cleanup;
    ret_code = 0;

cleanup:
    {
        double free_t0 = GetTime();
        MpiDestroySamFormatBatchMPI(fmt_cur);
        MpiDestroySamFormatBatchMPI(fmt_prev);
        if (fmt_cur) aligned_free_custom((unsigned char *)fmt_cur);
        if (fmt_prev) aligned_free_custom((unsigned char *)fmt_prev);
        MpiFreeBamToSamBlockSet(&input_blocks);
        MpiFreeBamToSamBlockSet(&un_blocks);
        MpiFreeBamToSamRecordSet(&record_set);
        stats->t_free_workspace += GetTime() - free_t0;
    }
    stats->t_fused_total += GetTime() - fused_t0;
    return ret_code;
}
