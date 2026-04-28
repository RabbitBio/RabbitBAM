#include "swbam_cgs.h"

#include <cstdio>
#include <cstring>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

extern "C" {
    void slave_cgs_bam_decodefunc();
    void slave_cgs_sam_format_tile();
}

namespace {

int InitSamFormatBatchCGS(CgsCrossAllocator &alloc, CgsSamFormatBatch *batch, bam1_t **records, sam_hdr_t *hdr) {
    batch->hdr = hdr;
    batch->records = records;
    batch->total_records = 0;
    char *out_storage = (char *)alloc.allocate(64, (size_t)CGS_NB * CGS_SAM_FORMAT_CORE_BUFFER_SIZE);
    if (!out_storage) return -1;

    for (int i = 0; i < CGS_NB; ++i) {
        batch->core_out_lines[i].l = 0;
        batch->core_out_lines[i].m = CGS_SAM_FORMAT_CORE_BUFFER_SIZE;
        batch->core_out_lines[i].s = out_storage + (size_t)i * CGS_SAM_FORMAT_CORE_BUFFER_SIZE;
        batch->core_line_bufs[i].l = 0;
        batch->core_line_bufs[i].m = 0;
        batch->core_line_bufs[i].s = nullptr;
        batch->status[i] = -1;
        batch->formatted_records[i] = 0;
    }
    return 0;
}

void ResetSamFormatBatchCGS(CgsSamFormatBatch *batch, int n_records) {
    batch->total_records = n_records;
    for (int i = 0; i < CGS_NB; ++i) {
        batch->core_out_lines[i].l = 0;
        batch->core_line_bufs[i].l = 0;
        batch->status[i] = 0;
        batch->formatted_records[i] = 0;
    }
}

} // namespace

int FusedBamToSamCGS(MemReader &reader, MemWriter &mem_writer, sam_hdr_t *hdr) {
    if (CheckCgsResources() != 0) return -1;

    double t0 = GetTime();
    double t_alloc_init = 0, t_setup_reset = 0, t_format_setup = 0, t_status_check = 0;
    double t_decomp = 0, t_format = 0, t_collect = 0, t_read = 0, t_write = 0;
    const int format_tile_records = CGS_NB * CGS_BAM2SAM_FORMAT_RECORDS_PER_CPE;

    double alloc_t0 = GetTime();
    CgsCrossAllocator alloc;
    CgsBamDecodePara *paras = (CgsBamDecodePara *)alloc.allocate(64, (size_t)CGS_NB * sizeof(CgsBamDecodePara));
    CgsBlockSet input_blocks = {};
    CgsBlockSet un_blocks = {};
    CgsRecordSet record_set = {};
    bam1_t **format_records = (bam1_t **)alloc.allocate(64, (size_t)format_tile_records * sizeof(bam1_t *));
    CgsSamFormatBatch *format_cur =
        (CgsSamFormatBatch *)alloc.allocate(64, sizeof(CgsSamFormatBatch));
    CgsSamFormatBatch *format_prev =
        (CgsSamFormatBatch *)alloc.allocate(64, sizeof(CgsSamFormatBatch));

    if (!paras || !format_records || !format_cur || !format_prev ||
        AllocBlockSet(alloc, &input_blocks, CGS_NB) != 0 ||
        AllocBlockSet(alloc, &un_blocks, CGS_NB) != 0 ||
        AllocRecordSet(alloc, &record_set, CGS_NB * (int)MAX_RECORDS_PER_BLOCK) != 0 ||
        InitSamFormatBatchCGS(alloc, format_cur, format_records, hdr) != 0 ||
        InitSamFormatBatchCGS(alloc, format_prev, format_records, hdr) != 0) {
        fprintf(stderr, "ERROR: failed to allocate CGS bam2sam cross memory. Increase -cross_size and retry.\n");
        return -1;
    }
    t_alloc_init += GetTime() - alloc_t0;

    long long input_blocks_count = 0;
    long long group_count = 0;
    long long total_records = 0;
    long long format_tiles = 0;
    long long overlapped_reads = 0;
    long long overlapped_writes = 0;
    bool has_pending_write = false;

    auto do_read_group = [&](int *n_blocks) -> int {
        double read_t0 = GetTime();
        int count = 0;
        for (int b = 0; b < CGS_NB; ++b) {
            bam_block *blk = &input_blocks.blocks[b];
            int ret = MemReadBlockCGS(reader.base, reader.size, reader.pos, blk);
            if (ret < 0 || blk->length == 28) break;
            blk->block_id = b;
            blk->pos = 0;
            count++;
        }
        *n_blocks = count;
        t_read += GetTime() - read_t0;
        return 0;
    };

    auto flush_pending_write = [&]() -> int {
        if (!has_pending_write) return 0;
        double write_t0 = GetTime();
        for (int i = 0; i < CGS_NB; ++i) {
            if (format_prev->status[i] != 0) {
                fprintf(stderr, "ERROR: CGS bam2sam format failed on CPE %d with status %d.\n",
                        i, format_prev->status[i]);
                return -1;
            }
            kstring_t *ks = &format_prev->core_out_lines[i];
            if (ks->l > 0 && WriteBytesToMemCGS(mem_writer, ks->s, ks->l) != 0) {
                fprintf(stderr, "ERROR: failed to append CGS SAM text to memory writer.\n");
                return -1;
            }
        }
        has_pending_write = false;
        t_write += GetTime() - write_t0;
        overlapped_writes++;
        return 0;
    };

    auto flush_format_tile = [&](int n_records, bool *next_read_done, int *next_n_blocks) -> int {
        if (n_records <= 0) return 0;
        double setup_t0 = GetTime();
        ResetSamFormatBatchCGS(format_cur, n_records);
        t_format_setup += GetTime() - setup_t0;
        double format_t0 = GetTime();
        __real_athread_spawn_cgs((void *)slave_cgs_sam_format_tile, format_cur, 1);
        int overlap_status = 0;
        if (flush_pending_write() != 0) overlap_status = -1;
        if (overlap_status == 0 && !*next_read_done) {
            if (do_read_group(next_n_blocks) != 0) overlap_status = -1;
            *next_read_done = true;
            overlapped_reads++;
        }
        athread_join_cgs();
        t_format += GetTime() - format_t0;


        if (overlap_status != 0) return -1;

        double check_t0 = GetTime();
        for (int i = 0; i < CGS_NB; ++i) {
            if (format_cur->status[i] != 0) {
                fprintf(stderr, "ERROR: CGS bam2sam format failed on CPE %d with status %d.\n",
                        i, format_cur->status[i]);
                return -1;
            }
        }
        t_status_check += GetTime() - check_t0;
        std::swap(format_cur, format_prev);
        has_pending_write = true;
        format_tiles++;
        return 0;
    };

    int n_blocks = 0;
    do_read_group(&n_blocks);
    while (n_blocks > 0) {
        input_blocks_count += n_blocks;
        double setup_t0 = GetTime();
        for (int b = 0; b < CGS_NB; ++b) {
            paras[b].block_id = b;
            paras[b].output_records = record_set.ptrs + (size_t)b * MAX_RECORDS_PER_BLOCK;
            paras[b].n_records = 0;
            if (b < n_blocks) {
                paras[b].input_block = &input_blocks.blocks[b];
                paras[b].un_comp_block = &un_blocks.blocks[b];
                paras[b].status = 0;
            } else {
                paras[b].input_block = nullptr;
                paras[b].un_comp_block = nullptr;
                paras[b].status = -1;
            }
        }
        t_setup_reset += GetTime() - setup_t0;

        double decomp_t0 = GetTime();
        __real_athread_spawn_cgs((void *)slave_cgs_bam_decodefunc, paras, 1);
        if (flush_pending_write() != 0) return -1;
        athread_join_cgs();
        t_decomp += GetTime() - decomp_t0;

        double check_t0 = GetTime();
        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].status != 0) {
                fprintf(stderr, "ERROR: CGS bam2sam decode failed on input block %d with status %d.\n",
                        b, paras[b].status);
                return -1;
            }
        }
        t_status_check += GetTime() - check_t0;

        double collect_t0 = GetTime();
        int format_count = 0;
        int next_n_blocks = 0;
        bool next_read_done = false;
        for (int b = 0; b < n_blocks; ++b) {
            int copied = 0;
            int n_records = paras[b].n_records;
            total_records += n_records;
            while (copied < n_records) {
                int room = format_tile_records - format_count;
                int take = n_records - copied;
                if (take > room) take = room;
                memcpy(format_records + format_count,
                       paras[b].output_records + copied,
                       (size_t)take * sizeof(bam1_t *));
                format_count += take;
                copied += take;
                if (format_count == format_tile_records) {
                    t_collect += GetTime() - collect_t0;
                    if (flush_format_tile(format_count, &next_read_done, &next_n_blocks) != 0) return -1;
                    collect_t0 = GetTime();
                    format_count = 0;
                }
            }
        }
        t_collect += GetTime() - collect_t0;
        if (flush_format_tile(format_count, &next_read_done, &next_n_blocks) != 0) return -1;
        if (!next_read_done) {
            do_read_group(&next_n_blocks);
        }

        group_count++;
        n_blocks = next_n_blocks;
    }

    if (flush_pending_write() != 0) return -1;

    printf("FusedBamToSamCGS finished. blocks=%lld groups=%lld records=%lld format_tiles=%lld cost=%.3f s\n",
           input_blocks_count, group_count, total_records, format_tiles, GetTime() - t0);
    printf("  decomp_slave=%.3f  format_slave=%.3f  collect=%.3f  read=%.3f  write=%.3f\n",
           t_decomp, t_format, t_collect, t_read, t_write);
    printf("  alloc_init=%.3f  setup_reset=%.3f  format_setup=%.3f  status_check=%.3f\n",
           t_alloc_init, t_setup_reset, t_format_setup, t_status_check);
    printf("  cgs_threads=%d  format_records_per_cpe=%d  format_tile_records=%d  overlapped_reads=%lld  overlapped_writes=%lld\n",
           CGS_NB, CGS_BAM2SAM_FORMAT_RECORDS_PER_CPE, format_tile_records,
           overlapped_reads, overlapped_writes);

    return 0;
}
