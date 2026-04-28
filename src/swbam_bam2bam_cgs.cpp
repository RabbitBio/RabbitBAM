#include "swbam_cgs.h"

#include <cstdio>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

extern "C" {
    void slave_cgs_decompress_filterfunc();
    void slave_cgs_decompress_bam2bam_passthrough();
    void slave_cgs_compressfunc();
}

int FusedBamToBamCGS(MemReader &reader, MemWriter &mem_writer, const BamFilterOptions &filter) {
    if (CheckCgsResources() != 0) return -1;

    double t0 = GetTime();
    double t_alloc_init = 0, t_setup_reset = 0, t_status_check = 0, t_compress_setup = 0;
    double t_decomp_filter = 0, t_pack = 0, t_compress = 0, t_read = 0, t_write = 0;
    const bool no_filter = bam_filter_is_noop(filter);

    double alloc_t0 = GetTime();
    CgsCrossAllocator alloc;
    Bam2BamPara *paras = (Bam2BamPara *)alloc.allocate(64, (size_t)CGS_NB * sizeof(Bam2BamPara));
    Comp_Para *comp_a = (Comp_Para *)alloc.allocate(64, (size_t)CGS_NB * sizeof(Comp_Para));
    Comp_Para *comp_b = (Comp_Para *)alloc.allocate(64, (size_t)CGS_NB * sizeof(Comp_Para));
    CgsBlockSet input_blocks = {};
    CgsBlockSet un_blocks = {};
    CgsBlockSet comp_un_a = {};
    CgsBlockSet comp_un_b = {};
    CgsBlockSet out_a = {};
    CgsBlockSet out_b = {};
    CgsRecordSet record_set = {};
    CgsPackWorkspace pack_workspace = {};

    if (!paras || !comp_a || !comp_b ||
        AllocBlockSet(alloc, &input_blocks, CGS_NB) != 0 ||
        AllocBlockSet(alloc, &un_blocks, CGS_NB) != 0 ||
        AllocBlockSet(alloc, &comp_un_a, CGS_NB) != 0 ||
        AllocBlockSet(alloc, &comp_un_b, CGS_NB) != 0 ||
        AllocBlockSet(alloc, &out_a, CGS_NB) != 0 ||
        AllocBlockSet(alloc, &out_b, CGS_NB) != 0 ||
        AllocRecordSet(alloc, &record_set, CGS_NB * (int)MAX_RECORDS_PER_BLOCK) != 0 ||
        AllocPackWorkspace(alloc, &pack_workspace, CGS_NB * (int)MAX_RECORDS_PER_BLOCK, CGS_NB) != 0) {
        fprintf(stderr, "ERROR: failed to allocate CGS cross memory. Increase -cross_size and retry.\n");
        return -1;
    }

    for (int i = 0; i < CGS_NB; ++i) {
        InitEmptyCompParaCGS(&comp_a[i], i);
        InitEmptyCompParaCGS(&comp_b[i], i);
    }
    t_alloc_init += GetTime() - alloc_t0;

    Comp_Para *comp_active = comp_a;
    Comp_Para *comp_pending = comp_b;
    CgsBlockSet *comp_un_active = &comp_un_a;
    CgsBlockSet *comp_un_pending = &comp_un_b;
    CgsBlockSet *out_active = &out_a;
    CgsBlockSet *out_pending = &out_b;
    bool has_pending = false;

    long long bgzf_blocks = 0;
    long long total_records = 0;
    long long kept_records = 0;
    long long dropped_records = 0;
    long long input_block_count = 0;
    long long group_count = 0;
    long long pack_records = 0;

    auto flush_pending = [&]() -> int {
        double flush_t0 = GetTime();
        if (!has_pending) return 0;
        for (int k = 0; k < CGS_NB; ++k) {
            if (comp_pending[k].status == 0 && comp_pending[k].output_block) {
                if (WriteBlockToMemCGS(mem_writer, comp_pending[k].output_block) != 0) {
                    fprintf(stderr, "ERROR: failed to append CGS compressed block to memory writer.\n");
                    return -1;
                }
                bgzf_blocks++;
            }
            InitEmptyCompParaCGS(&comp_pending[k], k);
        }
        has_pending = false;
        t_write += GetTime() - flush_t0;
        return 0;
    };

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

    auto do_compress = [&](int active_output_blocks, int *next_n_blocks) -> int {
        double setup_t0 = GetTime();
        for (int k = 0; k < active_output_blocks; ++k) {
            comp_active[k].block_id = k;
            comp_active[k].input_records = pack_workspace.plans[k].records;
            comp_active[k].n_records = pack_workspace.plans[k].n_records;
            comp_active[k].output_block = &out_active->blocks[k];
            comp_active[k].un_comp_block = &comp_un_active->blocks[k];
            comp_active[k].un_comp_size = (int)pack_workspace.plans[k].total_len;
            comp_active[k].output_size = 0;
            comp_active[k].status = 0;
        }
        for (int k = active_output_blocks; k < CGS_NB; ++k) {
            InitEmptyCompParaCGS(&comp_active[k], k);
        }
        t_compress_setup += GetTime() - setup_t0;

        double compress_t0 = GetTime();
        __real_athread_spawn_cgs((void *)slave_cgs_compressfunc, comp_active, 1);
        if (flush_pending() != 0) return -1;
        do_read_group(next_n_blocks);
        athread_join_cgs();
        t_compress += GetTime() - compress_t0;

        double check_t0 = GetTime();
        for (int k = 0; k < active_output_blocks; ++k) {
            if (comp_active[k].status != 0) {
                fprintf(stderr, "ERROR: CGS bam2bam compress failed on output block %d with status %d.\n",
                        k, comp_active[k].status);
                return -1;
            }
        }
        t_status_check += GetTime() - check_t0;

        has_pending = active_output_blocks > 0;
        std::swap(comp_active, comp_pending);
        std::swap(comp_un_active, comp_un_pending);
        std::swap(out_active, out_pending);
        ResetPackWorkspace(&pack_workspace);
        return 0;
    };

    int next_n_blocks = 0;
    do_read_group(&next_n_blocks);

    while (next_n_blocks > 0) {
        int n_blocks = next_n_blocks;
        next_n_blocks = 0;
        input_block_count += n_blocks;

        double setup_t0 = GetTime();
        for (int b = 0; b < CGS_NB; ++b) {
            paras[b].filter = filter;
            paras[b].block_id = b;
            paras[b].output_records = record_set.ptrs + (size_t)b * MAX_RECORDS_PER_BLOCK;
            paras[b].record_base = record_set.records + (size_t)b * MAX_RECORDS_PER_BLOCK;
            paras[b].bam_lens = record_set.bam_lens + (size_t)b * MAX_RECORDS_PER_BLOCK;
            paras[b].n_total_records = 0;
            paras[b].n_kept_records = 0;
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
        __real_athread_spawn_cgs((void *)(no_filter ? slave_cgs_decompress_bam2bam_passthrough
                                                    : slave_cgs_decompress_filterfunc),
                                 paras, 1);
        athread_join_cgs();
        t_decomp_filter += GetTime() - decomp_t0;

        double check_t0 = GetTime();
        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].status != 0) {
                fprintf(stderr, "ERROR: CGS bam2bam decompress/filter failed on input block %d with status %d.\n",
                        b, paras[b].status);
                return -1;
            }
        }
        t_status_check += GetTime() - check_t0;

        ResetPackWorkspace(&pack_workspace);
        double pack_t0 = GetTime();
        for (int b = 0; b < n_blocks; ++b) {
            total_records += paras[b].n_total_records;
            kept_records += paras[b].n_kept_records;
            dropped_records += paras[b].n_total_records - paras[b].n_kept_records;
            pack_records += paras[b].n_kept_records;

            for (int r = 0; r < paras[b].n_kept_records; ++r) {
                bam1_t *record = paras[b].output_records[r];
                uint32_t bam_len = paras[b].bam_lens[r];
                while (true) {
                    int ret = AppendRecordToPackWorkspace(&pack_workspace, record, bam_len);
                    if (ret == 0) break;
                    if (ret == -1) {
                        fprintf(stderr, "ERROR: CGS bam2bam pack encountered an oversized BAM record.\n");
                        return -1;
                    }
                    if (ret == -2) {
                        fprintf(stderr, "ERROR: CGS bam2bam pack workspace capacity exceeded.\n");
                        return -1;
                    }
                    if (SealCurrentPackBlock(&pack_workspace) != 0) {
                        fprintf(stderr, "ERROR: CGS bam2bam produced more than %d output blocks from one input group.\n",
                                CGS_NB);
                        return -1;
                    }
                }
            }
        }
        if (SealCurrentPackBlock(&pack_workspace) != 0) {
            fprintf(stderr, "ERROR: CGS bam2bam output block plan capacity exceeded.\n");
            return -1;
        }
        t_pack += GetTime() - pack_t0;

        group_count++;
        if (pack_workspace.active_blocks == 0) {
            if (flush_pending() != 0) return -1;
            do_read_group(&next_n_blocks);
            continue;
        }

        if (do_compress(pack_workspace.active_blocks, &next_n_blocks) != 0) return -1;
    }

    if (flush_pending() != 0) return -1;

    double keep_ratio = total_records > 0 ? (double)kept_records / (double)total_records : 1.0;
    printf("FusedBamToBamCGS finished. in_blocks=%lld groups=%lld total_records=%lld kept_records=%lld dropped_records=%lld bgzf_blocks=%lld cost=%.3f s\n",
           input_block_count, group_count, total_records, kept_records, dropped_records, bgzf_blocks, GetTime() - t0);
    printf("  decomp_filter_slave=%.3f  pack=%.3f  compress_slave=%.3f  read=%.3f  write=%.3f\n",
           t_decomp_filter, t_pack, t_compress, t_read, t_write);
    printf("  alloc_init=%.3f  setup_reset=%.3f  compress_setup=%.3f  status_check=%.3f\n",
           t_alloc_init, t_setup_reset, t_compress_setup, t_status_check);
    printf("  pack_records=%lld  output_bgzf_blocks=%lld  keep_ratio=%.6f  filter_mode=%s  cgs_threads=%d\n",
           pack_records, bgzf_blocks, keep_ratio, no_filter ? "passthrough" : "filtered", CGS_NB);

    return 0;
}
