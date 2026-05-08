#include "swbam_mpi.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

extern "C" {
    void slave_decompress_filterfunc();
    void slave_decompress_bam2bam_passthrough();
    void slave_compressfunc();
}

namespace {

struct MpiBlockSet {
    bam_block *blocks;
    unsigned char *data;
    int n;
};

struct MpiRecordSet {
    bam1_t *records;
    unsigned char *data;
    bam1_t **ptrs;
    uint32_t *bam_lens;
    int capacity;
};

struct MpiPackedBlockPlan {
    bam1_t **records;
    int n_records;
    uint32_t total_len;
};

struct MpiFlatPackWorkspace {
    std::vector<bam1_t *> records;
    MpiPackedBlockPlan plans[64];
    int active_blocks;
    int total_records;
    bam1_t **current_begin;
    int current_records;
    uint32_t current_len;
};

int MpiAllocateBlockSet(MpiBlockSet *set, int n) {
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

void MpiFreeBlockSet(MpiBlockSet *set) {
    if (set->blocks) aligned_free_custom((unsigned char *)set->blocks);
    if (set->data) aligned_free_custom(set->data);
    set->blocks = nullptr;
    set->data = nullptr;
    set->n = 0;
}

int MpiAllocateRecordSet(MpiRecordSet *set, int total_records) {
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

void MpiFreeRecordSet(MpiRecordSet *set) {
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

void MpiResetPackWorkspace(MpiFlatPackWorkspace *workspace) {
    workspace->active_blocks = 0;
    workspace->total_records = 0;
    workspace->current_begin = workspace->records.empty() ? nullptr : workspace->records.data();
    workspace->current_records = 0;
    workspace->current_len = 0;
}

void MpiInitPackWorkspace(MpiFlatPackWorkspace *workspace, size_t capacity) {
    workspace->records.resize(capacity);
    MpiResetPackWorkspace(workspace);
}

int MpiSealCurrentPackBlock(MpiFlatPackWorkspace *workspace) {
    if (workspace->current_records == 0) return 0;
    if (workspace->active_blocks >= 64) return -1;
    MpiPackedBlockPlan &plan = workspace->plans[workspace->active_blocks++];
    plan.records = workspace->current_begin;
    plan.n_records = workspace->current_records;
    plan.total_len = workspace->current_len;
    workspace->current_begin = workspace->records.data() + workspace->total_records;
    workspace->current_records = 0;
    workspace->current_len = 0;
    return 0;
}

int MpiAppendRecordToWorkspace(MpiFlatPackWorkspace *workspace, bam1_t *record, uint32_t bam_len) {
    const uint32_t packed_len = bam_len + 4;
    if (packed_len > BGZF_BLOCK_SIZE) return -1;
    if (workspace->current_records > 0 && workspace->current_len + packed_len > BGZF_BLOCK_SIZE) return 1;
    if (workspace->total_records >= (int)workspace->records.size()) return -2;
    if (workspace->current_records == 0) {
        workspace->current_begin = workspace->records.data() + workspace->total_records;
    }
    workspace->records[workspace->total_records++] = record;
    workspace->current_records++;
    workspace->current_len += packed_len;
    return 0;
}

void MpiInitEmptyCompPara(Comp_Para *para, int block_id) {
    para->block_id = block_id;
    para->input_records = nullptr;
    para->n_records = 0;
    para->un_comp_block = nullptr;
    para->un_comp_size = 0;
    para->output_block = nullptr;
    para->output_size = 0;
    para->status = -1;
}

int MpiMemReadBlock(char *base, size_t size, size_t &pos, bam_block *block) {
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

} // namespace

int FusedBamToBamMPI(MemReader &reader, MemWriter &mem_writer,
                     const BamFilterOptions &filter,
                     MpiBamToBamStats *stats) {
    const int NB = 64;
    const bool no_filter = bam_filter_is_noop(filter);
    double fused_t0 = GetTime();

    Bam2BamPara paras[NB];
    Comp_Para comp_a[NB], comp_b[NB];
    MpiBlockSet input_blocks = {};
    MpiBlockSet un_blocks = {};
    MpiBlockSet comp_un_a = {};
    MpiBlockSet comp_un_b = {};
    MpiBlockSet out_a = {};
    MpiBlockSet out_b = {};
    MpiRecordSet record_set = {};
    MpiFlatPackWorkspace pack_workspace;

    if (MpiAllocateBlockSet(&input_blocks, NB) != 0 ||
        MpiAllocateBlockSet(&un_blocks, NB) != 0 ||
        MpiAllocateBlockSet(&comp_un_a, NB) != 0 ||
        MpiAllocateBlockSet(&comp_un_b, NB) != 0 ||
        MpiAllocateBlockSet(&out_a, NB) != 0 ||
        MpiAllocateBlockSet(&out_b, NB) != 0 ||
        MpiAllocateRecordSet(&record_set, NB * (int)MAX_RECORDS_PER_BLOCK) != 0) {
        fprintf(stderr, "ERROR: failed to allocate MPI bam2bam 1CG workspace.\n");
        return -1;
    }
    MpiInitPackWorkspace(&pack_workspace, NB * MAX_RECORDS_PER_BLOCK);
    for (int i = 0; i < NB; ++i) {
        MpiInitEmptyCompPara(&comp_a[i], i);
        MpiInitEmptyCompPara(&comp_b[i], i);
    }
    // stats->t_alloc_init += GetTime() - fused_t0;

    Comp_Para *comp_active = comp_a;
    Comp_Para *comp_pending = comp_b;
    MpiBlockSet *comp_un_active = &comp_un_a;
    MpiBlockSet *comp_un_pending = &comp_un_b;
    MpiBlockSet *out_active = &out_a;
    MpiBlockSet *out_pending = &out_b;
    bool has_pending = false;

    auto flush_pending = [&]() -> int {
        if (!has_pending) return 0;
        double flush_t0 = GetTime();
        for (int k = 0; k < NB; ++k) {
            if (comp_pending[k].status == 0 && comp_pending[k].output_block) {
                if (MpiWriteBlockToMem(mem_writer, comp_pending[k].output_block) != 0) {
                    fprintf(stderr, "ERROR: MPI bam2bam failed to append compressed block.\n");
                    return -1;
                }
                stats->bgzf_blocks++;
            }
            MpiInitEmptyCompPara(&comp_pending[k], k);
        }
        has_pending = false;
        stats->t_write += GetTime() - flush_t0;
        return 0;
    };

    auto do_read_group = [&](int *n_blocks) -> int {
        double read_t0 = GetTime();
        int count = 0;
        for (int b = 0; b < NB; ++b) {
            bam_block *blk = &input_blocks.blocks[b];
            int ret = MpiMemReadBlock(reader.base, reader.size, reader.pos, blk);
            if (ret < 0 || blk->length == 28) break;
            blk->block_id = b;
            blk->pos = 0;
            count++;
        }
        *n_blocks = count;
        stats->t_read += GetTime() - read_t0;
        return 0;
    };

    auto do_compress = [&](int active_output_blocks, int *next_n_blocks) -> int {
        // double comp_prepare_t0 = GetTime();
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
        for (int k = active_output_blocks; k < NB; ++k) {
            MpiInitEmptyCompPara(&comp_active[k], k);
        }
        // stats->t_compress_prepare += GetTime() - comp_prepare_t0;

        double compress_t0 = GetTime();
        __real_athread_spawn((void *)slave_compressfunc, comp_active, 1);
        if (flush_pending() != 0) return -1;
        if (do_read_group(next_n_blocks) != 0) return -1;
        athread_join();
        stats->t_compress += GetTime() - compress_t0;

        // double comp_check_t0 = GetTime();
        for (int k = 0; k < active_output_blocks; ++k) {
            if (comp_active[k].status != 0) {
                fprintf(stderr, "ERROR: MPI bam2bam compress failed on output block %d with status %d.\n",
                        k, comp_active[k].status);
                return -1;
            }
        }

        has_pending = active_output_blocks > 0;
        std::swap(comp_active, comp_pending);
        std::swap(comp_un_active, comp_un_pending);
        std::swap(out_active, out_pending);
        MpiResetPackWorkspace(&pack_workspace);
        // stats->t_compress_check += GetTime() - comp_check_t0;
        return 0;
    };

    int next_n_blocks = 0;
    // double initial_read_t0 = GetTime();
    do_read_group(&next_n_blocks);
    // stats->t_initial_read += GetTime() - initial_read_t0;
    while (next_n_blocks > 0) {
        int n_blocks = next_n_blocks;
        next_n_blocks = 0;
        stats->input_blocks += n_blocks;

        // double prepare_t0 = GetTime();
        for (int b = 0; b < NB; ++b) {
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
        // stats->t_prepare_decomp += GetTime() - prepare_t0;

        double decomp_t0 = GetTime();
        __real_athread_spawn((void *)(no_filter ? slave_decompress_bam2bam_passthrough
                                                : slave_decompress_filterfunc),
                             paras, 1);
        athread_join();
        stats->t_decomp_filter += GetTime() - decomp_t0;

        // double decomp_check_t0 = GetTime();
        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].status != 0) {
                fprintf(stderr, "ERROR: MPI bam2bam decompress/filter failed on input block %d with status %d.\n",
                        b, paras[b].status);
                return -1;
            }
        }
        // stats->t_decomp_check += GetTime() - decomp_check_t0;

        MpiResetPackWorkspace(&pack_workspace);
        double pack_t0 = GetTime();
        for (int b = 0; b < n_blocks; ++b) {
            stats->total_records += paras[b].n_total_records;
            stats->kept_records += paras[b].n_kept_records;
            stats->dropped_records += paras[b].n_total_records - paras[b].n_kept_records;
            stats->pack_records += paras[b].n_kept_records;

            for (int r = 0; r < paras[b].n_kept_records; ++r) {
                bam1_t *record = paras[b].output_records[r];
                uint32_t bam_len = paras[b].bam_lens[r];
                while (true) {
                    int ret = MpiAppendRecordToWorkspace(&pack_workspace, record, bam_len);
                    if (ret == 0) break;
                    if (ret == -1) {
                        fprintf(stderr, "ERROR: MPI bam2bam pack encountered an oversized BAM record.\n");
                        return -1;
                    }
                    if (ret == -2) {
                        fprintf(stderr, "ERROR: MPI bam2bam pack workspace capacity exceeded.\n");
                        return -1;
                    }
                    if (MpiSealCurrentPackBlock(&pack_workspace) != 0) {
                        fprintf(stderr, "ERROR: MPI bam2bam produced more than %d output blocks from one input group.\n",
                                NB);
                        return -1;
                    }
                }
            }
        }
        if (MpiSealCurrentPackBlock(&pack_workspace) != 0) {
            fprintf(stderr, "ERROR: MPI bam2bam output block plan capacity exceeded.\n");
            return -1;
        }
        stats->t_pack += GetTime() - pack_t0;

        stats->group_count++;
        if (pack_workspace.active_blocks == 0) {
            // double empty_t0 = GetTime();
            if (flush_pending() != 0) return -1;
            if (do_read_group(&next_n_blocks) != 0) return -1;
            // stats->t_empty_group_path += GetTime() - empty_t0;
            continue;
        }

        if (do_compress(pack_workspace.active_blocks, &next_n_blocks) != 0) return -1;
    }

    if (flush_pending() != 0) return -1;

    // double free_t0 = GetTime();
    MpiFreeBlockSet(&input_blocks);
    MpiFreeBlockSet(&un_blocks);
    MpiFreeBlockSet(&comp_un_a);
    MpiFreeBlockSet(&comp_un_b);
    MpiFreeBlockSet(&out_a);
    MpiFreeBlockSet(&out_b);
    MpiFreeRecordSet(&record_set);
    // stats->t_free_workspace += GetTime() - free_t0;
    stats->t_fused_total += GetTime() - fused_t0;
    return 0;
}
