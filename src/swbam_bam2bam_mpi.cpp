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
    void slave_mpi_decompress_filterfunc();
    void slave_mpi_decompress_bam2bam_passthrough();
    void slave_mpi_compressfunc();
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
    int records_per_block;
    int n_blocks;
    size_t arena_stride;
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

int MpiAllocateRecordSet(MpiRecordSet *set,
                         int n_blocks,
                         int records_per_block,
                         size_t arena_stride) {
    int total_records = n_blocks * records_per_block;
    set->capacity = total_records;
    set->records_per_block = records_per_block;
    set->n_blocks = n_blocks;
    set->arena_stride = arena_stride;
    set->records = (bam1_t *)aligned_alloc_custom(64, (size_t)total_records * sizeof(bam1_t));
    set->data = aligned_alloc_custom(64, (size_t)n_blocks * arena_stride);
    set->ptrs = (bam1_t **)aligned_alloc_custom(64, (size_t)total_records * sizeof(bam1_t *));
    set->bam_lens = (uint32_t *)aligned_alloc_custom(64, (size_t)total_records * sizeof(uint32_t));
    if (!set->records || !set->data || !set->ptrs || !set->bam_lens) return -1;
    memset(set->records, 0, (size_t)total_records * sizeof(bam1_t));
    memset(set->data, 0, (size_t)n_blocks * arena_stride);
    memset(set->bam_lens, 0, (size_t)total_records * sizeof(uint32_t));
    for (int i = 0; i < total_records; ++i) {
        set->records[i].data = nullptr;
        set->records[i].m_data = 0;
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
    set->records_per_block = 0;
    set->n_blocks = 0;
    set->arena_stride = 0;
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

int MpiAppendRecordRangeToWorkspace(MpiFlatPackWorkspace *workspace,
                                    bam1_t **records,
                                    const uint32_t *bam_lens,
                                    int n_records,
                                    uint32_t total_len) {
    if (n_records <= 0) return 0;

    if (total_len > 0 &&
        total_len <= BGZF_BLOCK_SIZE &&
        (workspace->current_records == 0 ||
         workspace->current_len + total_len <= BGZF_BLOCK_SIZE)) {
        if (workspace->total_records + n_records > (int)workspace->records.size()) return -2;
        if (workspace->current_records == 0) {
            workspace->current_begin = workspace->records.data() + workspace->total_records;
        }
        memcpy(workspace->records.data() + workspace->total_records,
               records,
               (size_t)n_records * sizeof(bam1_t *));
        workspace->total_records += n_records;
        workspace->current_records += n_records;
        workspace->current_len += total_len;
        return 0;
    }

    int offset = 0;
    while (offset < n_records) {
        if (workspace->current_records == 0) {
            workspace->current_begin = workspace->records.data() + workspace->total_records;
        }

        uint32_t room = BGZF_BLOCK_SIZE - workspace->current_len;
        int take = 0;
        uint32_t take_len = 0;
        while (offset + take < n_records) {
            uint32_t packed_len = bam_lens[offset + take] + 4;
            if (packed_len > BGZF_BLOCK_SIZE) return -1;
            if (take_len + packed_len > room) break;
            take_len += packed_len;
            take++;
        }

        if (take == 0) {
            if (MpiSealCurrentPackBlock(workspace) != 0) return -3;
            continue;
        }
        if (workspace->total_records + take > (int)workspace->records.size()) return -2;

        memcpy(workspace->records.data() + workspace->total_records,
               records + offset,
               (size_t)take * sizeof(bam1_t *));
        workspace->total_records += take;
        workspace->current_records += take;
        workspace->current_len += take_len;
        offset += take;
    }

    return 0;
}

int MpiBuildPassthroughPlans(MpiFlatPackWorkspace *workspace,
                             Bam2BamPara *paras,
                             MpiBlockSet *un_blocks,
                             int n_blocks) {
    for (int b = 0; b < n_blocks; ++b) {
        int n_records = paras[b].n_kept_records;
        if (n_records <= 0) continue;

        uint32_t total_len = paras[b].kept_total_len;
        if (total_len == 0) total_len = un_blocks->blocks[b].length;
        if (total_len > BGZF_BLOCK_SIZE) {
            int ret = MpiAppendRecordRangeToWorkspace(workspace,
                                                      paras[b].output_records,
                                                      paras[b].bam_lens,
                                                      n_records,
                                                      total_len);
            if (ret != 0) return ret;
            continue;
        }

        if (workspace->current_records > 0 && MpiSealCurrentPackBlock(workspace) != 0) return -3;
        if (workspace->active_blocks >= 64) return -3;
        MpiPackedBlockPlan &plan = workspace->plans[workspace->active_blocks++];
        plan.records = paras[b].output_records;
        plan.n_records = n_records;
        plan.total_len = total_len;
        workspace->total_records += n_records;
    }

    if (workspace->current_records > 0 && MpiSealCurrentPackBlock(workspace) != 0) return -3;
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
    para->compress_level = 1;
    para->compress_serialize_cycles = 0;
    para->compress_alloc_cycles = 0;
    para->compress_deflate_cycles = 0;
    para->compress_footer_cycles = 0;
    para->compress_total_cycles = 0;
}

void MpiAccumulateCompressDetail(const Comp_Para *paras,
                                 int active_output_blocks,
                                 double compress_wall,
                                 MpiBamToBamStats *stats) {
    if (!stats || active_output_blocks <= 0 || compress_wall <= 0.0) return;

    const Comp_Para *critical = nullptr;
    uint64_t critical_total = 0;
    for (int k = 0; k < active_output_blocks; ++k) {
        if (paras[k].compress_total_cycles >= critical_total) {
            critical_total = paras[k].compress_total_cycles;
            critical = &paras[k];
        }
    }
    if (!critical || critical_total == 0) {
        stats->t_compress_other += compress_wall;
        return;
    }

    const double scale = compress_wall / (double)critical_total;
    const double serialize_time = scale * (double)critical->compress_serialize_cycles;
    const double alloc_time = scale * (double)critical->compress_alloc_cycles;
    const double deflate_time = scale * (double)critical->compress_deflate_cycles;
    const double footer_time = scale * (double)critical->compress_footer_cycles;
    double other_time = compress_wall - serialize_time - alloc_time - deflate_time - footer_time;
    if (other_time < 0.0) other_time = 0.0;

    stats->t_compress_serialize += serialize_time;
    stats->t_compress_alloc += alloc_time;
    stats->t_compress_deflate += deflate_time;
    stats->t_compress_footer += footer_time;
    stats->t_compress_other += other_time;
}

void MpiAccumulateDecompDetail(const Bam2BamPara *paras,
                               int active_blocks,
                               double decomp_wall,
                               MpiBamToBamStats *stats) {
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
                     int compress_level,
                     MpiBamToBamStats *stats) {
    const int NB = 64;
    const int records_per_block = (int)MPI_RECORDS_PER_BLOCK;
    const size_t block_arena_stride = MPI_BAM_BLOCK_ARENA_SIZE;
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
        MpiAllocateRecordSet(&record_set, NB, records_per_block, block_arena_stride) != 0) {
        fprintf(stderr, "ERROR: failed to allocate MPI bam2bam 1CG workspace.\n");
        return -1;
    }
    MpiInitPackWorkspace(&pack_workspace, (size_t)NB * records_per_block);
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
            comp_active[k].compress_level = compress_level;
            comp_active[k].compress_serialize_cycles = 0;
            comp_active[k].compress_alloc_cycles = 0;
            comp_active[k].compress_deflate_cycles = 0;
            comp_active[k].compress_footer_cycles = 0;
            comp_active[k].compress_total_cycles = 0;
        }
        for (int k = active_output_blocks; k < NB; ++k) {
            MpiInitEmptyCompPara(&comp_active[k], k);
        }
        // stats->t_compress_prepare += GetTime() - comp_prepare_t0;

        double compress_t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_compressfunc, comp_active, 1);
        #ifdef ENABLE_MASKING
        if (flush_pending() != 0) return -1;
        if (do_read_group(next_n_blocks) != 0) return -1;
        #endif
        athread_join();
        double compress_wall = GetTime() - compress_t0;
        stats->t_compress += compress_wall;
        MpiAccumulateCompressDetail(comp_active, active_output_blocks, compress_wall, stats);
        
        #ifndef ENABLE_MASKING
        if (flush_pending() != 0) return -1;
        if (do_read_group(next_n_blocks) != 0) return -1;
        #endif

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
            paras[b].output_records = record_set.ptrs + (size_t)b * records_per_block;
            paras[b].record_base = record_set.records + (size_t)b * records_per_block;
            paras[b].bam_lens = record_set.bam_lens + (size_t)b * records_per_block;
            paras[b].record_capacity = records_per_block;
            paras[b].data_arena = record_set.data + (size_t)b * block_arena_stride;
            paras[b].data_arena_capacity = block_arena_stride;
            paras[b].data_arena_used = 0;
            paras[b].record_index = 0;
            paras[b].actual_value = 0;
            paras[b].limit_value = 0;
            paras[b].limit_id = BOUNDS_LIMIT_NONE;
            paras[b].n_total_records = 0;
            paras[b].n_kept_records = 0;
            paras[b].kept_total_len = 0;
            paras[b].decomp_alloc_cycles = 0;
            paras[b].decomp_inflate_cycles = 0;
            paras[b].decomp_crc_cycles = 0;
            paras[b].decomp_parse_cycles = 0;
            paras[b].decomp_total_cycles = 0;
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
        __real_athread_spawn((void *)(no_filter ? slave_mpi_decompress_bam2bam_passthrough
                                                : slave_mpi_decompress_filterfunc),
                             paras, 1);
        athread_join();
        double decomp_wall = GetTime() - decomp_t0;
        stats->t_decomp_filter += decomp_wall;
        MpiAccumulateDecompDetail(paras, n_blocks, decomp_wall, stats);

        // double decomp_check_t0 = GetTime();
        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].status != 0) {
                if (paras[b].status == -3) {
                    fprintf(stderr,
                            "ERROR: MPI bam2bam capacity exceeded on input block %d. limit_id=%d limit=%lld actual=%lld record=%d.\n",
                            b, paras[b].limit_id, paras[b].limit_value,
                            paras[b].actual_value, paras[b].record_index);
                } else {
                    fprintf(stderr, "ERROR: MPI bam2bam decompress/filter failed on input block %d with status %d.\n",
                            b, paras[b].status);
                }
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
        }

        int pack_ret = 0;
        if (no_filter) {
            pack_ret = MpiBuildPassthroughPlans(&pack_workspace, paras, &un_blocks, n_blocks);
        } else {
            for (int b = 0; b < n_blocks; ++b) {
                pack_ret = MpiAppendRecordRangeToWorkspace(&pack_workspace,
                                                           paras[b].output_records,
                                                           paras[b].bam_lens,
                                                           paras[b].n_kept_records,
                                                           paras[b].kept_total_len);
                if (pack_ret != 0) break;
            }
            if (pack_ret == 0 && MpiSealCurrentPackBlock(&pack_workspace) != 0) {
                pack_ret = -3;
            }
        }
        if (pack_ret != 0) {
            if (pack_ret == -1) {
                fprintf(stderr, "ERROR: MPI bam2bam pack encountered an oversized BAM record.\n");
            } else if (pack_ret == -2) {
                fprintf(stderr, "ERROR: MPI bam2bam pack workspace capacity exceeded.\n");
            } else {
                fprintf(stderr, "ERROR: MPI bam2bam produced more than %d output blocks from one input group.\n",
                        NB);
            }
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
