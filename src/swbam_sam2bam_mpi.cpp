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
    void slave_mpi_copy_and_count();
    void slave_mpi_sam_parse_chunk();
    void slave_mpi_compressfunc();
}

namespace {

struct MpiSamToBamBlockSet {
    bam_block *blocks;
    unsigned char *data;
    int n;
};

struct MpiSamToBamRecordSet {
    bam1_t *records;
    unsigned char *data;
    int capacity;
};

struct MpiSamToBamPackedPlan {
    bam1_t **records;
    int n_records;
    uint32_t total_len;
};

struct MpiSamToBamPackWorkspace {
    bam1_t **records;
    MpiSamToBamPackedPlan *plans;
    int capacity;
    int max_plans;
    int active_blocks;
    int total_records;
    bam1_t **current_begin;
    int current_records;
    uint32_t current_len;
};

int MpiAllocateSamToBamBlockSet(MpiSamToBamBlockSet *set, int n) {
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

void MpiFreeSamToBamBlockSet(MpiSamToBamBlockSet *set) {
    if (set->blocks) aligned_free_custom((unsigned char *)set->blocks);
    if (set->data) aligned_free_custom(set->data);
    set->blocks = nullptr;
    set->data = nullptr;
    set->n = 0;
}

int MpiAllocateSamToBamRecordSet(MpiSamToBamRecordSet *set, int total_records) {
    set->capacity = total_records;
    set->records = (bam1_t *)aligned_alloc_custom(64, (size_t)total_records * sizeof(bam1_t));
    set->data = aligned_alloc_custom(64, (size_t)total_records * INIT_DATA_SIZE);
    if (!set->records || !set->data) return -1;
    return 0;
}

void MpiFreeSamToBamRecordSet(MpiSamToBamRecordSet *set) {
    if (set->records) aligned_free_custom((unsigned char *)set->records);
    if (set->data) aligned_free_custom(set->data);
    set->records = nullptr;
    set->data = nullptr;
    set->capacity = 0;
}

int MpiInitSamParseBatch(MpiSamParseBatch *batch,
                         char *text_storage,
                         uint32_t *bam_lens_storage) {
    memset(batch, 0, sizeof(MpiSamParseBatch));
    if (!text_storage || !bam_lens_storage) return -1;
    for (int i = 0; i < 64; ++i) {
        MpiSamParseChunk *chunk = &batch->chunks[i];
        chunk->text_buf = text_storage + (size_t)i * CHUNK_BUFFER_SIZE;
        chunk->bam_lens = bam_lens_storage + (size_t)i * MAX_BAMS_PER_CHUNK;
        chunk->src_ptr = nullptr;
        chunk->src_len = 0;
        chunk->text_len = 0;
        chunk->bams = nullptr;
        chunk->bam_data = nullptr;
        chunk->count = 0;
        chunk->parse_fast_records = 0;
        chunk->parse_fallback_records = 0;
        chunk->parse_total_cycles = 0;
        chunk->parse_core_cycles = 0;
        chunk->parse_aux_cycles = 0;
        chunk->parse_cg_cycles = 0;
        chunk->parse_fallback_cycles = 0;
    }
    return 0;
}

void MpiDestroySamParseBatch(MpiSamParseBatch *batch) {
    if (!batch) return;
    for (int i = 0; i < 64; ++i) {
        batch->chunks[i].text_buf = nullptr;
        batch->chunks[i].bams = nullptr;
        batch->chunks[i].bam_data = nullptr;
        batch->chunks[i].bam_lens = nullptr;
    }
}

void MpiResetPackWorkspace(MpiSamToBamPackWorkspace *workspace) {
    workspace->active_blocks = 0;
    workspace->total_records = 0;
    workspace->current_begin = workspace->records;
    workspace->current_records = 0;
    workspace->current_len = 0;
}

int MpiAllocatePackWorkspace(MpiSamToBamPackWorkspace *workspace,
                             int capacity,
                             int n_plans) {
    workspace->capacity = capacity;
    workspace->max_plans = n_plans;
    workspace->records = (bam1_t **)aligned_alloc_custom(64, (size_t)capacity * sizeof(bam1_t *));
    workspace->plans = (MpiSamToBamPackedPlan *)aligned_alloc_custom(
        64, (size_t)n_plans * sizeof(MpiSamToBamPackedPlan));
    if (!workspace->records || !workspace->plans) return -1;
    MpiResetPackWorkspace(workspace);
    return 0;
}

void MpiFreePackWorkspace(MpiSamToBamPackWorkspace *workspace) {
    if (workspace->records) aligned_free_custom((unsigned char *)workspace->records);
    if (workspace->plans) aligned_free_custom((unsigned char *)workspace->plans);
    workspace->records = nullptr;
    workspace->plans = nullptr;
    workspace->capacity = 0;
    workspace->max_plans = 0;
}

int MpiSealCurrentPackBlock(MpiSamToBamPackWorkspace *workspace) {
    if (workspace->current_records == 0) return 0;
    if (workspace->active_blocks >= workspace->max_plans) return -1;
    MpiSamToBamPackedPlan &plan = workspace->plans[workspace->active_blocks++];
    plan.records = workspace->current_begin;
    plan.n_records = workspace->current_records;
    plan.total_len = workspace->current_len;
    workspace->current_begin = workspace->records + workspace->total_records;
    workspace->current_records = 0;
    workspace->current_len = 0;
    return 0;
}

int MpiAppendRecordToPackWorkspace(MpiSamToBamPackWorkspace *workspace,
                                   bam1_t *record,
                                   uint32_t bam_len) {
    const uint32_t packed_len = bam_len + 4;
    if (packed_len > BGZF_BLOCK_SIZE) return -1;
    if (workspace->current_records > 0 &&
        workspace->current_len + packed_len > BGZF_BLOCK_SIZE) {
        return 1;
    }
    if (workspace->total_records >= workspace->capacity) return -2;
    if (workspace->current_records == 0) {
        workspace->current_begin = workspace->records + workspace->total_records;
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
                                 MpiSamToBamStats *stats) {
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

size_t MpiRoundUpSize(size_t value, size_t align) {
    if (align == 0) return value;
    size_t rem = value % align;
    if (rem == 0) return value;
    size_t add = align - rem;
    if (value > SIZE_MAX - add) return SIZE_MAX;
    return value + add;
}

size_t MpiChooseSamChunkSize(size_t local_input_size) {
    const size_t min_chunk = 128 * 1024;
    const size_t max_chunk = SAM_CHUNK_SIZE;
    const size_t align = 64 * 1024;

    if (local_input_size == 0) return min_chunk;
    size_t target = local_input_size / 64 + (local_input_size % 64 != 0);
    if (target < min_chunk) target = min_chunk;
    target = MpiRoundUpSize(target, align);
    if (target < min_chunk) target = min_chunk;
    if (target > max_chunk) target = max_chunk;
    return target;
}

int MpiSplitSamChunks(MemReader &reader,
                      MpiSamParseBatch *batch,
                      size_t chunk_size,
                      int *active_chunks,
                      int *max_chunk_len) {
    int count = 0;
    *max_chunk_len = 0;
    if (chunk_size == 0) chunk_size = SAM_CHUNK_SIZE;
    for (int i = 0; i < 64; ++i) {
        MpiSamParseChunk *chunk = &batch->chunks[i];
        chunk->src_ptr = nullptr;
        chunk->src_len = 0;
        chunk->text_len = 0;
        chunk->count = 0;
        chunk->parse_fast_records = 0;
        chunk->parse_fallback_records = 0;
        chunk->parse_total_cycles = 0;
        chunk->parse_core_cycles = 0;
        chunk->parse_aux_cycles = 0;
        chunk->parse_cg_cycles = 0;
        chunk->parse_fallback_cycles = 0;

        if (reader.pos >= reader.size) continue;

        size_t start_pos = reader.pos;
        size_t remaining = reader.size - start_pos;
        size_t end_pos = start_pos + (remaining > chunk_size ? chunk_size : remaining);
        if (end_pos < reader.size) {
            while (end_pos < reader.size && reader.base[end_pos] != '\n') end_pos++;
            if (end_pos < reader.size && reader.base[end_pos] == '\n') end_pos++;
        } else {
            end_pos = reader.size;
        }

        size_t read_len = end_pos - start_pos;
        if (read_len + 1 > CHUNK_BUFFER_SIZE) {
            fprintf(stderr,
                    "ERROR: MPI SAM chunk exceeds buffer. limit=%zu actual=%zu chunk=%d.\n",
                    (size_t)CHUNK_BUFFER_SIZE, read_len + 1, i);
            return -1;
        }
        chunk->src_ptr = reader.base + start_pos;
        chunk->src_len = read_len;
        reader.pos = end_pos;
        count++;
        if ((int)read_len > *max_chunk_len) *max_chunk_len = (int)read_len;
    }
    *active_chunks = count;
    return 0;
}

void MpiAccumulateSamParseDetail(MpiSamParseBatch *batch,
                                 int active_chunks,
                                 double parse_wall,
                                 MpiSamToBamStats *stats) {
    uint64_t max_cycles = 0;
    MpiSamParseChunk *critical = nullptr;
    for (int i = 0; i < active_chunks; ++i) {
        MpiSamParseChunk *chunk = &batch->chunks[i];
        if (chunk->parse_total_cycles > max_cycles) {
            max_cycles = chunk->parse_total_cycles;
            critical = chunk;
        }
    }
    if (!critical || max_cycles == 0) {
        stats->t_parse_other += parse_wall;
        return;
    }

    double scale = parse_wall / (double)max_cycles;
    double t_core = (double)critical->parse_core_cycles * scale;
    double t_aux = (double)critical->parse_aux_cycles * scale;
    double t_cg = (double)critical->parse_cg_cycles * scale;
    double t_fallback = (double)critical->parse_fallback_cycles * scale;
    double detail_sum = t_core + t_aux + t_cg + t_fallback;
    double t_other = parse_wall - detail_sum;
    if (t_other < 0.0) t_other = 0.0;

    stats->t_parse_core += t_core;
    stats->t_parse_aux += t_aux;
    stats->t_parse_cg += t_cg;
    stats->t_parse_fallback += t_fallback;
    stats->t_parse_other += t_other;
}

} // namespace

int FusedSamToBamMPI(MemReader &reader, MemWriter &mem_writer,
                     sam_hdr_t *hdr,
                     int compress_level,
                     MpiSamToBamStats *stats) {
    const int NB = 64;
    MpiSamToBamStats local_stats = {};
    if (!stats) stats = &local_stats;
    double fused_t0 = GetTime();
    const size_t sam_chunk_size = MpiChooseSamChunkSize(reader.size);

    // double alloc_t0 = GetTime();
    MpiSamParseBatch *batch = (MpiSamParseBatch *)aligned_alloc_custom(64, sizeof(MpiSamParseBatch));
    char *batch_text_storage = (char *)aligned_alloc_custom(64, (size_t)NB * CHUNK_BUFFER_SIZE);
    uint32_t *batch_bam_lens_storage = (uint32_t *)aligned_alloc_custom(
        64, (size_t)NB * MAX_BAMS_PER_CHUNK * sizeof(uint32_t));
    Comp_Para comp_a[NB], comp_b[NB];
    MpiSamToBamBlockSet comp_un_a = {};
    MpiSamToBamBlockSet comp_un_b = {};
    MpiSamToBamBlockSet out_a = {};
    MpiSamToBamBlockSet out_b = {};
    MpiSamToBamRecordSet record_set = {};
    MpiSamToBamPackWorkspace pack_workspace = {};

    if (!batch || !batch_text_storage || !batch_bam_lens_storage ||
        MpiInitSamParseBatch(batch, batch_text_storage, batch_bam_lens_storage) != 0 ||
        MpiAllocateSamToBamBlockSet(&comp_un_a, NB) != 0 ||
        MpiAllocateSamToBamBlockSet(&comp_un_b, NB) != 0 ||
        MpiAllocateSamToBamBlockSet(&out_a, NB) != 0 ||
        MpiAllocateSamToBamBlockSet(&out_b, NB) != 0 ||
        MpiAllocateSamToBamRecordSet(&record_set, FUSED_SAM2BAM_BAM_POOL_SIZE) != 0 ||
        MpiAllocatePackWorkspace(&pack_workspace, FUSED_SAM2BAM_BAM_POOL_SIZE, NB) != 0) {
        fprintf(stderr, "ERROR: failed to allocate MPI sam2bam 1CG workspace.\n");
        MpiDestroySamParseBatch(batch);
        if (batch) aligned_free_custom((unsigned char *)batch);
        if (batch_text_storage) aligned_free_custom((unsigned char *)batch_text_storage);
        if (batch_bam_lens_storage) aligned_free_custom((unsigned char *)batch_bam_lens_storage);
        MpiFreeSamToBamBlockSet(&comp_un_a);
        MpiFreeSamToBamBlockSet(&comp_un_b);
        MpiFreeSamToBamBlockSet(&out_a);
        MpiFreeSamToBamBlockSet(&out_b);
        MpiFreeSamToBamRecordSet(&record_set);
        MpiFreePackWorkspace(&pack_workspace);
        // stats->t_alloc_init += GetTime() - alloc_t0;
        stats->t_fused_total += GetTime() - fused_t0;
        return -1;
    }
    batch->hdr = hdr;

    for (int i = 0; i < NB; ++i) {
        MpiInitEmptyCompPara(&comp_a[i], i);
        MpiInitEmptyCompPara(&comp_b[i], i);
    }
    // stats->t_alloc_init += GetTime() - alloc_t0;

    Comp_Para *comp_active = comp_a;
    Comp_Para *comp_pending = comp_b;
    MpiSamToBamBlockSet *comp_un_active = &comp_un_a;
    MpiSamToBamBlockSet *comp_un_pending = &comp_un_b;
    MpiSamToBamBlockSet *out_active = &out_a;
    MpiSamToBamBlockSet *out_pending = &out_b;
    bool has_pending = false;
    int ret_code = -1;

    auto flush_pending = [&]() -> int {
        if (!has_pending) return 0;
        double write_t0 = GetTime();
        for (int k = 0; k < NB; ++k) {
            if (comp_pending[k].status == 0 && comp_pending[k].output_block) {
                if (MpiWriteBlockToMem(mem_writer, comp_pending[k].output_block) != 0) {
                    fprintf(stderr, "ERROR: failed to append MPI sam2bam compressed block.\n");
                    return -1;
                }
                stats->bgzf_blocks++;
            }
            MpiInitEmptyCompPara(&comp_pending[k], k);
        }
        has_pending = false;
        stats->t_write += GetTime() - write_t0;
        return 0;
    };

    auto do_compress = [&](int active_output_blocks) -> int {
        if (active_output_blocks <= 0) return 0;
        // double setup_t0 = GetTime();
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
        // stats->t_compress_setup += GetTime() - setup_t0;

        double compress_t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_compressfunc, comp_active, 1);
        #ifdef ENABLE_MASKING
        int flush_ret = flush_pending();
        #endif
        athread_join();
        double compress_wall = GetTime() - compress_t0;
        stats->t_compress += compress_wall;
        MpiAccumulateCompressDetail(comp_active, active_output_blocks, compress_wall, stats);

        #ifndef ENABLE_MASKING
        int flush_ret = flush_pending();
        #endif
        if (flush_ret != 0) return -1;

        // double check_t0 = GetTime();
        for (int k = 0; k < active_output_blocks; ++k) {
            if (comp_active[k].status != 0) {
                // stats->t_status_check += GetTime() - check_t0;
                fprintf(stderr, "ERROR: MPI sam2bam compress failed on output block %d with status %d.\n",
                        k, comp_active[k].status);
                return -1;
            }
        }
        // stats->t_status_check += GetTime() - check_t0;

        has_pending = true;
        std::swap(comp_active, comp_pending);
        std::swap(comp_un_active, comp_un_pending);
        std::swap(out_active, out_pending);
        MpiResetPackWorkspace(&pack_workspace);
        stats->compress_groups++;
        return 0;
    };

    while (reader.pos < reader.size) {
        int active_chunks = 0;
        int local_max_chunk_len = 0;
        // double split_t0 = GetTime();
        if (MpiSplitSamChunks(reader, batch, sam_chunk_size, &active_chunks, &local_max_chunk_len) != 0) goto cleanup;
        // stats->t_split += GetTime() - split_t0;
        if (active_chunks == 0) break;
        stats->input_chunks += active_chunks;

        double copy_t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_copy_and_count, batch, 1);
        #ifdef ENABLE_MASKING
        int flush_ret = flush_pending();
        #endif
        athread_join();
        stats->t_copy_count += GetTime() - copy_t0;

        #ifndef ENABLE_MASKING
        int flush_ret = flush_pending();
        #endif
        if (flush_ret != 0) goto cleanup;

        // double check_t0 = GetTime();
        int requested_records = 0;
        for (int i = 0; i < active_chunks; ++i) {
            MpiSamParseChunk *chunk = &batch->chunks[i];
            if (chunk->count > MAX_BAMS_PER_CHUNK) {
                // stats->t_status_check += GetTime() - check_t0;
                fprintf(stderr, "ERROR: MPI sam2bam chunk record limit exceeded. limit=%d actual=%d chunk=%d.\n",
                        MAX_BAMS_PER_CHUNK, chunk->count, i);
                goto cleanup;
            }
            requested_records += chunk->count;
            if (requested_records > FUSED_SAM2BAM_BAM_POOL_SIZE) {
                // stats->t_status_check += GetTime() - check_t0;
                fprintf(stderr, "ERROR: MPI sam2bam record pool exceeded. limit=%d actual=%d.\n",
                        FUSED_SAM2BAM_BAM_POOL_SIZE, requested_records);
                goto cleanup;
            }
        }
        // stats->t_status_check += GetTime() - check_t0;

        // double setup_t0 = GetTime();
        int offset = 0;
        for (int i = 0; i < active_chunks; ++i) {
            MpiSamParseChunk *chunk = &batch->chunks[i];
            chunk->bams = record_set.records + offset;
            chunk->bam_data = record_set.data + (size_t)offset * INIT_DATA_SIZE;
            offset += chunk->count;
        }
        // stats->t_setup_reset += GetTime() - setup_t0;

        double parse_t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_sam_parse_chunk, batch, 1);
        athread_join();
        double parse_wall = GetTime() - parse_t0;
        stats->t_parse += parse_wall;
        MpiAccumulateSamParseDetail(batch, active_chunks, parse_wall, stats);
        for (int i = 0; i < active_chunks; ++i) {
            MpiSamParseChunk *chunk = &batch->chunks[i];
            stats->parse_fast_records += chunk->parse_fast_records;
            stats->parse_fallback_records += chunk->parse_fallback_records;
        }

        double pack_t0 = GetTime();
        for (int i = 0; i < active_chunks; ++i) {
            MpiSamParseChunk *chunk = &batch->chunks[i];
            for (int j = 0; j < chunk->count; ++j) {
                bam1_t *b = chunk->bams + j;
                uint32_t bam_len = chunk->bam_lens[j];
                stats->total_records++;
                while (true) {
                    int ret = MpiAppendRecordToPackWorkspace(&pack_workspace, b, bam_len);
                    if (ret == 0) break;
                    if (ret == -1) {
                        fprintf(stderr, "ERROR: MPI sam2bam pack encountered an oversized BAM record.\n");
                        goto cleanup;
                    }
                    if (ret == -2) {
                        fprintf(stderr, "ERROR: MPI sam2bam pack workspace capacity exceeded.\n");
                        goto cleanup;
                    }
                    if (MpiSealCurrentPackBlock(&pack_workspace) != 0) {
                        fprintf(stderr, "ERROR: MPI sam2bam output block plan capacity exceeded.\n");
                        goto cleanup;
                    }
                    if (pack_workspace.active_blocks == NB) {
                        stats->t_pack += GetTime() - pack_t0;
                        if (do_compress(NB) != 0) goto cleanup;
                        pack_t0 = GetTime();
                    }
                }
            }
        }
        if (MpiSealCurrentPackBlock(&pack_workspace) != 0) {
            fprintf(stderr, "ERROR: MPI sam2bam final output block plan capacity exceeded.\n");
            goto cleanup;
        }
        stats->t_pack += GetTime() - pack_t0;
        if (pack_workspace.active_blocks > 0 && do_compress(pack_workspace.active_blocks) != 0) goto cleanup;
        stats->chunk_groups++;
    }

    if (flush_pending() != 0) goto cleanup;
    ret_code = 0;

cleanup:
    {
        // double free_t0 = GetTime();
        MpiDestroySamParseBatch(batch);
        if (batch) aligned_free_custom((unsigned char *)batch);
        if (batch_text_storage) aligned_free_custom((unsigned char *)batch_text_storage);
        if (batch_bam_lens_storage) aligned_free_custom((unsigned char *)batch_bam_lens_storage);
        MpiFreeSamToBamBlockSet(&comp_un_a);
        MpiFreeSamToBamBlockSet(&comp_un_b);
        MpiFreeSamToBamBlockSet(&out_a);
        MpiFreeSamToBamBlockSet(&out_b);
        MpiFreeSamToBamRecordSet(&record_set);
        MpiFreePackWorkspace(&pack_workspace);
        // stats->t_free_workspace += GetTime() - free_t0;
    }
    stats->t_fused_total += GetTime() - fused_t0;
    return ret_code;
}
