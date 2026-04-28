#include "swbam_cgs.h"

#include <algorithm>
#include <cstdio>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

extern "C" {
    void slave_cgs_copy_and_count();
    void slave_cgs_sam_parse_chunk();
    void slave_cgs_compressfunc();
}

namespace {

int InitSamParseBatchCGS(CgsCrossAllocator &alloc, CgsSamParseBatch *batch, sam_hdr_t *hdr) {
    batch->hdr = hdr;
    char *text_storage = (char *)alloc.allocate(64, (size_t)CGS_NB * CGS_SAM_CHUNK_BUFFER_SIZE);
    uint32_t *bam_lens = (uint32_t *)alloc.allocate(
        64, (size_t)CGS_NB * CGS_SAM2BAM_MAX_BAMS_PER_CHUNK * sizeof(uint32_t));
    if (!text_storage || !bam_lens) return -1;

    for (int i = 0; i < CGS_NB; ++i) {
        CgsSamParseChunk *chunk = &batch->chunks[i];
        chunk->src_ptr = nullptr;
        chunk->src_len = 0;
        chunk->text_buf = text_storage + (size_t)i * CGS_SAM_CHUNK_BUFFER_SIZE;
        chunk->text_len = 0;
        chunk->bams = nullptr;
        chunk->bam_lens = bam_lens + (size_t)i * CGS_SAM2BAM_MAX_BAMS_PER_CHUNK;
        chunk->count = 0;
        chunk->status = -1;
    }
    return 0;
}

int SplitSamChunksCGS(MemReader &reader, CgsSamParseBatch *batch, int *active_chunks, int *max_chunk_len) {
    int count = 0;
    *max_chunk_len = 0;
    for (int i = 0; i < CGS_NB; ++i) {
        CgsSamParseChunk *chunk = &batch->chunks[i];
        chunk->src_ptr = nullptr;
        chunk->src_len = 0;
        chunk->text_len = 0;
        chunk->count = 0;
        chunk->status = -1;

        if (reader.pos >= reader.size) continue;

        size_t start_pos = reader.pos;
        size_t end_pos = start_pos + CGS_SAM_CHUNK_SIZE;
        if (end_pos < reader.size) {
            while (end_pos < reader.size && reader.base[end_pos] != '\n') end_pos++;
            if (end_pos < reader.size && reader.base[end_pos] == '\n') end_pos++;
        } else {
            end_pos = reader.size;
        }

        size_t read_len = end_pos - start_pos;
        if (read_len + 1 > CGS_SAM_CHUNK_BUFFER_SIZE) {
            fprintf(stderr,
                    "ERROR: CGS SAM chunk exceeds buffer. limit=%zu actual=%zu chunk=%d. "
                    "Use shorter SAM lines or reduce CGS_SAM_CHUNK_SIZE.\n",
                    CGS_SAM_CHUNK_BUFFER_SIZE, read_len + 1, i);
            return -1;
        }
        chunk->src_ptr = reader.base + start_pos;
        chunk->src_len = read_len;
        chunk->status = 0;
        reader.pos = end_pos;
        count++;
        if ((int)read_len > *max_chunk_len) *max_chunk_len = (int)read_len;
    }
    *active_chunks = count;
    return 0;
}

} // namespace

int FusedSamToBamCGS(MemReader reader, MemWriter &mem_writer, sam_hdr_t *hdr) {
    if (CheckCgsResources() != 0) return -1;

    double t0 = GetTime();
    double t_alloc_init = 0, t_split = 0, t_setup_reset = 0, t_compress_setup = 0, t_status_check = 0;
    double t_copy_count = 0, t_parse = 0, t_pack = 0, t_compress = 0, t_write = 0;

    double alloc_t0 = GetTime();
    CgsCrossAllocator alloc;
    CgsSamParseBatch *batch = (CgsSamParseBatch *)alloc.allocate(64, sizeof(CgsSamParseBatch));
    Comp_Para *comp_a = (Comp_Para *)alloc.allocate(64, (size_t)CGS_NB * sizeof(Comp_Para));
    Comp_Para *comp_b = (Comp_Para *)alloc.allocate(64, (size_t)CGS_NB * sizeof(Comp_Para));
    CgsBlockSet comp_un_a = {};
    CgsBlockSet comp_un_b = {};
    CgsBlockSet out_a = {};
    CgsBlockSet out_b = {};
    CgsRecordSet record_set = {};
    CgsPackWorkspace pack_workspace = {};

    if (!batch || !comp_a || !comp_b ||
        InitSamParseBatchCGS(alloc, batch, hdr) != 0 ||
        AllocBlockSet(alloc, &comp_un_a, CGS_NB) != 0 ||
        AllocBlockSet(alloc, &comp_un_b, CGS_NB) != 0 ||
        AllocBlockSet(alloc, &out_a, CGS_NB) != 0 ||
        AllocBlockSet(alloc, &out_b, CGS_NB) != 0 ||
        AllocRecordSet(alloc, &record_set, CGS_SAM2BAM_RECORD_POOL_SIZE) != 0 ||
        AllocPackWorkspace(alloc, &pack_workspace, CGS_SAM2BAM_RECORD_POOL_SIZE, CGS_NB) != 0) {
        fprintf(stderr, "ERROR: failed to allocate CGS sam2bam cross memory. Increase -cross_size and retry.\n");
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

    long long bam1_t_nums = 0;
    long long group_nums = 0;
    long long bgzf_nums = 0;
    long long chunk_groups = 0;
    int max_lines_per_chunk = 0;
    int max_chunk_len = 0;

    auto flush_pending = [&]() -> int {
        double write_t0 = GetTime();
        if (!has_pending) return 0;
        for (int k = 0; k < CGS_NB; ++k) {
            if (comp_pending[k].status == 0 && comp_pending[k].output_block) {
                if (WriteBlockToMemCGS(mem_writer, comp_pending[k].output_block) != 0) {
                    fprintf(stderr, "ERROR: failed to append CGS compressed block to memory writer.\n");
                    return -1;
                }
                bgzf_nums++;
            }
            InitEmptyCompParaCGS(&comp_pending[k], k);
        }
        has_pending = false;
        t_write += GetTime() - write_t0;
        return 0;
    };

    auto do_compress = [&](int active_output_blocks) -> int {
        if (active_output_blocks <= 0) return 0;
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
        athread_join_cgs();
        t_compress += GetTime() - compress_t0;

        double check_t0 = GetTime();
        for (int k = 0; k < active_output_blocks; ++k) {
            if (comp_active[k].status != 0) {
                fprintf(stderr, "ERROR: CGS sam2bam compress failed on output block %d with status %d.\n",
                        k, comp_active[k].status);
                return -1;
            }
        }
        t_status_check += GetTime() - check_t0;

        has_pending = true;
        std::swap(comp_active, comp_pending);
        std::swap(comp_un_active, comp_un_pending);
        std::swap(out_active, out_pending);
        ResetPackWorkspace(&pack_workspace);
        group_nums++;
        return 0;
    };

    int active_chunks = 0;
    while (reader.pos < reader.size) {
        int local_max_chunk_len = 0;
        double split_t0 = GetTime();
        if (SplitSamChunksCGS(reader, batch, &active_chunks, &local_max_chunk_len) != 0) return -1;
        t_split += GetTime() - split_t0;
        if (active_chunks == 0) break;
        if (local_max_chunk_len > max_chunk_len) max_chunk_len = local_max_chunk_len;

        double copy_t0 = GetTime();
        __real_athread_spawn_cgs((void *)slave_cgs_copy_and_count, batch, 1);
        if (flush_pending() != 0) return -1;
        athread_join_cgs();
        t_copy_count += GetTime() - copy_t0;

        double check_t0 = GetTime();
        int requested_records = 0;
        for (int i = 0; i < active_chunks; ++i) {
            CgsSamParseChunk *chunk = &batch->chunks[i];
            if (chunk->status != 0) {
                fprintf(stderr, "ERROR: CGS sam2bam copy/count failed on chunk %d with status %d.\n",
                        i, chunk->status);
                return -1;
            }
            if (chunk->count > CGS_SAM2BAM_MAX_BAMS_PER_CHUNK) {
                fprintf(stderr, "ERROR: CGS sam2bam chunk record limit exceeded. limit=%d actual=%d chunk=%d.\n",
                        CGS_SAM2BAM_MAX_BAMS_PER_CHUNK, chunk->count, i);
                return -1;
            }
            if (chunk->count > max_lines_per_chunk) max_lines_per_chunk = chunk->count;
            requested_records += chunk->count;
            if (requested_records > CGS_SAM2BAM_RECORD_POOL_SIZE) {
                fprintf(stderr, "ERROR: CGS sam2bam record pool exceeded. limit=%d actual=%d.\n",
                        CGS_SAM2BAM_RECORD_POOL_SIZE, requested_records);
                return -1;
            }
        }
        t_status_check += GetTime() - check_t0;

        double setup_t0 = GetTime();
        int offset = 0;
        for (int i = 0; i < active_chunks; ++i) {
            CgsSamParseChunk *chunk = &batch->chunks[i];
            chunk->bams = record_set.ptrs + offset;
            offset += chunk->count;
        }
        t_setup_reset += GetTime() - setup_t0;

        double parse_t0 = GetTime();
        __real_athread_spawn_cgs((void *)slave_cgs_sam_parse_chunk, batch, 1);
        athread_join_cgs();
        t_parse += GetTime() - parse_t0;

        check_t0 = GetTime();
        for (int i = 0; i < active_chunks; ++i) {
            CgsSamParseChunk *chunk = &batch->chunks[i];
            if (chunk->status != 0) {
                fprintf(stderr, "ERROR: CGS sam2bam parse failed on chunk %d with status %d.\n",
                        i, chunk->status);
                return -1;
            }
        }
        t_status_check += GetTime() - check_t0;

        ResetPackWorkspace(&pack_workspace);
        double pack_t0 = GetTime();
        for (int i = 0; i < active_chunks; ++i) {
            CgsSamParseChunk *chunk = &batch->chunks[i];
            for (int j = 0; j < chunk->count; ++j) {
                bam1_t *b = chunk->bams[j];
                uint32_t bam_len = chunk->bam_lens[j];
                bam1_t_nums++;
                while (true) {
                    int ret = AppendRecordToPackWorkspace(&pack_workspace, b, bam_len);
                    if (ret == 0) break;
                    if (ret == -1) {
                        fprintf(stderr, "ERROR: CGS sam2bam pack encountered an oversized BAM record.\n");
                        return -1;
                    }
                    if (ret == -2) {
                        fprintf(stderr, "ERROR: CGS sam2bam pack workspace capacity exceeded.\n");
                        return -1;
                    }
                    if (SealCurrentPackBlock(&pack_workspace) != 0) {
                        fprintf(stderr, "ERROR: CGS sam2bam output block plan capacity exceeded.\n");
                        return -1;
                    }
                    if (pack_workspace.active_blocks == CGS_NB) {
                        t_pack += GetTime() - pack_t0;
                        if (do_compress(CGS_NB) != 0) return -1;
                        pack_t0 = GetTime();
                    }
                }
            }
        }
        if (SealCurrentPackBlock(&pack_workspace) != 0) {
            fprintf(stderr, "ERROR: CGS sam2bam final output block plan capacity exceeded.\n");
            return -1;
        }
        t_pack += GetTime() - pack_t0;
        if (pack_workspace.active_blocks > 0 && do_compress(pack_workspace.active_blocks) != 0) return -1;
        chunk_groups++;
    }

    if (flush_pending() != 0) return -1;

    printf("FusedSamToBamCGS finished. bam1_t=%lld, groups=%lld, chunk_groups=%lld, bgzf=%lld, cost %.6f\n",
           bam1_t_nums, group_nums, chunk_groups, bgzf_nums, GetTime() - t0);
    printf("  copy_count_slave=%lf  parse_slave=%lf  pack=%lf  compress_slave=%lf  write=%lf\n",
           t_copy_count, t_parse, t_pack, t_compress, t_write);
    printf("  alloc_init=%lf  split_chunks=%lf  setup_reset=%lf  compress_setup=%lf  status_check=%lf\n",
           t_alloc_init, t_split, t_setup_reset, t_compress_setup, t_status_check);
    printf("  cgs_threads=%d  sam_chunk_size=%zu  max_lines_per_chunk=%d  max_chunk_len=%d  output_bgzf_blocks=%lld\n",
           CGS_NB, CGS_SAM_CHUNK_SIZE, max_lines_per_chunk, max_chunk_len, bgzf_nums);

    return 0;
}
