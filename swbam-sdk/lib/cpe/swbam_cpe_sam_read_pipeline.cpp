#include "swbam/cpe_sam_read_pipeline.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <stdint.h>
#include <vector>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

extern "C" void slave_mpi_copy_and_count();
extern "C" void slave_mpi_sam_parse_chunk();

namespace swbam {
namespace cpe {
namespace {

const int kBatchChunks = 64;

struct SamChunkPiece {
    char *text_buf;
    size_t text_len;
    int count;
    size_t estimated_bam_data;
    size_t max_line_len;
};

size_t RoundUp(size_t value, size_t align) {
    if (!align || value % align == 0) return value;
    const size_t add = align - value % align;
    return value > SIZE_MAX - add ? SIZE_MAX : value + add;
}

size_t ChooseChunkSize(size_t input_size) {
    const size_t minimum = 128 * 1024;
    if (!input_size) return minimum;
    size_t target = input_size / 64 + (input_size % 64 != 0);
    target = RoundUp(std::max(target, minimum), 64 * 1024);
    return std::min(target, static_cast<size_t>(SAM_CHUNK_SIZE));
}

class SamReadWorkspace {
public:
    SamReadWorkspace()
        : count_batch(nullptr), parse_batch(nullptr), text_storage(nullptr),
          bam_lens(nullptr), bam_offsets(nullptr), records(nullptr),
          arenas(kBatchChunks, nullptr), arena_capacity(kBatchChunks, 0) {}

    ~SamReadWorkspace() {
        if (count_batch) {
            for (int i = 0; i < kBatchChunks; ++i) {
                MpiSamParseChunk &chunk = count_batch->chunks[i];
                if (chunk.owns_text_buf && chunk.text_buf)
                    aligned_free_custom(reinterpret_cast<unsigned char *>(chunk.text_buf));
            }
            aligned_free_custom(reinterpret_cast<unsigned char *>(count_batch));
        }
        if (parse_batch) aligned_free_custom(reinterpret_cast<unsigned char *>(parse_batch));
        if (text_storage) aligned_free_custom(reinterpret_cast<unsigned char *>(text_storage));
        if (bam_lens) aligned_free_custom(reinterpret_cast<unsigned char *>(bam_lens));
        if (bam_offsets) aligned_free_custom(reinterpret_cast<unsigned char *>(bam_offsets));
        if (records) aligned_free_custom(reinterpret_cast<unsigned char *>(records));
        for (size_t i = 0; i < arenas.size(); ++i)
            if (arenas[i]) aligned_free_custom(arenas[i]);
    }

    int Initialize(const sam_hdr_t *header) {
        count_batch = reinterpret_cast<MpiSamParseBatch *>(
            aligned_alloc_custom(64, sizeof(MpiSamParseBatch)));
        parse_batch = reinterpret_cast<MpiSamParseBatch *>(
            aligned_alloc_custom(64, sizeof(MpiSamParseBatch)));
        if (count_batch) memset(count_batch, 0, sizeof(*count_batch));
        if (parse_batch) memset(parse_batch, 0, sizeof(*parse_batch));
        text_storage = reinterpret_cast<char *>(aligned_alloc_custom(
            64, static_cast<size_t>(kBatchChunks) * CHUNK_BUFFER_SIZE));
        bam_lens = reinterpret_cast<uint32_t *>(aligned_alloc_custom(
            64, static_cast<size_t>(kBatchChunks) * MAX_BAMS_PER_CHUNK * sizeof(uint32_t)));
        bam_offsets = reinterpret_cast<uint32_t *>(aligned_alloc_custom(
            64, static_cast<size_t>(kBatchChunks) * MAX_BAMS_PER_CHUNK * sizeof(uint32_t)));
        records = reinterpret_cast<bam1_t *>(aligned_alloc_custom(
            64, static_cast<size_t>(FUSED_SAM2BAM_BAM_POOL_SIZE) * sizeof(bam1_t)));
        if (!count_batch || !parse_batch || !text_storage || !bam_lens ||
            !bam_offsets || !records) return -1;
        memset(records, 0, static_cast<size_t>(FUSED_SAM2BAM_BAM_POOL_SIZE) * sizeof(bam1_t));
        count_batch->hdr = header;
        parse_batch->hdr = header;
        for (int i = 0; i < kBatchChunks; ++i) {
            MpiSamParseChunk &chunk = count_batch->chunks[i];
            chunk.text_buf = text_storage + static_cast<size_t>(i) * CHUNK_BUFFER_SIZE;
            chunk.text_buf_capacity = CHUNK_BUFFER_SIZE;
            chunk.bam_lens = bam_lens + static_cast<size_t>(i) * MAX_BAMS_PER_CHUNK;
            chunk.bam_offsets = bam_offsets + static_cast<size_t>(i) * MAX_BAMS_PER_CHUNK;
            MpiSamParseChunk &parsed = parse_batch->chunks[i];
            parsed.bam_lens = chunk.bam_lens;
            parsed.bam_offsets = chunk.bam_offsets;
        }
        return 0;
    }

    int EnsureArena(int slot, MpiSamParseChunk *chunk, size_t required) {
        required = RoundUp(std::max(required, static_cast<size_t>(INIT_DATA_SIZE)), 64);
        if (arena_capacity[slot] < required) {
            unsigned char *next = aligned_alloc_custom(64, required);
            if (!next) return -1;
            if (arenas[slot]) aligned_free_custom(arenas[slot]);
            arenas[slot] = next;
            arena_capacity[slot] = required;
        }
        chunk->bam_data = arenas[slot];
        chunk->bam_data_capacity = arena_capacity[slot];
        return 0;
    }

    MpiSamParseBatch *count_batch;
    MpiSamParseBatch *parse_batch;
    char *text_storage;
    uint32_t *bam_lens;
    uint32_t *bam_offsets;
    bam1_t *records;
    std::vector<unsigned char *> arenas;
    std::vector<size_t> arena_capacity;
};

int SplitChunks(MemReader *reader, SamReadWorkspace *workspace,
                size_t chunk_size, int *active) {
    *active = 0;
    for (int i = 0; i < kBatchChunks; ++i) {
        MpiSamParseChunk &chunk = workspace->count_batch->chunks[i];
        if (chunk.owns_text_buf && chunk.text_buf)
            aligned_free_custom(reinterpret_cast<unsigned char *>(chunk.text_buf));
        chunk.text_buf = workspace->text_storage + static_cast<size_t>(i) * CHUNK_BUFFER_SIZE;
        chunk.text_buf_capacity = CHUNK_BUFFER_SIZE;
        chunk.owns_text_buf = 0;
        chunk.src_ptr = nullptr;
        chunk.src_len = chunk.text_len = 0;
        chunk.count = chunk.status = 0;
        chunk.record_index = 0;
        chunk.actual_value = chunk.limit_value = 0;
        chunk.limit_id = BOUNDS_LIMIT_NONE;
        chunk.max_line_len = chunk.estimated_bam_data = chunk.bam_data_used = 0;
        chunk.parse_fast_records = chunk.parse_fallback_records = 0;
        chunk.parse_total_cycles = chunk.parse_core_cycles = chunk.parse_aux_cycles = 0;
        chunk.parse_cg_cycles = chunk.parse_fallback_cycles = 0;

        if (reader->pos >= reader->size) continue;
        const size_t start = reader->pos;
        const size_t remaining = reader->size - start;
        size_t end = start + std::min(remaining, chunk_size);
        if (end < reader->size) {
            while (end < reader->size && reader->base[end] != '\n') ++end;
            if (end < reader->size) ++end;
        }
        const size_t length = end - start;
        if (length + 1 > chunk.text_buf_capacity) {
            const size_t capacity = RoundUp(length + 1, 64);
            chunk.text_buf = reinterpret_cast<char *>(aligned_alloc_custom(64, capacity));
            if (!chunk.text_buf) {
                fprintf(stderr, "ERROR: MPI SAM chunk side buffer allocation failed. actual=%zu chunk=%d.\n",
                        length + 1, i);
                return -1;
            }
            chunk.text_buf_capacity = capacity;
            chunk.owns_text_buf = 1;
        }
        chunk.src_ptr = reader->base + start;
        chunk.src_len = length;
        reader->pos = end;
        ++*active;
    }
    return 0;
}

int BuildPieces(MpiSamParseBatch *batch, int active,
                std::vector<SamChunkPiece> *pieces) {
    pieces->clear();
    for (int i = 0; i < active; ++i) {
        MpiSamParseChunk &chunk = batch->chunks[i];
        if (chunk.count <= 0 || chunk.text_len == 0) continue;
        if (chunk.count <= MAX_BAMS_PER_CHUNK) {
            pieces->push_back({chunk.text_buf, chunk.text_len, chunk.count,
                               chunk.estimated_bam_data, chunk.max_line_len});
            continue;
        }
        char *const end = chunk.text_buf + chunk.text_len;
        char *begin = chunk.text_buf;
        int piece_count = 0;
        int total_count = 0;
        size_t estimated = 0;
        size_t max_line = 0;
        for (char *ptr = begin; ptr < end;) {
            char *line_end = ptr;
            while (line_end < end && *line_end != '\n' && *line_end != '\0') ++line_end;
            int line_len = static_cast<int>(line_end - ptr);
            if (line_len > 0 && ptr[line_len - 1] == '\r') --line_len;
            if (line_len > 0 && piece_count == MAX_BAMS_PER_CHUNK) {
                pieces->push_back({begin, static_cast<size_t>(ptr - begin),
                                   piece_count, estimated, max_line});
                begin = ptr;
                piece_count = 0;
                estimated = max_line = 0;
            }
            if (line_len > 0) {
                ++piece_count;
                ++total_count;
                estimated += static_cast<size_t>(line_len) + 64u;
                max_line = std::max(max_line, static_cast<size_t>(line_len));
            }
            ptr = line_end < end && *line_end == '\n' ? line_end + 1 : line_end;
        }
        if (piece_count > 0) {
            pieces->push_back({begin, static_cast<size_t>(end - begin),
                               piece_count, estimated, max_line});
        }
        if (total_count != chunk.count) {
            fprintf(stderr,
                    "ERROR: MPI sam2bam chunk split count mismatch. chunk=%d counted=%d split=%d.\n",
                    i, chunk.count, total_count);
            return -1;
        }
    }
    return 0;
}

void AccumulateParseDetail(const MpiSamParseBatch *batch, int count,
                           double wall, CpeSamReadTiming *timing) {
    const MpiSamParseChunk *critical = nullptr;
    uint64_t cycles = 0;
    for (int i = 0; i < count; ++i) {
        if (batch->chunks[i].parse_total_cycles > cycles) {
            cycles = batch->chunks[i].parse_total_cycles;
            critical = &batch->chunks[i];
        }
    }
    if (!critical || !cycles) {
        timing->parse_other += wall;
        return;
    }
    const double scale = wall / static_cast<double>(cycles);
    const double core = critical->parse_core_cycles * scale;
    const double aux = critical->parse_aux_cycles * scale;
    const double cg = critical->parse_cg_cycles * scale;
    const double fallback = critical->parse_fallback_cycles * scale;
    timing->parse_core += core;
    timing->parse_aux += aux;
    timing->parse_cg += cg;
    timing->parse_fallback += fallback;
    timing->parse_other += std::max(0.0, wall - core - aux - cg - fallback);
}

} // namespace

CpeSamReadTiming::CpeSamReadTiming()
    : input_chunks(0), chunk_groups(0), parse_fast_records(0),
      parse_fallback_records(0), split(0.0), copy_count(0.0), parse(0.0),
      parse_core(0.0), parse_aux(0.0), parse_cg(0.0),
      parse_fallback(0.0), parse_other(0.0) {}

SamReadKernelSpec DefaultSamReadKernelSpec() {
    return {reinterpret_cast<void *>(slave_mpi_copy_and_count),
            reinterpret_cast<void *>(slave_mpi_sam_parse_chunk)};
}

int RunCpeSamReadPipeline(MemReader *reader, const sam_hdr_t *header,
                          SamParsedBatchPostProcessor *post_processor,
                          CpeSamReadTiming *timing,
                          const SamReadKernelSpec &kernels) {
    if (!reader || !header || !post_processor || !timing ||
        !kernels.copy_count_entry || !kernels.parse_entry) return -1;
    *timing = CpeSamReadTiming();
    SamReadWorkspace workspace;
    if (workspace.Initialize(header) != 0) {
        fprintf(stderr, "ERROR: failed to allocate SAM read pipeline workspace.\n");
        return -1;
    }
    const size_t chunk_size = ChooseChunkSize(reader->size);
    std::vector<SamChunkPiece> pieces;
    while (reader->pos < reader->size) {
        int active = 0;
        const double split_t0 = GetTime();
        const int split_ret = SplitChunks(reader, &workspace, chunk_size, &active);
        timing->split += GetTime() - split_t0;
        if (split_ret != 0) return -1;
        if (!active) break;
        timing->input_chunks += active;

        const double copy_t0 = GetTime();
        int flush_ret = 0;
#ifdef PLATFORM_SUNWAY
        __real_athread_spawn(kernels.copy_count_entry, workspace.count_batch, 1);
        flush_ret = post_processor->FlushPendingOutput();
        athread_join();
#else
        return -1;
#endif
        timing->copy_count += GetTime() - copy_t0;
        if (flush_ret != 0) return -1;
        for (int i = 0; i < active; ++i) {
            const MpiSamParseChunk &chunk = workspace.count_batch->chunks[i];
            if (chunk.status != 0) {
                fprintf(stderr,
                        "ERROR: MPI sam2bam copy/count failed. chunk=%d status=%d limit=%lld actual=%lld.\n",
                        i, chunk.status, chunk.limit_value, chunk.actual_value);
                return -1;
            }
        }
        if (BuildPieces(workspace.count_batch, active, &pieces) != 0) return -1;

        for (size_t start = 0; start < pieces.size();) {
            size_t end = start;
            int wave_records = 0;
            while (end < pieces.size() && end - start < kBatchChunks) {
                if (wave_records > 0 &&
                    wave_records + pieces[end].count > FUSED_SAM2BAM_BAM_POOL_SIZE) break;
                wave_records += pieces[end].count;
                ++end;
            }
            if (end == start) return -1;
            const int wave_chunks = static_cast<int>(end - start);
            int offset = 0;
            for (int i = 0; i < kBatchChunks; ++i) {
                MpiSamParseChunk &chunk = workspace.parse_batch->chunks[i];
                chunk.status = chunk.parse_fast_records = chunk.parse_fallback_records = 0;
                chunk.parse_total_cycles = chunk.parse_core_cycles = 0;
                chunk.parse_aux_cycles = chunk.parse_cg_cycles = chunk.parse_fallback_cycles = 0;
                chunk.bam_data_used = 0;
                chunk.record_index = 0;
                chunk.actual_value = chunk.limit_value = 0;
                chunk.limit_id = BOUNDS_LIMIT_NONE;
                if (i < wave_chunks) {
                    const SamChunkPiece &piece = pieces[start + static_cast<size_t>(i)];
                    chunk.text_buf = piece.text_buf;
                    chunk.text_buf_capacity = chunk.text_len = piece.text_len;
                    chunk.owns_text_buf = 0;
                    chunk.count = piece.count;
                    chunk.estimated_bam_data = piece.estimated_bam_data;
                    chunk.max_line_len = piece.max_line_len;
                    chunk.bams = workspace.records + offset;
                    size_t need = piece.text_len + static_cast<size_t>(piece.count) * 96u + 4096u;
                    if (piece.estimated_bam_data > need) need = piece.estimated_bam_data + 4096u;
                    if (workspace.EnsureArena(i, &chunk, need) != 0) {
                        fprintf(stderr, "ERROR: MPI sam2bam failed to allocate BAM arena for parse slot %d.\n", i);
                        return -1;
                    }
                    offset += piece.count;
                } else {
                    chunk.count = 0;
                    chunk.text_len = 0;
                    chunk.text_buf = nullptr;
                    chunk.text_buf_capacity = 0;
                    chunk.owns_text_buf = 0;
                    chunk.bams = nullptr;
                    chunk.bam_data = nullptr;
                    chunk.bam_data_capacity = 0;
                    chunk.estimated_bam_data = chunk.max_line_len = 0;
                }
            }

            for (int attempt = 0; attempt < 2; ++attempt) {
                const double parse_t0 = GetTime();
#ifdef PLATFORM_SUNWAY
                __real_athread_spawn(kernels.parse_entry, workspace.parse_batch, 1);
                athread_join();
#else
                return -1;
#endif
                const double wall = GetTime() - parse_t0;
                timing->parse += wall;
                AccumulateParseDetail(workspace.parse_batch, wave_chunks, wall, timing);
                bool retry = false;
                for (int i = 0; i < wave_chunks; ++i) {
                    MpiSamParseChunk &chunk = workspace.parse_batch->chunks[i];
                    if (chunk.status == -3 &&
                        chunk.limit_id == BOUNDS_LIMIT_INIT_DATA_SIZE && attempt == 0) {
                        const size_t need = std::max(chunk.bam_data_capacity * 2u,
                            static_cast<size_t>(chunk.actual_value) + 4096u);
                        if (workspace.EnsureArena(i, &chunk, need) != 0) return -1;
                        retry = true;
                    } else if (chunk.status != 0) {
                        fprintf(stderr,
                                "ERROR: MPI sam2bam parse failed. chunk=%d status=%d limit=%lld actual=%lld.\n",
                                i, chunk.status, chunk.limit_value, chunk.actual_value);
                        return -1;
                    }
                }
                if (!retry) break;
                for (int i = 0; i < wave_chunks; ++i) {
                    workspace.parse_batch->chunks[i].count = pieces[start + static_cast<size_t>(i)].count;
                    workspace.parse_batch->chunks[i].text_len = pieces[start + static_cast<size_t>(i)].text_len;
                    workspace.parse_batch->chunks[i].bam_data_used = 0;
                }
            }

            for (int i = 0; i < wave_chunks; ++i) {
                timing->parse_fast_records += workspace.parse_batch->chunks[i].parse_fast_records;
                timing->parse_fallback_records += workspace.parse_batch->chunks[i].parse_fallback_records;
            }
            if (post_processor->PostProcessParsedBatch(
                    workspace.parse_batch->chunks, static_cast<size_t>(wave_chunks)) != 0) return -1;
            start = end;
        }
        ++timing->chunk_groups;
    }
    return 0;
}

int RunCpeSamReadPipeline(MemReader *reader, const sam_hdr_t *header,
                          SamParsedBatchPostProcessor *post_processor,
                          CpeSamReadTiming *timing) {
    return RunCpeSamReadPipeline(reader, header, post_processor, timing,
                                 DefaultSamReadKernelSpec());
}

} // namespace cpe
} // namespace swbam
