#include "swbam_mpi.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

#include <mpi.h>

extern "C" {
    void slave_mpi_decompress_bam2bam_passthrough();
    void slave_mpi_markdup_extract();
    void slave_mpi_markdup_rewrite();
    void slave_mpi_compressfunc();
    void slave_mpi_sort_compress_payload();
}

int MpiCommonLoadFileToMemory(const std::string &path,
                              char **data, size_t *size);
int MpiCommonScanBgzfBlocksInMemory(
    const char *base, size_t size, long long body_start,
    std::vector<long long> *offsets,
    std::vector<long long> *lengths);
int MpiCommonSelectBlockRangeFromMemory(
    char *base, size_t input_size,
    const std::vector<long long> &offsets,
    const std::vector<long long> &lengths,
    long long begin, long long end,
    char **data, size_t *size);
int MpiCommonDumpMemoryToFile(const std::string &path,
                              const char *data, size_t size);
int MpiCommonInitMemWriter(MemWriter &w, size_t cap);
int MpiCommonBuildBamHeaderMemory(
    sam_hdr_t *hdr, int compress_level,
    char **data, size_t *size);
int MpiMarkdupFindDuplicatesStreamingHash(
    std::vector<MpiMarkdupCandidateShared> *owner_candidates,
    const std::vector<unsigned char> &owner_qnames,
    int comm_size,
    std::vector<std::vector<uint64_t> > *duplicates_by_source,
    MpiMarkdupStats *stats);

namespace {

const int kMarkdupNB = 64;
const int kMarkdupTagCount = 8;
const int kMarkdupExchangeChunk = INT_MAX / 2;

const unsigned char kMarkdupBgzfEofBlock[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
};

static int MdEnvFlagEnabled(const char *name) {
    const char *value = getenv(name);
    return value && value[0] && value[0] != '0';
}

struct MdBlockSet {
    bam_block *blocks;
    unsigned char *data;
    int n;
};

struct MdRecordSet {
    bam1_t *records;
    unsigned char *data;
    bam1_t **ptrs;
    uint32_t *bam_lens;
    int records_per_block;
    int n_blocks;
    size_t arena_stride;
};

struct MdPackPlan {
    bam1_t **records;
    int n_records;
    uint32_t total_len;
};

struct MdPackWorkspace {
    std::vector<bam1_t *> records;
    MdPackPlan plans[kMarkdupNB];
    int active_blocks;
    int total_records;
    bam1_t **current_begin;
    int current_records;
    uint32_t current_len;
};

struct MdCandidateWorkspace {
    MpiMarkdupCandidateShared *candidates;
    unsigned char *qnames;
    int candidates_per_block;
    size_t qname_stride;
};

struct MdMarkdupCandidateBufferDeleter {
    void operator()(MpiMarkdupCandidateShared *p) const {
        if (p) aligned_free_custom((unsigned char *)p);
    }
};

struct MdMarkdupKeyBufferDeleter {
    void operator()(MpiMarkdupKeyShared *p) const {
        if (p) aligned_free_custom((unsigned char *)p);
    }
};

struct MdByteBufferDeleter {
    void operator()(unsigned char *p) const {
        if (p) aligned_free_custom(p);
    }
};

struct MdDuplicateId {
    int32_t source_rank;
    uint32_t pad;
    uint64_t ordinal;
};

struct MdCoordBoundary {
    int has_records;
    int32_t first_ref;
    int64_t first_coord;
    int32_t last_ref;
    int64_t last_coord;
};

struct MdRemoteCandidateRef {
    size_t index;
    int owner;
};

struct MdBucketRange {
    uint64_t code;
    size_t begin;
    size_t end;
};

static int MdAllRanksOk(int local_ok) {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    return global_ok;
}

static int MdAllRanksOkTimed(int local_ok, MpiMarkdupStats *stats) {
    double t0 = GetTime();
    const int global_ok = MdAllRanksOk(local_ok);
    stats->t_rank_sync += GetTime() - t0;
    return global_ok;
}

static double MdReduceMax(double local) {
    double result = 0.0;
    MPI_Reduce(&local, &result, 1, MPI_DOUBLE, MPI_MAX, 0,
               MPI_COMM_WORLD);
    return result;
}

static int MdParseMemory(const std::string &text, size_t *value) {
    *value = 0;
    if (text.empty()) return 0;
    errno = 0;
    char *end = nullptr;
    unsigned long long base = strtoull(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str()) return -1;
    while (*end == ' ' || *end == '\t') ++end;
    unsigned long long multiplier = 1;
    if (*end) {
        char suffix = *end++;
        if (suffix == 'K' || suffix == 'k') multiplier = 1024ull;
        else if (suffix == 'M' || suffix == 'm') multiplier = 1024ull * 1024ull;
        else if (suffix == 'G' || suffix == 'g') multiplier = 1024ull * 1024ull * 1024ull;
        else if (suffix == 'T' || suffix == 't') multiplier = 1024ull * 1024ull * 1024ull * 1024ull;
        else return -1;
        if (*end == 'B' || *end == 'b') ++end;
        while (*end == ' ' || *end == '\t') ++end;
        if (*end) return -1;
    }
    if (base > (unsigned long long)SIZE_MAX / multiplier) return -1;
    *value = (size_t)(base * multiplier);
    return 0;
}

static size_t MdVectorBytes(
        const std::vector<MpiMarkdupCandidateShared> &candidates,
        const std::vector<unsigned char> &qnames,
        const std::vector<uint8_t> &bitmap) {
    size_t total = candidates.capacity() * sizeof(MpiMarkdupCandidateShared);
    if (qnames.capacity() > SIZE_MAX - total) return SIZE_MAX;
    total += qnames.capacity();
    if (bitmap.capacity() > SIZE_MAX - total) return SIZE_MAX;
    return total + bitmap.capacity();
}

static size_t MdFlatHashCapacityForEstimate(size_t expected) {
    if (expected == 0) return 0;
    if (expected > (SIZE_MAX / 2) - 16) return 0;
    size_t wanted = expected * 2 + 16;
    size_t cap = 16;
    while (cap < wanted) {
        if (cap > SIZE_MAX / 2) return 0;
        cap <<= 1;
    }
    return cap;
}

static size_t MdFlatHashBytesForEstimate(size_t owner_count) {
    const size_t cap = MdFlatHashCapacityForEstimate(owner_count);
    if (owner_count != 0 && cap == 0) return SIZE_MAX;
    const size_t best_slot =
        sizeof(MpiMarkdupKeyShared) + sizeof(uint64_t) +
        sizeof(size_t) + sizeof(unsigned char);
    if (cap > SIZE_MAX / best_slot) return SIZE_MAX;
    size_t total = cap * best_slot;
    if (cap > SIZE_MAX / best_slot ||
        total > SIZE_MAX - cap * best_slot) {
        return SIZE_MAX;
    }
    total += cap * best_slot;
    return total;
}

static int MdCheckMemory(size_t bytes, size_t limit,
                         MpiMarkdupStats *stats, int rank,
                         const char *stage) {
    if ((long long)bytes > stats->tracked_peak_bytes) {
        stats->tracked_peak_bytes =
            bytes > (size_t)LLONG_MAX ? LLONG_MAX : (long long)bytes;
    }
    if (limit != 0 && bytes > limit) {
        fprintf(stderr,
                "[rank %d] ERROR: markdup memory limit exceeded at %s: "
                "required=%zu limit=%zu. Increase -m.\n",
                rank, stage, bytes, limit);
        return -1;
    }
    return 0;
}

static int MdAllocateBlockSet(MdBlockSet *set, int n) {
    memset(set, 0, sizeof(*set));
    set->n = n;
    set->blocks = (bam_block *)aligned_alloc_custom(
        64, (size_t)n * sizeof(bam_block));
    set->data = aligned_alloc_custom(
        64, (size_t)n * BGZF_MAX_BLOCK_SIZE);
    if (!set->blocks || !set->data) return -1;
    memset(set->blocks, 0, (size_t)n * sizeof(bam_block));
    for (int i = 0; i < n; ++i) {
        set->blocks[i].data =
            set->data + (size_t)i * BGZF_MAX_BLOCK_SIZE;
        set->blocks[i].block_id = i;
    }
    return 0;
}

static void MdFreeBlockSet(MdBlockSet *set) {
    if (set->blocks) aligned_free_custom((unsigned char *)set->blocks);
    if (set->data) aligned_free_custom(set->data);
    memset(set, 0, sizeof(*set));
}

static int MdAllocateRecordSet(MdRecordSet *set, int n_blocks,
                               int records_per_block,
                               size_t arena_stride) {
    memset(set, 0, sizeof(*set));
    const size_t total = (size_t)n_blocks * records_per_block;
    set->records_per_block = records_per_block;
    set->n_blocks = n_blocks;
    set->arena_stride = arena_stride;
    set->records = (bam1_t *)aligned_alloc_custom(
        64, total * sizeof(bam1_t));
    set->data = aligned_alloc_custom(
        64, (size_t)n_blocks * arena_stride);
    set->ptrs = (bam1_t **)aligned_alloc_custom(
        64, total * sizeof(bam1_t *));
    set->bam_lens = (uint32_t *)aligned_alloc_custom(
        64, total * sizeof(uint32_t));
    if (!set->records || !set->data || !set->ptrs || !set->bam_lens) {
        return -1;
    }
    memset(set->records, 0, total * sizeof(bam1_t));
    memset(set->bam_lens, 0, total * sizeof(uint32_t));
    for (size_t i = 0; i < total; ++i) {
        set->records[i].mempolicy = BAM_USER_OWNS_DATA;
        set->ptrs[i] = set->records + i;
    }
    return 0;
}

static void MdFreeRecordSet(MdRecordSet *set) {
    if (set->records) aligned_free_custom((unsigned char *)set->records);
    if (set->data) aligned_free_custom(set->data);
    if (set->ptrs) aligned_free_custom((unsigned char *)set->ptrs);
    if (set->bam_lens) aligned_free_custom(
        (unsigned char *)set->bam_lens);
    memset(set, 0, sizeof(*set));
}

static int MdAllocateCandidateWorkspace(MdCandidateWorkspace *workspace,
                                        int records_per_block) {
    memset(workspace, 0, sizeof(*workspace));
    workspace->candidates_per_block = records_per_block * 2;
    workspace->qname_stride = BGZF_MAX_BLOCK_SIZE;
    workspace->candidates =
        (MpiMarkdupCandidateShared *)aligned_alloc_custom(
            64, (size_t)kMarkdupNB *
                    workspace->candidates_per_block *
                    sizeof(MpiMarkdupCandidateShared));
    workspace->qnames = aligned_alloc_custom(
        64, (size_t)kMarkdupNB * workspace->qname_stride);
    return workspace->candidates && workspace->qnames ? 0 : -1;
}

static void MdFreeCandidateWorkspace(MdCandidateWorkspace *workspace) {
    if (workspace->candidates) {
        aligned_free_custom((unsigned char *)workspace->candidates);
    }
    if (workspace->qnames) aligned_free_custom(workspace->qnames);
    memset(workspace, 0, sizeof(*workspace));
}

static int MdMemReadBlock(MemReader &reader, bam_block *block) {
    if (reader.pos >= reader.size ||
        reader.size - reader.pos < BLOCK_HEADER_LENGTH) {
        return -1;
    }
    const unsigned char *base =
        (const unsigned char *)reader.base + reader.pos;
    int length = (int)base[16] | ((int)base[17] << 8);
    length++;
    if (length <= 0 || (size_t)length > reader.size - reader.pos ||
        length > BGZF_MAX_BLOCK_SIZE) {
        return -1;
    }
    memcpy(block->data, base, (size_t)length);
    block->length = (unsigned int)length;
    block->pos = 0;
    block->errcode = 0;
    block->block_address = (int64_t)reader.pos;
    reader.pos += (size_t)length;
    return 0;
}

static void MdInitDecompPara(Bam2BamPara *para, int block_id,
                             MdBlockSet *input, MdBlockSet *uncompressed,
                             MdRecordSet *records, int active) {
    const int records_per_block = records->records_per_block;
    memset(para, 0, sizeof(*para));
    para->block_id = block_id;
    para->output_records =
        records->ptrs + (size_t)block_id * records_per_block;
    para->record_base =
        records->records + (size_t)block_id * records_per_block;
    para->bam_lens =
        records->bam_lens + (size_t)block_id * records_per_block;
    para->record_capacity = records_per_block;
    para->data_arena =
        records->data + (size_t)block_id * records->arena_stride;
    para->data_arena_capacity = records->arena_stride;
    para->filter.min_mapq = -1;
    para->filter.max_mapq = -1;
    para->filter.require_flag = 0;
    para->filter.exclude_flag = 0;
    para->filter.ref_tid = -1;
    para->filter.min_read_len = -1;
    para->filter.max_read_len = -1;
    if (active) {
        para->input_block = input->blocks + block_id;
        para->un_comp_block = uncompressed->blocks + block_id;
        para->status = 0;
    } else {
        para->input_block = nullptr;
        para->un_comp_block = nullptr;
        para->status = -1;
    }
}

static int MdReadGroup(MemReader &reader, MdBlockSet *input,
                       int *n_blocks, MpiMarkdupStats *stats) {
    double t0 = GetTime();
    int count = 0;
    for (int i = 0; i < kMarkdupNB; ++i) {
        if (MdMemReadBlock(reader, input->blocks + i) != 0) break;
        input->blocks[i].block_id = i;
        count++;
    }
    *n_blocks = count;
    stats->t_read += GetTime() - t0;
    return 0;
}

static const char *MdExtractErrorText(int status) {
    switch (status) {
    case -10:
        return "missing MC:Z tag; run fixmate -m first";
    case -11:
        return "MC tag has wrong type; run fixmate -m first";
    case -12:
        return "missing ms integer tag; run fixmate -m first";
    case -13:
        return "ms tag has wrong type; run fixmate -m first";
    case -20:
        return "input is not coordinate sorted";
    default:
        return "candidate extraction failed";
    }
}

static int MdExtractCandidates(
        MemReader &reader, long long global_block_begin,
        int rank, int include_fails, size_t memory_limit,
        std::vector<MpiMarkdupCandidateShared> *local_candidates,
        std::vector<unsigned char> *local_qnames,
        uint64_t *local_records,
        int *range_has_records,
        int *range_first_tid, int *range_first_pos,
        int *range_last_tid, int *range_last_pos,
        MpiMarkdupStats *stats) {

    // 工作区分配
    const int records_per_block = (int)MPI_RECORDS_PER_BLOCK;
    struct MdExtractSlot {
        MdBlockSet input;
        MdBlockSet uncompressed;
        MdRecordSet records;
        MdCandidateWorkspace candidate_workspace;
        Bam2BamPara decomp[kMarkdupNB];
        MpiMarkdupExtractPara extract[kMarkdupNB];
        int n_blocks;
        long long block_group_base;
        uint64_t ordinal_base;
        uint64_t ordinal_end;
    };

    MdExtractSlot slots[2];
    memset(slots, 0, sizeof(slots));

    auto free_slot = [](MdExtractSlot *slot) {
        MdFreeBlockSet(&slot->input);
        MdFreeBlockSet(&slot->uncompressed);
        MdFreeRecordSet(&slot->records);
        MdFreeCandidateWorkspace(&slot->candidate_workspace);
        slot->n_blocks = 0;
    };
    auto free_slots = [&]() {
        free_slot(slots + 0);
        free_slot(slots + 1);
    };

    for (int s = 0; s < 2; ++s) {
        if (MdAllocateBlockSet(&slots[s].input, kMarkdupNB) != 0 ||
            MdAllocateBlockSet(&slots[s].uncompressed, kMarkdupNB) != 0 ||
            MdAllocateRecordSet(&slots[s].records, kMarkdupNB,
                                records_per_block,
                                MPI_BAM_BLOCK_ARENA_SIZE) != 0 ||
            MdAllocateCandidateWorkspace(
                &slots[s].candidate_workspace,
                records_per_block) != 0) {
            fprintf(stderr,
                    "[rank %d] ERROR: failed to allocate markdup "
                    "candidate pipeline workspace.\n",
                    rank);
            free_slots();
            return -1;
        }
    }

    // 固定内存估算：candidate extract 使用双槽流水线。
    const MdCandidateWorkspace &cw0 = slots[0].candidate_workspace;
    const size_t fixed_workspace_per_slot =
        (size_t)kMarkdupNB * BGZF_MAX_BLOCK_SIZE * 2 +
        (size_t)kMarkdupNB * MPI_BAM_BLOCK_ARENA_SIZE +
        (size_t)kMarkdupNB * records_per_block *
            (sizeof(bam1_t) + sizeof(bam1_t *) + sizeof(uint32_t)) +
        (size_t)kMarkdupNB * cw0.candidates_per_block *
            sizeof(MpiMarkdupCandidateShared) +
        (size_t)kMarkdupNB * cw0.qname_stride;
    size_t fixed_workspace = fixed_workspace_per_slot;
    if (fixed_workspace > SIZE_MAX - fixed_workspace_per_slot) {
        free_slots();
        return -1;
    }
    fixed_workspace += fixed_workspace_per_slot;
    if (MdCheckMemory(fixed_workspace, memory_limit, stats, rank,
                      "candidate pipeline workspace") != 0) {
        free_slots();
        return -1;
    }

    // ordinal 是本 rank 内 record 编号，作为索引
    uint64_t ordinal = 0;
    long long block_group_base = global_block_begin;
    int previous_tid = -1;
    int previous_pos = -1;
    int have_previous = 0;
    *range_has_records = 0;
    *range_first_tid = -1;
    *range_first_pos = -1;
    *range_last_tid = -1;
    *range_last_pos = -1;

    auto read_slot = [&](MdExtractSlot *slot, long long block_base,
                         uint64_t ordinal_base) -> int {
        slot->block_group_base = block_base;
        slot->ordinal_base = ordinal_base;
        slot->ordinal_end = ordinal_base;
        if (MdReadGroup(reader, &slot->input, &slot->n_blocks,
                        stats) != 0) {
            return -1;
        }
        if (slot->n_blocks > 0) {
            stats->input_blocks += slot->n_blocks;
        }
        return 0;
    };

    auto prepare_decomp = [&](MdExtractSlot *slot) {
        for (int b = 0; b < kMarkdupNB; ++b) {
            MdInitDecompPara(slot->decomp + b, b, &slot->input,
                             &slot->uncompressed, &slot->records,
                             b < slot->n_blocks);
        }
    };

    auto launch_decomp = [&](MdExtractSlot *slot) {
        prepare_decomp(slot);
        __real_athread_spawn(
            (void *)slave_mpi_decompress_bam2bam_passthrough,
            slot->decomp, 1);
    };

    auto join_decomp = [&](MdExtractSlot *slot) -> int {
        double cpe_sync_t0 = GetTime();
        athread_join();
        const double sync_dt = GetTime() - cpe_sync_t0;
        stats->t_cpe_sync += sync_dt;
        stats->t_candidate_decomp += sync_dt;
        for (int b = 0; b < slot->n_blocks; ++b) {
            if (slot->decomp[b].status != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: markdup candidate BAM decode failed "
                        "at global block %lld status=%d record=%d.\n",
                        rank, slot->block_group_base + b,
                        slot->decomp[b].status,
                        slot->decomp[b].record_index);
                return -1;
            }
        }
        return 0;
    };

    auto prepare_extract = [&](MdExtractSlot *slot) {
        uint64_t block_ordinal = slot->ordinal_base;
        for (int b = 0; b < kMarkdupNB; ++b) {
            memset(slot->extract + b, 0, sizeof(slot->extract[b]));
            slot->extract[b].block_id = b;
            slot->extract[b].status = b < slot->n_blocks ? 0 : -1;
            if (b >= slot->n_blocks) continue;
            slot->extract[b].records = slot->decomp[b].output_records;
            slot->extract[b].n_records = slot->decomp[b].n_total_records;
            slot->extract[b].ordinal_base = block_ordinal;
            slot->extract[b].global_block_index =
                (uint64_t)(slot->block_group_base + b);
            slot->extract[b].source_rank = rank;
            slot->extract[b].include_fails = include_fails;
            slot->extract[b].candidates =
                slot->candidate_workspace.candidates +
                (size_t)b *
                    slot->candidate_workspace.candidates_per_block;
            slot->extract[b].candidate_capacity =
                slot->candidate_workspace.candidates_per_block;
            slot->extract[b].qname_arena =
                slot->candidate_workspace.qnames +
                (size_t)b * slot->candidate_workspace.qname_stride;
            slot->extract[b].qname_capacity =
                slot->candidate_workspace.qname_stride;
            block_ordinal +=
                (uint64_t)slot->decomp[b].n_total_records;
        }
        slot->ordinal_end = block_ordinal;
    };

    auto launch_extract = [&](MdExtractSlot *slot) {
        __real_athread_spawn((void *)slave_mpi_markdup_extract,
                             slot->extract, 1);
    };

    auto join_extract = [&](MdExtractSlot *slot) {
        (void)slot;
        double cpe_sync_t0 = GetTime();
        athread_join();
        const double sync_dt = GetTime() - cpe_sync_t0;
        stats->t_cpe_sync += sync_dt;
        stats->t_candidate_extract += sync_dt;
    };

    auto run_slot_cpe = [&](MdExtractSlot *slot) -> int {
        launch_decomp(slot);
        if (join_decomp(slot) != 0) return -1;
        prepare_extract(slot);
        launch_extract(slot);
        join_extract(slot);
        return 0;
    };

    auto process_slot = [&](MdExtractSlot *slot) -> int {
        double status_total = 0.0;
        double merge_total = 0.0;
        for (int b = 0; b < slot->n_blocks; ++b) {
            double status_t0 = GetTime();
            // 检查从核提取是否成功
            if (slot->extract[b].status != 0) {
                stats->t_extract_status +=
                    status_total + (GetTime() - status_t0);
                fprintf(stderr,
                        "[rank %d] ERROR: %s at global block %lld "
                        "record=%d status=%d.\n",
                        rank,
                        MdExtractErrorText(slot->extract[b].status),
                        slot->block_group_base + b,
                        slot->extract[b].record_index,
                        slot->extract[b].status);
                return -1;
            }

            // 检查 coordinate sorted
            // MdExtractCandidates 内部检查 rank 内 block 顺序
            // MdValidateRankBoundaries 检查 rank 与 rank 之间的坐标顺序
            if (slot->extract[b].has_records) {
                if (!*range_has_records) {
                    *range_first_tid = slot->extract[b].first_tid;
                    *range_first_pos = slot->extract[b].first_pos;
                    *range_has_records = 1;
                }
                if (slot->extract[b].first_tid >= 0 && have_previous &&
                    (slot->extract[b].first_tid < previous_tid ||
                     (slot->extract[b].first_tid == previous_tid &&
                      slot->extract[b].first_pos < previous_pos))) {
                    fprintf(stderr,
                            "[rank %d] ERROR: input is not coordinate sorted "
                            "near global block %lld.\n",
                            rank, slot->block_group_base + b);
                    stats->t_extract_status +=
                        status_total + (GetTime() - status_t0);
                    return -1;
                }
                previous_tid = slot->extract[b].last_tid;
                previous_pos = slot->extract[b].last_pos;
                *range_last_tid = slot->extract[b].last_tid;
                *range_last_pos = slot->extract[b].last_pos;
                have_previous = 1;
            }
            status_total += GetTime() - status_t0;

            // 合并 qname 和 candidates 到 local_candidates 中
            double merge_t0 = GetTime();
            const uint64_t qname_base =
                (uint64_t)local_qnames->size();
            if (slot->extract[b].qname_used > 0) {
                local_qnames->insert(
                    local_qnames->end(), slot->extract[b].qname_arena,
                    slot->extract[b].qname_arena +
                        slot->extract[b].qname_used);
            }
            for (int i = 0; i < slot->extract[b].n_candidates; ++i) {
                MpiMarkdupCandidateShared candidate =
                    slot->extract[b].candidates[i];
                if (!candidate.key.single) {
                    candidate.qname_offset += qname_base;
                }
                local_candidates->push_back(candidate);
            }
            stats->total_records += slot->decomp[b].n_total_records;
            stats->examined_records += slot->extract[b].examined;
            stats->excluded_records += slot->extract[b].excluded;
            stats->pair_candidates +=
                slot->extract[b].pair_candidates;
            stats->single_candidates +=
                slot->extract[b].single_candidates;
            merge_total += GetTime() - merge_t0;
        }
        stats->t_extract_status += status_total;
        stats->t_extract_merge += merge_total;
        ordinal = slot->ordinal_end;
        stats->group_count++;

        // 动态内存检查
        double memcheck_t0 = GetTime();
        size_t dynamic = MdVectorBytes(
            *local_candidates, *local_qnames,
            std::vector<uint8_t>());
        if (dynamic == SIZE_MAX ||
            fixed_workspace > SIZE_MAX - dynamic ||
            MdCheckMemory(fixed_workspace + dynamic, memory_limit,
                          stats, rank,
                          "local candidate accumulation") != 0) {
            stats->t_extract_memcheck += GetTime() - memcheck_t0;
            return -1;
        }
        stats->t_extract_memcheck += GetTime() - memcheck_t0;
        return 0;
    };

    int current = 0;
    int next = 1;
    if (read_slot(slots + current, global_block_begin, ordinal) != 0) {
        free_slots();
        return -1;
    }
    if (slots[current].n_blocks == 0) {
        *local_records = 0;
        free_slots();
        return 0;
    }
    if (run_slot_cpe(slots + current) != 0) {
        free_slots();
        return -1;
    }

    for (;;) {
        const long long next_block_base =
            slots[current].block_group_base + slots[current].n_blocks;
        const uint64_t next_ordinal_base = slots[current].ordinal_end;
        if (read_slot(slots + next, next_block_base,
                      next_ordinal_base) != 0) {
            free_slots();
            return -1;
        }
        const int have_next = slots[next].n_blocks > 0;
        if (have_next) {
            launch_decomp(slots + next);
        }

        int process_status = process_slot(slots + current);
        int decomp_status = 0;
        if (have_next) {
            decomp_status = join_decomp(slots + next);
        }
        if (process_status != 0 || decomp_status != 0) {
            free_slots();
            return -1;
        }
        if (!have_next) break;

        prepare_extract(slots + next);
        launch_extract(slots + next);
        join_extract(slots + next);
        std::swap(current, next);
    }

    *local_records = ordinal;
    free_slots();
    return 0;
}

static bool MdKeyLess(const MpiMarkdupKeyShared &a,
                      const MpiMarkdupKeyShared &b) {
    if (a.single != b.single) return a.single < b.single;
    if (a.this_ref != b.this_ref) return a.this_ref < b.this_ref;
    if (a.this_coord != b.this_coord) {
        return a.this_coord < b.this_coord;
    }
    if (a.orientation != b.orientation) {
        return a.orientation < b.orientation;
    }
    if (!a.single) {
        if (a.other_ref != b.other_ref) {
            return a.other_ref < b.other_ref;
        }
        if (a.other_coord != b.other_coord) {
            return a.other_coord < b.other_coord;
        }
        if (a.leftmost != b.leftmost) {
            return a.leftmost < b.leftmost;
        }
    }
    return false;
}

static bool MdKeyEqual(const MpiMarkdupKeyShared &a,
                       const MpiMarkdupKeyShared &b) {
    return !MdKeyLess(a, b) && !MdKeyLess(b, a);
}

static uint64_t MdHashKey(const MpiMarkdupKeyShared &key) {
    const unsigned char *p =
        (const unsigned char *)&key;
    size_t n = sizeof(key);
    uint64_t hash = 1469598103934665603ull;
    for (size_t i = 0; i < n; ++i) {
        hash ^= p[i];
        hash *= 1099511628211ull;
    }
    return hash;
}

static int MdSendrecvBytes(const void *send_data,
                           unsigned long long send_bytes, int dest,
                           void *recv_data,
                           unsigned long long recv_bytes, int src,
                           int tag) {
    unsigned long long sent = 0;
    unsigned long long received = 0;
    while (sent < send_bytes || received < recv_bytes) {
        int send_count = (int)std::min(
            (unsigned long long)kMarkdupExchangeChunk,
            send_bytes - sent);
        int recv_count = (int)std::min(
            (unsigned long long)kMarkdupExchangeChunk,
            recv_bytes - received);
        const unsigned char *send_ptr =
            send_count ? (const unsigned char *)send_data + sent
                       : nullptr;
        unsigned char *recv_ptr =
            recv_count ? (unsigned char *)recv_data + received
                       : nullptr;
        if (MPI_Sendrecv((void *)send_ptr, send_count, MPI_BYTE,
                         dest, tag, recv_ptr, recv_count, MPI_BYTE,
                         src, tag, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE) != MPI_SUCCESS) {
            return -1;
        }
        sent += (unsigned long long)send_count;
        received += (unsigned long long)recv_count;
    }
    return 0;
}

/*
把本 rank 提取出来的 local_candidates
按照 duplicate key 的 hash 分配给 owner rank
通过 MPI 交换
让每个 owner rank 拿到自己负责判断 duplicate 的候选集合
*/
static int MdExchangeCandidates(
        const std::vector<MpiMarkdupCandidateShared> &local_candidates,
        const std::vector<unsigned char> &local_qnames,
        int rank, int comm_size, size_t memory_limit,
        std::vector<MpiMarkdupCandidateShared> *owner_candidates,
        std::vector<unsigned char> *owner_qnames,
        MpiMarkdupStats *stats) {
    double exchange_total_t0 = GetTime();
    std::vector<unsigned long long> send_candidate_counts(
        (size_t)comm_size, 0);
    std::vector<unsigned long long> recv_candidate_counts(
        (size_t)comm_size, 0);
    std::vector<unsigned long long> send_qname_counts(
        (size_t)comm_size, 0);
    std::vector<unsigned long long> recv_qname_counts(
        (size_t)comm_size, 0);
    std::vector<size_t> send_candidate_offsets((size_t)comm_size + 1, 0);
    std::vector<size_t> send_qname_offsets((size_t)comm_size + 1, 0);
    std::unique_ptr<MpiMarkdupCandidateShared,
                    MdMarkdupCandidateBufferDeleter> send_candidates_flat;
    std::unique_ptr<unsigned char,
                    MdByteBufferDeleter> send_qnames_flat;
    int local_ok = 1;

    // 第一阶段：两遍式按 owner 连续打包 send_candidates / send_qnames
    double pack_t0 = GetTime();
    for (size_t i = 0; i < local_candidates.size(); ++i) {
        const MpiMarkdupCandidateShared &candidate =
            local_candidates[i];
        const int owner = (int)(MdHashKey(candidate.key) %
                          (uint64_t)comm_size);
        if (send_candidate_counts[(size_t)owner] ==
            ULLONG_MAX) {
            local_ok = 0;
            break;
        }
        send_candidate_counts[(size_t)owner]++;
        if (!candidate.key.single) {
            if (candidate.qname_offset > local_qnames.size() ||
                candidate.qname_len >
                    local_qnames.size() -
                    (size_t)candidate.qname_offset) {
                fprintf(stderr,
                        "[rank %d] ERROR: invalid local markdup QNAME range.\n",
                        rank);
                local_ok = 0;
                break;
            }
            const unsigned long long qname_len =
                (unsigned long long)candidate.qname_len;
            if (send_qname_counts[(size_t)owner] >
                ULLONG_MAX - qname_len) {
                local_ok = 0;
                break;
            }
            send_qname_counts[(size_t)owner] += qname_len;
        }
    }
    if (!MdAllRanksOk(local_ok)) return -1;

    size_t total_send_candidates = 0;
    size_t total_send_qnames = 0;
    for (int i = 0; i < comm_size; ++i) {
        send_candidate_offsets[(size_t)i] = total_send_candidates;
        send_qname_offsets[(size_t)i] = total_send_qnames;
        if (send_candidate_counts[(size_t)i] >
                (unsigned long long)(SIZE_MAX - total_send_candidates) ||
            send_qname_counts[(size_t)i] >
                (unsigned long long)(SIZE_MAX - total_send_qnames)) {
            local_ok = 0;
            break;
        }
        total_send_candidates +=
            (size_t)send_candidate_counts[(size_t)i];
        total_send_qnames +=
            (size_t)send_qname_counts[(size_t)i];
    }
    send_candidate_offsets[(size_t)comm_size] = total_send_candidates;
    send_qname_offsets[(size_t)comm_size] = total_send_qnames;
    if (!MdAllRanksOk(local_ok)) return -1;

    if (total_send_candidates > 0) {
        if (total_send_candidates >
            SIZE_MAX / sizeof(MpiMarkdupCandidateShared)) {
            local_ok = 0;
        } else {
            send_candidates_flat.reset(
                (MpiMarkdupCandidateShared *)aligned_alloc_custom(
                    64, total_send_candidates *
                            sizeof(MpiMarkdupCandidateShared)));
            if (!send_candidates_flat) local_ok = 0;
        }
    }
    if (total_send_qnames > 0) {
        send_qnames_flat.reset(
            aligned_alloc_custom(64, total_send_qnames));
        if (!send_qnames_flat) local_ok = 0;
    }
    if (!MdAllRanksOk(local_ok)) return -1;
    std::vector<size_t> candidate_write = send_candidate_offsets;
    std::vector<size_t> qname_write = send_qname_offsets;

    for (size_t i = 0; i < local_candidates.size(); ++i) {
        MpiMarkdupCandidateShared candidate = local_candidates[i];
        const int owner = (int)(MdHashKey(candidate.key) %
                          (uint64_t)comm_size);
        if (!candidate.key.single) {
            const size_t qname_src = (size_t)candidate.qname_offset;
            const size_t qname_len = (size_t)candidate.qname_len;
            const size_t qname_dst = qname_write[(size_t)owner];
            candidate.qname_offset =
                (uint64_t)(qname_dst -
                           send_qname_offsets[(size_t)owner]);
            if (qname_len > 0) {
                memcpy(send_qnames_flat.get() + qname_dst,
                       local_qnames.data() + qname_src,
                       qname_len);
            }
            qname_write[(size_t)owner] += qname_len;
        }
        send_candidates_flat.get()[candidate_write[(size_t)owner]++] =
            candidate;
    }
    stats->t_candidate_pack += GetTime() - pack_t0;

    // 第二阶段：统计将要收多少数据，并做内存检查
    double setup_t0 = GetTime();
    size_t send_candidate_bytes_total = 0;
    size_t exchange_memory = 0;
    if (total_send_candidates >
        SIZE_MAX / sizeof(MpiMarkdupCandidateShared)) {
        local_ok = 0;
    } else {
        send_candidate_bytes_total =
            total_send_candidates * sizeof(MpiMarkdupCandidateShared);
        exchange_memory = send_candidate_bytes_total;
    }
    if (total_send_qnames > SIZE_MAX - exchange_memory) {
        local_ok = 0;
    } else {
        exchange_memory += total_send_qnames;
    }
    size_t local_memory =
        local_candidates.capacity() *
            sizeof(MpiMarkdupCandidateShared) +
        local_qnames.capacity();

    MPI_Alltoall(send_candidate_counts.data(), 1,
                 MPI_UNSIGNED_LONG_LONG,
                 recv_candidate_counts.data(), 1,
                 MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
    MPI_Alltoall(send_qname_counts.data(), 1,
                 MPI_UNSIGNED_LONG_LONG,
                 recv_qname_counts.data(), 1,
                 MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

    unsigned long long owner_candidate_count = 0;
    unsigned long long owner_qname_count = 0;
    for (int i = 0; i < comm_size; ++i) {
        if (owner_candidate_count >
                ULLONG_MAX - recv_candidate_counts[(size_t)i] ||
            owner_qname_count >
                ULLONG_MAX - recv_qname_counts[(size_t)i]) {
            local_ok = 0;
            break;
        }
        owner_candidate_count +=
            recv_candidate_counts[(size_t)i];
        owner_qname_count += recv_qname_counts[(size_t)i];
    }
    if (owner_candidate_count >
            (unsigned long long)(SIZE_MAX /
                sizeof(MpiMarkdupCandidateShared)) ||
        owner_qname_count > (unsigned long long)SIZE_MAX) {
        local_ok = 0;
    } else {
        size_t owner_candidate_bytes =
            (size_t)owner_candidate_count *
                sizeof(MpiMarkdupCandidateShared);
        size_t owner_bytes = 0;
        if ((size_t)owner_qname_count >
            SIZE_MAX - owner_candidate_bytes) {
            local_ok = 0;
        } else {
            owner_bytes =
                owner_candidate_bytes + (size_t)owner_qname_count;
        }
        if (local_ok &&
            (local_memory > SIZE_MAX - exchange_memory ||
             local_memory + exchange_memory >
                 SIZE_MAX - owner_bytes ||
             MdCheckMemory(local_memory + exchange_memory +
                               owner_bytes,
                           memory_limit, stats, rank,
                           "candidate owner exchange") != 0)) {
            local_ok = 0;
        }
        if (local_ok) {
            owner_candidates->reserve((size_t)owner_candidate_count);
            owner_qnames->reserve((size_t)owner_qname_count);
        }
    }
    if (!MdAllRanksOk(local_ok)) {
        return -1;
    }
    stats->t_candidate_exchange_setup += GetTime() - setup_t0;

    // 第三阶段：真正 MPI 交换
    double mpi_t0 = GetTime();
    for (int step = 0; step < comm_size; ++step) {
        int dest = (rank + step) % comm_size;
        int src = (rank - step + comm_size) % comm_size;
        if (step == 0) {
            const uint64_t qbase = (uint64_t)owner_qnames->size();
            size_t candidate_base = owner_candidates->size();
            const size_t qname_begin = send_qname_offsets[(size_t)rank];
            const size_t qname_count =
                (size_t)send_qname_counts[(size_t)rank];
            const size_t candidate_begin =
                send_candidate_offsets[(size_t)rank];
            const size_t candidate_count =
                (size_t)send_candidate_counts[(size_t)rank];
            if (qname_count > 0) {
                owner_qnames->insert(
                    owner_qnames->end(),
                    send_qnames_flat.get() + qname_begin,
                    send_qnames_flat.get() + qname_begin + qname_count);
            }
            if (candidate_count > 0) {
                owner_candidates->insert(
                    owner_candidates->end(),
                    send_candidates_flat.get() + candidate_begin,
                    send_candidates_flat.get() + candidate_begin +
                        candidate_count);
            }
            for (size_t i = candidate_base;
                 i < owner_candidates->size(); ++i) {
                if (!(*owner_candidates)[i].key.single) {
                    (*owner_candidates)[i].qname_offset += qbase;
                }
            }
            continue;
        }

        unsigned long long send_candidate_count =
            send_candidate_counts[(size_t)dest];
        unsigned long long send_qname_count =
            send_qname_counts[(size_t)dest];
        unsigned long long recv_candidate_count =
            recv_candidate_counts[(size_t)src];
        unsigned long long recv_qname_count =
            recv_qname_counts[(size_t)src];
        if (recv_candidate_count >
            (unsigned long long)(SIZE_MAX /
                sizeof(MpiMarkdupCandidateShared)) ||
            recv_qname_count > (unsigned long long)SIZE_MAX) {
            return -1;
        }
        size_t candidate_base = owner_candidates->size();
        size_t qname_base = owner_qnames->size();
        owner_candidates->resize(
            candidate_base + (size_t)recv_candidate_count);
        owner_qnames->resize(
            qname_base + (size_t)recv_qname_count);

        unsigned long long send_candidate_bytes =
            send_candidate_count *
            sizeof(MpiMarkdupCandidateShared);
        unsigned long long recv_candidate_bytes =
            recv_candidate_count *
            sizeof(MpiMarkdupCandidateShared);
        if (MdSendrecvBytes(
                send_candidate_count
                    ? send_candidates_flat.get() +
                          send_candidate_offsets[(size_t)dest]
                    : nullptr,
                send_candidate_bytes, dest,
                recv_candidate_count
                    ? owner_candidates->data() + candidate_base
                    : nullptr,
                recv_candidate_bytes, src, 4101) != 0 ||
            MdSendrecvBytes(
                send_qname_count
                    ? send_qnames_flat.get() +
                          send_qname_offsets[(size_t)dest]
                    : nullptr,
                send_qname_count, dest,
                recv_qname_count
                    ? owner_qnames->data() + qname_base
                    : nullptr,
                recv_qname_count, src, 4102) != 0) {
            return -1;
        }
        for (size_t i = candidate_base;
             i < owner_candidates->size(); ++i) {
            if (!(*owner_candidates)[i].key.single) {
                (*owner_candidates)[i].qname_offset +=
                    (uint64_t)qname_base;
            }
        }
        stats->mpi_candidate_bytes +=
            (long long)send_candidate_bytes +
            (long long)send_qname_count;

    }

    // 最后更新统计信息
    stats->t_candidate_exchange_mpi += GetTime() - mpi_t0;
    stats->t_candidate_exchange += GetTime() - exchange_total_t0;
    stats->owner_candidates =
        (long long)owner_candidates->size();
    stats->qname_bytes = (long long)owner_qnames->size();
    return 0;
}

static int MdCompareCoord(int32_t a_ref, int64_t a_coord,
                          int32_t b_ref, int64_t b_coord) {
    if (a_ref != b_ref) return a_ref < b_ref ? -1 : 1;
    if (a_coord != b_coord) return a_coord < b_coord ? -1 : 1;
    return 0;
}

static int MdBuildCoordBoundaries(
        int range_has_records, int range_first_tid, int range_first_pos,
        int range_last_tid, int range_last_pos, int comm_size,
        std::vector<MdCoordBoundary> *boundaries) {
    long long local[5] = {
        (long long)range_has_records,
        range_has_records ? (long long)range_first_tid + 1 : 0,
        range_has_records ? (long long)range_first_pos + 1 : 0,
        range_has_records ? (long long)range_last_tid + 1 : 0,
        range_has_records ? (long long)range_last_pos + 1 : 0
    };
    std::vector<long long> gathered((size_t)comm_size * 5, 0);
    MPI_Allgather(local, 5, MPI_LONG_LONG, gathered.data(), 5,
                  MPI_LONG_LONG, MPI_COMM_WORLD);
    boundaries->assign((size_t)comm_size, MdCoordBoundary());
    for (int r = 0; r < comm_size; ++r) {
        const long long *src = gathered.data() + (size_t)r * 5;
        MdCoordBoundary &b = (*boundaries)[(size_t)r];
        b.has_records = src[0] != 0;
        b.first_ref = (int32_t)src[1];
        b.first_coord = (int64_t)src[2];
        b.last_ref = (int32_t)src[3];
        b.last_coord = (int64_t)src[4];
    }
    return 0;
}

static int MdCoordinateOwner(const MpiMarkdupKeyShared &key,
                             const std::vector<MdCoordBoundary> &boundaries,
                             int fallback_rank) {
    int last_nonempty = fallback_rank;
    const int32_t ref = key.this_ref;
    const int64_t coord = key.this_coord;
    for (size_t r = 0; r < boundaries.size(); ++r) {
        const MdCoordBoundary &b = boundaries[r];
        if (!b.has_records) continue;
        last_nonempty = (int)r;
        if (MdCompareCoord(ref, coord, b.last_ref, b.last_coord) <= 0) {
            return (int)r;
        }
    }
    return last_nonempty;
}

static int MdCoordinateOwnerFast(const MpiMarkdupKeyShared &key,
                                 const std::vector<MdCoordBoundary> &boundaries,
                                 int rank) {
    if (rank < 0 || rank >= (int)boundaries.size()) {
        return MdCoordinateOwner(key, boundaries, rank);
    }

    const MdCoordBoundary &local = boundaries[(size_t)rank];
    const int32_t ref = key.this_ref;
    const int64_t coord = key.this_coord;
    if (local.has_records) {
        const int cmp_first =
            MdCompareCoord(ref, coord, local.first_ref,
                           local.first_coord);
        const int cmp_last =
            MdCompareCoord(ref, coord, local.last_ref,
                           local.last_coord);
        if (cmp_first >= 0 && cmp_last <= 0) {
            return rank;
        }

        if (cmp_last > 0) {
            int last_nonempty = rank;
            for (size_t r = (size_t)rank + 1;
                 r < boundaries.size(); ++r) {
                const MdCoordBoundary &b = boundaries[r];
                if (!b.has_records) continue;
                last_nonempty = (int)r;
                if (MdCompareCoord(ref, coord, b.last_ref,
                                   b.last_coord) <= 0) {
                    return (int)r;
                }
            }
            return last_nonempty;
        }
    }

    return MdCoordinateOwner(key, boundaries, rank);
}

static int MdAppendOwnerCandidate(
        const MpiMarkdupCandidateShared &input,
        const std::vector<unsigned char> &local_qnames,
        std::vector<MpiMarkdupCandidateShared> *owner_candidates,
        std::vector<unsigned char> *owner_qnames,
        int rank) {
    MpiMarkdupCandidateShared candidate = input;
    if (!candidate.key.single) {
        const size_t qname_src = (size_t)candidate.qname_offset;
        const size_t qname_len = (size_t)candidate.qname_len;
        if (qname_src > local_qnames.size() ||
            qname_len > local_qnames.size() - qname_src) {
            fprintf(stderr,
                    "[rank %d] ERROR: invalid local markdup QNAME range.\n",
                    rank);
            return -1;
        }
        candidate.qname_offset = (uint64_t)owner_qnames->size();
        if (qname_len > 0) {
            owner_qnames->insert(owner_qnames->end(),
                                 local_qnames.data() + qname_src,
                                 local_qnames.data() + qname_src + qname_len);
        }
    }
    owner_candidates->push_back(candidate);
    return 0;
}

static int MdExchangeCandidatesByCoordinate(
        std::vector<MpiMarkdupCandidateShared> *local_candidates_in,
        std::vector<unsigned char> *local_qnames_in,
        int range_has_records, int range_first_tid, int range_first_pos,
        int range_last_tid, int range_last_pos,
        int rank, int comm_size, size_t memory_limit,
        std::vector<MpiMarkdupCandidateShared> *owner_candidates,
        std::vector<unsigned char> *owner_qnames,
        MpiMarkdupStats *stats) {
    double exchange_total_t0 = GetTime();
    std::vector<MpiMarkdupCandidateShared> &local_candidates =
        *local_candidates_in;
    std::vector<unsigned char> &local_qnames = *local_qnames_in;
    std::vector<MdCoordBoundary> boundaries;
    MdBuildCoordBoundaries(range_has_records, range_first_tid,
                           range_first_pos, range_last_tid,
                           range_last_pos, comm_size, &boundaries);

    std::vector<unsigned long long> send_candidate_counts(
        (size_t)comm_size, 0);
    std::vector<unsigned long long> recv_candidate_counts(
        (size_t)comm_size, 0);
    std::vector<unsigned long long> send_qname_counts(
        (size_t)comm_size, 0);
    std::vector<unsigned long long> recv_qname_counts(
        (size_t)comm_size, 0);
    std::vector<size_t> send_candidate_offsets((size_t)comm_size + 1, 0);
    std::vector<size_t> send_qname_offsets((size_t)comm_size + 1, 0);
    std::unique_ptr<MpiMarkdupCandidateShared,
                    MdMarkdupCandidateBufferDeleter> send_candidates_flat;
    std::unique_ptr<unsigned char, MdByteBufferDeleter> send_qnames_flat;
    std::vector<MdRemoteCandidateRef> remote_candidates;
    size_t self_candidate_count = 0;
    size_t self_qname_count = 0;
    int local_ok = 1;

    double pack_t0 = GetTime();
    for (size_t i = 0; i < local_candidates.size(); ++i) {
        const MpiMarkdupCandidateShared &candidate = local_candidates[i];
        const int owner =
            MdCoordinateOwnerFast(candidate.key, boundaries, rank);
        if (owner < 0 || owner >= comm_size) {
            local_ok = 0;
            break;
        }
        if (owner == rank) {
            self_candidate_count++;
            if (!candidate.key.single) {
                if (candidate.qname_offset > local_qnames.size() ||
                    candidate.qname_len >
                        local_qnames.size() -
                            (size_t)candidate.qname_offset) {
                    fprintf(stderr,
                            "[rank %d] ERROR: invalid local markdup QNAME range.\n",
                            rank);
                    local_ok = 0;
                    break;
                }
                self_qname_count += (size_t)candidate.qname_len;
            }
            continue;
        }
        MdRemoteCandidateRef ref;
        ref.index = i;
        ref.owner = owner;
        remote_candidates.push_back(ref);
        if (send_candidate_counts[(size_t)owner] == ULLONG_MAX) {
            local_ok = 0;
            break;
        }
        send_candidate_counts[(size_t)owner]++;
        if (!candidate.key.single) {
            if (candidate.qname_offset > local_qnames.size() ||
                candidate.qname_len >
                    local_qnames.size() -
                        (size_t)candidate.qname_offset) {
                fprintf(stderr,
                        "[rank %d] ERROR: invalid local markdup QNAME range.\n",
                        rank);
                local_ok = 0;
                break;
            }
            const unsigned long long qname_len =
                (unsigned long long)candidate.qname_len;
            if (send_qname_counts[(size_t)owner] >
                ULLONG_MAX - qname_len) {
                local_ok = 0;
                break;
            }
            send_qname_counts[(size_t)owner] += qname_len;
        }
    }
    if (!MdAllRanksOk(local_ok)) return -1;

    size_t total_send_candidates = 0;
    size_t total_send_qnames = 0;
    for (int i = 0; i < comm_size; ++i) {
        send_candidate_offsets[(size_t)i] = total_send_candidates;
        send_qname_offsets[(size_t)i] = total_send_qnames;
        if (send_candidate_counts[(size_t)i] >
                (unsigned long long)(SIZE_MAX - total_send_candidates) ||
            send_qname_counts[(size_t)i] >
                (unsigned long long)(SIZE_MAX - total_send_qnames)) {
            local_ok = 0;
            break;
        }
        total_send_candidates +=
            (size_t)send_candidate_counts[(size_t)i];
        total_send_qnames += (size_t)send_qname_counts[(size_t)i];
    }
    send_candidate_offsets[(size_t)comm_size] = total_send_candidates;
    send_qname_offsets[(size_t)comm_size] = total_send_qnames;
    if (!MdAllRanksOk(local_ok)) return -1;

    if (total_send_candidates > 0) {
        if (total_send_candidates >
            SIZE_MAX / sizeof(MpiMarkdupCandidateShared)) {
            local_ok = 0;
        } else {
            send_candidates_flat.reset(
                (MpiMarkdupCandidateShared *)aligned_alloc_custom(
                    64, total_send_candidates *
                            sizeof(MpiMarkdupCandidateShared)));
            if (!send_candidates_flat) local_ok = 0;
        }
    }
    if (total_send_qnames > 0) {
        send_qnames_flat.reset(
            aligned_alloc_custom(64, total_send_qnames));
        if (!send_qnames_flat) local_ok = 0;
    }
    if (!MdAllRanksOk(local_ok)) return -1;

    double setup_t0 = GetTime();
    MPI_Alltoall(send_candidate_counts.data(), 1,
                 MPI_UNSIGNED_LONG_LONG,
                 recv_candidate_counts.data(), 1,
                 MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
    MPI_Alltoall(send_qname_counts.data(), 1,
                 MPI_UNSIGNED_LONG_LONG,
                 recv_qname_counts.data(), 1,
                 MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

    unsigned long long recv_candidate_total = 0;
    unsigned long long recv_qname_total = 0;
    for (int i = 0; i < comm_size; ++i) {
        if (recv_candidate_total >
                ULLONG_MAX - recv_candidate_counts[(size_t)i] ||
            recv_qname_total >
                ULLONG_MAX - recv_qname_counts[(size_t)i]) {
            local_ok = 0;
            break;
        }
        recv_candidate_total += recv_candidate_counts[(size_t)i];
        recv_qname_total += recv_qname_counts[(size_t)i];
    }
    if (!MdAllRanksOk(local_ok)) return -1;

    if (recv_candidate_total >
            (unsigned long long)(SIZE_MAX /
                sizeof(MpiMarkdupCandidateShared)) ||
        recv_qname_total > (unsigned long long)SIZE_MAX ||
        self_candidate_count >
            SIZE_MAX - (size_t)recv_candidate_total ||
        self_qname_count > SIZE_MAX - (size_t)recv_qname_total) {
        local_ok = 0;
    }

    size_t send_candidate_bytes_total = 0;
    size_t exchange_memory = 0;
    if (local_ok) {
        send_candidate_bytes_total =
            total_send_candidates * sizeof(MpiMarkdupCandidateShared);
        exchange_memory = send_candidate_bytes_total;
        if (total_send_qnames > SIZE_MAX - exchange_memory) {
            local_ok = 0;
        } else {
            exchange_memory += total_send_qnames;
        }
    }
    size_t local_memory =
        local_candidates.capacity() *
            sizeof(MpiMarkdupCandidateShared) +
        local_qnames.capacity();
    if (remote_candidates.capacity() >
        SIZE_MAX / sizeof(remote_candidates[0])) {
        local_ok = 0;
    } else if (local_ok) {
        const size_t remote_ref_memory =
            remote_candidates.capacity() *
            sizeof(remote_candidates[0]);
        if (local_memory > SIZE_MAX - remote_ref_memory) {
            local_ok = 0;
        } else {
            local_memory += remote_ref_memory;
        }
    }
    if (local_ok) {
        size_t recv_owner_bytes = 0;
        if ((size_t)recv_candidate_total >
            SIZE_MAX / sizeof(MpiMarkdupCandidateShared)) {
            local_ok = 0;
        } else {
            recv_owner_bytes =
                (size_t)recv_candidate_total *
                sizeof(MpiMarkdupCandidateShared);
        }
        if (local_ok &&
            (size_t)recv_qname_total > SIZE_MAX - recv_owner_bytes) {
            local_ok = 0;
        } else if (local_ok) {
            recv_owner_bytes += (size_t)recv_qname_total;
        }
        if (local_ok &&
            (local_memory > SIZE_MAX - exchange_memory ||
             local_memory + exchange_memory >
                 SIZE_MAX - recv_owner_bytes ||
             MdCheckMemory(local_memory + exchange_memory +
                               recv_owner_bytes,
                           memory_limit, stats, rank,
                           "coordinate candidate owner exchange") != 0)) {
            local_ok = 0;
        }
    }
    if (!MdAllRanksOk(local_ok)) return -1;
    stats->t_candidate_exchange_setup += GetTime() - setup_t0;

    std::vector<size_t> candidate_write = send_candidate_offsets;
    std::vector<size_t> qname_write = send_qname_offsets;
    for (size_t i = 0; i < remote_candidates.size(); ++i) {
        const MdRemoteCandidateRef &ref = remote_candidates[i];
        MpiMarkdupCandidateShared candidate =
            local_candidates[ref.index];
        const int owner = ref.owner;
        if (!candidate.key.single) {
            const size_t qname_src = (size_t)candidate.qname_offset;
            const size_t qname_len = (size_t)candidate.qname_len;
            const size_t qname_dst = qname_write[(size_t)owner];
            candidate.qname_offset =
                (uint64_t)(qname_dst -
                           send_qname_offsets[(size_t)owner]);
            if (qname_len > 0) {
                memcpy(send_qnames_flat.get() + qname_dst,
                       local_qnames.data() + qname_src, qname_len);
            }
            qname_write[(size_t)owner] += qname_len;
        }
        send_candidates_flat.get()[candidate_write[(size_t)owner]++] =
            candidate;
    }
    for (size_t i = remote_candidates.size(); i > 0; --i) {
        const size_t remove_index = remote_candidates[i - 1].index;
        const size_t last_index = local_candidates.size() - 1;
        if (remove_index != last_index) {
            local_candidates[remove_index] =
                local_candidates[last_index];
        }
        local_candidates.pop_back();
    }
    remote_candidates.clear();
    owner_candidates->swap(local_candidates);
    owner_qnames->swap(local_qnames);
    const size_t needed_candidates =
        owner_candidates->size() + (size_t)recv_candidate_total;
    const size_t needed_qnames =
        owner_qnames->size() + (size_t)recv_qname_total;
    if (owner_candidates->capacity() < needed_candidates) {
        owner_candidates->reserve(needed_candidates);
    }
    if (owner_qnames->capacity() < needed_qnames) {
        owner_qnames->reserve(needed_qnames);
    }
    stats->t_candidate_pack += GetTime() - pack_t0;
    if (!MdAllRanksOk(local_ok)) return -1;

    double mpi_t0 = GetTime();
    for (int step = 1; step < comm_size; ++step) {
        int dest = (rank + step) % comm_size;
        int src = (rank - step + comm_size) % comm_size;
        unsigned long long send_candidate_count =
            send_candidate_counts[(size_t)dest];
        unsigned long long send_qname_count =
            send_qname_counts[(size_t)dest];
        unsigned long long recv_candidate_count =
            recv_candidate_counts[(size_t)src];
        unsigned long long recv_qname_count =
            recv_qname_counts[(size_t)src];
        if (recv_candidate_count >
            (unsigned long long)(SIZE_MAX /
                sizeof(MpiMarkdupCandidateShared)) ||
            recv_qname_count > (unsigned long long)SIZE_MAX) {
            return -1;
        }
        size_t candidate_base = owner_candidates->size();
        size_t qname_base = owner_qnames->size();
        owner_candidates->resize(
            candidate_base + (size_t)recv_candidate_count);
        owner_qnames->resize(qname_base + (size_t)recv_qname_count);

        unsigned long long send_candidate_bytes =
            send_candidate_count * sizeof(MpiMarkdupCandidateShared);
        unsigned long long recv_candidate_bytes =
            recv_candidate_count * sizeof(MpiMarkdupCandidateShared);
        if (MdSendrecvBytes(
                send_candidate_count
                    ? send_candidates_flat.get() +
                          send_candidate_offsets[(size_t)dest]
                    : nullptr,
                send_candidate_bytes, dest,
                recv_candidate_count
                    ? owner_candidates->data() + candidate_base
                    : nullptr,
                recv_candidate_bytes, src, 4111) != 0 ||
            MdSendrecvBytes(
                send_qname_count
                    ? send_qnames_flat.get() +
                          send_qname_offsets[(size_t)dest]
                    : nullptr,
                send_qname_count, dest,
                recv_qname_count
                    ? owner_qnames->data() + qname_base
                    : nullptr,
                recv_qname_count, src, 4112) != 0) {
            return -1;
        }
        for (size_t i = candidate_base;
             i < owner_candidates->size(); ++i) {
            if (!(*owner_candidates)[i].key.single) {
                (*owner_candidates)[i].qname_offset +=
                    (uint64_t)qname_base;
            }
        }
        stats->mpi_candidate_bytes +=
            (long long)send_candidate_bytes +
            (long long)send_qname_count;
    }

    stats->t_candidate_exchange_mpi += GetTime() - mpi_t0;
    stats->t_candidate_exchange += GetTime() - exchange_total_t0;
    stats->owner_candidates = (long long)owner_candidates->size();
    stats->qname_bytes = (long long)owner_qnames->size();
    return 0;
}

static int MdCompareQname(
        const MpiMarkdupCandidateShared &a,
        const MpiMarkdupCandidateShared &b,
        const std::vector<unsigned char> &qnames) {
    if (a.qname_offset > qnames.size() ||
        a.qname_len > qnames.size() - (size_t)a.qname_offset ||
        b.qname_offset > qnames.size() ||
        b.qname_len > qnames.size() - (size_t)b.qname_offset) {
        return 0;
    }
    size_t common = std::min((size_t)a.qname_len,
                             (size_t)b.qname_len);
    int cmp = common ? memcmp(
        qnames.data() + (size_t)a.qname_offset,
        qnames.data() + (size_t)b.qname_offset,
        common) : 0;
    if (cmp != 0) return cmp;
    if (a.qname_len < b.qname_len) return -1;
    if (a.qname_len > b.qname_len) return 1;
    return 0;
}

static bool MdPairBetter(
        const MpiMarkdupCandidateShared &a,
        const MpiMarkdupCandidateShared &b,
        const std::vector<unsigned char> &qnames) {
    if (a.qc_fail != b.qc_fail) return a.qc_fail < b.qc_fail;
    if (a.score != b.score) return a.score > b.score;
    int qcmp = MdCompareQname(a, b, qnames);
    if (qcmp != 0) return qcmp < 0;
    return a.global_order < b.global_order;
}

static bool MdCandidateGroupLess(
        const MpiMarkdupCandidateShared &a,
        const MpiMarkdupCandidateShared &b) {
    if (MdKeyLess(a.key, b.key)) return true;
    if (MdKeyLess(b.key, a.key)) return false;
    if (a.global_order != b.global_order) {
        return a.global_order < b.global_order;
    }
    if (a.source_rank != b.source_rank) {
        return a.source_rank < b.source_rank;
    }
    return a.ordinal < b.ordinal;
}

static bool MdKeyVectorContains(
        const std::vector<MpiMarkdupKeyShared> &keys,
        const MpiMarkdupKeyShared &key) {
    std::vector<MpiMarkdupKeyShared>::const_iterator it =
        std::lower_bound(
            keys.begin(), keys.end(), key,
            [](const MpiMarkdupKeyShared &a,
               const MpiMarkdupKeyShared &b) {
                return MdKeyLess(a, b);
            });
    return it != keys.end() && MdKeyEqual(*it, key);
}

static bool MdKeyRangeContains(
        const std::vector<MpiMarkdupKeyShared> &keys,
        size_t begin, size_t end,
        const MpiMarkdupKeyShared &key) {
    std::vector<MpiMarkdupKeyShared>::const_iterator first =
        keys.begin() + (ptrdiff_t)begin;
    std::vector<MpiMarkdupKeyShared>::const_iterator last =
        keys.begin() + (ptrdiff_t)end;
    std::vector<MpiMarkdupKeyShared>::const_iterator it =
        std::lower_bound(
            first, last, key,
            [](const MpiMarkdupKeyShared &a,
               const MpiMarkdupKeyShared &b) {
                return MdKeyLess(a, b);
            });
    return it != last && MdKeyEqual(*it, key);
}

static uint64_t MdCoordBucketCode(const MpiMarkdupKeyShared &key) {
    uint64_t hash = 1469598103934665603ull;
    auto mix = [&](uint64_t value) {
        for (int i = 0; i < 8; ++i) {
            hash ^= (unsigned char)(value & 0xffu);
            hash *= 1099511628211ull;
            value >>= 8;
        }
    };
    mix((uint64_t)(uint32_t)key.this_ref);
    mix((uint64_t)key.this_coord);
    return hash;
}

static size_t MdFindBucketIndex(
        const std::vector<uint64_t> &unique_codes, uint64_t code) {
    std::vector<uint64_t>::const_iterator it =
        std::lower_bound(unique_codes.begin(), unique_codes.end(), code);
    if (it == unique_codes.end() || *it != code) return SIZE_MAX;
    return (size_t)(it - unique_codes.begin());
}

static int MdBucketCandidateRangeByCoord(
        std::vector<MpiMarkdupCandidateShared> *candidates,
        size_t begin, size_t end,
        std::vector<MdBucketRange> *ranges) {
    ranges->clear();
    if (end <= begin) return 0;
    const size_t n = end - begin;
    const MpiMarkdupCandidateShared *src = candidates->data() + begin;

    if (n == 1) {
        MdBucketRange range;
        range.code = MdCoordBucketCode(src[0].key);
        range.begin = begin;
        range.end = end;
        ranges->push_back(range);
        return 0;
    }

    std::vector<uint64_t> unique_codes;
    unique_codes.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        unique_codes.push_back(MdCoordBucketCode(src[i].key));
    }
    std::sort(unique_codes.begin(), unique_codes.end());
    unique_codes.erase(std::unique(unique_codes.begin(),
                                   unique_codes.end()),
                       unique_codes.end());
    if (unique_codes.size() > (size_t)UINT_MAX) return -1;

    const size_t nb = unique_codes.size();
    std::vector<unsigned int> bucket_ids(n);
    std::vector<size_t> offsets(nb + 1, 0);
    for (size_t i = 0; i < n; ++i) {
        const uint64_t code = MdCoordBucketCode(src[i].key);
        const size_t bucket = MdFindBucketIndex(unique_codes, code);
        if (bucket == SIZE_MAX) return -1;
        bucket_ids[i] = (unsigned int)bucket;
        offsets[bucket + 1]++;
    }
    for (size_t i = 1; i <= nb; ++i) {
        offsets[i] += offsets[i - 1];
    }

    std::unique_ptr<MpiMarkdupCandidateShared,
                    MdMarkdupCandidateBufferDeleter> scratch(
        (MpiMarkdupCandidateShared *)aligned_alloc_custom(
            64, n * sizeof(MpiMarkdupCandidateShared)));
    if (!scratch) return -1;

    std::vector<size_t> write = offsets;
    for (size_t i = 0; i < n; ++i) {
        const size_t bucket = (size_t)bucket_ids[i];
        scratch.get()[write[bucket]++] = src[i];
    }
    memcpy(candidates->data() + begin, scratch.get(),
           n * sizeof(MpiMarkdupCandidateShared));

    ranges->reserve(nb);
    for (size_t i = 0; i < nb; ++i) {
        if (offsets[i] == offsets[i + 1]) continue;
        MdBucketRange range;
        range.code = unique_codes[i];
        range.begin = begin + offsets[i];
        range.end = begin + offsets[i + 1];
        ranges->push_back(range);
    }
    return 0;
}

static int MdBucketKeyRangeByCoord(
        std::vector<MpiMarkdupKeyShared> *keys,
        std::vector<MdBucketRange> *ranges) {
    ranges->clear();
    const size_t n = keys->size();
    if (n == 0) return 0;
    if (n == 1) {
        MdBucketRange range;
        range.code = MdCoordBucketCode((*keys)[0]);
        range.begin = 0;
        range.end = 1;
        ranges->push_back(range);
        return 0;
    }

    std::vector<uint64_t> unique_codes;
    unique_codes.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        unique_codes.push_back(MdCoordBucketCode((*keys)[i]));
    }
    std::sort(unique_codes.begin(), unique_codes.end());
    unique_codes.erase(std::unique(unique_codes.begin(),
                                   unique_codes.end()),
                       unique_codes.end());
    if (unique_codes.size() > (size_t)UINT_MAX) return -1;

    const size_t nb = unique_codes.size();
    std::vector<unsigned int> bucket_ids(n);
    std::vector<size_t> offsets(nb + 1, 0);
    for (size_t i = 0; i < n; ++i) {
        const uint64_t code = MdCoordBucketCode((*keys)[i]);
        const size_t bucket = MdFindBucketIndex(unique_codes, code);
        if (bucket == SIZE_MAX) return -1;
        bucket_ids[i] = (unsigned int)bucket;
        offsets[bucket + 1]++;
    }
    for (size_t i = 1; i <= nb; ++i) {
        offsets[i] += offsets[i - 1];
    }

    std::unique_ptr<MpiMarkdupKeyShared, MdMarkdupKeyBufferDeleter> scratch(
        (MpiMarkdupKeyShared *)aligned_alloc_custom(
            64, n * sizeof(MpiMarkdupKeyShared)));
    if (!scratch) return -1;

    std::vector<size_t> write = offsets;
    for (size_t i = 0; i < n; ++i) {
        const size_t bucket = (size_t)bucket_ids[i];
        scratch.get()[write[bucket]++] = (*keys)[i];
    }
    memcpy(keys->data(), scratch.get(),
           n * sizeof(MpiMarkdupKeyShared));

    ranges->reserve(nb);
    for (size_t i = 0; i < nb; ++i) {
        if (offsets[i] == offsets[i + 1]) continue;
        MdBucketRange range;
        range.code = unique_codes[i];
        range.begin = offsets[i];
        range.end = offsets[i + 1];
        ranges->push_back(range);
    }
    return 0;
}

static int MdPushDuplicate(
        const MpiMarkdupCandidateShared &candidate, int comm_size,
        std::vector<std::vector<uint64_t> > *duplicates_by_source) {
    const int source = candidate.source_rank;
    if (source < 0 || source >= comm_size) return -1;
    (*duplicates_by_source)[(size_t)source].push_back(candidate.ordinal);
    return 0;
}

static int MdScanPairDuplicateRange(
        const std::vector<MpiMarkdupCandidateShared> &candidates,
        const std::vector<unsigned char> &qnames,
        size_t begin, size_t end, int comm_size,
        std::vector<std::vector<uint64_t> > *duplicates_by_source,
        MpiMarkdupStats *stats) {
    while (begin < end) {
        size_t group_end = begin + 1;
        while (group_end < end &&
               MdKeyEqual(candidates[begin].key,
                          candidates[group_end].key)) {
            ++group_end;
        }

        size_t winner = begin;
        for (size_t i = begin + 1; i < group_end; ++i) {
            if (MdPairBetter(candidates[i], candidates[winner], qnames)) {
                winner = i;
            }
        }
        for (size_t i = begin; i < group_end; ++i) {
            if (i == winner) continue;
            if (MdPushDuplicate(candidates[i], comm_size,
                                duplicates_by_source) != 0) {
                return -1;
            }
            stats->pair_duplicates++;
        }
        begin = group_end;
    }
    return 0;
}

static int MdScanSingleDuplicateRange(
        const std::vector<MpiMarkdupCandidateShared> &candidates,
        const std::vector<MpiMarkdupKeyShared> &paired_marker_keys,
        size_t begin, size_t end,
        size_t marker_begin, size_t marker_end,
        int comm_size,
        std::vector<std::vector<uint64_t> > *duplicates_by_source,
        MpiMarkdupStats *stats) {
    while (begin < end) {
        size_t group_end = begin + 1;
        while (group_end < end &&
               MdKeyEqual(candidates[begin].key,
                          candidates[group_end].key)) {
            ++group_end;
        }

        if (MdKeyRangeContains(paired_marker_keys, marker_begin,
                               marker_end, candidates[begin].key)) {
            for (size_t i = begin; i < group_end; ++i) {
                if (MdPushDuplicate(candidates[i], comm_size,
                                    duplicates_by_source) != 0) {
                    return -1;
                }
                stats->single_duplicates++;
            }
        } else if (group_end - begin > 1) {
            size_t winner = begin;
            for (size_t i = begin + 1; i < group_end; ++i) {
                const MpiMarkdupCandidateShared &candidate = candidates[i];
                const MpiMarkdupCandidateShared &best = candidates[winner];
                if (candidate.score > best.score ||
                    (candidate.score == best.score &&
                     candidate.global_order < best.global_order)) {
                    winner = i;
                }
            }
            for (size_t i = begin; i < group_end; ++i) {
                if (i == winner) continue;
                if (MdPushDuplicate(candidates[i], comm_size,
                                    duplicates_by_source) != 0) {
                    return -1;
                }
                stats->single_duplicates++;
            }
        }
        begin = group_end;
    }
    return 0;
}

static int MdFindDuplicatesCoordinateBuckets(
        std::vector<MpiMarkdupCandidateShared> *owner_candidates,
        const std::vector<unsigned char> &owner_qnames,
        int comm_size,
        std::vector<std::vector<uint64_t> > *duplicates_by_source,
        MpiMarkdupStats *stats) {
    double group_t0 = GetTime();
    std::vector<MpiMarkdupCandidateShared>::iterator pair_end =
        std::partition(
            owner_candidates->begin(), owner_candidates->end(),
            [](const MpiMarkdupCandidateShared &candidate) {
                return !candidate.key.single;
            });
    const size_t pair_count =
        (size_t)(pair_end - owner_candidates->begin());

    std::vector<MpiMarkdupKeyShared> paired_marker_keys;
    paired_marker_keys.reserve(owner_candidates->size() - pair_count);
    size_t real_single_write = pair_count;
    for (size_t i = pair_count; i < owner_candidates->size(); ++i) {
        const MpiMarkdupCandidateShared candidate =
            (*owner_candidates)[i];
        if (candidate.paired_marker) {
            paired_marker_keys.push_back(candidate.key);
        } else {
            (*owner_candidates)[real_single_write++] = candidate;
        }
    }
    owner_candidates->resize(real_single_write);
    const size_t real_single_end = real_single_write;

    double group_sort_t0 = GetTime();
    std::vector<MdBucketRange> pair_ranges;
    std::vector<MdBucketRange> marker_ranges;
    std::vector<MdBucketRange> single_ranges;
    if (MdBucketCandidateRangeByCoord(owner_candidates, 0, pair_count,
                                      &pair_ranges) != 0 ||
        MdBucketKeyRangeByCoord(&paired_marker_keys,
                                &marker_ranges) != 0 ||
        MdBucketCandidateRangeByCoord(owner_candidates, pair_count,
                                      real_single_end,
                                      &single_ranges) != 0) {
        return -1;
    }
    for (size_t i = 0; i < pair_ranges.size(); ++i) {
        const MdBucketRange &range = pair_ranges[i];
        if (range.end - range.begin > 1) {
            std::sort(owner_candidates->begin() +
                          (ptrdiff_t)range.begin,
                      owner_candidates->begin() +
                          (ptrdiff_t)range.end,
                      MdCandidateGroupLess);
        }
    }
    for (size_t i = 0; i < marker_ranges.size(); ++i) {
        const MdBucketRange &range = marker_ranges[i];
        if (range.end - range.begin > 1) {
            std::sort(paired_marker_keys.begin() +
                          (ptrdiff_t)range.begin,
                      paired_marker_keys.begin() +
                          (ptrdiff_t)range.end,
                      [](const MpiMarkdupKeyShared &a,
                         const MpiMarkdupKeyShared &b) {
                          return MdKeyLess(a, b);
                      });
        }
    }
    for (size_t i = 0; i < single_ranges.size(); ++i) {
        const MdBucketRange &range = single_ranges[i];
        if (range.end - range.begin > 1) {
            std::sort(owner_candidates->begin() +
                          (ptrdiff_t)range.begin,
                      owner_candidates->begin() +
                          (ptrdiff_t)range.end,
                      MdCandidateGroupLess);
        }
    }
    stats->t_group_sort += GetTime() - group_sort_t0;

    double group_scan_t0 = GetTime();
    duplicates_by_source->assign(
        (size_t)comm_size, std::vector<uint64_t>());

    for (size_t i = 0; i < pair_ranges.size(); ++i) {
        const MdBucketRange &range = pair_ranges[i];
        if (MdScanPairDuplicateRange(
                *owner_candidates, owner_qnames, range.begin,
                range.end, comm_size, duplicates_by_source,
                stats) != 0) {
            return -1;
        }
    }

    size_t marker_idx = 0;
    for (size_t i = 0; i < single_ranges.size(); ++i) {
        const MdBucketRange &range = single_ranges[i];
        while (marker_idx < marker_ranges.size() &&
               marker_ranges[marker_idx].code < range.code) {
            ++marker_idx;
        }
        size_t marker_begin = 0;
        size_t marker_end = 0;
        if (marker_idx < marker_ranges.size() &&
            marker_ranges[marker_idx].code == range.code) {
            marker_begin = marker_ranges[marker_idx].begin;
            marker_end = marker_ranges[marker_idx].end;
        }
        if (MdScanSingleDuplicateRange(
                *owner_candidates, paired_marker_keys,
                range.begin, range.end,
                marker_begin, marker_end, comm_size,
                duplicates_by_source, stats) != 0) {
            return -1;
        }
    }
    stats->t_group_scan += GetTime() - group_scan_t0;
    stats->t_group += GetTime() - group_t0;
    return 0;
}

static int MdFindDuplicates(
        std::vector<MpiMarkdupCandidateShared> *owner_candidates,
        const std::vector<unsigned char> &owner_qnames,
        int comm_size,
        std::vector<std::vector<uint64_t> > *duplicates_by_source,
        MpiMarkdupStats *stats) {

    double group_t0 = GetTime();
    std::vector<MpiMarkdupCandidateShared>::iterator pair_end =
        std::partition(
            owner_candidates->begin(), owner_candidates->end(),
            [](const MpiMarkdupCandidateShared &candidate) {
                return !candidate.key.single;
            });
    const size_t pair_count =
        (size_t)(pair_end - owner_candidates->begin());

    std::vector<MpiMarkdupKeyShared> paired_marker_keys;
    paired_marker_keys.reserve(owner_candidates->size() - pair_count);
    size_t real_single_write = pair_count;
    for (size_t i = pair_count; i < owner_candidates->size(); ++i) {
        const MpiMarkdupCandidateShared candidate =
            (*owner_candidates)[i];
        if (candidate.paired_marker) {
            paired_marker_keys.push_back(candidate.key);
        } else {
            (*owner_candidates)[real_single_write++] = candidate;
        }
    }
    owner_candidates->resize(real_single_write);
    const size_t real_single_count = real_single_write - pair_count;

    double group_sort_t0 = GetTime();
    std::sort(owner_candidates->begin(),
              owner_candidates->begin() + (ptrdiff_t)pair_count,
              MdCandidateGroupLess);
    std::sort(paired_marker_keys.begin(), paired_marker_keys.end(),
              [](const MpiMarkdupKeyShared &a,
                 const MpiMarkdupKeyShared &b) {
                  return MdKeyLess(a, b);
              });
    paired_marker_keys.erase(
        std::unique(
            paired_marker_keys.begin(), paired_marker_keys.end(),
            [](const MpiMarkdupKeyShared &a,
               const MpiMarkdupKeyShared &b) {
                return MdKeyEqual(a, b);
            }),
        paired_marker_keys.end());
    std::sort(owner_candidates->begin() + (ptrdiff_t)pair_count,
              owner_candidates->end(), MdCandidateGroupLess);
    stats->t_group_sort += GetTime() - group_sort_t0;

    double group_scan_t0 = GetTime();
    duplicates_by_source->assign(
        (size_t)comm_size, std::vector<uint64_t>());

    size_t begin = 0;
    while (begin < pair_count) {
        size_t end = begin + 1;
        while (end < pair_count &&
               MdKeyEqual((*owner_candidates)[begin].key,
                          (*owner_candidates)[end].key)) {
            ++end;
        }

        size_t winner = begin;
        for (size_t i = begin + 1; i < end; ++i) {
            if (MdPairBetter((*owner_candidates)[i],
                             (*owner_candidates)[winner],
                             owner_qnames)) {
                winner = i;
            }
        }
        for (size_t i = begin; i < end; ++i) {
            if (i == winner) continue;
            int source = (*owner_candidates)[i].source_rank;
            if (source < 0 || source >= comm_size) return -1;
            (*duplicates_by_source)[(size_t)source].push_back(
                (*owner_candidates)[i].ordinal);
            stats->pair_duplicates++;
        }
        begin = end;
    }

    begin = pair_count;
    const size_t real_single_end = pair_count + real_single_count;
    while (begin < real_single_end) {
        size_t end = begin + 1;
        while (end < real_single_end &&
               MdKeyEqual((*owner_candidates)[begin].key,
                          (*owner_candidates)[end].key)) {
            ++end;
        }

        if (MdKeyVectorContains(
                paired_marker_keys, (*owner_candidates)[begin].key)) {
            for (size_t i = begin; i < end; ++i) {
                int source = (*owner_candidates)[i].source_rank;
                if (source < 0 || source >= comm_size) return -1;
                (*duplicates_by_source)[(size_t)source]
                    .push_back((*owner_candidates)[i].ordinal);
                stats->single_duplicates++;
            }
        } else if (end - begin > 1) {
            size_t winner = begin;
            for (size_t i = begin + 1; i < end; ++i) {
                const MpiMarkdupCandidateShared &candidate =
                    (*owner_candidates)[i];
                const MpiMarkdupCandidateShared &best =
                    (*owner_candidates)[winner];
                if (candidate.score > best.score ||
                    (candidate.score == best.score &&
                     candidate.global_order < best.global_order)) {
                    winner = i;
                }
            }
            for (size_t i = begin; i < end; ++i) {
                if (i == winner) continue;
                int source = (*owner_candidates)[i].source_rank;
                if (source < 0 || source >= comm_size) return -1;
                (*duplicates_by_source)[(size_t)source]
                    .push_back((*owner_candidates)[i].ordinal);
                stats->single_duplicates++;
            }
        }
        begin = end;
    }
    stats->t_group_scan += GetTime() - group_scan_t0;
    stats->t_group += GetTime() - group_t0;
    return 0;
}

static int MdExchangeDuplicateResults(
        std::vector<std::vector<uint64_t> > *duplicates_by_source,
        int rank, int comm_size, uint64_t local_records,
        size_t memory_limit, std::vector<uint8_t> *bitmap,
        MpiMarkdupStats *stats) {
    double exchange_t0 = GetTime();
    // 先对 duplicate ordinal 去重
    for (int i = 0; i < comm_size; ++i) {
        std::vector<uint64_t> &values =
            (*duplicates_by_source)[(size_t)i];
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()),
                     values.end());
    }
    // 创建本 rank 的 duplicate bitmap
    int local_ok =
        local_records <= (uint64_t)SIZE_MAX * 8ull ? 1 : 0;
    if (!MdAllRanksOk(local_ok)) return -1;
    bitmap->assign((size_t)((local_records + 7) / 8), 0);

    // 先用 MPI_Alltoall 交换每个 rank 会收多少结果
    std::vector<unsigned long long> send_counts((size_t)comm_size, 0);
    std::vector<unsigned long long> recv_counts((size_t)comm_size, 0);
    for (int i = 0; i < comm_size; ++i) {
        send_counts[(size_t)i] =
            (unsigned long long)
                (*duplicates_by_source)[(size_t)i].size();
    }
    MPI_Alltoall(send_counts.data(), 1, MPI_UNSIGNED_LONG_LONG,
                 recv_counts.data(), 1, MPI_UNSIGNED_LONG_LONG,
                 MPI_COMM_WORLD);

    // 内存检查
    unsigned long long max_receive = 0;
    size_t result_memory = bitmap->capacity();
    for (int i = 0; i < comm_size; ++i) {
        max_receive = std::max(max_receive,
                               recv_counts[(size_t)i]);
        size_t bytes =
            (*duplicates_by_source)[(size_t)i].capacity() *
            sizeof(uint64_t);
        if (bytes > SIZE_MAX - result_memory) {
            local_ok = 0;
            break;
        }
        result_memory += bytes;
    }
    if (max_receive >
        (unsigned long long)(SIZE_MAX / sizeof(uint64_t))) {
        local_ok = 0;
    } else if (local_ok) {
        size_t receive_bytes =
            (size_t)max_receive * sizeof(uint64_t);
        if (receive_bytes > SIZE_MAX - result_memory ||
            MdCheckMemory(result_memory + receive_bytes,
                          memory_limit, stats, rank,
                          "duplicate result exchange") != 0) {
            local_ok = 0;
        }
    }
    if (!MdAllRanksOk(local_ok)) return -1;

    // 使用 ring MPI_Sendrecv 回传 duplicate 结果
    for (int step = 0; step < comm_size; ++step) {
        int dest = (rank + step) % comm_size;
        int src = (rank - step + comm_size) % comm_size;
        if (step == 0) {
            const std::vector<uint64_t> &values =
                (*duplicates_by_source)[(size_t)rank];
            for (size_t i = 0; i < values.size(); ++i) {
                if (values[i] >= local_records) {
                    local_ok = 0;
                    continue;
                }
                (*bitmap)[(size_t)(values[i] >> 3)] |=
                    (uint8_t)(1u << (values[i] & 7u));
            }
            continue;
        }
        unsigned long long send_count =
            (unsigned long long)
                (*duplicates_by_source)[(size_t)dest].size();
        unsigned long long recv_count = 0;
        if (MPI_Sendrecv(&send_count, 1, MPI_UNSIGNED_LONG_LONG,
                         dest, 4200,
                         &recv_count, 1, MPI_UNSIGNED_LONG_LONG,
                         src, 4200, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE) != MPI_SUCCESS) {
            return -1;
        }
        std::vector<uint64_t> received((size_t)recv_count);
        if (MdSendrecvBytes(
                (*duplicates_by_source)[(size_t)dest].data(),
                send_count * sizeof(uint64_t), dest,
                received.data(),
                recv_count * sizeof(uint64_t), src, 4201) != 0) {
            return -1;
        }
        stats->mpi_result_bytes +=
            (long long)(send_count * sizeof(uint64_t));
        for (size_t i = 0; i < received.size(); ++i) {
            if (received[i] >= local_records) {
                fprintf(stderr,
                        "[rank %d] ERROR: received invalid duplicate ordinal "
                        "%llu, local_records=%llu.\n",
                        rank,
                        (unsigned long long)received[i],
                        (unsigned long long)local_records);
                local_ok = 0;
                continue;
            }
            (*bitmap)[(size_t)(received[i] >> 3)] |=
                (uint8_t)(1u << (received[i] & 7u));
        }
    }
    stats->t_result_exchange += GetTime() - exchange_t0;
    return MdAllRanksOk(local_ok) ? 0 : -1;
}

static void MdResetPackWorkspace(MdPackWorkspace *workspace) {
    workspace->active_blocks = 0;
    workspace->total_records = 0;
    workspace->current_begin =
        workspace->records.empty() ? nullptr
                                   : workspace->records.data();
    workspace->current_records = 0;
    workspace->current_len = 0;
}

static void MdInitPackWorkspace(MdPackWorkspace *workspace,
                                size_t capacity) {
    workspace->records.resize(capacity);
    MdResetPackWorkspace(workspace);
}

static int MdSealPackBlock(MdPackWorkspace *workspace) {
    if (workspace->current_records == 0) return 0;
    if (workspace->active_blocks >= kMarkdupNB) return -1;
    MdPackPlan &plan =
        workspace->plans[workspace->active_blocks++];
    plan.records = workspace->current_begin;
    plan.n_records = workspace->current_records;
    plan.total_len = workspace->current_len;
    workspace->current_begin =
        workspace->records.data() + workspace->total_records;
    workspace->current_records = 0;
    workspace->current_len = 0;
    return 0;
}

static int MdAppendRecords(MdPackWorkspace *workspace,
                           bam1_t **records,
                           const uint32_t *bam_lens,
                           int n_records,
                           uint32_t total_len) {
    if (n_records <= 0) return 0;
    if (total_len > 0 && total_len <= BGZF_BLOCK_SIZE &&
        (workspace->current_records == 0 ||
         workspace->current_len + total_len <=
             BGZF_BLOCK_SIZE)) {
        if (workspace->total_records + n_records >
            (int)workspace->records.size()) {
            return -2;
        }
        if (workspace->current_records == 0) {
            workspace->current_begin =
                workspace->records.data() +
                workspace->total_records;
        }
        memcpy(workspace->records.data() +
                   workspace->total_records,
               records, (size_t)n_records *
                            sizeof(bam1_t *));
        workspace->total_records += n_records;
        workspace->current_records += n_records;
        workspace->current_len += total_len;
        return 0;
    }

    int offset = 0;
    while (offset < n_records) {
        if (workspace->current_records == 0) {
            workspace->current_begin =
                workspace->records.data() +
                workspace->total_records;
        }
        uint32_t room =
            BGZF_BLOCK_SIZE - workspace->current_len;
        int take = 0;
        uint32_t take_len = 0;
        while (offset + take < n_records) {
            uint32_t packed_len =
                bam_lens[offset + take] + 4;
            if (packed_len > BGZF_BLOCK_SIZE) return -1;
            if (take_len + packed_len > room) break;
            take_len += packed_len;
            ++take;
        }
        if (take == 0) {
            if (MdSealPackBlock(workspace) != 0) return -3;
            continue;
        }
        if (workspace->total_records + take >
            (int)workspace->records.size()) {
            return -2;
        }
        memcpy(workspace->records.data() +
                   workspace->total_records,
               records + offset,
               (size_t)take * sizeof(bam1_t *));
        workspace->total_records += take;
        workspace->current_records += take;
        workspace->current_len += take_len;
        offset += take;
        if (workspace->current_len == BGZF_BLOCK_SIZE &&
            MdSealPackBlock(workspace) != 0) {
            return -3;
        }
    }
    return 0;
}

static void MdInitEmptyRawCompress(MpiSortRawCompressPara *para,
                                   int block_id) {
    memset(para, 0, sizeof(*para));
    para->block_id = block_id;
    para->status = -1;
    para->compress_level = 1;
}

static void MdSetupRawCompress(MpiSortRawCompressPara *para,
                               int block_id, bam_block *payload,
                               bam_block *output,
                               int compress_level) {
    memset(para, 0, sizeof(*para));
    para->block_id = block_id;
    para->un_comp_block = payload;
    para->un_comp_size = payload ? payload->pos : 0;
    para->output_block = output;
    para->output_size = 0;
    para->status = para->un_comp_size > 0 ? 0 : -1;
    para->compress_level = compress_level;
}

static void MdResetPayloadBlock(bam_block *block, int block_id) {
    block->pos = 0;
    block->length = 0;
    block->errcode = 0;
    block->block_id = block_id;
    block->block_address = 0;
}

static uint32_t MdLoadLe32(const unsigned char *p) {
    uint32_t v = 0;
    memcpy(&v, p, sizeof(v));
    return v;
}

static int MdAppendRawPayload(bam_block *payload,
                              const unsigned char *record,
                              uint32_t raw_len) {
    if (!payload || !record || raw_len == 0 ||
        raw_len > BGZF_BLOCK_SIZE) {
        return -1;
    }
    if (payload->pos > 0 &&
        payload->pos + (int)raw_len > BGZF_BLOCK_SIZE) {
        return 1;
    }
    if (payload->pos + (int)raw_len > BGZF_BLOCK_SIZE) {
        return -1;
    }
    memcpy(payload->data + payload->pos, record, raw_len);
    payload->pos += (int)raw_len;
    payload->length = payload->pos;
    return 0;
}

static int MdAppendRawPayloadPatchedFlag(bam_block *payload,
                                         const unsigned char *record,
                                         uint32_t raw_len,
                                         uint16_t final_flag) {
    if (!payload || !record || raw_len == 0 ||
        raw_len > BGZF_BLOCK_SIZE) {
        return -1;
    }
    if (payload->pos > 0 &&
        payload->pos + (int)raw_len > BGZF_BLOCK_SIZE) {
        return 1;
    }
    if (payload->pos + (int)raw_len > BGZF_BLOCK_SIZE) {
        return -1;
    }

    unsigned char *dst =
        (unsigned char *)payload->data + payload->pos;
    memcpy(dst, record, raw_len);
    uint32_t flag_nc = MdLoadLe32(dst + 16);
    flag_nc = (flag_nc & 0xffffu) |
              ((uint32_t)final_flag << 16);
    memcpy(dst + 16, &flag_nc, sizeof(flag_nc));
    payload->pos += (int)raw_len;
    payload->length = payload->pos;
    return 0;
}

static const unsigned char *MdRawAuxPayloadEnd(
        const unsigned char *tag, const unsigned char *end) {
    if (!tag || tag + 3 > end) return nullptr;
    const unsigned char type = tag[2];
    const unsigned char *p = tag + 3;
    switch (type) {
    case 'A':
    case 'c':
    case 'C':
        return p + 1 <= end ? p + 1 : nullptr;
    case 's':
    case 'S':
        return p + 2 <= end ? p + 2 : nullptr;
    case 'i':
    case 'I':
    case 'f':
        return p + 4 <= end ? p + 4 : nullptr;
    case 'd':
        return p + 8 <= end ? p + 8 : nullptr;
    case 'Z':
    case 'H': {
        while (p < end) {
            if (*p == '\0') return p + 1;
            ++p;
        }
        return nullptr;
    }
    case 'B': {
        if (p + 5 > end) return nullptr;
        size_t elem_size = 0;
        switch (p[0]) {
        case 'c':
        case 'C':
        case 'A':
            elem_size = 1;
            break;
        case 's':
        case 'S':
            elem_size = 2;
            break;
        case 'i':
        case 'I':
        case 'f':
            elem_size = 4;
            break;
        default:
            return nullptr;
        }
        const uint32_t count = MdLoadLe32(p + 1);
        if (count > (uint32_t)((size_t)(end - (p + 5)) /
                               elem_size)) {
            return nullptr;
        }
        return p + 5 + (size_t)count * elem_size;
    }
    default:
        return nullptr;
    }
}

static int MdAppendRawPayloadClearOld(bam_block *payload,
                                      const unsigned char *record,
                                      uint32_t block_len,
                                      uint16_t final_flag) {
    if (!payload || !record || block_len < 32) return -1;
    const uint32_t raw_len = block_len + 4u;
    if (raw_len > BGZF_BLOCK_SIZE) return -1;

    const uint32_t bin_mq_nl = MdLoadLe32(record + 12);
    const uint32_t flag_nc = MdLoadLe32(record + 16);
    const uint32_t raw_l_qname = bin_mq_nl & 0xffu;
    const uint32_t n_cigar = flag_nc & 0xffffu;
    const uint32_t l_qseq = MdLoadLe32(record + 20);
    const uint64_t aux_start64 =
        36ull + (uint64_t)raw_l_qname +
        4ull * (uint64_t)n_cigar +
        (((uint64_t)l_qseq + 1ull) >> 1) +
        (uint64_t)l_qseq;
    if (raw_l_qname == 0 || aux_start64 > raw_len) return -1;

    const unsigned char *aux_start =
        record + (size_t)aux_start64;
    const unsigned char *record_end = record + raw_len;
    const size_t aux_len = (size_t)(record_end - aux_start);
    if (aux_len == 0 || memchr(aux_start, 'd', aux_len) == nullptr) {
        return MdAppendRawPayloadPatchedFlag(
            payload, record, raw_len, final_flag);
    }

    const unsigned char *scan = aux_start;
    uint32_t removed = 0;
    while (scan + 3 <= record_end) {
        const unsigned char *next =
            MdRawAuxPayloadEnd(scan, record_end);
        if (!next) break;
        if (scan[0] == 'd' &&
            (scan[1] == 't' || scan[1] == 'o')) {
            removed += (uint32_t)(next - scan);
        }
        scan = next;
    }

    if (removed == 0) {
        return MdAppendRawPayloadPatchedFlag(
            payload, record, raw_len, final_flag);
    }

    const uint32_t output_block_len = block_len - removed;
    const uint32_t output_raw_len = output_block_len + 4u;
    if (payload->pos > 0 &&
        payload->pos + (int)output_raw_len > BGZF_BLOCK_SIZE) {
        return 1;
    }
    if (payload->pos + (int)output_raw_len > BGZF_BLOCK_SIZE) {
        return -1;
    }

    unsigned char *dst =
        (unsigned char *)payload->data + payload->pos;
    memcpy(dst, &output_block_len, 4);
    memcpy(dst + 4, record + 4, (size_t)aux_start64 - 4u);
    uint32_t patched_flag_nc = MdLoadLe32(dst + 16);
    patched_flag_nc = (patched_flag_nc & 0xffffu) |
                      ((uint32_t)final_flag << 16);
    memcpy(dst + 16, &patched_flag_nc, sizeof(patched_flag_nc));

    unsigned char *out = dst + (size_t)aux_start64;
    scan = aux_start;
    while (scan < record_end) {
        const unsigned char *next =
            MdRawAuxPayloadEnd(scan, record_end);
        if (!next) {
            memcpy(out, scan, (size_t)(record_end - scan));
            out += record_end - scan;
            break;
        }
        if (!(scan[0] == 'd' &&
              (scan[1] == 't' || scan[1] == 'o'))) {
            memcpy(out, scan, (size_t)(next - scan));
            out += next - scan;
        }
        scan = next;
    }
    if ((uint32_t)(out - dst) != output_raw_len) return -2;
    payload->pos += (int)output_raw_len;
    payload->length = payload->pos;
    return 0;
}

static int MdRawRecordCanPatchFlagOnly(const unsigned char *raw,
                                       uint32_t block_len,
                                       const bam1_t *record,
                                       uint32_t final_bam_len) {
    if (!raw || !record || !record->data) return 0;
    const bam1_core_t *core = &record->core;
    if (final_bam_len != block_len) return 0;
    if (core->n_cigar > 0xffff) return 0;
    const int raw_l_qname =
        (int)core->l_qname - (int)core->l_extranul;
    if (raw_l_qname <= 0 || raw_l_qname > 255) return 0;

    const uint32_t raw_tid = MdLoadLe32(raw + 4);
    const uint32_t raw_pos = MdLoadLe32(raw + 8);
    const uint32_t raw_bin_mq_nl = MdLoadLe32(raw + 12);
    const uint32_t raw_flag_nc = MdLoadLe32(raw + 16);
    const uint32_t raw_l_qseq = MdLoadLe32(raw + 20);
    const uint32_t raw_mtid = MdLoadLe32(raw + 24);
    const uint32_t raw_mpos = MdLoadLe32(raw + 28);
    const uint32_t raw_isize = MdLoadLe32(raw + 32);

    if ((int32_t)raw_tid != core->tid ||
        (int32_t)raw_pos != core->pos ||
        (int32_t)raw_mtid != core->mtid ||
        (int32_t)raw_mpos != core->mpos ||
        (int32_t)raw_isize != core->isize) {
        return 0;
    }
    if ((int32_t)raw_l_qseq != core->l_qseq) return 0;
    if ((raw_bin_mq_nl & 0xffu) != (uint32_t)raw_l_qname) {
        return 0;
    }
    if (((raw_bin_mq_nl >> 8) & 0xffu) !=
        (uint32_t)core->qual) {
        return 0;
    }
    if ((raw_bin_mq_nl >> 16) != (uint32_t)core->bin) {
        return 0;
    }
    if ((raw_flag_nc & 0xffffu) !=
        ((uint32_t)core->n_cigar & 0xffffu)) {
        return 0;
    }
    return 1;
}

static int MdAppendBamRecordPayload(bam_block *payload,
                                    const bam1_t *record) {
    if (!payload || !record || !record->data) return -1;
    const bam1_core_t *core = &record->core;
    const int raw_l_qname =
        (int)core->l_qname - (int)core->l_extranul;
    if (raw_l_qname <= 0 || raw_l_qname > 255 ||
        core->n_cigar > 0xffff ||
        record->l_data < (int)core->l_qname) {
        return -1;
    }
    const uint32_t rest_len =
        (uint32_t)(record->l_data - core->l_qname);
    const uint32_t block_len =
        (uint32_t)(record->l_data - core->l_extranul + 32);
    const uint32_t raw_len = block_len + 4u;
    if (raw_len > BGZF_BLOCK_SIZE) return -1;
    if (payload->pos > 0 &&
        payload->pos + (int)raw_len > BGZF_BLOCK_SIZE) {
        return 1;
    }
    if (payload->pos + (int)raw_len > BGZF_BLOCK_SIZE) {
        return -1;
    }

    uint32_t fields[8];
    fields[0] = (uint32_t)core->tid;
    fields[1] = (uint32_t)core->pos;
    fields[2] = ((uint32_t)core->bin << 16) |
                ((uint32_t)core->qual << 8) |
                (uint32_t)raw_l_qname;
    fields[3] = ((uint32_t)core->flag << 16) |
                ((uint32_t)core->n_cigar & 0xffffu);
    fields[4] = (uint32_t)core->l_qseq;
    fields[5] = (uint32_t)core->mtid;
    fields[6] = (uint32_t)core->mpos;
    fields[7] = (uint32_t)core->isize;

    unsigned char *dst =
        (unsigned char *)payload->data + payload->pos;
    memcpy(dst, &block_len, 4);
    memcpy(dst + 4, fields, sizeof(fields));
    memcpy(dst + 36, record->data, (size_t)raw_l_qname);
    memcpy(dst + 36 + raw_l_qname,
           record->data + core->l_qname, rest_len);
    payload->pos += (int)raw_len;
    payload->length = payload->pos;
    return 0;
}

static int MdPrepareRemoveDupsPayloadFromRaw(
        Bam2BamPara *decomp, MdBlockSet *uncompressed,
        int n_blocks, const std::vector<uint8_t> &bitmap,
        uint64_t ordinal_base, MdBlockSet *payload_blocks,
        MdBlockSet *output_blocks, MpiSortRawCompressPara *paras,
        int compress_level, int *active_blocks,
        uint64_t *ordinal_end, MpiMarkdupStats *stats) {
    int active = 0;
    uint64_t ordinal = ordinal_base;
    if (!decomp || !uncompressed || !payload_blocks ||
        !output_blocks || !paras || !active_blocks ||
        !ordinal_end || payload_blocks->n < kMarkdupNB ||
        output_blocks->n < kMarkdupNB) {
        return -1;
    }

    MdResetPayloadBlock(payload_blocks->blocks + active, active);
    for (int b = 0; b < n_blocks; ++b) {
        bam_block *src = uncompressed->blocks + b;
        size_t pos = 0;
        for (int r = 0; r < decomp[b].n_total_records; ++r) {
            if (pos + 36u > src->length) return -11;
            const unsigned char *raw = src->data + pos;
            const uint32_t block_len = MdLoadLe32(raw);
            const uint64_t raw_len64 = (uint64_t)block_len + 4u;
            if (block_len < 32 ||
                raw_len64 > BGZF_BLOCK_SIZE ||
                raw_len64 > (uint64_t)src->length - pos) {
                return -12;
            }
            const uint32_t flag_nc = MdLoadLe32(raw + 16);
            const uint16_t flag = (uint16_t)(flag_nc >> 16);
            const size_t byte = (size_t)(ordinal >> 3);
            const uint8_t mask = (uint8_t)(1u << (ordinal & 7u));
            const int duplicate =
                byte < bitmap.size() && (bitmap[byte] & mask);

            if (duplicate && stats) stats->marked_records++;
            if (duplicate || (flag & BAM_FDUP)) {
                if (stats) stats->removed_records++;
                ++ordinal;
                pos += (size_t)raw_len64;
                continue;
            }

            while (true) {
                int ret = MdAppendRawPayload(
                    payload_blocks->blocks + active,
                    raw, (uint32_t)raw_len64);
                if (ret == 0) break;
                if (ret < 0) return ret;

                MdSetupRawCompress(paras + active, active,
                                   payload_blocks->blocks + active,
                                   output_blocks->blocks + active,
                                   compress_level);
                ++active;
                if (active >= kMarkdupNB) return -13;
                MdResetPayloadBlock(payload_blocks->blocks + active,
                                    active);
            }
            ++ordinal;
            pos += (size_t)raw_len64;
        }
        if (pos != src->length) return -14;
    }

    if (payload_blocks->blocks[active].pos > 0) {
        MdSetupRawCompress(paras + active, active,
                           payload_blocks->blocks + active,
                           output_blocks->blocks + active,
                           compress_level);
        ++active;
    }
    for (int i = active; i < kMarkdupNB; ++i) {
        MdInitEmptyRawCompress(paras + i, i);
    }
    *active_blocks = active;
    *ordinal_end = ordinal;
    return 0;
}

static int MdPrepareClearRemoveDupsPayloadFromRaw(
        Bam2BamPara *decomp, MdBlockSet *uncompressed,
        int n_blocks, const std::vector<uint8_t> &bitmap,
        uint64_t ordinal_base, MdBlockSet *payload_blocks,
        MdBlockSet *output_blocks, MpiSortRawCompressPara *paras,
        int compress_level, int *active_blocks,
        uint64_t *ordinal_end, MpiMarkdupStats *stats) {
    int active = 0;
    uint64_t ordinal = ordinal_base;
    if (!decomp || !uncompressed || !payload_blocks ||
        !output_blocks || !paras || !active_blocks ||
        !ordinal_end || payload_blocks->n < kMarkdupNB ||
        output_blocks->n < kMarkdupNB) {
        return -1;
    }

    MdResetPayloadBlock(payload_blocks->blocks + active, active);
    for (int b = 0; b < n_blocks; ++b) {
        bam_block *src = uncompressed->blocks + b;
        size_t pos = 0;
        for (int r = 0; r < decomp[b].n_total_records; ++r) {
            if (pos + 36u > src->length) return -11;
            const unsigned char *raw = src->data + pos;
            const uint32_t block_len = MdLoadLe32(raw);
            const uint64_t raw_len64 = (uint64_t)block_len + 4u;
            if (block_len < 32 ||
                raw_len64 > BGZF_BLOCK_SIZE ||
                raw_len64 > (uint64_t)src->length - pos) {
                return -12;
            }
            const uint32_t flag_nc = MdLoadLe32(raw + 16);
            const uint16_t old_flag = (uint16_t)(flag_nc >> 16);
            uint16_t final_flag =
                (uint16_t)(old_flag & (uint16_t)~BAM_FDUP);
            const size_t byte = (size_t)(ordinal >> 3);
            const uint8_t mask = (uint8_t)(1u << (ordinal & 7u));
            const int duplicate =
                byte < bitmap.size() && (bitmap[byte] & mask);

            if (old_flag & BAM_FDUP) {
                if (stats) stats->cleared_records++;
            }
            if (duplicate) {
                if (stats) {
                    stats->marked_records++;
                    stats->removed_records++;
                }
                ++ordinal;
                pos += (size_t)raw_len64;
                continue;
            }

            while (true) {
                int ret = MdAppendRawPayloadClearOld(
                    payload_blocks->blocks + active,
                    raw, block_len, final_flag);
                if (ret == 0) break;
                if (ret < 0) return ret;

                MdSetupRawCompress(paras + active, active,
                                   payload_blocks->blocks + active,
                                   output_blocks->blocks + active,
                                   compress_level);
                ++active;
                if (active >= kMarkdupNB) return -13;
                MdResetPayloadBlock(payload_blocks->blocks + active,
                                    active);
            }
            ++ordinal;
            pos += (size_t)raw_len64;
        }
        if (pos != src->length) return -14;
    }

    if (payload_blocks->blocks[active].pos > 0) {
        MdSetupRawCompress(paras + active, active,
                           payload_blocks->blocks + active,
                           output_blocks->blocks + active,
                           compress_level);
        ++active;
    }
    for (int i = active; i < kMarkdupNB; ++i) {
        MdInitEmptyRawCompress(paras + i, i);
    }
    *active_blocks = active;
    *ordinal_end = ordinal;
    return 0;
}

static int MdPrepareRewrittenPayloadFromRecords(
        Bam2BamPara *decomp, MdBlockSet *uncompressed,
        MpiMarkdupRewritePara *rewrite, int n_blocks,
        MdBlockSet *payload_blocks, MdBlockSet *output_blocks,
        MpiSortRawCompressPara *paras, int compress_level,
        int *active_blocks, MpiMarkdupStats *stats) {
    if (!decomp || !uncompressed || !rewrite ||
        !payload_blocks || !output_blocks || !paras ||
        !active_blocks || uncompressed->n < n_blocks ||
        payload_blocks->n < kMarkdupNB ||
        output_blocks->n < kMarkdupNB) {
        return -1;
    }
    int active = 0;
    MdResetPayloadBlock(payload_blocks->blocks + active, active);
    for (int b = 0; b < n_blocks; ++b) {
        bam_block *src = uncompressed->blocks + b;
        size_t pos = 0;
        if (rewrite[b].n_kept_records != decomp[b].n_total_records ||
            rewrite[b].n_records != decomp[b].n_total_records) {
            return -2;
        }
        for (int i = 0; i < decomp[b].n_total_records; ++i) {
            if (pos + 36u > src->length) return -11;
            const unsigned char *raw = src->data + pos;
            const uint32_t block_len = MdLoadLe32(raw);
            const uint64_t raw_len64 = (uint64_t)block_len + 4u;
            if (block_len < 32 ||
                raw_len64 > BGZF_BLOCK_SIZE ||
                raw_len64 > (uint64_t)src->length - pos) {
                return -12;
            }
            bam1_t *record = rewrite[b].records[i];
            if (!record) return -1;
            if (record->core.flag & BAM_FDUP) {
                if (stats) stats->removed_records++;
                pos += (size_t)raw_len64;
                continue;
            }

            while (true) {
                int ret = 0;
                if (MdRawRecordCanPatchFlagOnly(
                        raw, block_len, record,
                        rewrite[b].bam_lens[i])) {
                    ret = MdAppendRawPayloadPatchedFlag(
                        payload_blocks->blocks + active,
                        raw, (uint32_t)raw_len64,
                        record->core.flag);
                } else {
                    ret = MdAppendBamRecordPayload(
                        payload_blocks->blocks + active, record);
                }
                if (ret == 0) break;
                if (ret < 0) return ret;

                MdSetupRawCompress(paras + active, active,
                                   payload_blocks->blocks + active,
                                   output_blocks->blocks + active,
                                   compress_level);
                ++active;
                if (active >= kMarkdupNB) return -13;
                MdResetPayloadBlock(payload_blocks->blocks + active,
                                    active);
            }
            pos += (size_t)raw_len64;
        }
        if (pos != src->length) return -14;
    }
    if (payload_blocks->blocks[active].pos > 0) {
        MdSetupRawCompress(paras + active, active,
                           payload_blocks->blocks + active,
                           output_blocks->blocks + active,
                           compress_level);
        ++active;
    }
    for (int i = active; i < kMarkdupNB; ++i) {
        MdInitEmptyRawCompress(paras + i, i);
    }
    *active_blocks = active;
    return 0;
}

static void MdInitEmptyComp(Comp_Para *para, int block_id) {
    memset(para, 0, sizeof(*para));
    para->block_id = block_id;
    para->status = -1;
}

static int MdRewriteOutput(
        MemReader &reader, const std::vector<uint8_t> &bitmap,
        uint64_t expected_records, int clear_old,
        int remove_dups, int compress_level,
        MemWriter &writer, int rank,
        MpiMarkdupStats *stats) {
    const int records_per_block = (int)MPI_RECORDS_PER_BLOCK;
    MdBlockSet input = {};
    MdBlockSet uncompressed = {};
    MdBlockSet compress_un_a = {};
    MdBlockSet compress_un_b = {};
    MdBlockSet output_a = {};
    MdBlockSet output_b = {};
    MdRecordSet records = {};
    MdPackWorkspace pack;
    Bam2BamPara decomp[kMarkdupNB];
    MpiMarkdupRewritePara rewrite[kMarkdupNB];
    Comp_Para comp_a[kMarkdupNB];
    Comp_Para comp_b[kMarkdupNB];
    MpiSortRawCompressPara raw_comp[kMarkdupNB];

    if (MdAllocateBlockSet(&input, kMarkdupNB) != 0 ||
        MdAllocateBlockSet(&uncompressed, kMarkdupNB) != 0 ||
        MdAllocateBlockSet(&compress_un_a, kMarkdupNB) != 0 ||
        MdAllocateBlockSet(&compress_un_b, kMarkdupNB) != 0 ||
        MdAllocateBlockSet(&output_a, kMarkdupNB) != 0 ||
        MdAllocateBlockSet(&output_b, kMarkdupNB) != 0 ||
        MdAllocateRecordSet(&records, kMarkdupNB,
                            records_per_block,
                            MPI_BAM_BLOCK_ARENA_SIZE) != 0) {
        fprintf(stderr,
                "[rank %d] ERROR: failed to allocate markdup rewrite workspace.\n",
                rank);
        MdFreeBlockSet(&input);
        MdFreeBlockSet(&uncompressed);
        MdFreeBlockSet(&compress_un_a);
        MdFreeBlockSet(&compress_un_b);
        MdFreeBlockSet(&output_a);
        MdFreeBlockSet(&output_b);
        MdFreeRecordSet(&records);
        return -1;
    }
    MdInitPackWorkspace(
        &pack, (size_t)kMarkdupNB * records_per_block);
    for (int i = 0; i < kMarkdupNB; ++i) {
        MdInitEmptyComp(comp_a + i, i);
        MdInitEmptyComp(comp_b + i, i);
        MdInitEmptyRawCompress(raw_comp + i, i);
    }

    Comp_Para *comp_active = comp_a;
    Comp_Para *comp_pending = comp_b;
    MdBlockSet *compress_un_active = &compress_un_a;
    MdBlockSet *compress_un_pending = &compress_un_b;
    MdBlockSet *output_active = &output_a;
    MdBlockSet *output_pending = &output_b;
    bool has_pending = false;
    const int fused_remove_path =
        clear_old || MdEnvFlagEnabled("RABBITBAM_MARKDUP_FUSED_REMOVE");

    auto flush_pending = [&]() -> int {
        if (!has_pending) return 0;
        double t0 = GetTime();
        for (int i = 0; i < kMarkdupNB; ++i) {
            if (comp_pending[i].status == 0 &&
                comp_pending[i].output_block) {
                if (MpiWriteBlockToMem(
                        writer,
                        comp_pending[i].output_block) != 0) {
                    return -1;
                }
                stats->bgzf_blocks++;
            }
            MdInitEmptyComp(comp_pending + i, i);
        }
        has_pending = false;
        stats->t_write += GetTime() - t0;
        return 0;
    };

    auto read_group = [&](int *count) -> int {
        return MdReadGroup(reader, &input, count, stats);
    };

    auto compress_group = [&](int active_blocks,
                              int *next_blocks) -> int {
        for (int i = 0; i < active_blocks; ++i) {
            memset(comp_active + i, 0,
                   sizeof(comp_active[i]));
            comp_active[i].block_id = i;
            comp_active[i].input_records =
                pack.plans[i].records;
            comp_active[i].n_records =
                pack.plans[i].n_records;
            comp_active[i].un_comp_block =
                compress_un_active->blocks + i;
            comp_active[i].output_block =
                output_active->blocks + i;
            comp_active[i].compress_level =
                compress_level;
            comp_active[i].status = 0;
        }
        for (int i = active_blocks; i < kMarkdupNB; ++i) {
            MdInitEmptyComp(comp_active + i, i);
        }
        double t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_compressfunc,
                             comp_active, 1);
        if (flush_pending() != 0 ||
            read_group(next_blocks) != 0) {
            return -1;
        }
        double cpe_sync_t0 = GetTime();
        athread_join();
        stats->t_cpe_sync += GetTime() - cpe_sync_t0;
        stats->t_compress += GetTime() - t0;
        for (int i = 0; i < active_blocks; ++i) {
            if (comp_active[i].status != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: markdup compress failed "
                        "at output block %d status=%d.\n",
                        rank, i, comp_active[i].status);
                return -1;
            }
        }
        has_pending = active_blocks > 0;
        std::swap(comp_active, comp_pending);
        std::swap(compress_un_active,
                  compress_un_pending);
        std::swap(output_active, output_pending);
        MdResetPackWorkspace(&pack);
        return 0;
    };

    uint64_t ordinal = 0;
    int next_blocks = 0;
    read_group(&next_blocks);
    while (next_blocks > 0) {
        int n_blocks = next_blocks;
        next_blocks = 0;
        for (int b = 0; b < kMarkdupNB; ++b) {
            MdInitDecompPara(decomp + b, b, &input,
                             &uncompressed, &records,
                             b < n_blocks);
        }
        double decomp_t0 = GetTime();
        __real_athread_spawn(
            (void *)slave_mpi_decompress_bam2bam_passthrough,
            decomp, 1);
        double cpe_sync_t0 = GetTime();
        athread_join();
        stats->t_cpe_sync += GetTime() - cpe_sync_t0;
        stats->t_rewrite_decomp += GetTime() - decomp_t0;
        for (int b = 0; b < n_blocks; ++b) {
            if (decomp[b].status != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: markdup rewrite decode "
                        "failed at local ordinal %llu status=%d.\n",
                        rank, (unsigned long long)ordinal,
                        decomp[b].status);
                return -1;
            }
        }

        // 普通 -r 默认走已验证的 raw payload 删除路径；需要测试
        // host-filter 融合路径时设置 RABBITBAM_MARKDUP_FUSED_REMOVE=1。
        if (remove_dups && !clear_old && !fused_remove_path) {
            int active_payload_blocks = 0;
            uint64_t next_ordinal = ordinal;
            double payload_t0 = GetTime();
            int payload_status = MdPrepareRemoveDupsPayloadFromRaw(
                decomp, &uncompressed, n_blocks, bitmap,
                ordinal, compress_un_active, output_active,
                raw_comp, compress_level, &active_payload_blocks,
                &next_ordinal, stats);
            stats->t_pack += GetTime() - payload_t0;
            if (payload_status != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: markdup remove-dups raw "
                        "payload pack failed status=%d ordinal=%llu.\n",
                        rank, payload_status,
                        (unsigned long long)next_ordinal);
                return -1;
            }
            ordinal = next_ordinal;
            if (active_payload_blocks == 0) {
                if (read_group(&next_blocks) != 0) return -1;
                continue;
            }

            double compress_t0 = GetTime();
            __real_athread_spawn(
                (void *)slave_mpi_sort_compress_payload,
                raw_comp, 1);
            cpe_sync_t0 = GetTime();
            athread_join();
            stats->t_cpe_sync += GetTime() - cpe_sync_t0;
            stats->t_compress += GetTime() - compress_t0;
            double write_t0 = GetTime();
            for (int i = 0; i < active_payload_blocks; ++i) {
                if (raw_comp[i].status != 0) {
                    fprintf(stderr,
                            "[rank %d] ERROR: markdup remove-dups "
                            "raw payload compress failed block=%d status=%d.\n",
                            rank, i, raw_comp[i].status);
                    return -1;
                }
                if (MpiWriteBlockToMem(writer,
                        output_active->blocks + i) != 0) {
                    return -1;
                }
                stats->bgzf_blocks++;
            }
            stats->t_write += GetTime() - write_t0;
            if (read_group(&next_blocks) != 0) return -1;
            continue;
        }

        if (remove_dups && clear_old && fused_remove_path) {
            int active_payload_blocks = 0;
            uint64_t next_ordinal = ordinal;
            double payload_t0 = GetTime();
            int payload_status = MdPrepareClearRemoveDupsPayloadFromRaw(
                decomp, &uncompressed, n_blocks, bitmap,
                ordinal, compress_un_active, output_active,
                raw_comp, compress_level, &active_payload_blocks,
                &next_ordinal, stats);
            stats->t_pack += GetTime() - payload_t0;
            if (payload_status != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: markdup clear/remove raw "
                        "payload pack failed status=%d ordinal=%llu.\n",
                        rank, payload_status,
                        (unsigned long long)next_ordinal);
                return -1;
            }
            ordinal = next_ordinal;
            if (active_payload_blocks == 0) {
                if (read_group(&next_blocks) != 0) return -1;
                continue;
            }

            double compress_t0 = GetTime();
            __real_athread_spawn(
                (void *)slave_mpi_sort_compress_payload,
                raw_comp, 1);
            cpe_sync_t0 = GetTime();
            athread_join();
            stats->t_cpe_sync += GetTime() - cpe_sync_t0;
            stats->t_compress += GetTime() - compress_t0;
            double write_t0 = GetTime();
            for (int i = 0; i < active_payload_blocks; ++i) {
                if (raw_comp[i].status != 0) {
                    fprintf(stderr,
                            "[rank %d] ERROR: markdup clear/remove "
                            "raw payload compress failed block=%d status=%d.\n",
                            rank, i, raw_comp[i].status);
                    return -1;
                }
                if (MpiWriteBlockToMem(writer,
                        output_active->blocks + i) != 0) {
                    return -1;
                }
                stats->bgzf_blocks++;
            }
            stats->t_write += GetTime() - write_t0;
            if (read_group(&next_blocks) != 0) return -1;
            continue;
        }

        // 普通路径：修改 FLAG / 清旧标记 / 重新打包
        uint64_t block_ordinal = ordinal;
        for (int b = 0; b < kMarkdupNB; ++b) {
            memset(rewrite + b, 0, sizeof(rewrite[b]));
            rewrite[b].block_id = b;
            rewrite[b].status = b < n_blocks ? 0 : -1;
            if (b >= n_blocks) continue;
            rewrite[b].records =
                decomp[b].output_records;
            rewrite[b].bam_lens =
                decomp[b].bam_lens;
            rewrite[b].n_records =
                decomp[b].n_total_records;
            rewrite[b].ordinal_base = block_ordinal;
            rewrite[b].duplicate_bitmap =
                bitmap.empty() ? nullptr : bitmap.data();
            rewrite[b].duplicate_bitmap_bytes =
                bitmap.size();
            rewrite[b].clear_old = clear_old;
            // 删除由 MPE 在 pack 前过滤，CPE 只负责清旧/标记。
            rewrite[b].remove_dups = 0;
            block_ordinal +=
                (uint64_t)decomp[b].n_total_records;
        }
        double rewrite_t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_markdup_rewrite,
                             rewrite, 1);
        cpe_sync_t0 = GetTime();
        athread_join();
        stats->t_cpe_sync += GetTime() - cpe_sync_t0;
        stats->t_rewrite += GetTime() - rewrite_t0;
        ordinal = block_ordinal;

        for (int b = 0; b < n_blocks; ++b) {
            if (rewrite[b].status != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: markdup rewrite failed "
                        "at block %d record=%d status=%d.\n",
                        rank, b, rewrite[b].record_index,
                        rewrite[b].status);
                return -1;
            }
            stats->marked_records +=
                rewrite[b].marked_records;
            stats->cleared_records +=
                rewrite[b].cleared_records;
            stats->removed_records +=
                rewrite[b].removed_records;
        }

        if (remove_dups && fused_remove_path) {
            int active_payload_blocks = 0;
            double payload_t0 = GetTime();
            int payload_status = MdPrepareRewrittenPayloadFromRecords(
                decomp, &uncompressed, rewrite, n_blocks,
                compress_un_active, output_active, raw_comp,
                compress_level, &active_payload_blocks, stats);
            stats->t_pack += GetTime() - payload_t0;
            if (payload_status != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: markdup rewritten payload "
                        "pack failed status=%d.\n",
                        rank, payload_status);
                return -1;
            }
            if (active_payload_blocks == 0) {
                if (read_group(&next_blocks) != 0) return -1;
                continue;
            }

            double compress_t0 = GetTime();
            __real_athread_spawn(
                (void *)slave_mpi_sort_compress_payload,
                raw_comp, 1);
            cpe_sync_t0 = GetTime();
            athread_join();
            stats->t_cpe_sync += GetTime() - cpe_sync_t0;
            stats->t_compress += GetTime() - compress_t0;
            double write_t0 = GetTime();
            for (int i = 0; i < active_payload_blocks; ++i) {
                if (raw_comp[i].status != 0) {
                    fprintf(stderr,
                            "[rank %d] ERROR: markdup rewritten "
                            "payload compress failed block=%d status=%d.\n",
                            rank, i, raw_comp[i].status);
                    return -1;
                }
                if (MpiWriteBlockToMem(writer,
                        output_active->blocks + i) != 0) {
                    return -1;
                }
                stats->bgzf_blocks++;
            }
            stats->t_write += GetTime() - write_t0;
            if (read_group(&next_blocks) != 0) return -1;
            continue;
        }

        // MdAppendRecords：把保留 records 重新打包进行压缩
        MdResetPackWorkspace(&pack);
        double pack_t0 = GetTime();
        int pack_status = 0;
        for (int b = 0; b < n_blocks; ++b) {
            pack_status = MdAppendRecords(
                &pack, rewrite[b].records,
                rewrite[b].bam_lens,
                rewrite[b].n_kept_records,
                rewrite[b].kept_total_len);
            if (pack_status != 0) break;
        }
        if (pack_status == 0) {
            pack_status = MdSealPackBlock(&pack);
        }
        stats->t_pack += GetTime() - pack_t0;
        if (pack_status != 0) {
            fprintf(stderr,
                    "[rank %d] ERROR: markdup output pack failed "
                    "status=%d.\n",
                    rank, pack_status);
            return -1;
        }
        if (pack.active_blocks == 0) {
            if (flush_pending() != 0 ||
                read_group(&next_blocks) != 0) {
                return -1;
            }
            continue;
        }
        if (compress_group(pack.active_blocks,
                           &next_blocks) != 0) {
            return -1;
        }
    }
    if (flush_pending() != 0) return -1;
    if (ordinal != expected_records) {
        fprintf(stderr,
                "[rank %d] ERROR: markdup pass record count changed: "
                "candidate=%llu rewrite=%llu.\n",
                rank, (unsigned long long)expected_records,
                (unsigned long long)ordinal);
        return -1;
    }

    MdFreeBlockSet(&input);
    MdFreeBlockSet(&uncompressed);
    MdFreeBlockSet(&compress_un_a);
    MdFreeBlockSet(&compress_un_b);
    MdFreeBlockSet(&output_a);
    MdFreeBlockSet(&output_b);
    MdFreeRecordSet(&records);
    return 0;
}

static int MdValidateRankBoundaries(
        int local_has_records, int first_tid, int first_pos,
        int last_tid, int last_pos, int rank, int comm_size) {
    long long local[5] = {
        local_has_records, first_tid, first_pos, last_tid, last_pos
    };
    std::vector<long long> all((size_t)comm_size * 5);
    MPI_Allgather(local, 5, MPI_LONG_LONG,
                  all.data(), 5, MPI_LONG_LONG,
                  MPI_COMM_WORLD);
    int previous_tid = -1;
    int previous_pos = -1;
    int have_previous = 0;
    for (int r = 0; r < comm_size; ++r) {
        const long long *entry = all.data() + (size_t)r * 5;
        if (!entry[0]) continue;
        int current_first_tid = (int)entry[1];
        int current_first_pos = (int)entry[2];
        if (current_first_tid >= 0 && have_previous &&
            (current_first_tid < previous_tid ||
             (current_first_tid == previous_tid &&
              current_first_pos < previous_pos))) {
            if (rank == 0) {
                fprintf(stderr,
                        "ERROR: markdup input is not coordinate sorted "
                        "across rank boundary before rank %d.\n", r);
            }
            return -1;
        }
        previous_tid = (int)entry[3];
        previous_pos = (int)entry[4];
        have_previous = 1;
    }
    return 0;
}

struct MdOutputMemory {
    char *data;
    size_t size;
    char *header_data;
    size_t header_size;
    long long total_body_size;
    std::vector<long long> body_sizes;
    std::vector<long long> body_offsets;
};

struct MdOutputStageCosts {
    double stage44;
    double stage45_header;
    double stage45_malloc;
    double stage46;
};

static int MdPrepareOutputMemory(
        const MemWriter &local_writer, sam_hdr_t *header,
        int compress_level, int rank, int comm_size,
        MpiMarkdupStats *stats, MdOutputMemory *output,
        MdOutputStageCosts *costs) {
    memset(costs, 0, sizeof(*costs));
    output->data = nullptr;
    output->size = 0;
    output->header_data = nullptr;
    output->header_size = 0;
    output->total_body_size = 0;
    output->body_sizes.clear();
    output->body_offsets.clear();

    // 4.4 统计各 rank 输出大小 / gather layout
    double stage44_t0 = GetTime();
    long long local_size = (long long)local_writer.size;
    std::vector<long long> sizes((size_t)comm_size, 0);
    MPI_Allgather(&local_size, 1, MPI_LONG_LONG,
                  sizes.data(), 1, MPI_LONG_LONG,
                  MPI_COMM_WORLD);
    std::vector<long long> offsets((size_t)comm_size, 0);
    long long total_body = 0;
    for (int i = 0; i < comm_size; ++i) {
        offsets[(size_t)i] = total_body;
        if (sizes[(size_t)i] < 0 ||
            total_body > LLONG_MAX - sizes[(size_t)i]) {
            return -1;
        }
        total_body += sizes[(size_t)i];
    }
    output->body_sizes = sizes;
    output->body_offsets = offsets;
    output->total_body_size = total_body;
    costs->stage44 = MdReduceMax(GetTime() - stage44_t0);
    if (rank == 0) {
        printf("Complete the 4.4 gather body sizes cost %lf\n",
               costs->stage44);
    }

    char *header_memory = nullptr;
    size_t header_size = 0;
    int local_ok = 1;
    double stage45_header_t0 = 0.0;
    double stage45_malloc_t0 = 0.0;
    double sim_write_t0 = 0.0;
    char *simulated_write_mem = nullptr;
    size_t simulated_write_size = 0;
    size_t simulated_pos = 0;
    volatile unsigned long long simulated_write_guard = 0;

    // 4.5a header/layout
    stage45_header_t0 = GetTime();
    if (rank == 0) {
        if (sam_hdr_add_pg(header, "RabbitBAM-MPI",
                           "PN", "RabbitBAM-MPI",
                           "CL", "markdup",
                           NULL) != 0) {
            fprintf(stderr,
                    "ERROR: failed to add markdup @PG header line.\n");
            local_ok = 0;
        }
        if (local_ok &&
            MpiCommonBuildBamHeaderMemory(
                header, compress_level,
                &header_memory, &header_size) != 0) {
            fprintf(stderr,
                    "ERROR: failed to build markdup BAM header.\n");
            local_ok = 0;
        }
        unsigned long long final_size =
            (unsigned long long)header_size +
            (unsigned long long)total_body +
            sizeof(kMarkdupBgzfEofBlock);
        if (local_ok &&
            (final_size > (unsigned long long)SIZE_MAX)) {
            fprintf(stderr,
                    "ERROR: markdup output exceeds addressable memory.\n");
            local_ok = 0;
        }
        if (local_ok) {
            output->size = (size_t)final_size;
            output->header_data = header_memory;
            output->header_size = header_size;
            header_memory = nullptr;
        }
    }
    local_ok = MdAllRanksOk(local_ok);
    costs->stage45_header =
        MdReduceMax(GetTime() - stage45_header_t0);
    if (rank == 0 && local_ok) {
        printf("Complete the 4.5a header/layout cost %lf\n",
               costs->stage45_header);
    }
    if (!local_ok) goto fail;

    // 4.5b 分配各 rank 的本地模拟写内存，这部分不计入核心处理时间
    stage45_malloc_t0 = GetTime();
    simulated_write_size = local_writer.size;
    if (rank == 0) {
        if (output->header_size > SIZE_MAX - simulated_write_size) {
            fprintf(stderr,
                    "[rank %d] ERROR: simulated markdup output memory "
                    "estimate overflow.\n", rank);
            local_ok = 0;
        } else {
            simulated_write_size += output->header_size;
        }
    }
    if (local_ok && simulated_write_size > 0) {
        simulated_write_mem = (char *)malloc(simulated_write_size);
        if (!simulated_write_mem) {
            fprintf(stderr,
                    "[rank %d] ERROR: failed to allocate markdup "
                    "simulated output memory. size=%zu\n",
                    rank, simulated_write_size);
            local_ok = 0;
        }
    }
    local_ok = MdAllRanksOk(local_ok);
    costs->stage45_malloc =
        MdReduceMax(GetTime() - stage45_malloc_t0);
    if (rank == 0 && local_ok) {
        printf("Complete the 4.5b malloc simulated write memory cost %lf\n",
               costs->stage45_malloc);
    }
    if (!local_ok) goto fail;

    // 4.6 各 rank 根据逻辑偏移模拟写入高速磁盘内存，只统计 memcpy 时间
    sim_write_t0 = GetTime();
    simulated_pos = 0;
    if (rank == 0 && output->header_size > 0) {
        memcpy(simulated_write_mem + simulated_pos,
               output->header_data, output->header_size);
        simulated_pos += output->header_size;
    }
    if (local_writer.size > 0) {
        memcpy(simulated_write_mem + simulated_pos,
               local_writer.data, local_writer.size);
        simulated_pos += local_writer.size;
    }
    if (simulated_pos > 0) {
        unsigned char *guard_ptr =
            (unsigned char *)simulated_write_mem;
        simulated_write_guard += guard_ptr[0];
        simulated_write_guard += guard_ptr[simulated_pos - 1];
    }
    {
        double sim_write_cost = GetTime() - sim_write_t0;
        costs->stage46 = MdReduceMax(sim_write_cost);
    }
    if (rank == 0 && local_ok) {
        printf("Complete the 4.6 simulated header/body distributed write cost %lf\n",
               costs->stage46);
    }
    if (simulated_write_guard == (unsigned long long)-1 && rank < 0) {
        fprintf(stderr, "unused markdup simulated write guard %llu\n",
                simulated_write_guard);
    }
    if (simulated_write_mem) free(simulated_write_mem);
    if (header_memory) free(header_memory);
    return MdAllRanksOk(local_ok) ? 0 : -1;

fail:
    if (simulated_write_mem) free(simulated_write_mem);
    if (header_memory) free(header_memory);
    if (output->data) {
        free(output->data);
        output->data = nullptr;
        output->size = 0;
    }
    if (output->header_data) {
        free(output->header_data);
        output->header_data = nullptr;
        output->header_size = 0;
    }
    return -1;
}

static double MdAccountedFusedTime(const MpiMarkdupStats &stats) {
    return stats.t_candidate_decomp +
           stats.t_candidate_extract +
           stats.t_extract_status +
           stats.t_extract_merge +
           stats.t_extract_memcheck +
           stats.t_candidate_exchange +
           stats.t_boundary_check +
           stats.t_candidate_cleanup +
           stats.t_group_prepare +
           stats.t_group +
           stats.t_owner_cleanup +
           stats.t_result_exchange +
           stats.t_result_cleanup +
           stats.t_rewrite_decomp +
           stats.t_rewrite +
           stats.t_pack +
           stats.t_compress +
           stats.t_rank_sync;
}

static void MdPrintStats(const MpiMarkdupStats &stats,
                         int rank, int comm_size) {
    const int count = 16;
    long long local_counts[count] = {
        stats.input_blocks,
        stats.group_count,
        stats.total_records,
        stats.examined_records,
        stats.excluded_records,
        stats.pair_candidates,
        stats.single_candidates,
        stats.owner_candidates,
        stats.pair_duplicates,
        stats.single_duplicates,
        stats.marked_records,
        stats.cleared_records,
        stats.removed_records,
        stats.qname_bytes,
        stats.mpi_candidate_bytes,
        stats.mpi_result_bytes
    };
    const double accounted = MdAccountedFusedTime(stats);
    const double unaccounted = stats.t_fused_total - accounted;

    const int kLineBytes = 8192;
    char local_lines[kLineBytes];
    memset(local_lines, 0, sizeof(local_lines));
    int used = snprintf(
        local_lines, sizeof(local_lines),
        "[rank %d] markdup_stats blocks=%lld groups=%lld "
        "records=%lld examined=%lld excluded=%lld "
        "pair_candidates=%lld single_candidates=%lld "
        "owner_candidates=%lld qname_bytes=%lld "
        "pair_duplicates=%lld single_duplicates=%lld "
        "marked=%lld cleared=%lld removed=%lld "
        "mpi_candidate_bytes=%lld result_bytes=%lld "
        "bgzf=%lld tracked_peak=%lld\n",
        rank, stats.input_blocks, stats.group_count,
        stats.total_records, stats.examined_records,
        stats.excluded_records, stats.pair_candidates,
        stats.single_candidates, stats.owner_candidates,
        stats.qname_bytes, stats.pair_duplicates,
        stats.single_duplicates, stats.marked_records,
        stats.cleared_records, stats.removed_records,
        stats.mpi_candidate_bytes, stats.mpi_result_bytes,
        stats.bgzf_blocks, stats.tracked_peak_bytes);
    if (used < 0) used = 0;
    if (used >= kLineBytes) used = kLineBytes - 1;
    snprintf(
        local_lines + used, (size_t)(kLineBytes - used),
        "[rank %d] markdup_timing candidate_decomp=%.6f "
        "candidate_extract=%.6f extract_status=%.6f "
        "extract_merge=%.6f extract_memcheck=%.6f "
        "candidate_pack=%.6f candidate_mpi=%.6f "
        "candidate_exchange=%.6f boundary=%.6f "
        "group_sort=%.6f group_scan=%.6f group=%.6f "
        "flat_init=%.6f flat_probe=%.6f result=%.6f "
        "rewrite_decomp=%.6f "
        "rewrite=%.6f pack=%.6f compress=%.6f "
        "read=%.6f write=%.6f rank_sync=%.6f "
        "cpe_sync=%.6f accounted=%.6f "
        "unaccounted=%.6f fused=%.6f\n",
        rank, stats.t_candidate_decomp,
        stats.t_candidate_extract, stats.t_extract_status,
        stats.t_extract_merge, stats.t_extract_memcheck,
        stats.t_candidate_pack,
        stats.t_candidate_exchange_mpi,
        stats.t_candidate_exchange,
        stats.t_boundary_check, stats.t_group_sort,
        stats.t_group_scan, stats.t_group,
        stats.t_flat_init, stats.t_flat_probe,
        stats.t_result_exchange,
        stats.t_rewrite_decomp, stats.t_rewrite, stats.t_pack,
        stats.t_compress, stats.t_read, stats.t_write,
        stats.t_rank_sync, stats.t_cpe_sync, accounted,
        unaccounted, stats.t_fused_total);
    std::vector<char> gathered_lines;
    if (rank == 0) {
        gathered_lines.resize((size_t)comm_size * kLineBytes);
    }
    MPI_Gather(local_lines, kLineBytes, MPI_CHAR,
               rank == 0 ? gathered_lines.data() : nullptr,
               kLineBytes, MPI_CHAR, 0, MPI_COMM_WORLD);

    long long sums[count] = {};
    MPI_Reduce(local_counts, sums, count, MPI_LONG_LONG,
               MPI_SUM, 0, MPI_COMM_WORLD);
    long long peak = stats.tracked_peak_bytes;
    long long peak_max = 0;
    MPI_Reduce(&peak, &peak_max, 1, MPI_LONG_LONG,
               MPI_MAX, 0, MPI_COMM_WORLD);

    const int time_count = 26;
    double local_times[time_count] = {
        stats.t_candidate_decomp,
        stats.t_candidate_extract,
        stats.t_extract_status,
        stats.t_extract_merge,
        stats.t_extract_memcheck,
        stats.t_candidate_pack,
        stats.t_candidate_exchange_mpi,
        stats.t_candidate_exchange,
        stats.t_boundary_check,
        stats.t_group_sort,
        stats.t_group_scan,
        stats.t_group,
        stats.t_flat_init,
        stats.t_flat_probe,
        stats.t_result_exchange,
        stats.t_rewrite_decomp,
        stats.t_rewrite,
        stats.t_pack,
        stats.t_compress,
        stats.t_read,
        stats.t_write,
        stats.t_rank_sync,
        stats.t_cpe_sync,
        accounted,
        unaccounted,
        stats.t_fused_total
    };
    double time_sums[time_count] = {};
    double time_max[time_count] = {};
    MPI_Reduce(local_times, time_sums, time_count, MPI_DOUBLE,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_times, time_max, time_count, MPI_DOUBLE,
               MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int r = 0; r < comm_size; ++r) {
            fputs(gathered_lines.data() + (size_t)r * kLineBytes,
                  stdout);
        }
        printf("FusedMarkdupMPI finished. ranks=%d blocks=%lld "
               "records=%lld examined=%lld excluded=%lld\n",
               comm_size, sums[0], sums[2], sums[3], sums[4]);
        printf("  candidates pair=%lld single=%lld owner=%lld "
               "qname_bytes=%lld\n",
               sums[5], sums[6], sums[7], sums[13]);
        printf("  duplicates pair=%lld single=%lld marked=%lld "
               "cleared=%lld removed=%lld\n",
               sums[8], sums[9], sums[10], sums[11], sums[12]);
        printf("  mpi candidate_bytes=%lld result_bytes=%lld "
               "tracked_peak_max=%lld\n",
               sums[14], sums[15], peak_max);
        printf("  timing_sum candidate_decomp=%.3f "
               "candidate_extract=%.3f extract_status=%.3f "
               "extract_merge=%.3f extract_memcheck=%.3f "
               "candidate_pack=%.3f candidate_mpi=%.3f "
               "candidate_exchange=%.3f boundary=%.3f "
               "group_sort=%.3f group_scan=%.3f group=%.3f "
               "flat_init=%.3f flat_probe=%.3f result=%.3f "
               "rewrite_decomp=%.3f "
               "rewrite=%.3f pack=%.3f compress=%.3f "
               "read=%.3f write=%.3f rank_sync=%.3f "
               "cpe_sync=%.3f accounted=%.3f "
               "unaccounted=%.3f fused=%.3f\n",
               time_sums[0], time_sums[1], time_sums[2],
               time_sums[3], time_sums[4], time_sums[5],
               time_sums[6], time_sums[7], time_sums[8],
               time_sums[9], time_sums[10], time_sums[11],
               time_sums[12], time_sums[13], time_sums[14],
               time_sums[15], time_sums[16], time_sums[17],
               time_sums[18], time_sums[19], time_sums[20],
               time_sums[21], time_sums[22], time_sums[23],
               time_sums[24], time_sums[25]);
        printf("  timing_max candidate_decomp=%.3f "
               "candidate_extract=%.3f extract_status=%.3f "
               "extract_merge=%.3f extract_memcheck=%.3f "
               "candidate_pack=%.3f candidate_mpi=%.3f "
               "candidate_exchange=%.3f boundary=%.3f "
               "group_sort=%.3f group_scan=%.3f group=%.3f "
               "flat_init=%.3f flat_probe=%.3f result=%.3f "
               "rewrite_decomp=%.3f "
               "rewrite=%.3f pack=%.3f compress=%.3f "
               "read=%.3f write=%.3f rank_sync=%.3f "
               "cpe_sync=%.3f accounted=%.3f "
               "unaccounted=%.3f fused=%.3f\n",
               time_max[0], time_max[1], time_max[2],
               time_max[3], time_max[4], time_max[5],
               time_max[6], time_max[7], time_max[8],
               time_max[9], time_max[10], time_max[11],
               time_max[12], time_max[13], time_max[14],
               time_max[15], time_max[16], time_max[17],
               time_max[18], time_max[19], time_max[20],
               time_max[21], time_max[22], time_max[23],
               time_max[24], time_max[25]);
        fflush(stdout);
    }
}

} // namespace

int ProcessMarkdupMPI(CmdInfo *cmd_info) {

    // 1. 初始化阶段：拿 rank、解析内存参数
    double total_t0 = GetTime();
    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int exit_code = 1;
    int local_ok = 1;
    char *input_memory = nullptr;
    size_t input_size = 0;
    hFILE *input_hfile = nullptr;
    samFile *input = nullptr;
    sam_hdr_t *header = nullptr;
    long long body_start = 0;
    std::vector<long long> block_offsets;
    std::vector<long long> block_lengths;
    long long n_blocks = 0;
    long long local_block_begin = 0;
    long long local_block_end = 0;
    char *rank_input = nullptr;
    size_t rank_input_size = 0;
    MemReader reader = {};
    MemWriter writer = {};
    size_t memory_limit = 0;
    MpiMarkdupStats stats = {};
    MdOutputMemory output_memory = {};
    MdOutputStageCosts output_stage_costs = {};
    double stage41_cost = 0.0;
    double stage42_cost = 0.0;
    double stage43_cost = 0.0;
    double body_total_t0 = 0.0;

    std::vector<MpiMarkdupCandidateShared> local_candidates;
    std::vector<unsigned char> local_qnames;
    std::vector<MpiMarkdupCandidateShared> owner_candidates;
    std::vector<unsigned char> owner_qnames;
    std::vector<std::vector<uint64_t> > duplicates_by_source;
    std::vector<uint8_t> duplicate_bitmap;
    uint64_t local_records = 0;
    int range_has_records = 0;
    int range_first_tid = -1;
    int range_first_pos = -1;
    int range_last_tid = -1;
    int range_last_pos = -1;

    if (MdParseMemory(cmd_info->markdup_memory_,
                      &memory_limit) != 0) {
        if (rank == 0) {
            fprintf(stderr,
                    "ERROR: invalid markdup memory limit '%s'.\n",
                    cmd_info->markdup_memory_.c_str());
        }
        local_ok = 0;
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;
    if (rank == 0) {
        printf("111Complete the initialization cost %lf-----\n",
               MdReduceMax(GetTime() - total_t0));
    } else {
        MdReduceMax(GetTime() - total_t0);
    }


    // 2. 222Complete the memory cost：每个 rank 预读整个 BAM
    {
        double t0 = GetTime();
        if (MpiCommonLoadFileToMemory(
                cmd_info->in_file_name_,
                &input_memory, &input_size) != 0) {
            fprintf(stderr,
                    "[rank %d] ERROR: cannot preload markdup input %s.\n",
                    rank, cmd_info->in_file_name_.c_str());
            local_ok = 0;
        }
        double cost = MdReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("222Complete the memory cost %lf--\n", cost);
        }
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;


    // 3. 333Complete the head cost：从内存打开 BAM，读 header，检查格式
    {
        double t0 = GetTime();
        input_hfile = hopen("mem:", "rb:",
                            input_memory, input_size);
        if (input_hfile) {
            input = (samFile *)hts_hopen(
                input_hfile, "data", "rb");
            if (input) input_hfile = nullptr;
        }
        if (!input) {
            fprintf(stderr,
                    "[rank %d] ERROR: cannot open preloaded markdup BAM.\n",
                    rank);
            local_ok = 0;
        }
        if (local_ok) header = sam_hdr_read(input);
        if (!header) local_ok = 0;
        if (local_ok &&
            input->format.format != bam &&
            input->format.format != binary_format) {
            if (rank == 0) {
                fprintf(stderr,
                        "ERROR: RabbitBAM-MPI markdup only supports BAM input.\n");
            }
            local_ok = 0;
        }
        if (local_ok) {
            kstring_t sort_order = {0, 0, nullptr};
            if (sam_hdr_find_tag_hd(header, "SO",
                                    &sort_order) == 0 &&
                sort_order.s &&
                strcmp(sort_order.s, "queryname") == 0) {
                if (rank == 0) {
                    fprintf(stderr,
                            "ERROR: markdup input is queryname sorted; "
                            "coordinate sort it first.\n");
                }
                local_ok = 0;
            }
            free(sort_order.s);
        }
        if (local_ok) {
            body_start =
                (long long)input->fp.bgzf->block_address;
            if (body_start < 0 ||
                (unsigned long long)body_start >
                    (unsigned long long)input_size) {
                local_ok = 0;
            }
        }
        double cost = MdReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("333Complete the head cost %lf---\n", cost);
        }
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;

    body_total_t0 = GetTime();
    {
        // 4.1：扫描 BGZF block，并按 block 分给各 rank
        double t0 = GetTime();
        if (rank == 0) {
            printf("Enable MPI BAM MARKDUP mode (%d MPE + %d CPEs)!!!\n",
                   comm_size, comm_size * 64);
            printf("MPI BAM output compression level=%d\n",
                   cmd_info->compress_level_);
            if (MpiCommonScanBgzfBlocksInMemory(
                    input_memory, input_size, body_start,
                    &block_offsets, &block_lengths) != 0) {
                fprintf(stderr,
                        "ERROR: failed to scan markdup input BGZF blocks.\n");
                local_ok = 0;
            }
            n_blocks = (long long)block_offsets.size();
            if (n_blocks > INT_MAX) local_ok = 0;
        }
        MPI_Bcast(&local_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (!local_ok) goto cleanup;
        MPI_Bcast(&n_blocks, 1, MPI_LONG_LONG, 0,
                  MPI_COMM_WORLD);
        if (rank != 0) {
            block_offsets.resize((size_t)n_blocks);
            block_lengths.resize((size_t)n_blocks);
        }
        if (n_blocks > 0) {
            MPI_Bcast(block_offsets.data(), (int)n_blocks,
                      MPI_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(block_lengths.data(), (int)n_blocks,
                      MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        }
        local_block_begin = n_blocks * rank / comm_size;
        local_block_end =
            n_blocks * (rank + 1) / comm_size;
        if (MpiCommonSelectBlockRangeFromMemory(
                input_memory, input_size,
                block_offsets, block_lengths,
                local_block_begin, local_block_end,
                &rank_input, &rank_input_size) != 0) {
            local_ok = 0;
        }
        stage41_cost = MdReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("MPI BAM scan complete. data_blocks=%lld "
                   "body_start=%lld\n",
                   n_blocks, body_start);
            printf("Complete the 4.1 scan/split/broadcast/select cost %lf\n",
                   stage41_cost);
        }
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;

    // 4.2：初始化本 rank 的内存 reader / writer
    {
        double t0 = GetTime();
        reader.base = rank_input;
        reader.size = rank_input_size;
        reader.pos = 0;
        if (MpiCommonInitMemWriter(
                writer, rank_input_size ?
                            rank_input_size :
                            64 * 1024 * 1024) != 0) {
            local_ok = 0;
        }
        stage42_cost = MdReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.2 init reader/writer cost %lf\n",
                   stage42_cost);
        }
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;

    // 4.3 核心：提取候选、交换、找 duplicate、重写 BAM
    {
        double fused_t0 = GetTime();
        // 1. MdExtractCandidates：提取本 rank 的候选记录
        if (MdExtractCandidates(
                reader, local_block_begin, rank,
                cmd_info->markdup_include_fails_,
                memory_limit, &local_candidates,
                &local_qnames, &local_records,
                &range_has_records,
                &range_first_tid, &range_first_pos,
                &range_last_tid, &range_last_pos,
                &stats) != 0) {
            local_ok = 0;
        }
        if (!MdAllRanksOkTimed(local_ok, &stats)) goto cleanup;

        // 2. MdValidateRankBoundaries：检查 rank 边界
        {
            double boundary_t0 = GetTime();
            if (MdValidateRankBoundaries(
                range_has_records, range_first_tid, range_first_pos,
                range_last_tid, range_last_pos, rank,
                comm_size) != 0) {
                local_ok = 0;
            }
            stats.t_boundary_check += GetTime() - boundary_t0;
        }
        if (!MdAllRanksOkTimed(local_ok, &stats)) goto cleanup;

        // 3. MdExchangeCandidatesByCoordinate：
        // coordinate-sorted markdup fast path.  Most candidates are owned by
        // the local rank; only boundary/unclipped candidates are sent.
        // RABBITBAM_MARKDUP_HASH_OWNER=1 keeps the old hash-owner path for
        // quick correctness bisecting.
        const char *force_hash_owner =
            getenv("RABBITBAM_MARKDUP_HASH_OWNER");
        if (force_hash_owner && force_hash_owner[0] == '1') {
            if (MdExchangeCandidates(
                    local_candidates, local_qnames,
                    rank, comm_size, memory_limit,
                    &owner_candidates, &owner_qnames,
                    &stats) != 0) {
                local_ok = 0;
            }
        } else if (MdExchangeCandidatesByCoordinate(
                    &local_candidates, &local_qnames,
                    range_has_records, range_first_tid, range_first_pos,
                    range_last_tid, range_last_pos,
                    rank, comm_size, memory_limit,
                    &owner_candidates, &owner_qnames,
                    &stats) != 0) {
                local_ok = 0;
        }
        if (!MdAllRanksOkTimed(local_ok, &stats)) goto cleanup;

        {
            double cleanup_t0 = GetTime();
            local_candidates.clear();
            local_candidates.shrink_to_fit();
            local_qnames.clear();
            local_qnames.shrink_to_fit();
            stats.t_candidate_cleanup += GetTime() - cleanup_t0;
        }

        const char *use_coord_bucket =
            getenv("RABBITBAM_MARKDUP_COORD_BUCKET");
        const char *use_full_sort =
            getenv("RABBITBAM_MARKDUP_FULL_SORT");
        const char *use_global_sort =
            getenv("RABBITBAM_MARKDUP_GLOBAL_SORT");
        const int enable_coord_bucket =
            use_coord_bucket && use_coord_bucket[0] == '1';
        const int enable_full_sort =
            (use_full_sort && use_full_sort[0] == '1') ||
            (use_global_sort && use_global_sort[0] == '1');
        const int enable_flat_hash =
            !enable_full_sort && !enable_coord_bucket;

        {
            double group_prepare_t0 = GetTime();
            size_t grouping_bytes =
                owner_candidates.capacity() *
                    sizeof(MpiMarkdupCandidateShared) +
                owner_qnames.capacity();
            size_t marker_key_bytes =
                owner_candidates.size() >
                        SIZE_MAX / sizeof(MpiMarkdupKeyShared)
                    ? SIZE_MAX
                    : owner_candidates.size() *
                        sizeof(MpiMarkdupKeyShared);
            size_t worst_results =
                owner_candidates.size() >
                        SIZE_MAX / sizeof(uint64_t)
                    ? SIZE_MAX
                    : owner_candidates.size() *
                        sizeof(uint64_t);
            size_t flat_hash_bytes = enable_flat_hash
                ? MdFlatHashBytesForEstimate(owner_candidates.size())
                : 0;
            size_t bucket_scratch_bytes = 0;
            if (enable_coord_bucket && owner_candidates.size() >
                    SIZE_MAX / (sizeof(MpiMarkdupCandidateShared) +
                                sizeof(uint64_t) +
                                sizeof(unsigned int) +
                                sizeof(size_t))) {
                bucket_scratch_bytes = SIZE_MAX;
            } else if (enable_coord_bucket) {
                bucket_scratch_bytes =
                    owner_candidates.size() *
                    (sizeof(MpiMarkdupCandidateShared) +
                     sizeof(uint64_t) + sizeof(unsigned int) +
                     sizeof(size_t));
            }
            size_t group_total = grouping_bytes;
            if (marker_key_bytes == SIZE_MAX ||
                worst_results == SIZE_MAX ||
                flat_hash_bytes == SIZE_MAX ||
                bucket_scratch_bytes == SIZE_MAX ||
                marker_key_bytes > SIZE_MAX - group_total) {
                local_ok = 0;
            } else {
                group_total += marker_key_bytes;
            }
            if (local_ok && worst_results > SIZE_MAX - group_total) {
                local_ok = 0;
            } else if (local_ok) {
                group_total += worst_results;
            }
            if (local_ok && flat_hash_bytes > SIZE_MAX - group_total) {
                local_ok = 0;
            } else if (local_ok) {
                group_total += flat_hash_bytes;
            }
            if (local_ok && bucket_scratch_bytes > SIZE_MAX - group_total) {
                local_ok = 0;
            } else if (local_ok) {
                group_total += bucket_scratch_bytes;
            }
            if (local_ok &&
                MdCheckMemory(group_total, memory_limit, &stats, rank,
                              "duplicate grouping") != 0) {
                local_ok = 0;
            }
            stats.t_group_prepare += GetTime() - group_prepare_t0;
        }
        if (!MdAllRanksOkTimed(local_ok, &stats)) goto cleanup;

        // 4. 在 owner rank 上判断重复。默认使用已验证最快的
        // contiguous flat-hash grouping。RABBITBAM_MARKDUP_FULL_SORT=1
        // 可回退 full sort/group，RABBITBAM_MARKDUP_COORD_BUCKET=1
        // 可启用 coordinate bucket 实验路径。
        if (enable_full_sort) {
            if (MdFindDuplicates(
                    &owner_candidates, owner_qnames,
                    comm_size, &duplicates_by_source,
                    &stats) != 0) {
                local_ok = 0;
            }
        } else if (enable_coord_bucket) {
            if (MdFindDuplicatesCoordinateBuckets(
                    &owner_candidates, owner_qnames,
                    comm_size, &duplicates_by_source,
                    &stats) != 0) {
                local_ok = 0;
            }
        } else {
            if (MpiMarkdupFindDuplicatesStreamingHash(
                    &owner_candidates, owner_qnames,
                    comm_size, &duplicates_by_source,
                    &stats) != 0) {
                local_ok = 0;
            }
        }
        if (!MdAllRanksOkTimed(local_ok, &stats)) goto cleanup;
        {
            double cleanup_t0 = GetTime();
            owner_candidates.clear();
            owner_candidates.shrink_to_fit();
            owner_qnames.clear();
            owner_qnames.shrink_to_fit();
            stats.t_owner_cleanup += GetTime() - cleanup_t0;
        }

        // 5. MdExchangeDuplicateResults：把 duplicate 结果发回原 rank
        if (MdExchangeDuplicateResults(
                &duplicates_by_source, rank,
                comm_size, local_records,
                memory_limit, &duplicate_bitmap,
                &stats) != 0) {
            local_ok = 0;
        }
        if (!MdAllRanksOkTimed(local_ok, &stats)) goto cleanup;
        {
            double cleanup_t0 = GetTime();
            duplicates_by_source.clear();
            stats.t_result_cleanup += GetTime() - cleanup_t0;
        }

        // 6. MdRewriteOutput：重写本 rank 的 BAM body
        const int remove_dups = cmd_info->markdup_remove_dups_;
        const int legacy_remove_path =
            MdEnvFlagEnabled("RABBITBAM_MARKDUP_LEGACY_REMOVE");
        const int fused_remove_dups =
            remove_dups &&
            (!cmd_info->markdup_clear_ || !legacy_remove_path);
        reader.pos = 0;
        if (MdRewriteOutput(
                reader, duplicate_bitmap, local_records,
                cmd_info->markdup_clear_,
                fused_remove_dups,
                cmd_info->compress_level_,
                writer, rank, &stats) != 0) {
            local_ok = 0;
        }
        if (!MdAllRanksOkTimed(local_ok, &stats)) goto cleanup;

        // 兼容回退路径：设置 RABBITBAM_MARKDUP_LEGACY_REMOVE=1 时，
        // -c -r 仍使用两段式稳定实现，便于和融合路径对照。
        if (remove_dups && !fused_remove_dups) {
            MemReader filter_reader = {};
            MemWriter filtered_writer = {};
            MpiBamToBamStats filter_stats = {};
            BamFilterOptions filter = {};
            filter.min_mapq = -1;
            filter.max_mapq = -1;
            filter.require_flag = 0;
            filter.exclude_flag = BAM_FDUP;
            filter.ref_tid = -2;
            filter.min_read_len = -1;
            filter.max_read_len = -1;

            filter_reader.base = writer.data;
            filter_reader.size = writer.size;
            filter_reader.pos = 0;
            if (MpiCommonInitMemWriter(
                    filtered_writer,
                    writer.size ? writer.size : 64 * 1024 * 1024) != 0) {
                local_ok = 0;
            }
            if (!MdAllRanksOkTimed(local_ok, &stats)) {
                if (filtered_writer.data) free(filtered_writer.data);
                goto cleanup;
            }
            if (FusedBamToBamMPI(filter_reader, filtered_writer,
                                 filter, cmd_info->compress_level_,
                                 &filter_stats) != 0) {
                local_ok = 0;
            }
            if (!MdAllRanksOkTimed(local_ok, &stats)) {
                if (filtered_writer.data) free(filtered_writer.data);
                goto cleanup;
            }
            free(writer.data);
            writer = filtered_writer;
            stats.removed_records += filter_stats.dropped_records;
            stats.bgzf_blocks = filter_stats.bgzf_blocks;
            stats.t_rewrite_decomp += filter_stats.t_decomp_filter;
            stats.t_pack += filter_stats.t_pack;
            stats.t_compress += filter_stats.t_compress;
            stats.t_write += filter_stats.t_write;
        }
        stats.t_fused_total = GetTime() - fused_t0;
        stage43_cost = MdReduceMax(stats.t_fused_total);
        if (rank == 0) {
            printf("Complete the 4.3 FusedMarkdupMPI cost %lf\n",
                   stage43_cost);
        }
    }

    // 4.4~4.6 MdPrepareOutputMemory：准备最终输出文件内存
    if (MdPrepareOutputMemory(writer, header,
                              cmd_info->compress_level_,
                              rank, comm_size, &stats,
                              &output_memory,
                              &output_stage_costs) != 0) {
        local_ok = 0;
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;

    if (rank == 0) {
        double core_total = stage41_cost + stage42_cost +
                            stage43_cost + output_stage_costs.stage44 +
                            output_stage_costs.stage45_header +
                            output_stage_costs.stage46;
        printf("Complete the total (4.1~4.6) cost %lf-----\n",
               core_total);
    }

    // 4.7 MdPrintStats打印统计信息
    {
        double stage47_t0 = GetTime();
        MdPrintStats(stats, rank, comm_size);
        double stage47_cost =
            MdReduceMax(GetTime() - stage47_t0);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.7 rank stats reduce/print cost %lf\n",
                   stage47_cost);
        }
    }

    {
        double body_total_cost =
            MdReduceMax(GetTime() - body_total_t0);
        if (rank == 0 && local_ok) {
            printf("444Complete the total body cost %lf\n",
                   body_total_cost);
        }
    }

    // 5. rank 0 为验证结果汇总输出内存并写入文件，不计入核心处理时间
    {
        double verify_alloc_t0 = GetTime();
        if (rank == 0) {
            output_memory.data =
                output_memory.size ? (char *)malloc(output_memory.size)
                                   : nullptr;
            if (output_memory.size && !output_memory.data) {
                fprintf(stderr,
                        "ERROR: failed to allocate markdup final output memory.\n");
                local_ok = 0;
            }
            if (local_ok && output_memory.header_size > 0) {
                memcpy(output_memory.data,
                       output_memory.header_data,
                       output_memory.header_size);
            }
        }
        double verify_alloc_cost =
            MdReduceMax(GetTime() - verify_alloc_t0);
        if (rank == 0 && local_ok) {
            printf("555Prepare verification output memory cost %lf--\n",
                   verify_alloc_cost);
        }
        if (!MdAllRanksOk(local_ok)) goto cleanup;

        double verify_gather_t0 = GetTime();
        if (rank == 0) {
            if (writer.size > 0) {
                memcpy(output_memory.data + output_memory.header_size +
                           output_memory.body_offsets[(size_t)rank],
                       writer.data, writer.size);
            }
            for (int src = 1; src < comm_size; ++src) {
                unsigned long long received = 0;
                unsigned long long bytes =
                    (unsigned long long)
                        output_memory.body_sizes[(size_t)src];
                while (received < bytes) {
                    int chunk = (int)std::min(
                        (unsigned long long)kMarkdupExchangeChunk,
                        bytes - received);
                    MPI_Recv(
                        output_memory.data + output_memory.header_size +
                            output_memory.body_offsets[(size_t)src] +
                            received,
                        chunk, MPI_BYTE, src, 4300,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    received += (unsigned long long)chunk;
                }
            }
            if (local_ok) {
                memcpy(output_memory.data + output_memory.header_size +
                           output_memory.total_body_size,
                       kMarkdupBgzfEofBlock,
                       sizeof(kMarkdupBgzfEofBlock));
            }
        } else {
            unsigned long long sent = 0;
            unsigned long long bytes =
                (unsigned long long)writer.size;
            while (sent < bytes) {
                int chunk = (int)std::min(
                    (unsigned long long)kMarkdupExchangeChunk,
                    bytes - sent);
                MPI_Send(writer.data + sent, chunk, MPI_BYTE,
                         0, 4300, MPI_COMM_WORLD);
                sent += (unsigned long long)chunk;
            }
        }
        double verify_gather_cost =
            MdReduceMax(GetTime() - verify_gather_t0);
        if (rank == 0 && local_ok) {
            printf("555Gather verification output memory cost %lf--\n",
                   verify_gather_cost);
        }
        if (!MdAllRanksOk(local_ok)) goto cleanup;

        double dump_t0 = 0.0;
        double dump_cost = 0.0;
        if (rank == 0) {
            dump_t0 = GetTime();
            if (MpiCommonDumpMemoryToFile(
                    cmd_info->out_file_name_,
                    output_memory.data, output_memory.size) != 0) {
                fprintf(stderr,
                        "ERROR: failed to write markdup output %s.\n",
                        cmd_info->out_file_name_.c_str());
                local_ok = 0;
            }
            dump_cost = GetTime() - dump_t0;
        }
        dump_cost = MdReduceMax(dump_cost);
        if (rank == 0 && local_ok) {
            printf("555Dump memory to output file cost %lf--\n",
                   dump_cost);
        }
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;
    exit_code = 0;

cleanup:
    if (output_memory.data) free(output_memory.data);
    if (output_memory.header_data) free(output_memory.header_data);
    if (writer.data) free(writer.data);
    if (header) sam_hdr_destroy(header);
    if (input) {
        if (hts_close(input) < 0) {
            fprintf(stderr,
                    "[rank %d] ERROR: closing markdup input failed.\n",
                    rank);
        }
        input_memory = nullptr;
    } else if (input_hfile) {
        if (hclose(input_hfile) != 0) {
            fprintf(stderr,
                    "[rank %d] ERROR: closing markdup memory hFILE failed.\n",
                    rank);
        }
        input_memory = nullptr;
    } else if (input_memory) {
        free(input_memory);
        input_memory = nullptr;
    }
    double close_cost = MdReduceMax(GetTime() - total_t0);
    if (rank == 0) {
        printf("666markdup total process cost %lf-----\n",
               close_cost);
    }
    return exit_code;
}
