#include "swbam_mpi.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

struct MdDuplicateId {
    int32_t source_rank;
    uint32_t pad;
    uint64_t ordinal;
};

static int MdAllRanksOk(int local_ok) {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
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
    const int records_per_block = (int)MPI_RECORDS_PER_BLOCK;
    MdBlockSet input = {};
    MdBlockSet uncompressed = {};
    MdRecordSet records = {};
    MdCandidateWorkspace candidate_workspace = {};
    Bam2BamPara decomp[kMarkdupNB];
    MpiMarkdupExtractPara extract[kMarkdupNB];

    if (MdAllocateBlockSet(&input, kMarkdupNB) != 0 ||
        MdAllocateBlockSet(&uncompressed, kMarkdupNB) != 0 ||
        MdAllocateRecordSet(&records, kMarkdupNB, records_per_block,
                            MPI_BAM_BLOCK_ARENA_SIZE) != 0 ||
        MdAllocateCandidateWorkspace(&candidate_workspace,
                                     records_per_block) != 0) {
        fprintf(stderr,
                "[rank %d] ERROR: failed to allocate markdup candidate workspace.\n",
                rank);
        MdFreeBlockSet(&input);
        MdFreeBlockSet(&uncompressed);
        MdFreeRecordSet(&records);
        MdFreeCandidateWorkspace(&candidate_workspace);
        return -1;
    }

    const size_t fixed_workspace =
        (size_t)kMarkdupNB * BGZF_MAX_BLOCK_SIZE * 2 +
        (size_t)kMarkdupNB * MPI_BAM_BLOCK_ARENA_SIZE +
        (size_t)kMarkdupNB * records_per_block *
            (sizeof(bam1_t) + sizeof(bam1_t *) + sizeof(uint32_t)) +
        (size_t)kMarkdupNB *
            candidate_workspace.candidates_per_block *
            sizeof(MpiMarkdupCandidateShared) +
        (size_t)kMarkdupNB * candidate_workspace.qname_stride;
    if (MdCheckMemory(fixed_workspace, memory_limit, stats, rank,
                      "candidate workspace") != 0) {
        MdFreeBlockSet(&input);
        MdFreeBlockSet(&uncompressed);
        MdFreeRecordSet(&records);
        MdFreeCandidateWorkspace(&candidate_workspace);
        return -1;
    }

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
    int n_blocks = 0;
    MdReadGroup(reader, &input, &n_blocks, stats);
    while (n_blocks > 0) {
        stats->input_blocks += n_blocks;
        for (int b = 0; b < kMarkdupNB; ++b) {
            MdInitDecompPara(decomp + b, b, &input, &uncompressed,
                             &records, b < n_blocks);
        }

        double decomp_t0 = GetTime();
        __real_athread_spawn(
            (void *)slave_mpi_decompress_bam2bam_passthrough,
            decomp, 1);
        athread_join();
        stats->t_candidate_decomp += GetTime() - decomp_t0;
        for (int b = 0; b < n_blocks; ++b) {
            if (decomp[b].status != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: markdup candidate BAM decode failed "
                        "at global block %lld status=%d record=%d.\n",
                        rank, block_group_base + b, decomp[b].status,
                        decomp[b].record_index);
                MdFreeBlockSet(&input);
                MdFreeBlockSet(&uncompressed);
                MdFreeRecordSet(&records);
                MdFreeCandidateWorkspace(&candidate_workspace);
                return -1;
            }
        }

        uint64_t block_ordinal = ordinal;
        for (int b = 0; b < kMarkdupNB; ++b) {
            memset(extract + b, 0, sizeof(extract[b]));
            extract[b].block_id = b;
            extract[b].status = b < n_blocks ? 0 : -1;
            if (b >= n_blocks) continue;
            extract[b].records = decomp[b].output_records;
            extract[b].n_records = decomp[b].n_total_records;
            extract[b].ordinal_base = block_ordinal;
            extract[b].global_block_index =
                (uint64_t)(block_group_base + b);
            extract[b].source_rank = rank;
            extract[b].include_fails = include_fails;
            extract[b].candidates =
                candidate_workspace.candidates +
                (size_t)b *
                    candidate_workspace.candidates_per_block;
            extract[b].candidate_capacity =
                candidate_workspace.candidates_per_block;
            extract[b].qname_arena =
                candidate_workspace.qnames +
                (size_t)b * candidate_workspace.qname_stride;
            extract[b].qname_capacity =
                candidate_workspace.qname_stride;
            block_ordinal +=
                (uint64_t)decomp[b].n_total_records;
        }

        double extract_t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_markdup_extract,
                             extract, 1);
        athread_join();
        stats->t_candidate_extract += GetTime() - extract_t0;

        for (int b = 0; b < n_blocks; ++b) {
            if (extract[b].status != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: %s at global block %lld "
                        "record=%d status=%d.\n",
                        rank, MdExtractErrorText(extract[b].status),
                        block_group_base + b, extract[b].record_index,
                        extract[b].status);
                MdFreeBlockSet(&input);
                MdFreeBlockSet(&uncompressed);
                MdFreeRecordSet(&records);
                MdFreeCandidateWorkspace(&candidate_workspace);
                return -1;
            }
            if (extract[b].has_records) {
                if (!*range_has_records) {
                    *range_first_tid = extract[b].first_tid;
                    *range_first_pos = extract[b].first_pos;
                    *range_has_records = 1;
                }
                if (extract[b].first_tid >= 0 && have_previous &&
                    (extract[b].first_tid < previous_tid ||
                     (extract[b].first_tid == previous_tid &&
                      extract[b].first_pos < previous_pos))) {
                    fprintf(stderr,
                            "[rank %d] ERROR: input is not coordinate sorted "
                            "near global block %lld.\n",
                            rank, block_group_base + b);
                    MdFreeBlockSet(&input);
                    MdFreeBlockSet(&uncompressed);
                    MdFreeRecordSet(&records);
                    MdFreeCandidateWorkspace(&candidate_workspace);
                    return -1;
                }
                previous_tid = extract[b].last_tid;
                previous_pos = extract[b].last_pos;
                *range_last_tid = extract[b].last_tid;
                *range_last_pos = extract[b].last_pos;
                have_previous = 1;
            }

            const uint64_t qname_base =
                (uint64_t)local_qnames->size();
            if (extract[b].qname_used > 0) {
                local_qnames->insert(
                    local_qnames->end(), extract[b].qname_arena,
                    extract[b].qname_arena + extract[b].qname_used);
            }
            for (int i = 0; i < extract[b].n_candidates; ++i) {
                MpiMarkdupCandidateShared candidate =
                    extract[b].candidates[i];
                if (!candidate.key.single) {
                    candidate.qname_offset += qname_base;
                }
                local_candidates->push_back(candidate);
            }
            stats->total_records += decomp[b].n_total_records;
            stats->examined_records += extract[b].examined;
            stats->excluded_records += extract[b].excluded;
            stats->pair_candidates +=
                extract[b].pair_candidates;
            stats->single_candidates +=
                extract[b].single_candidates;
        }
        ordinal = block_ordinal;
        block_group_base += n_blocks;
        stats->group_count++;

        size_t dynamic = MdVectorBytes(
            *local_candidates, *local_qnames,
            std::vector<uint8_t>());
        if (dynamic == SIZE_MAX ||
            fixed_workspace > SIZE_MAX - dynamic ||
            MdCheckMemory(fixed_workspace + dynamic, memory_limit,
                          stats, rank,
                          "local candidate accumulation") != 0) {
            MdFreeBlockSet(&input);
            MdFreeBlockSet(&uncompressed);
            MdFreeRecordSet(&records);
            MdFreeCandidateWorkspace(&candidate_workspace);
            return -1;
        }
        MdReadGroup(reader, &input, &n_blocks, stats);
    }

    *local_records = ordinal;
    MdFreeBlockSet(&input);
    MdFreeBlockSet(&uncompressed);
    MdFreeRecordSet(&records);
    MdFreeCandidateWorkspace(&candidate_workspace);
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

static int MdExchangeCandidates(
        const std::vector<MpiMarkdupCandidateShared> &local_candidates,
        const std::vector<unsigned char> &local_qnames,
        int rank, int comm_size, size_t memory_limit,
        std::vector<MpiMarkdupCandidateShared> *owner_candidates,
        std::vector<unsigned char> *owner_qnames,
        MpiMarkdupStats *stats) {
    std::vector<std::vector<MpiMarkdupCandidateShared> > send_candidates(
        (size_t)comm_size);
    std::vector<std::vector<unsigned char> > send_qnames(
        (size_t)comm_size);
    int local_ok = 1;

    for (size_t i = 0; i < local_candidates.size(); ++i) {
        MpiMarkdupCandidateShared candidate = local_candidates[i];
        int owner = (int)(MdHashKey(candidate.key) %
                          (uint64_t)comm_size);
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
            candidate.qname_offset =
                (uint64_t)send_qnames[(size_t)owner].size();
            send_qnames[(size_t)owner].insert(
                send_qnames[(size_t)owner].end(),
                local_qnames.begin() +
                    (size_t)local_candidates[i].qname_offset,
                local_qnames.begin() +
                    (size_t)local_candidates[i].qname_offset +
                    local_candidates[i].qname_len);
        }
        send_candidates[(size_t)owner].push_back(candidate);
    }
    if (!MdAllRanksOk(local_ok)) return -1;

    size_t exchange_memory = 0;
    for (int i = 0; i < comm_size; ++i) {
        size_t bytes =
            send_candidates[(size_t)i].capacity() *
                sizeof(MpiMarkdupCandidateShared) +
            send_qnames[(size_t)i].capacity();
        if (bytes > SIZE_MAX - exchange_memory) {
            local_ok = 0;
        } else {
            exchange_memory += bytes;
        }
    }
    size_t local_memory =
        local_candidates.capacity() *
            sizeof(MpiMarkdupCandidateShared) +
        local_qnames.capacity();
    unsigned long long send_candidate_counts[kMarkdupTagCount] = {};
    unsigned long long recv_candidate_counts[kMarkdupTagCount] = {};
    unsigned long long send_qname_counts[kMarkdupTagCount] = {};
    unsigned long long recv_qname_counts[kMarkdupTagCount] = {};
    if (comm_size > kMarkdupTagCount) {
        std::vector<unsigned long long> send_candidates_dynamic(
            (size_t)comm_size);
        std::vector<unsigned long long> recv_candidates_dynamic(
            (size_t)comm_size);
        std::vector<unsigned long long> send_qnames_dynamic(
            (size_t)comm_size);
        std::vector<unsigned long long> recv_qnames_dynamic(
            (size_t)comm_size);
        for (int i = 0; i < comm_size; ++i) {
            send_candidates_dynamic[(size_t)i] =
                (unsigned long long)send_candidates[(size_t)i].size();
            send_qnames_dynamic[(size_t)i] =
                (unsigned long long)send_qnames[(size_t)i].size();
        }
        MPI_Alltoall(send_candidates_dynamic.data(), 1,
                     MPI_UNSIGNED_LONG_LONG,
                     recv_candidates_dynamic.data(), 1,
                     MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
        MPI_Alltoall(send_qnames_dynamic.data(), 1,
                     MPI_UNSIGNED_LONG_LONG,
                     recv_qnames_dynamic.data(), 1,
                     MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
        unsigned long long owner_candidate_count = 0;
        unsigned long long owner_qname_count = 0;
        for (int i = 0; i < comm_size; ++i) {
            owner_candidate_count +=
                recv_candidates_dynamic[(size_t)i];
            owner_qname_count += recv_qnames_dynamic[(size_t)i];
        }
        if (owner_candidate_count >
                (unsigned long long)(SIZE_MAX /
                    sizeof(MpiMarkdupCandidateShared)) ||
            owner_qname_count > (unsigned long long)SIZE_MAX) {
            local_ok = 0;
        } else {
            size_t owner_bytes =
                (size_t)owner_candidate_count *
                    sizeof(MpiMarkdupCandidateShared) +
                (size_t)owner_qname_count;
            if (local_memory > SIZE_MAX - exchange_memory ||
                local_memory + exchange_memory >
                    SIZE_MAX - owner_bytes ||
                MdCheckMemory(local_memory + exchange_memory +
                                  owner_bytes,
                              memory_limit, stats, rank,
                              "candidate owner exchange") != 0) {
                local_ok = 0;
            } else {
                owner_candidates->reserve(
                    (size_t)owner_candidate_count);
                owner_qnames->reserve((size_t)owner_qname_count);
            }
        }
    } else {
        for (int i = 0; i < comm_size; ++i) {
            send_candidate_counts[i] =
                (unsigned long long)send_candidates[(size_t)i].size();
            send_qname_counts[i] =
                (unsigned long long)send_qnames[(size_t)i].size();
        }
        MPI_Alltoall(send_candidate_counts, 1,
                     MPI_UNSIGNED_LONG_LONG,
                     recv_candidate_counts, 1,
                     MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
        MPI_Alltoall(send_qname_counts, 1,
                     MPI_UNSIGNED_LONG_LONG,
                     recv_qname_counts, 1,
                     MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
        unsigned long long owner_candidate_count = 0;
        unsigned long long owner_qname_count = 0;
        for (int i = 0; i < comm_size; ++i) {
            owner_candidate_count += recv_candidate_counts[i];
            owner_qname_count += recv_qname_counts[i];
        }
        if (owner_candidate_count >
                (unsigned long long)(SIZE_MAX /
                    sizeof(MpiMarkdupCandidateShared)) ||
            owner_qname_count > (unsigned long long)SIZE_MAX) {
            local_ok = 0;
        } else {
            size_t owner_bytes =
                (size_t)owner_candidate_count *
                    sizeof(MpiMarkdupCandidateShared) +
                (size_t)owner_qname_count;
            if (local_memory > SIZE_MAX - exchange_memory ||
                local_memory + exchange_memory >
                    SIZE_MAX - owner_bytes ||
                MdCheckMemory(local_memory + exchange_memory +
                                  owner_bytes,
                              memory_limit, stats, rank,
                              "candidate owner exchange") != 0) {
                local_ok = 0;
            } else {
                owner_candidates->reserve(
                    (size_t)owner_candidate_count);
                owner_qnames->reserve((size_t)owner_qname_count);
            }
        }
    }
    if (!MdAllRanksOk(local_ok)) {
        return -1;
    }

    double exchange_t0 = GetTime();
    for (int step = 0; step < comm_size; ++step) {
        int dest = (rank + step) % comm_size;
        int src = (rank - step + comm_size) % comm_size;
        if (step == 0) {
            const uint64_t qbase = (uint64_t)owner_qnames->size();
            size_t candidate_base = owner_candidates->size();
            owner_qnames->insert(
                owner_qnames->end(),
                send_qnames[(size_t)rank].begin(),
                send_qnames[(size_t)rank].end());
            owner_candidates->insert(
                owner_candidates->end(),
                send_candidates[(size_t)rank].begin(),
                send_candidates[(size_t)rank].end());
            for (size_t i = candidate_base;
                 i < owner_candidates->size(); ++i) {
                if (!(*owner_candidates)[i].key.single) {
                    (*owner_candidates)[i].qname_offset += qbase;
                }
            }
            continue;
        }

        unsigned long long send_counts[2] = {
            (unsigned long long)send_candidates[(size_t)dest].size(),
            (unsigned long long)send_qnames[(size_t)dest].size()
        };
        unsigned long long recv_counts[2] = {0, 0};
        if (MPI_Sendrecv(send_counts, 2, MPI_UNSIGNED_LONG_LONG,
                         dest, 4100,
                         recv_counts, 2, MPI_UNSIGNED_LONG_LONG,
                         src, 4100, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE) != MPI_SUCCESS) {
            return -1;
        }
        if (recv_counts[0] >
            (unsigned long long)(SIZE_MAX /
                sizeof(MpiMarkdupCandidateShared)) ||
            recv_counts[1] > (unsigned long long)SIZE_MAX) {
            return -1;
        }
        size_t candidate_base = owner_candidates->size();
        size_t qname_base = owner_qnames->size();
        owner_candidates->resize(
            candidate_base + (size_t)recv_counts[0]);
        owner_qnames->resize(
            qname_base + (size_t)recv_counts[1]);

        unsigned long long send_candidate_bytes =
            send_counts[0] *
            sizeof(MpiMarkdupCandidateShared);
        unsigned long long recv_candidate_bytes =
            recv_counts[0] *
            sizeof(MpiMarkdupCandidateShared);
        if (MdSendrecvBytes(
                send_candidates[(size_t)dest].data(),
                send_candidate_bytes, dest,
                recv_counts[0]
                    ? owner_candidates->data() + candidate_base
                    : nullptr,
                recv_candidate_bytes, src, 4101) != 0 ||
            MdSendrecvBytes(
                send_qnames[(size_t)dest].data(),
                send_counts[1], dest,
                recv_counts[1]
                    ? owner_qnames->data() + qname_base
                    : nullptr,
                recv_counts[1], src, 4102) != 0) {
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
            (long long)send_counts[1];

    }
    stats->t_candidate_exchange += GetTime() - exchange_t0;
    stats->owner_candidates =
        (long long)owner_candidates->size();
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

static int MdFindDuplicates(
        std::vector<MpiMarkdupCandidateShared> *owner_candidates,
        const std::vector<unsigned char> &owner_qnames,
        int comm_size,
        std::vector<std::vector<uint64_t> > *duplicates_by_source,
        MpiMarkdupStats *stats) {
    double group_t0 = GetTime();
    std::sort(owner_candidates->begin(), owner_candidates->end(),
              [](const MpiMarkdupCandidateShared &a,
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
    });

    duplicates_by_source->assign(
        (size_t)comm_size, std::vector<uint64_t>());
    size_t begin = 0;
    while (begin < owner_candidates->size()) {
        size_t end = begin + 1;
        while (end < owner_candidates->size() &&
               MdKeyEqual((*owner_candidates)[begin].key,
                          (*owner_candidates)[end].key)) {
            ++end;
        }

        if (!(*owner_candidates)[begin].key.single) {
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
        } else {
            bool has_pair = false;
            for (size_t i = begin; i < end; ++i) {
                if ((*owner_candidates)[i].paired_marker) {
                    has_pair = true;
                    break;
                }
            }
            if (has_pair) {
                for (size_t i = begin; i < end; ++i) {
                    if ((*owner_candidates)[i].paired_marker) continue;
                    int source =
                        (*owner_candidates)[i].source_rank;
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
                         candidate.global_order <
                             best.global_order)) {
                        winner = i;
                    }
                }
                for (size_t i = begin; i < end; ++i) {
                    if (i == winner) continue;
                    int source =
                        (*owner_candidates)[i].source_rank;
                    if (source < 0 || source >= comm_size) return -1;
                    (*duplicates_by_source)[(size_t)source]
                        .push_back((*owner_candidates)[i].ordinal);
                    stats->single_duplicates++;
                }
            }
        }
        begin = end;
    }
    stats->t_group += GetTime() - group_t0;
    return 0;
}

static int MdExchangeDuplicateResults(
        std::vector<std::vector<uint64_t> > *duplicates_by_source,
        int rank, int comm_size, uint64_t local_records,
        size_t memory_limit, std::vector<uint8_t> *bitmap,
        MpiMarkdupStats *stats) {
    for (int i = 0; i < comm_size; ++i) {
        std::vector<uint64_t> &values =
            (*duplicates_by_source)[(size_t)i];
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()),
                     values.end());
    }
    int local_ok =
        local_records <= (uint64_t)SIZE_MAX * 8ull ? 1 : 0;
    if (!MdAllRanksOk(local_ok)) return -1;
    bitmap->assign((size_t)((local_records + 7) / 8), 0);

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

    double exchange_t0 = GetTime();
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

static int MdDeleteAuxTagHost(bam1_t *record, const char tag[2]) {
    uint8_t *aux = bam_aux_get(record, tag);
    if (!aux) return 0;
    return bam_aux_del(record, aux);
}

static int MdRewriteBlockHost(MpiMarkdupRewritePara *para) {
    if (!para || para->status != 0 || !para->records || !para->bam_lens) {
        return 0;
    }
    int kept = 0;
    uint32_t kept_len = 0;
    para->marked_records = 0;
    para->cleared_records = 0;
    para->removed_records = 0;

    for (int i = 0; i < para->n_records; ++i) {
        bam1_t *record = para->records[i];
        if (!record) {
            para->status = -2;
            para->record_index = i;
            return -1;
        }
        const uint64_t ordinal = para->ordinal_base + (uint64_t)i;
        const size_t byte = (size_t)(ordinal >> 3);
        const uint8_t mask = (uint8_t)(1u << (ordinal & 7u));
        const int duplicate =
            byte < para->duplicate_bitmap_bytes &&
            (para->duplicate_bitmap[byte] & mask);

        if (para->clear_old) {
            if (record->core.flag & BAM_FDUP) para->cleared_records++;
            record->core.flag &= (uint16_t)~BAM_FDUP;
            if (MdDeleteAuxTagHost(record, "dt") < 0 ||
                MdDeleteAuxTagHost(record, "do") < 0) {
                para->status = -2;
                para->record_index = i;
                return -1;
            }
        }
        if (duplicate) {
            record->core.flag |= BAM_FDUP;
            para->marked_records++;
        }
        if (para->remove_dups && (record->core.flag & BAM_FDUP)) {
            para->removed_records++;
            continue;
        }

        const uint32_t bam_len =
            (uint32_t)(record->l_data - record->core.l_extranul + 32);
        para->records[kept] = record;
        para->bam_lens[kept] = bam_len;
        kept_len += bam_len + 4;
        kept++;
    }

    para->n_kept_records = kept;
    para->kept_total_len = kept_len;
    para->status = 0;
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

static int MdAppendRecordPayload(bam_block *payload, const bam1_t *record) {
    if (!payload || !record || !record->data) return -1;
    const bam1_core_t *core = &record->core;
    const int raw_l_qname = core->l_qname - core->l_extranul;
    if (raw_l_qname <= 0 || raw_l_qname > 255 ||
        core->l_extranul > core->l_qname ||
        record->l_data < core->l_qname ||
        core->n_cigar > 0xffff ||
        core->pos > INT_MAX ||
        core->mpos > INT_MAX) {
        return -1;
    }

    const uint32_t block_len =
        (uint32_t)(record->l_data - core->l_extranul + 32);
    const uint32_t packed_len = block_len + 4;
    if (packed_len > BGZF_BLOCK_SIZE) return -1;
    if (payload->pos > 0 &&
        payload->pos + (int)packed_len > BGZF_BLOCK_SIZE) {
        return 1;
    }
    if (payload->pos + (int)packed_len > BGZF_BLOCK_SIZE) return -1;

    uint8_t *dst = payload->data + payload->pos;
    uint32_t fields[8];
    fields[0] = (uint32_t)core->tid;
    fields[1] = (uint32_t)core->pos;
    fields[2] = (uint32_t)core->bin << 16 |
                (uint32_t)core->qual << 8 |
                (uint32_t)raw_l_qname;
    fields[3] = (uint32_t)core->flag << 16 |
                (uint32_t)(core->n_cigar & 0xffff);
    fields[4] = (uint32_t)core->l_qseq;
    fields[5] = (uint32_t)core->mtid;
    fields[6] = (uint32_t)core->mpos;
    fields[7] = (uint32_t)core->isize;

    memcpy(dst, &block_len, 4);
    memcpy(dst + 4, fields, sizeof(fields));
    memcpy(dst + 36, record->data, (size_t)raw_l_qname);
    const uint32_t rest_len =
        (uint32_t)(record->l_data - core->l_qname);
    memcpy(dst + 36 + raw_l_qname,
           record->data + core->l_qname,
           (size_t)rest_len);
    payload->pos += (int)packed_len;
    payload->length = payload->pos;
    return 0;
}

static int MdPreparePayloadFromRewrite(
        MpiMarkdupRewritePara *rewrite, int n_blocks,
        MdBlockSet *payload_blocks, MdBlockSet *output_blocks,
        MpiSortRawCompressPara *paras, int compress_level,
        int *active_blocks) {
    int active = 0;
    if (payload_blocks->n < kMarkdupNB ||
        output_blocks->n < kMarkdupNB) {
        return -1;
    }
    MdResetPayloadBlock(payload_blocks->blocks + active, active);

    for (int b = 0; b < n_blocks; ++b) {
        for (int r = 0; r < rewrite[b].n_kept_records; ++r) {
            bam1_t *record = rewrite[b].records[r];
            while (true) {
                int ret = MdAppendRecordPayload(
                    payload_blocks->blocks + active, record);
                if (ret == 0) break;
                if (ret < 0) return ret;

                MdSetupRawCompress(paras + active, active,
                                   payload_blocks->blocks + active,
                                   output_blocks->blocks + active,
                                   compress_level);
                ++active;
                if (active >= kMarkdupNB) return -3;
                MdResetPayloadBlock(payload_blocks->blocks + active,
                                    active);
            }
        }
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
        athread_join();
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
        athread_join();
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
            rewrite[b].remove_dups = remove_dups;
            block_ordinal +=
                (uint64_t)decomp[b].n_total_records;
        }
        double rewrite_t0 = GetTime();
        if (remove_dups) {
            for (int b = 0; b < n_blocks; ++b) {
                if (MdRewriteBlockHost(rewrite + b) != 0) break;
            }
        } else {
            __real_athread_spawn((void *)slave_mpi_markdup_rewrite,
                                 rewrite, 1);
            athread_join();
        }
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
        if (remove_dups) {
            int active_payload_blocks = 0;
            double payload_t0 = GetTime();
            int payload_status = MdPreparePayloadFromRewrite(
                rewrite, n_blocks, compress_un_active,
                output_active, raw_comp, compress_level,
                &active_payload_blocks);
            stats->t_pack += GetTime() - payload_t0;
            if (payload_status != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: markdup remove-dups payload "
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
            athread_join();
            stats->t_compress += GetTime() - compress_t0;
            double write_t0 = GetTime();
            for (int i = 0; i < active_payload_blocks; ++i) {
                if (raw_comp[i].status != 0) {
                    fprintf(stderr,
                            "[rank %d] ERROR: markdup remove-dups "
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

static int MdGatherOutput(
        const MemWriter &local_writer, sam_hdr_t *header,
        int compress_level, const std::string &output_path,
        int rank, int comm_size, MpiMarkdupStats *stats) {
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

    char *header_memory = nullptr;
    size_t header_size = 0;
    char *output_memory = nullptr;
    size_t output_size = 0;
    int local_ok = 1;
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
            output_size = (size_t)final_size;
            output_memory =
                output_size ? (char *)malloc(output_size) : nullptr;
            if (output_size && !output_memory) local_ok = 0;
        }
    }
    if (!MdAllRanksOk(local_ok)) {
        if (header_memory) free(header_memory);
        if (output_memory) free(output_memory);
        return -1;
    }

    double gather_t0 = GetTime();
    if (rank == 0) {
        memcpy(output_memory, header_memory, header_size);
        if (local_writer.size > 0) {
            memcpy(output_memory + header_size,
                   local_writer.data, local_writer.size);
        }
        for (int src = 1; src < comm_size; ++src) {
            unsigned long long received = 0;
            unsigned long long bytes =
                (unsigned long long)sizes[(size_t)src];
            while (received < bytes) {
                int chunk = (int)std::min(
                    (unsigned long long)kMarkdupExchangeChunk,
                    bytes - received);
                MPI_Recv(output_memory + header_size +
                             offsets[(size_t)src] + received,
                         chunk, MPI_BYTE, src, 4300,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                received += (unsigned long long)chunk;
            }
        }
        memcpy(output_memory + header_size + total_body,
               kMarkdupBgzfEofBlock,
               sizeof(kMarkdupBgzfEofBlock));
    } else {
        unsigned long long sent = 0;
        unsigned long long bytes =
            (unsigned long long)local_writer.size;
        while (sent < bytes) {
            int chunk = (int)std::min(
                (unsigned long long)kMarkdupExchangeChunk,
                bytes - sent);
            MPI_Send(local_writer.data + sent, chunk, MPI_BYTE,
                     0, 4300, MPI_COMM_WORLD);
            sent += (unsigned long long)chunk;
        }
    }
    stats->t_write += GetTime() - gather_t0;

    double dump_t0 = GetTime();
    if (rank == 0 &&
        MpiCommonDumpMemoryToFile(
            output_path, output_memory, output_size) != 0) {
        fprintf(stderr,
                "ERROR: failed to write markdup output %s.\n",
                output_path.c_str());
        local_ok = 0;
    }
    double dump_cost = MdReduceMax(GetTime() - dump_t0);
    if (rank == 0 && local_ok) {
        printf("555Dump memory to output file cost %lf--\n",
               dump_cost);
    }
    if (header_memory) free(header_memory);
    if (output_memory) free(output_memory);
    return MdAllRanksOk(local_ok) ? 0 : -1;
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
    long long sums[count] = {};
    MPI_Reduce(local_counts, sums, count, MPI_LONG_LONG,
               MPI_SUM, 0, MPI_COMM_WORLD);
    long long peak = stats.tracked_peak_bytes;
    long long peak_max = 0;
    MPI_Reduce(&peak, &peak_max, 1, MPI_LONG_LONG,
               MPI_MAX, 0, MPI_COMM_WORLD);

    double local_times[10] = {
        stats.t_candidate_decomp,
        stats.t_candidate_extract,
        stats.t_candidate_exchange,
        stats.t_group,
        stats.t_result_exchange,
        stats.t_rewrite_decomp,
        stats.t_rewrite,
        stats.t_pack,
        stats.t_compress,
        stats.t_fused_total
    };
    double time_sums[10] = {};
    MPI_Reduce(local_times, time_sums, 10, MPI_DOUBLE,
               MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
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
               "candidate_extract=%.3f exchange=%.3f group=%.3f "
               "result=%.3f rewrite_decomp=%.3f rewrite=%.3f "
               "pack=%.3f compress=%.3f fused=%.3f\n",
               time_sums[0], time_sums[1], time_sums[2],
               time_sums[3], time_sums[4], time_sums[5],
               time_sums[6], time_sums[7], time_sums[8],
               time_sums[9]);
    }
}

} // namespace

int ProcessMarkdupMPI(CmdInfo *cmd_info) {
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

    {
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
        double cost = MdReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("MPI BAM scan complete. data_blocks=%lld "
                   "body_start=%lld\n",
                   n_blocks, body_start);
            printf("Complete the 4.1 scan/split/broadcast/select cost %lf\n",
                   cost);
        }
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;

    reader.base = rank_input;
    reader.size = rank_input_size;
    reader.pos = 0;
    if (MpiCommonInitMemWriter(
            writer, rank_input_size ?
                        rank_input_size :
                        64 * 1024 * 1024) != 0) {
        local_ok = 0;
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;

    {
        double fused_t0 = GetTime();
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
        if (!MdAllRanksOk(local_ok)) goto cleanup;

        if (MdValidateRankBoundaries(
                range_has_records, range_first_tid, range_first_pos,
                range_last_tid, range_last_pos, rank,
                comm_size) != 0) {
            local_ok = 0;
        }
        if (!MdAllRanksOk(local_ok)) goto cleanup;

        if (MdExchangeCandidates(
                local_candidates, local_qnames,
                rank, comm_size, memory_limit,
                &owner_candidates, &owner_qnames,
                &stats) != 0) {
            local_ok = 0;
        }
        if (!MdAllRanksOk(local_ok)) goto cleanup;

        local_candidates.clear();
        local_candidates.shrink_to_fit();
        local_qnames.clear();
        local_qnames.shrink_to_fit();

        {
            size_t grouping_bytes =
                owner_candidates.capacity() *
                    sizeof(MpiMarkdupCandidateShared) +
                owner_qnames.capacity();
            size_t worst_results =
                owner_candidates.size() >
                        SIZE_MAX / sizeof(uint64_t)
                    ? SIZE_MAX
                    : owner_candidates.size() *
                        sizeof(uint64_t);
            if (worst_results == SIZE_MAX ||
                grouping_bytes > SIZE_MAX - worst_results ||
                MdCheckMemory(grouping_bytes + worst_results,
                              memory_limit, &stats, rank,
                              "duplicate grouping") != 0) {
                local_ok = 0;
            }
        }
        if (!MdAllRanksOk(local_ok)) goto cleanup;

        if (MdFindDuplicates(
                &owner_candidates, owner_qnames,
                comm_size, &duplicates_by_source,
                &stats) != 0) {
            local_ok = 0;
        }
        if (!MdAllRanksOk(local_ok)) goto cleanup;
        owner_candidates.clear();
        owner_candidates.shrink_to_fit();
        owner_qnames.clear();
        owner_qnames.shrink_to_fit();

        if (MdExchangeDuplicateResults(
                &duplicates_by_source, rank,
                comm_size, local_records,
                memory_limit, &duplicate_bitmap,
                &stats) != 0) {
            local_ok = 0;
        }
        if (!MdAllRanksOk(local_ok)) goto cleanup;
        duplicates_by_source.clear();

        const int remove_dups = cmd_info->markdup_remove_dups_;
        reader.pos = 0;
        if (MdRewriteOutput(
                reader, duplicate_bitmap, local_records,
                cmd_info->markdup_clear_,
                0,
                cmd_info->compress_level_,
                writer, rank, &stats) != 0) {
            local_ok = 0;
        }
        if (!MdAllRanksOk(local_ok)) goto cleanup;
        if (remove_dups) {
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
            if (!MdAllRanksOk(local_ok)) {
                if (filtered_writer.data) free(filtered_writer.data);
                goto cleanup;
            }
            if (FusedBamToBamMPI(filter_reader, filtered_writer,
                                 filter, cmd_info->compress_level_,
                                 &filter_stats) != 0) {
                local_ok = 0;
            }
            if (!MdAllRanksOk(local_ok)) {
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
        double fused_cost = MdReduceMax(stats.t_fused_total);
        if (rank == 0) {
            printf("Complete the 4.3 FusedMarkdupMPI cost %lf\n",
                   fused_cost);
        }
    }

    MdPrintStats(stats, rank, comm_size);

    if (MdGatherOutput(writer, header,
                       cmd_info->compress_level_,
                       cmd_info->out_file_name_,
                       rank, comm_size, &stats) != 0) {
        local_ok = 0;
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;
    exit_code = 0;

cleanup:
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
