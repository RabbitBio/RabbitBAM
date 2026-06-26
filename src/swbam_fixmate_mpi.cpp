#include "swbam_mpi.h"

#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

#include <mpi.h>

extern "C" {
    void slave_mpi_decompress_bam2bam_passthrough();
    void slave_mpi_fixmate_plan();
    void slave_mpi_fixmate_rewrite();
    void slave_mpi_compressfunc();
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
int MpiCommonReadBamHeaderFromMemory(
    const char *data, size_t size,
    sam_hdr_t **header, long long *body_start);

namespace {

const int kFixmateNB = 64;
const int kFixmateExchangeChunk = INT_MAX / 2;
const int kFixmateQnameBytes = 256;

const unsigned char kFixmateBgzfEofBlock[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
};

struct FmBlockSet {
    bam_block *blocks;
    unsigned char *data;
    int n;
};

struct FmRecordSet {
    bam1_t *records;
    unsigned char *data;
    bam1_t **ptrs;
    uint32_t *bam_lens;
    int records_per_block;
    int n_blocks;
    size_t arena_stride;
};

struct FmStoredGroup {
    std::string qname;
    std::vector<bam1_core_t> cores;
    std::vector<uint64_t> ids;
    std::vector<uint32_t> offsets;
    std::vector<uint32_t> lengths;
    std::vector<unsigned char> data;

    void clear() {
        qname.clear();
        cores.clear();
        ids.clear();
        offsets.clear();
        lengths.clear();
        data.clear();
    }

    bool empty() const {
        return cores.empty();
    }
};

struct FmBoundaryInfo {
    int has_records;
    int single_group;
    char first_qname[kFixmateQnameBytes];
    char last_qname[kFixmateQnameBytes];
};

struct FmWireRecordHeader {
    bam1_core_t core;
    uint64_t id;
    uint32_t l_data;
    uint32_t reserved;
};

struct FmPackPlan {
    bam1_t **records;
    int n_records;
    uint32_t total_len;
};

static int FmAllRanksOk(int local_ok) {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    return global_ok;
}

static int FmAllRanksOkTimed(int local_ok, MpiFixmateStats *stats) {
    double t0 = GetTime();
    const int global_ok = FmAllRanksOk(local_ok);
    stats->t_rank_sync += GetTime() - t0;
    return global_ok;
}

static double FmReduceMax(double local) {
    double result = 0.0;
    MPI_Reduce(&local, &result, 1, MPI_DOUBLE, MPI_MAX, 0,
               MPI_COMM_WORLD);
    return result;
}

struct FmOutputStageCosts {
    double stage44;
    double stage45_header;
    double stage45_malloc;
    double stage46;
};

static int FmAllocateBlockSet(FmBlockSet *set, int n) {
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

static void FmFreeBlockSet(FmBlockSet *set) {
    if (set->blocks) aligned_free_custom((unsigned char *)set->blocks);
    if (set->data) aligned_free_custom(set->data);
    memset(set, 0, sizeof(*set));
}

static int FmAllocateRecordSet(FmRecordSet *set, int n_blocks,
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
    if (!set->records || !set->data || !set->ptrs ||
        !set->bam_lens) {
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

static void FmFreeRecordSet(FmRecordSet *set) {
    if (set->records) aligned_free_custom((unsigned char *)set->records);
    if (set->data) aligned_free_custom(set->data);
    if (set->ptrs) aligned_free_custom((unsigned char *)set->ptrs);
    if (set->bam_lens) {
        aligned_free_custom((unsigned char *)set->bam_lens);
    }
    memset(set, 0, sizeof(*set));
}

static int FmMemReadBlock(MemReader &reader, bam_block *block) {
    if (reader.pos >= reader.size ||
        reader.size - reader.pos < BLOCK_HEADER_LENGTH) {
        return -1;
    }
    const unsigned char *base =
        (const unsigned char *)reader.base + reader.pos;
    int length = (int)base[16] | ((int)base[17] << 8);
    length++;
    if (length <= 0 || length > BGZF_MAX_BLOCK_SIZE ||
        (size_t)length > reader.size - reader.pos) {
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

static int FmReadGroup(MemReader &reader, FmBlockSet *input,
                       int *n_blocks, MpiFixmateStats *stats) {
    double t0 = GetTime();
    int count = 0;
    for (int i = 0; i < kFixmateNB; ++i) {
        if (FmMemReadBlock(reader, input->blocks + i) != 0) break;
        input->blocks[i].block_id = i;
        ++count;
    }
    *n_blocks = count;
    stats->t_read += GetTime() - t0;
    return 0;
}

static void FmInitDecompPara(
        Bam2BamPara *para, int block_id,
        FmBlockSet *input, FmBlockSet *uncompressed,
        FmRecordSet *records, int active) {
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
        para->status = -1;
    }
}

static int FmAppendStoredRecord(FmStoredGroup *group,
                                const bam1_t *record) {
    const char *qname = bam_get_qname(record);
    if (group->empty()) {
        group->qname = qname;
    } else if (group->qname != qname) {
        return -1;
    }
    if (record->l_data < 0 ||
        group->data.size() >
            (size_t)UINT32_MAX - (size_t)record->l_data) {
        return -1;
    }
    group->cores.push_back(record->core);
    group->ids.push_back(record->id);
    group->offsets.push_back((uint32_t)group->data.size());
    group->lengths.push_back((uint32_t)record->l_data);
    group->data.insert(group->data.end(), record->data,
                       record->data + record->l_data);
    return 0;
}

static int FmAppendStoredRange(FmStoredGroup *group,
                               const std::vector<bam1_t *> &records,
                               size_t begin, size_t end) {
    for (size_t i = begin; i < end; ++i) {
        if (FmAppendStoredRecord(group, records[i]) != 0) return -1;
    }
    return 0;
}

static int FmBuildStoredViews(
        FmStoredGroup *group, std::vector<bam1_t> *records,
        std::vector<bam1_t *> *ptrs) {
    records->resize(group->cores.size());
    ptrs->resize(group->cores.size());
    for (size_t i = 0; i < group->cores.size(); ++i) {
        bam1_t &record = (*records)[i];
        memset(&record, 0, sizeof(record));
        record.core = group->cores[i];
        record.id = group->ids[i];
        record.data = group->data.data() + group->offsets[i];
        record.l_data = (int)group->lengths[i];
        record.m_data = group->lengths[i];
        record.mempolicy = BAM_USER_OWNS_DATA;
        (*ptrs)[i] = &record;
    }
    return 0;
}

static int FmSerializeStoredGroup(
        const FmStoredGroup &group,
        std::vector<unsigned char> *wire) {
    const uint64_t count = group.cores.size();
    size_t total = sizeof(count);
    for (size_t i = 0; i < group.cores.size(); ++i) {
        if (total > SIZE_MAX - sizeof(FmWireRecordHeader) ||
            total + sizeof(FmWireRecordHeader) >
                SIZE_MAX - group.lengths[i]) {
            return -1;
        }
        total += sizeof(FmWireRecordHeader) + group.lengths[i];
    }
    wire->resize(total);
    unsigned char *out = wire->data();
    memcpy(out, &count, sizeof(count));
    out += sizeof(count);
    for (size_t i = 0; i < group.cores.size(); ++i) {
        FmWireRecordHeader header;
        memset(&header, 0, sizeof(header));
        header.core = group.cores[i];
        header.id = group.ids[i];
        header.l_data = group.lengths[i];
        memcpy(out, &header, sizeof(header));
        out += sizeof(header);
        memcpy(out, group.data.data() + group.offsets[i],
               group.lengths[i]);
        out += group.lengths[i];
    }
    return 0;
}

static int FmAppendSerializedGroup(
        FmStoredGroup *group,
        const unsigned char *wire, size_t wire_size) {
    if (!wire || wire_size < sizeof(uint64_t)) return -1;
    uint64_t count = 0;
    memcpy(&count, wire, sizeof(count));
    size_t pos = sizeof(count);
    for (uint64_t i = 0; i < count; ++i) {
        if (wire_size - pos < sizeof(FmWireRecordHeader)) return -1;
        FmWireRecordHeader header;
        memcpy(&header, wire + pos, sizeof(header));
        pos += sizeof(header);
        if (header.l_data > wire_size - pos) return -1;
        bam1_t record;
        memset(&record, 0, sizeof(record));
        record.core = header.core;
        record.id = header.id;
        record.data = (uint8_t *)(wire + pos);
        record.l_data = (int)header.l_data;
        record.m_data = header.l_data;
        record.mempolicy = BAM_USER_OWNS_DATA;
        if (FmAppendStoredRecord(group, &record) != 0) return -1;
        pos += header.l_data;
    }
    return pos == wire_size ? 0 : -1;
}

static void FmInitEmptyComp(Comp_Para *para, int block_id) {
    memset(para, 0, sizeof(*para));
    para->block_id = block_id;
    para->status = -1;
}

static int FmCompressPlans(
        const std::vector<bam1_t *> &records,
        const std::vector<uint32_t> &bam_lens,
        int compress_level, MemWriter *writer,
        int rank, MpiFixmateStats *stats) {
    double workspace_t0 = GetTime();
    FmBlockSet un_a = {};
    FmBlockSet un_b = {};
    FmBlockSet out_a = {};
    FmBlockSet out_b = {};
    if (FmAllocateBlockSet(&un_a, kFixmateNB) != 0 ||
        FmAllocateBlockSet(&un_b, kFixmateNB) != 0 ||
        FmAllocateBlockSet(&out_a, kFixmateNB) != 0 ||
        FmAllocateBlockSet(&out_b, kFixmateNB) != 0) {
        FmFreeBlockSet(&un_a);
        FmFreeBlockSet(&un_b);
        FmFreeBlockSet(&out_a);
        FmFreeBlockSet(&out_b);
        stats->t_compress_setup += GetTime() - workspace_t0;
        return -1;
    }
    stats->t_compress_setup += GetTime() - workspace_t0;

    double pack_t0 = GetTime();
    std::vector<FmPackPlan> plans;
    plans.reserve((records.size() + 1) / 2);
    size_t begin = 0;
    while (begin < records.size()) {
        size_t end = begin;
        uint32_t total = 0;
        while (end < records.size()) {
            const uint32_t packed = bam_lens[end] + 4;
            if (packed > BGZF_BLOCK_SIZE) {
                fprintf(stderr,
                        "[rank %d] ERROR: fixmate record exceeds BGZF payload.\n",
                        rank);
                FmFreeBlockSet(&un_a);
                FmFreeBlockSet(&un_b);
                FmFreeBlockSet(&out_a);
                FmFreeBlockSet(&out_b);
                return -1;
            }
            if (total + packed > BGZF_BLOCK_SIZE) break;
            total += packed;
            ++end;
        }
        FmPackPlan plan;
        plan.records =
            const_cast<bam1_t **>(records.data() + begin);
        plan.n_records = (int)(end - begin);
        plan.total_len = total;
        plans.push_back(plan);
        begin = end;
    }
    stats->t_pack += GetTime() - pack_t0;

    double setup_t0 = GetTime();
    Comp_Para comp_a[kFixmateNB];
    Comp_Para comp_b[kFixmateNB];
    for (int i = 0; i < kFixmateNB; ++i) {
        FmInitEmptyComp(comp_a + i, i);
        FmInitEmptyComp(comp_b + i, i);
    }
    stats->t_compress_setup += GetTime() - setup_t0;
    Comp_Para *active = comp_a;
    Comp_Para *pending = comp_b;
    FmBlockSet *active_un = &un_a;
    FmBlockSet *pending_un = &un_b;
    FmBlockSet *active_out = &out_a;
    FmBlockSet *pending_out = &out_b;
    int pending_count = 0;

    auto flush_pending = [&]() -> int {
        double t0 = GetTime();
        for (int i = 0; i < pending_count; ++i) {
            if (pending[i].status != 0 ||
                !pending[i].output_block ||
                MpiWriteBlockToMem(
                    *writer, pending[i].output_block) != 0) {
                return -1;
            }
            stats->bgzf_blocks++;
        }
        stats->t_write += GetTime() - t0;
        pending_count = 0;
        return 0;
    };

    size_t plan_pos = 0;
    while (plan_pos < plans.size()) {
        const int active_count = (int)std::min(
            (size_t)kFixmateNB, plans.size() - plan_pos);
        for (int i = 0; i < active_count; ++i) {
            memset(active + i, 0, sizeof(active[i]));
            active[i].block_id = i;
            active[i].input_records = plans[plan_pos + i].records;
            active[i].n_records = plans[plan_pos + i].n_records;
            active[i].un_comp_block = active_un->blocks + i;
            active[i].output_block = active_out->blocks + i;
            active[i].compress_level = compress_level;
            active[i].status = 0;
        }
        for (int i = active_count; i < kFixmateNB; ++i) {
            FmInitEmptyComp(active + i, i);
        }
        double t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_compressfunc,
                             active, 1);
        if (flush_pending() != 0) {
            athread_join();
            double free_t0 = GetTime();
            FmFreeBlockSet(&un_a);
            FmFreeBlockSet(&un_b);
            FmFreeBlockSet(&out_a);
            FmFreeBlockSet(&out_b);
            stats->t_workspace_free += GetTime() - free_t0;
            return -1;
        }
        athread_join();
        stats->t_compress += GetTime() - t0;
        for (int i = 0; i < active_count; ++i) {
            if (active[i].status != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: fixmate compression failed "
                        "at block %d status=%d.\n",
                        rank, i, active[i].status);
                double free_t0 = GetTime();
                FmFreeBlockSet(&un_a);
                FmFreeBlockSet(&un_b);
                FmFreeBlockSet(&out_a);
                FmFreeBlockSet(&out_b);
                stats->t_workspace_free += GetTime() - free_t0;
                return -1;
            }
        }
        pending_count = active_count;
        std::swap(active, pending);
        std::swap(active_un, pending_un);
        std::swap(active_out, pending_out);
        plan_pos += active_count;
    }
    const int result = flush_pending();
    double free_t0 = GetTime();
    FmFreeBlockSet(&un_a);
    FmFreeBlockSet(&un_b);
    FmFreeBlockSet(&out_a);
    FmFreeBlockSet(&out_b);
    stats->t_workspace_free += GetTime() - free_t0;
    return result;
}

static int FmProcessGroups(
        std::vector<bam1_t *> records,
        const std::vector<MpiFixmateGroupShared> &groups,
        int compress_level, MemWriter *writer,
        int rank, MpiFixmateStats *stats) {
    if (records.empty()) return 0;
    double alloc_t0 = GetTime();
    std::vector<MpiFixmateRecordPlanShared> plans(records.size());
    stats->t_process_alloc += GetTime() - alloc_t0;

    double setup_t0 = GetTime();
    MpiFixmatePlanPara plan_paras[kFixmateNB];
    for (int c = 0; c < kFixmateNB; ++c) {
        memset(plan_paras + c, 0, sizeof(plan_paras[c]));
        plan_paras[c].records = records.data();
        plan_paras[c].plans = plans.data();
        plan_paras[c].groups = groups.data();
        plan_paras[c].group_begin =
            (int)((long long)groups.size() * c / kFixmateNB);
        plan_paras[c].group_end =
            (int)((long long)groups.size() * (c + 1) / kFixmateNB);
        plan_paras[c].status = 0;
    }
    stats->t_process_setup += GetTime() - setup_t0;

    double t0 = GetTime();
    __real_athread_spawn((void *)slave_mpi_fixmate_plan,
                         plan_paras, 1);
    athread_join();
    stats->t_plan += GetTime() - t0;
    for (int c = 0; c < kFixmateNB; ++c) {
        if (plan_paras[c].status != 0) {
            fprintf(stderr,
                    "[rank %d] ERROR: fixmate plan failed "
                    "status=%d record=%d.\n",
                    rank, plan_paras[c].status,
                    plan_paras[c].record_index);
            return -1;
        }
        stats->paired_groups += plan_paras[c].paired_groups;
        stats->singleton_groups +=
            plan_paras[c].singleton_groups;
        stats->secondary_records +=
            plan_paras[c].secondary_records;
        stats->supplementary_records +=
            plan_paras[c].supplementary_records;
        stats->mq_updates += plan_paras[c].mq_updates;
        stats->mc_updates += plan_paras[c].mc_updates;
        stats->ms_updates += plan_paras[c].ms_updates;
    }
    stats->group_count += (long long)groups.size();

    double output_index_t0 = GetTime();
    std::vector<uint64_t> offsets(records.size() + 1, 0);
    for (size_t i = 0; i < records.size(); ++i) {
        if (plans[i].output_data_len >
            UINT64_MAX - offsets[i]) {
            return -1;
        }
        offsets[i + 1] =
            offsets[i] + plans[i].output_data_len;
    }
    stats->t_output_index += GetTime() - output_index_t0;

    if (offsets.back() > SIZE_MAX) return -1;
    alloc_t0 = GetTime();
    unsigned char *output_data = aligned_alloc_custom(
        64, offsets.back() ? (size_t)offsets.back() : 1);
    bam1_t *output_records =
        (bam1_t *)aligned_alloc_custom(
            64, records.size() * sizeof(bam1_t));
    stats->t_process_alloc += GetTime() - alloc_t0;
    if (!output_data || !output_records) {
        if (output_data) aligned_free_custom(output_data);
        if (output_records) {
            aligned_free_custom(
                (unsigned char *)output_records);
        }
        return -1;
    }

    setup_t0 = GetTime();
    MpiFixmateRewritePara rewrite[kFixmateNB];
    for (int c = 0; c < kFixmateNB; ++c) {
        memset(rewrite + c, 0, sizeof(rewrite[c]));
        rewrite[c].records = records.data();
        rewrite[c].output_records = output_records;
        rewrite[c].plans = plans.data();
        rewrite[c].output_offsets = offsets.data();
        rewrite[c].output_data = output_data;
        rewrite[c].output_capacity = (size_t)offsets.back();
        rewrite[c].record_begin =
            (int)((long long)records.size() * c / kFixmateNB);
        rewrite[c].record_end =
            (int)((long long)records.size() * (c + 1) / kFixmateNB);
        rewrite[c].status = 0;
    }
    stats->t_process_setup += GetTime() - setup_t0;

    t0 = GetTime();
    __real_athread_spawn((void *)slave_mpi_fixmate_rewrite,
                         rewrite, 1);
    athread_join();
    stats->t_rewrite += GetTime() - t0;
    for (int c = 0; c < kFixmateNB; ++c) {
        if (rewrite[c].status != 0) {
            fprintf(stderr,
                    "[rank %d] ERROR: fixmate rewrite failed "
                    "status=%d record=%d.\n",
                    rank, rewrite[c].status,
                    rewrite[c].record_index);
            double free_t0 = GetTime();
            aligned_free_custom(output_data);
            aligned_free_custom(
                (unsigned char *)output_records);
            stats->t_workspace_free += GetTime() - free_t0;
            return -1;
        }
    }

    output_index_t0 = GetTime();
    std::vector<bam1_t *> output_ptrs(records.size());
    std::vector<uint32_t> bam_lens(records.size());
    for (size_t i = 0; i < records.size(); ++i) {
        output_ptrs[i] = output_records + i;
        const uint64_t length =
            (uint64_t)output_records[i].l_data -
            output_records[i].core.l_extranul + 32;
        if (length > UINT32_MAX) {
            double free_t0 = GetTime();
            aligned_free_custom(output_data);
            aligned_free_custom(
                (unsigned char *)output_records);
            stats->t_workspace_free += GetTime() - free_t0;
            return -1;
        }
        bam_lens[i] = (uint32_t)length;
    }
    stats->t_output_index += GetTime() - output_index_t0;

    const int ret = FmCompressPlans(
        output_ptrs, bam_lens, compress_level,
        writer, rank, stats);
    double free_t0 = GetTime();
    aligned_free_custom(output_data);
    aligned_free_custom((unsigned char *)output_records);
    stats->t_workspace_free += GetTime() - free_t0;
    return ret;
}

static int FmProcessStoredGroup(
        FmStoredGroup *group, int compress_level,
        MemWriter *writer, int rank,
        MpiFixmateStats *stats) {
    double setup_t0 = GetTime();
    std::vector<bam1_t> records;
    std::vector<bam1_t *> ptrs;
    if (FmBuildStoredViews(group, &records, &ptrs) != 0) {
        stats->t_process_setup += GetTime() - setup_t0;
        return -1;
    }
    MpiFixmateGroupShared descriptor;
    descriptor.begin = 0;
    descriptor.count = (int)ptrs.size();
    std::vector<MpiFixmateGroupShared> groups(1, descriptor);
    stats->t_process_setup += GetTime() - setup_t0;
    return FmProcessGroups(
        std::move(ptrs), groups, compress_level, writer, rank, stats);
}

static int FmProcessDirectRange(
        const std::vector<bam1_t *> &all_records,
        const std::vector<std::pair<size_t, size_t> > &ranges,
        size_t range_begin, size_t range_end,
        int compress_level, MemWriter *writer,
        int rank, MpiFixmateStats *stats) {
    if (range_begin >= range_end) return 0;
    double setup_t0 = GetTime();
    const size_t record_begin = ranges[range_begin].first;
    const size_t record_end = ranges[range_end - 1].second;
    std::vector<bam1_t *> records(
        all_records.begin() + record_begin,
        all_records.begin() + record_end);
    std::vector<MpiFixmateGroupShared> groups;
    groups.reserve(range_end - range_begin);
    for (size_t i = range_begin; i < range_end; ++i) {
        MpiFixmateGroupShared group;
        group.begin =
            (int)(ranges[i].first - record_begin);
        group.count =
            (int)(ranges[i].second - ranges[i].first);
        groups.push_back(group);
    }
    stats->t_process_setup += GetTime() - setup_t0;
    return FmProcessGroups(
        std::move(records), groups, compress_level, writer, rank, stats);
}

static void FmBuildRanges(
        const std::vector<bam1_t *> &records, size_t begin,
        std::vector<std::pair<size_t, size_t> > *ranges) {
    ranges->clear();
    size_t pos = begin;
    while (pos < records.size()) {
        size_t end = pos + 1;
        const char *qname = bam_get_qname(records[pos]);
        while (end < records.size() &&
               strcmp(qname, bam_get_qname(records[end])) == 0) {
            ++end;
        }
        ranges->push_back(std::make_pair(pos, end));
        pos = end;
    }
}

static int FmStreamLocalGroups(
        MemReader &reader, int compress_level,
        MemWriter *middle_writer,
        FmStoredGroup *leading,
        FmStoredGroup *trailing,
        int *single_group,
        int rank, MpiFixmateStats *stats) {
    const int records_per_block = (int)MPI_RECORDS_PER_BLOCK;
    FmBlockSet input = {};
    FmBlockSet uncompressed = {};
    FmRecordSet record_set = {};
    Bam2BamPara decomp[kFixmateNB];
    if (FmAllocateBlockSet(&input, kFixmateNB) != 0 ||
        FmAllocateBlockSet(&uncompressed, kFixmateNB) != 0 ||
        FmAllocateRecordSet(
            &record_set, kFixmateNB, records_per_block,
            MPI_BAM_BLOCK_ARENA_SIZE) != 0) {
        FmFreeBlockSet(&input);
        FmFreeBlockSet(&uncompressed);
        FmFreeRecordSet(&record_set);
        return -1;
    }

    FmStoredGroup pending;
    bool leading_ready = false;
    int n_blocks = 0;
    FmReadGroup(reader, &input, &n_blocks, stats);
    while (n_blocks > 0) {
        for (int b = 0; b < kFixmateNB; ++b) {
            FmInitDecompPara(
                decomp + b, b, &input, &uncompressed,
                &record_set, b < n_blocks);
        }
        double t0 = GetTime();
        __real_athread_spawn(
            (void *)slave_mpi_decompress_bam2bam_passthrough,
            decomp, 1);
        athread_join();
        stats->t_decompress += GetTime() - t0;

        double collect_t0 = GetTime();
        std::vector<bam1_t *> records;
        for (int b = 0; b < n_blocks; ++b) {
            if (decomp[b].status != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: fixmate decode failed "
                        "block=%d status=%d record=%d "
                        "actual=%lld limit=%lld limit_id=%d.\n",
                        rank, b, decomp[b].status,
                        decomp[b].record_index,
                        decomp[b].actual_value,
                        decomp[b].limit_value,
                        decomp[b].limit_id);
                FmFreeBlockSet(&input);
                FmFreeBlockSet(&uncompressed);
                FmFreeRecordSet(&record_set);
                return -1;
            }
            stats->input_blocks++;
            stats->total_records += decomp[b].n_total_records;
            records.insert(
                records.end(), decomp[b].output_records,
                decomp[b].output_records +
                    decomp[b].n_total_records);
        }
        stats->t_record_collect += GetTime() - collect_t0;

        t0 = GetTime();
        size_t cursor = 0;
        if (!pending.empty()) {
            while (cursor < records.size() &&
                   pending.qname ==
                       bam_get_qname(records[cursor])) {
                if (FmAppendStoredRecord(
                        &pending, records[cursor]) != 0) {
                    FmFreeBlockSet(&input);
                    FmFreeBlockSet(&uncompressed);
                    FmFreeRecordSet(&record_set);
                    return -1;
                }
                ++cursor;
            }
            if (cursor < records.size()) {
                if (!leading_ready) {
                    *leading = std::move(pending);
                    leading_ready = true;
                } else {
                    stats->t_group_scan += GetTime() - t0;
                    if (FmProcessStoredGroup(
                            &pending, compress_level,
                            middle_writer, rank, stats) != 0) {
                        FmFreeBlockSet(&input);
                        FmFreeBlockSet(&uncompressed);
                        FmFreeRecordSet(&record_set);
                        return -1;
                    }
                    t0 = GetTime();
                }
                pending.clear();
            }
        }

        if (cursor < records.size()) {
            std::vector<std::pair<size_t, size_t> > ranges;
            FmBuildRanges(records, cursor, &ranges);
            const size_t last = ranges.size() - 1;
            size_t direct_begin = 0;
            if (!leading_ready && last > 0) {
                if (FmAppendStoredRange(
                        leading, records,
                        ranges[0].first,
                        ranges[0].second) != 0) {
                    FmFreeBlockSet(&input);
                    FmFreeBlockSet(&uncompressed);
                    FmFreeRecordSet(&record_set);
                    return -1;
                }
                leading_ready = true;
                direct_begin = 1;
            }
            if (direct_begin < last) {
                stats->t_group_scan += GetTime() - t0;
                if (FmProcessDirectRange(
                        records, ranges, direct_begin, last,
                        compress_level, middle_writer,
                        rank, stats) != 0) {
                    FmFreeBlockSet(&input);
                    FmFreeBlockSet(&uncompressed);
                    FmFreeRecordSet(&record_set);
                    return -1;
                }
                t0 = GetTime();
            }
            if (FmAppendStoredRange(
                    &pending, records,
                    ranges[last].first,
                    ranges[last].second) != 0) {
                FmFreeBlockSet(&input);
                FmFreeBlockSet(&uncompressed);
                FmFreeRecordSet(&record_set);
                return -1;
            }
        }
        stats->t_group_scan += GetTime() - t0;
        FmReadGroup(reader, &input, &n_blocks, stats);
    }

    if (!pending.empty()) {
        if (!leading_ready) {
            *leading = std::move(pending);
            *single_group = 1;
        } else {
            *trailing = std::move(pending);
            *single_group = 0;
        }
    } else {
        *single_group = leading_ready ? 1 : 0;
    }
    FmFreeBlockSet(&input);
    FmFreeBlockSet(&uncompressed);
    FmFreeRecordSet(&record_set);
    return 0;
}

static int FmPreviousNonempty(
        const std::vector<FmBoundaryInfo> &infos, int rank) {
    for (int r = rank - 1; r >= 0; --r) {
        if (infos[(size_t)r].has_records) return r;
    }
    return -1;
}

static int FmBoundaryOwner(
        const std::vector<FmBoundaryInfo> &infos, int rank) {
    if (!infos[(size_t)rank].has_records) return rank;
    const char *qname = infos[(size_t)rank].first_qname;
    int previous = FmPreviousNonempty(infos, rank);
    if (previous < 0 ||
        strcmp(infos[(size_t)previous].last_qname, qname) != 0) {
        return rank;
    }
    int owner = previous;
    while (infos[(size_t)owner].single_group) {
        previous = FmPreviousNonempty(infos, owner);
        if (previous < 0 ||
            strcmp(infos[(size_t)previous].last_qname, qname) != 0) {
            break;
        }
        owner = previous;
    }
    return owner;
}

static int FmExchangeBoundaryGroups(
        FmStoredGroup *leading, FmStoredGroup *trailing,
        int single_group, int rank, int comm_size,
        int *continues_previous,
        MpiFixmateStats *stats) {
    FmBoundaryInfo local;
    memset(&local, 0, sizeof(local));
    local.has_records = !leading->empty();
    local.single_group = single_group;
    if (local.has_records) {
        snprintf(local.first_qname, sizeof(local.first_qname),
                 "%s", leading->qname.c_str());
        const std::string &last =
            single_group ? leading->qname : trailing->qname;
        snprintf(local.last_qname, sizeof(local.last_qname),
                 "%s", last.c_str());
    }
    std::vector<FmBoundaryInfo> infos((size_t)comm_size);
    MPI_Allgather(&local, sizeof(local), MPI_BYTE,
                  infos.data(), sizeof(local), MPI_BYTE,
                  MPI_COMM_WORLD);
    std::vector<int> owners((size_t)comm_size);
    for (int r = 0; r < comm_size; ++r) {
        owners[(size_t)r] = FmBoundaryOwner(infos, r);
    }
    *continues_previous =
        local.has_records && owners[(size_t)rank] != rank;

    int local_ok = 1;
    std::vector<unsigned char> local_wire;
    if (*continues_previous &&
        FmSerializeStoredGroup(*leading, &local_wire) != 0) {
        local_ok = 0;
    }
    if (!FmAllRanksOk(local_ok)) return -1;
    unsigned long long local_size = local_wire.size();
    std::vector<unsigned long long> sizes((size_t)comm_size, 0);
    MPI_Allgather(&local_size, 1, MPI_UNSIGNED_LONG_LONG,
                  sizes.data(), 1, MPI_UNSIGNED_LONG_LONG,
                  MPI_COMM_WORLD);

    double t0 = GetTime();
    int received_fragment = 0;
    for (int source = 0; source < comm_size; ++source) {
        const int owner = owners[(size_t)source];
        if (owner == source) continue;
        const unsigned long long size = sizes[(size_t)source];
        if (rank == source) {
            unsigned long long sent = 0;
            while (sent < size) {
                const int chunk = (int)std::min(
                    (unsigned long long)kFixmateExchangeChunk,
                    size - sent);
                MPI_Send(local_wire.data() + sent, chunk, MPI_BYTE,
                         owner, 5100 + source,
                         MPI_COMM_WORLD);
                sent += (unsigned long long)chunk;
            }
            stats->boundary_bytes += (long long)size;
        } else if (rank == owner) {
            std::vector<unsigned char> wire((size_t)size);
            unsigned long long received = 0;
            while (received < size) {
                const int chunk = (int)std::min(
                    (unsigned long long)kFixmateExchangeChunk,
                    size - received);
                MPI_Recv(wire.data() + received, chunk, MPI_BYTE,
                         source, 5100 + source,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                received += (unsigned long long)chunk;
            }
            FmStoredGroup *target =
                single_group ? leading : trailing;
            if (FmAppendSerializedGroup(
                    target, wire.data(), wire.size()) != 0) {
                local_ok = 0;
            }
            received_fragment = 1;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (received_fragment) stats->boundary_groups++;
    stats->t_boundary_exchange += GetTime() - t0;
    return FmAllRanksOk(local_ok) ? 0 : -1;
}

static int FmPrependAppend(
        MemWriter *middle, const MemWriter &prefix,
        const MemWriter &suffix) {
    if (prefix.size > SIZE_MAX - middle->size ||
        prefix.size + middle->size > SIZE_MAX - suffix.size) {
        return -1;
    }
    const size_t total = prefix.size + middle->size + suffix.size;
    if (total > middle->capacity) {
        char *new_data = (char *)realloc(middle->data, total);
        if (!new_data && total > 0) return -1;
        middle->data = new_data;
        middle->capacity = total;
    }
    if (prefix.size > 0 && middle->size > 0) {
        memmove(middle->data + prefix.size,
                middle->data, middle->size);
    }
    if (prefix.size > 0) {
        memcpy(middle->data, prefix.data, prefix.size);
    }
    if (suffix.size > 0) {
        memcpy(middle->data + prefix.size + middle->size,
               suffix.data, suffix.size);
    }
    middle->size = total;
    return 0;
}

static void FmPrintStats(
        const MpiFixmateStats &stats, int rank, int comm_size);

static int FmGatherOutput(
        const MemWriter &local_writer, sam_hdr_t *header,
        int compress_level, const std::string &output_path,
        MpiMemoryBam *output_bam,
        int rank, int comm_size, MpiFixmateStats *stats,
        FmOutputStageCosts *costs,
        double stage41_cost, double stage42_cost,
        double stage43_cost, double body_total_t0) {
    memset(costs, 0, sizeof(*costs));

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
    costs->stage44 = FmReduceMax(GetTime() - stage44_t0);
    if (rank == 0) {
        printf("Complete the 4.4 gather body sizes cost %lf\n",
               costs->stage44);
    }

    char *header_memory = nullptr;
    size_t header_size = 0;
    char *output_memory = nullptr;
    size_t output_size = 0;
    int local_ok = 1;

    double stage45_header_t0 = GetTime();
    if (rank == 0) {
        if (sam_hdr_add_pg(header, "RabbitBAM-MPI",
                           "PN", "RabbitBAM-MPI",
                           "CL", "fixmate -m",
                           NULL) != 0) {
            fprintf(stderr,
                    "ERROR: failed to add fixmate @PG header line.\n");
            local_ok = 0;
        }
        if (local_ok &&
            MpiCommonBuildBamHeaderMemory(
                header, compress_level,
                &header_memory, &header_size) != 0) {
            local_ok = 0;
        }
        const unsigned long long final_size =
            (unsigned long long)header_size +
            (unsigned long long)total_body +
            sizeof(kFixmateBgzfEofBlock);
        if (local_ok && final_size > SIZE_MAX) local_ok = 0;
        if (local_ok) {
            output_size = (size_t)final_size;
        }
    }
    costs->stage45_header =
        FmReduceMax(GetTime() - stage45_header_t0);
    if (rank == 0 && local_ok) {
        printf("Complete the 4.5a header/layout cost %lf\n",
               costs->stage45_header);
    }
    if (!FmAllRanksOk(local_ok)) {
        free(header_memory);
        return -1;
    }

    double stage45_malloc_t0 = GetTime();
    const size_t simulated_write_size =
        local_writer.size +
        (rank == 0 ? header_size + sizeof(kFixmateBgzfEofBlock) : 0);
    char *simulated_write_mem =
        simulated_write_size ? (char *)malloc(simulated_write_size)
                             : nullptr;
    if (simulated_write_size && !simulated_write_mem) {
        local_ok = 0;
    }
    costs->stage45_malloc =
        FmReduceMax(GetTime() - stage45_malloc_t0);
    if (rank == 0 && local_ok) {
        printf("Complete the 4.5b malloc simulated write memory cost %lf\n",
               costs->stage45_malloc);
    }
    if (!FmAllRanksOk(local_ok)) {
        free(header_memory);
        free(simulated_write_mem);
        return -1;
    }

    double stage46_t0 = GetTime();
    volatile unsigned long long simulated_guard = 0;
    size_t simulated_pos = 0;
    if (rank == 0 && header_size > 0) {
        memcpy(simulated_write_mem + simulated_pos,
               header_memory, header_size);
        simulated_pos += header_size;
    }
    if (local_writer.size > 0) {
        memcpy(simulated_write_mem + simulated_pos,
               local_writer.data, local_writer.size);
        simulated_pos += local_writer.size;
    }
    if (rank == 0) {
        memcpy(simulated_write_mem + simulated_pos,
               kFixmateBgzfEofBlock,
               sizeof(kFixmateBgzfEofBlock));
        simulated_pos += sizeof(kFixmateBgzfEofBlock);
    }
    if (simulated_write_size > 0) {
        simulated_guard +=
            (unsigned char)simulated_write_mem[0];
        simulated_guard +=
            (unsigned char)simulated_write_mem[simulated_write_size - 1];
    }
    (void)simulated_guard;
    costs->stage46 = FmReduceMax(GetTime() - stage46_t0);
    if (rank == 0) {
        printf("Complete the 4.6 simulated header/body distributed write cost %lf\n",
               costs->stage46);
    }
    free(simulated_write_mem);

    if (rank == 0) {
        const double core_total =
            stage41_cost + stage42_cost + stage43_cost +
            costs->stage44 + costs->stage45_header +
            costs->stage46;
        printf("Complete the total (4.1~4.6) cost %lf-----\n",
               core_total);
    }

    {
        double stage47_t0 = GetTime();
        FmPrintStats(*stats, rank, comm_size);
        const double stage47_cost =
            FmReduceMax(GetTime() - stage47_t0);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.7 rank stats reduce/print cost %lf\n",
                   stage47_cost);
        }
    }

    {
        const double body_total_cost =
            FmReduceMax(GetTime() - body_total_t0);
        if (rank == 0 && local_ok) {
            printf("444Complete the total body cost %lf\n",
                   body_total_cost);
        }
    }

    double verify_alloc_t0 = GetTime();
    if (rank == 0) {
        output_memory =
            output_size ? (char *)malloc(output_size) : nullptr;
        if (output_size && !output_memory) local_ok = 0;
    }
    double verify_alloc_cost =
        FmReduceMax(GetTime() - verify_alloc_t0);
    if (rank == 0 && local_ok) {
        printf("555Prepare verification output memory cost %lf--\n",
               verify_alloc_cost);
    }
    if (!FmAllRanksOk(local_ok)) {
        free(header_memory);
        free(output_memory);
        return -1;
    }

    double gather_t0 = GetTime();
    if (rank == 0) {
        memcpy(output_memory, header_memory, header_size);
        if (local_writer.size) {
            memcpy(output_memory + header_size,
                   local_writer.data, local_writer.size);
        }
        for (int source = 1; source < comm_size; ++source) {
            unsigned long long received = 0;
            const unsigned long long bytes =
                (unsigned long long)sizes[(size_t)source];
            while (received < bytes) {
                const int chunk = (int)std::min(
                    (unsigned long long)kFixmateExchangeChunk,
                    bytes - received);
                MPI_Recv(output_memory + header_size +
                             offsets[(size_t)source] + received,
                         chunk, MPI_BYTE, source, 5200,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                received += (unsigned long long)chunk;
            }
        }
        memcpy(output_memory + header_size + total_body,
               kFixmateBgzfEofBlock,
               sizeof(kFixmateBgzfEofBlock));
    } else {
        unsigned long long sent = 0;
        const unsigned long long bytes = local_writer.size;
        while (sent < bytes) {
            const int chunk = (int)std::min(
                (unsigned long long)kFixmateExchangeChunk,
                bytes - sent);
            MPI_Send(local_writer.data + sent, chunk, MPI_BYTE,
                     0, 5200, MPI_COMM_WORLD);
            sent += (unsigned long long)chunk;
        }
    }
    const double gather_cost = FmReduceMax(GetTime() - gather_t0);
    if (rank == 0 && local_ok) {
        printf("555Gather verification output memory cost %lf--\n",
               gather_cost);
    }

    double dump_t0 = GetTime();
    if (output_bam) {
        if (rank == 0) {
            output_bam->data = output_memory;
            output_bam->size = output_size;
            output_memory = nullptr;
        }
    } else if (rank == 0 &&
        MpiCommonDumpMemoryToFile(
            output_path, output_memory, output_size) != 0) {
        fprintf(stderr,
                "ERROR: failed to write fixmate output %s.\n",
                output_path.c_str());
        local_ok = 0;
    }
    const double dump_cost = FmReduceMax(GetTime() - dump_t0);
    if (rank == 0 && local_ok) {
        printf(output_bam
                   ? "555Keep output memory cost %lf--\n"
                   : "555Dump memory to output file cost %lf--\n",
               dump_cost);
    }
    free(header_memory);
    free(output_memory);
    return FmAllRanksOk(local_ok) ? 0 : -1;
}

static double FmAccountedFusedTime(const MpiFixmateStats &stats) {
    return stats.t_read +
           stats.t_decompress +
           stats.t_record_collect +
           stats.t_group_scan +
           stats.t_process_setup +
           stats.t_process_alloc +
           stats.t_plan +
           stats.t_output_index +
           stats.t_boundary_exchange +
           stats.t_rewrite +
           stats.t_pack +
           stats.t_compress_setup +
           stats.t_compress +
           stats.t_rank_sync +
           stats.t_workspace_free;
}

static void FmPrintStats(
        const MpiFixmateStats &stats, int rank, int comm_size) {
    const int count = 13;
    long long local[count] = {
        stats.input_blocks,
        stats.group_count,
        stats.total_records,
        stats.paired_groups,
        stats.singleton_groups,
        stats.secondary_records,
        stats.supplementary_records,
        stats.boundary_groups,
        stats.boundary_bytes,
        stats.mq_updates,
        stats.mc_updates,
        stats.ms_updates,
        stats.bgzf_blocks
    };
    const double accounted = FmAccountedFusedTime(stats);
    const double unaccounted = stats.t_fused_total - accounted;

    const int kLineBytes = 8192;
    char local_lines[kLineBytes];
    memset(local_lines, 0, sizeof(local_lines));
    int used = snprintf(
        local_lines, sizeof(local_lines),
        "[rank %d] fixmate_stats blocks=%lld groups=%lld "
        "records=%lld paired=%lld singleton=%lld "
        "secondary=%lld supplementary=%lld boundary_groups=%lld "
        "boundary_bytes=%lld MQ=%lld MC=%lld ms=%lld bgzf=%lld\n",
        rank, stats.input_blocks, stats.group_count,
        stats.total_records, stats.paired_groups,
        stats.singleton_groups, stats.secondary_records,
        stats.supplementary_records, stats.boundary_groups,
        stats.boundary_bytes, stats.mq_updates,
        stats.mc_updates, stats.ms_updates,
        stats.bgzf_blocks);
    if (used < 0) used = 0;
    if (used >= kLineBytes) used = kLineBytes - 1;
    snprintf(
        local_lines + used, (size_t)(kLineBytes - used),
        "[rank %d] fixmate_timing read=%.6f decompress=%.6f "
        "record_collect=%.6f group_scan=%.6f "
        "process_setup=%.6f process_alloc=%.6f plan=%.6f "
        "output_index=%.6f boundary=%.6f rewrite=%.6f "
        "pack=%.6f compress_setup=%.6f compress=%.6f "
        "write=%.6f rank_sync=%.6f workspace_free=%.6f "
        "accounted=%.6f unaccounted=%.6f fused=%.6f\n",
        rank, stats.t_read, stats.t_decompress,
        stats.t_record_collect, stats.t_group_scan,
        stats.t_process_setup, stats.t_process_alloc,
        stats.t_plan, stats.t_output_index,
        stats.t_boundary_exchange, stats.t_rewrite,
        stats.t_pack, stats.t_compress_setup,
        stats.t_compress, stats.t_write,
        stats.t_rank_sync, stats.t_workspace_free,
        accounted, unaccounted, stats.t_fused_total);

    std::vector<char> gathered_lines;
    if (rank == 0) {
        gathered_lines.resize((size_t)comm_size * kLineBytes);
    }
    MPI_Gather(local_lines, kLineBytes, MPI_CHAR,
               rank == 0 ? gathered_lines.data() : nullptr,
               kLineBytes, MPI_CHAR, 0, MPI_COMM_WORLD);

    long long sums[count] = {};
    MPI_Reduce(local, sums, count, MPI_LONG_LONG,
               MPI_SUM, 0, MPI_COMM_WORLD);
    const int time_count = 19;
    double times[time_count] = {
        stats.t_read,
        stats.t_decompress,
        stats.t_record_collect,
        stats.t_group_scan,
        stats.t_process_setup,
        stats.t_process_alloc,
        stats.t_plan,
        stats.t_output_index,
        stats.t_boundary_exchange,
        stats.t_rewrite,
        stats.t_pack,
        stats.t_compress_setup,
        stats.t_compress,
        stats.t_write,
        stats.t_rank_sync,
        stats.t_workspace_free,
        accounted,
        unaccounted,
        stats.t_fused_total
    };
    double time_sums[time_count] = {};
    double time_max[time_count] = {};
    MPI_Reduce(times, time_sums, time_count, MPI_DOUBLE,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(times, time_max, time_count, MPI_DOUBLE,
               MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int r = 0; r < comm_size; ++r) {
            fputs(gathered_lines.data() + (size_t)r * kLineBytes,
                  stdout);
        }
        printf("FusedFixmateMPI finished. ranks=%d blocks=%lld "
               "records=%lld groups=%lld paired=%lld singleton=%lld\n",
               comm_size, sums[0], sums[2], sums[1],
               sums[3], sums[4]);
        printf("  secondary=%lld supplementary=%lld "
               "boundary_groups=%lld boundary_bytes=%lld\n",
               sums[5], sums[6], sums[7], sums[8]);
        printf("  tag_updates MQ=%lld MC=%lld ms=%lld "
               "bgzf_blocks=%lld\n",
               sums[9], sums[10], sums[11], sums[12]);
        printf("  timing_sum read=%.3f decompress=%.3f "
               "record_collect=%.3f group_scan=%.3f "
               "process_setup=%.3f process_alloc=%.3f "
               "plan=%.3f output_index=%.3f boundary=%.3f "
               "rewrite=%.3f pack=%.3f compress_setup=%.3f "
               "compress=%.3f write=%.3f rank_sync=%.3f "
               "workspace_free=%.3f accounted=%.3f "
               "unaccounted=%.3f fused=%.3f\n",
               time_sums[0], time_sums[1], time_sums[2],
               time_sums[3], time_sums[4], time_sums[5],
               time_sums[6], time_sums[7], time_sums[8],
               time_sums[9], time_sums[10], time_sums[11],
               time_sums[12], time_sums[13], time_sums[14],
               time_sums[15], time_sums[16], time_sums[17],
               time_sums[18]);
        printf("  timing_max read=%.3f decompress=%.3f "
               "record_collect=%.3f group_scan=%.3f "
               "process_setup=%.3f process_alloc=%.3f "
               "plan=%.3f output_index=%.3f boundary=%.3f "
               "rewrite=%.3f pack=%.3f compress_setup=%.3f "
               "compress=%.3f write=%.3f rank_sync=%.3f "
               "workspace_free=%.3f accounted=%.3f "
               "unaccounted=%.3f fused=%.3f\n",
               time_max[0], time_max[1], time_max[2],
               time_max[3], time_max[4], time_max[5],
               time_max[6], time_max[7], time_max[8],
               time_max[9], time_max[10], time_max[11],
               time_max[12], time_max[13], time_max[14],
               time_max[15], time_max[16], time_max[17],
               time_max[18]);
        fflush(stdout);
    }
}

} // namespace

int MpiFixmateMemoryToMemory(CmdInfo *cmd_info,
                             const char *input_memory_const,
                             size_t input_size,
                             MpiMemoryBam *output_bam,
                             double *core_cost) {
    double total_t0 = GetTime();
    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int exit_code = 1;
    int local_ok = 1;
    char *input_memory = const_cast<char *>(input_memory_const);
    sam_hdr_t *header = nullptr;
    long long body_start = 0;
    std::vector<long long> block_offsets;
    std::vector<long long> block_lengths;
    long long n_blocks = 0;
    long long local_begin = 0;
    long long local_end = 0;
    char *rank_input = nullptr;
    size_t rank_input_size = 0;
    MemReader reader = {};
    MemWriter middle = {};
    MemWriter prefix = {};
    MemWriter suffix = {};
    FmStoredGroup leading;
    FmStoredGroup trailing;
    int single_group = 0;
    int continues_previous = 0;
    MpiFixmateStats stats = {};
    FmOutputStageCosts output_stage_costs = {};
    double stage41_cost = 0.0;
    double stage42_cost = 0.0;
    double stage43_cost = 0.0;
    double body_total_t0 = 0.0;

    if (output_bam) {
        output_bam->data = nullptr;
        output_bam->size = 0;
    }
    if (core_cost) *core_cost = 0.0;

    if (!input_memory || input_size == 0 || !output_bam ||
        !cmd_info->fixmate_mate_score_) {
        if (rank == 0) {
            fprintf(stderr,
                    "ERROR: dedup-pipeline fixmate requires non-empty BAM input and -m mode.\n");
        }
        local_ok = 0;
    }
    if (!FmAllRanksOk(local_ok)) goto cleanup;
    {
        const double cost = FmReduceMax(GetTime() - total_t0);
        if (rank == 0) {
            printf("111Complete the initialization cost %lf-----\n",
                   cost);
        }
    }

    {
        double t0 = GetTime();
        if (rank == 0 &&
            MpiCommonReadBamHeaderFromMemory(
                input_memory, input_size, &header,
                &body_start) != 0) {
            fprintf(stderr,
                    "ERROR: failed to parse fixmate pipeline BAM header.\n");
            local_ok = 0;
        }
        if (rank == 0 && local_ok) {
            kstring_t order = {0, 0, nullptr};
            if (sam_hdr_find_tag_hd(header, "SO", &order) == 0 &&
                order.s && strcmp(order.s, "coordinate") == 0) {
                if (rank == 0) {
                    fprintf(stderr,
                            "ERROR: fixmate requires name-collated or "
                            "queryname-sorted input, not coordinate order.\n");
                }
                local_ok = 0;
            }
            free(order.s);
        }
        MPI_Bcast(&local_ok, 1, MPI_INT, 0,
                  MPI_COMM_WORLD);
        MPI_Bcast(&body_start, 1, MPI_LONG_LONG, 0,
                  MPI_COMM_WORLD);
        const double cost = FmReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("333Complete the head cost %lf---\n", cost);
        }
    }
    if (!FmAllRanksOk(local_ok)) goto cleanup;

    body_total_t0 = GetTime();
    {
        double t0 = GetTime();
        if (rank == 0) {
            printf("Enable MPI BAM FIXMATE -m mode "
                   "(%d MPE + %d CPEs)!!!\n",
                   comm_size, comm_size * 64);
            printf("MPI BAM output compression level=%d\n",
                   cmd_info->compress_level_);
            if (MpiCommonScanBgzfBlocksInMemory(
                    input_memory, input_size, body_start,
                    &block_offsets, &block_lengths) != 0) {
                local_ok = 0;
            }
            n_blocks = (long long)block_offsets.size();
            if (n_blocks > INT_MAX) local_ok = 0;
        }
        MPI_Bcast(&local_ok, 1, MPI_INT, 0,
                  MPI_COMM_WORLD);
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
        local_begin = n_blocks * rank / comm_size;
        local_end = n_blocks * (rank + 1) / comm_size;
        if (MpiCommonSelectBlockRangeFromMemory(
                input_memory, input_size,
                block_offsets, block_lengths,
                local_begin, local_end,
                &rank_input, &rank_input_size) != 0) {
            local_ok = 0;
        }
        stage41_cost = FmReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("MPI BAM scan complete. data_blocks=%lld "
                   "body_start=%lld\n", n_blocks, body_start);
            printf("Complete the 4.1 scan/split/broadcast/select cost %lf\n",
                   stage41_cost);
        }
    }
    if (!FmAllRanksOk(local_ok)) goto cleanup;

    {
        double t0 = GetTime();
        reader.base = rank_input;
        reader.size = rank_input_size;
        reader.pos = 0;
        if (MpiCommonInitMemWriter(
                middle, rank_input_size ?
                    rank_input_size : 1024 * 1024) != 0 ||
            MpiCommonInitMemWriter(prefix, 1024 * 1024) != 0 ||
            MpiCommonInitMemWriter(suffix, 1024 * 1024) != 0) {
            local_ok = 0;
        }
        stage42_cost = FmReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.2 init reader/writer cost %lf\n",
                   stage42_cost);
        }
    }
    if (!FmAllRanksOk(local_ok)) goto cleanup;

    {
        double fused_t0 = GetTime();
        if (FmStreamLocalGroups(
                reader, cmd_info->compress_level_,
                &middle, &leading, &trailing,
                &single_group, rank, &stats) != 0) {
            local_ok = 0;
        }
        if (!FmAllRanksOkTimed(local_ok, &stats)) goto cleanup;

        if (FmExchangeBoundaryGroups(
                &leading, &trailing, single_group,
                rank, comm_size, &continues_previous,
                &stats) != 0) {
            local_ok = 0;
        }
        if (!FmAllRanksOkTimed(local_ok, &stats)) goto cleanup;

        if (!leading.empty()) {
            if (single_group) {
                if (!continues_previous &&
                    FmProcessStoredGroup(
                        &leading, cmd_info->compress_level_,
                        &prefix, rank, &stats) != 0) {
                    local_ok = 0;
                }
            } else {
                if (!continues_previous &&
                    FmProcessStoredGroup(
                        &leading, cmd_info->compress_level_,
                        &prefix, rank, &stats) != 0) {
                    local_ok = 0;
                }
                if (local_ok &&
                    FmProcessStoredGroup(
                        &trailing, cmd_info->compress_level_,
                        &suffix, rank, &stats) != 0) {
                    local_ok = 0;
                }
            }
        }
        if (local_ok &&
            FmPrependAppend(&middle, prefix, suffix) != 0) {
            local_ok = 0;
        }
        if (!FmAllRanksOkTimed(local_ok, &stats)) goto cleanup;
        stats.t_fused_total = GetTime() - fused_t0;
        stage43_cost = FmReduceMax(stats.t_fused_total);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.3 FusedFixmateMPI cost %lf\n",
                   stage43_cost);
        }
    }
    if (!FmAllRanksOk(local_ok)) goto cleanup;

    if (FmGatherOutput(
            middle, header, cmd_info->compress_level_,
            cmd_info->out_file_name_, output_bam,
            rank, comm_size, &stats, &output_stage_costs,
            stage41_cost, stage42_cost, stage43_cost,
            body_total_t0) != 0) {
        local_ok = 0;
    }
    if (!FmAllRanksOk(local_ok)) goto cleanup;
    if (core_cost) {
        *core_cost = stage41_cost + stage42_cost + stage43_cost +
                     output_stage_costs.stage44 +
                     output_stage_costs.stage45_header +
                     output_stage_costs.stage46;
    }
    exit_code = 0;

cleanup:
    if (header) sam_hdr_destroy(header);
    free(middle.data);
    free(prefix.data);
    free(suffix.data);
    if (rank == 0) {
        printf("666fixmate total process cost %lf-----\n",
               GetTime() - total_t0);
    }
    return exit_code;
}

int ProcessFixmateMPI(CmdInfo *cmd_info) {
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
    long long local_begin = 0;
    long long local_end = 0;
    char *rank_input = nullptr;
    size_t rank_input_size = 0;
    MemReader reader = {};
    MemWriter middle = {};
    MemWriter prefix = {};
    MemWriter suffix = {};
    FmStoredGroup leading;
    FmStoredGroup trailing;
    int single_group = 0;
    int continues_previous = 0;
    MpiFixmateStats stats = {};
    FmOutputStageCosts output_stage_costs = {};
    double stage41_cost = 0.0;
    double stage42_cost = 0.0;
    double stage43_cost = 0.0;
    double body_total_t0 = 0.0;

    if (!cmd_info->fixmate_mate_score_) {
        if (rank == 0) {
            fprintf(stderr,
                    "ERROR: RabbitBAM-MPI fixmate v1 requires -m.\n");
        }
        local_ok = 0;
    }
    if (!FmAllRanksOk(local_ok)) goto cleanup;
    {
        const double cost = FmReduceMax(GetTime() - total_t0);
        if (rank == 0) {
            printf("111Complete the initialization cost %lf-----\n",
                   cost);
        }
    }

    {
        double t0 = GetTime();
        if (MpiCommonLoadFileToMemory(
                cmd_info->in_file_name_,
                &input_memory, &input_size) != 0) {
            fprintf(stderr,
                    "[rank %d] ERROR: cannot preload fixmate input %s.\n",
                    rank, cmd_info->in_file_name_.c_str());
            local_ok = 0;
        }
        const double cost = FmReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("222Complete the memory cost %lf--\n", cost);
        }
    }
    if (!FmAllRanksOk(local_ok)) goto cleanup;

    {
        double t0 = GetTime();
        input_hfile = hopen("mem:", "rb:", input_memory, input_size);
        if (input_hfile) {
            input = (samFile *)hts_hopen(
                input_hfile, "data", "rb");
            if (input) input_hfile = nullptr;
        }
        if (!input) local_ok = 0;
        if (local_ok) header = sam_hdr_read(input);
        if (!header) local_ok = 0;
        if (local_ok && input->format.format != bam &&
            input->format.format != binary_format) {
            if (rank == 0) {
                fprintf(stderr,
                        "ERROR: RabbitBAM-MPI fixmate only supports BAM input.\n");
            }
            local_ok = 0;
        }
        if (local_ok) {
            kstring_t order = {0, 0, nullptr};
            if (sam_hdr_find_tag_hd(header, "SO", &order) == 0 &&
                order.s && strcmp(order.s, "coordinate") == 0) {
                if (rank == 0) {
                    fprintf(stderr,
                            "ERROR: fixmate requires name-collated or "
                            "queryname-sorted input, not coordinate order.\n");
                }
                local_ok = 0;
            }
            free(order.s);
        }
        if (local_ok) {
            body_start =
                (long long)input->fp.bgzf->block_address;
            if (body_start < 0 ||
                (unsigned long long)body_start > input_size) {
                local_ok = 0;
            }
        }
        const double cost = FmReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("333Complete the head cost %lf---\n", cost);
        }
    }
    if (!FmAllRanksOk(local_ok)) goto cleanup;

    body_total_t0 = GetTime();
    {
        double t0 = GetTime();
        if (rank == 0) {
            printf("Enable MPI BAM FIXMATE -m mode "
                   "(%d MPE + %d CPEs)!!!\n",
                   comm_size, comm_size * 64);
            printf("MPI BAM output compression level=%d\n",
                   cmd_info->compress_level_);
            if (MpiCommonScanBgzfBlocksInMemory(
                    input_memory, input_size, body_start,
                    &block_offsets, &block_lengths) != 0) {
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
        local_begin = n_blocks * rank / comm_size;
        local_end = n_blocks * (rank + 1) / comm_size;
        if (MpiCommonSelectBlockRangeFromMemory(
                input_memory, input_size,
                block_offsets, block_lengths,
                local_begin, local_end,
                &rank_input, &rank_input_size) != 0) {
            local_ok = 0;
        }
        stage41_cost = FmReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("MPI BAM scan complete. data_blocks=%lld "
                   "body_start=%lld\n", n_blocks, body_start);
            printf("Complete the 4.1 scan/split/broadcast/select cost %lf\n",
                   stage41_cost);
        }
    }
    if (!FmAllRanksOk(local_ok)) goto cleanup;

    {
        double t0 = GetTime();
        reader.base = rank_input;
        reader.size = rank_input_size;
        reader.pos = 0;
        if (MpiCommonInitMemWriter(
                middle, rank_input_size ?
                    rank_input_size : 1024 * 1024) != 0 ||
            MpiCommonInitMemWriter(prefix, 1024 * 1024) != 0 ||
            MpiCommonInitMemWriter(suffix, 1024 * 1024) != 0) {
            local_ok = 0;
        }
        stage42_cost = FmReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.2 init reader/writer cost %lf\n",
                   stage42_cost);
        }
    }
    if (!FmAllRanksOk(local_ok)) goto cleanup;

    //4.3 核心计算部分
    /* 先流式处理 rank 内部完整 QNAME group，
    再交换 rank 边界上的 leading/trailing group，
    判断当前 leading 是否属于前一个 rank，
    最后只由 group 的负责 rank 处理边界 group，
    并把 prefix/middle/suffix 按顺序拼成当前 rank 输出。*/
    {
        double fused_t0 = GetTime();
        // 1. FmStreamLocalGroups：先处理本 rank 内部完整 group
        if (FmStreamLocalGroups(
                reader, cmd_info->compress_level_,
                &middle, &leading, &trailing,
                &single_group, rank, &stats) != 0) {
            local_ok = 0;
        }
        if (!FmAllRanksOkTimed(local_ok, &stats)) goto cleanup;

        // 2. FmExchangeBoundaryGroups：跨 rank 交换 leading / trailing
        if (FmExchangeBoundaryGroups(
                &leading, &trailing, single_group,
                rank, comm_size, &continues_previous,
                &stats) != 0) {
            local_ok = 0;
        }
        if (!FmAllRanksOkTimed(local_ok, &stats)) goto cleanup;

        // 3. 处理 leading / trailing 边界 group
        if (!leading.empty()) {
            if (single_group) {
                if (!continues_previous &&
                    // FmProcessStoredGroup：真正对一个 QNAME group 做 fixmate
                    FmProcessStoredGroup(
                        &leading, cmd_info->compress_level_,
                        &prefix, rank, &stats) != 0) {
                    local_ok = 0;
                }
            } else {
                if (!continues_previous &&
                    FmProcessStoredGroup(
                        &leading, cmd_info->compress_level_,
                        &prefix, rank, &stats) != 0) {
                    local_ok = 0;
                }
                if (local_ok &&
                    FmProcessStoredGroup(
                        &trailing, cmd_info->compress_level_,
                        &suffix, rank, &stats) != 0) {
                    local_ok = 0;
                }
            }
        }
        if (local_ok &&
            // 4. FmPrependAppend：把 prefix / suffix 拼回 middle
            FmPrependAppend(&middle, prefix, suffix) != 0) {
            local_ok = 0;
        }
        if (!FmAllRanksOkTimed(local_ok, &stats)) goto cleanup;
        stats.t_fused_total = GetTime() - fused_t0;
        stage43_cost = FmReduceMax(stats.t_fused_total);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.3 FusedFixmateMPI cost %lf\n",
                   stage43_cost);
        }
    }
    if (!FmAllRanksOk(local_ok)) goto cleanup;

    if (FmGatherOutput(
            middle, header, cmd_info->compress_level_,
            cmd_info->out_file_name_, nullptr, rank, comm_size,
            &stats, &output_stage_costs,
            stage41_cost, stage42_cost, stage43_cost,
            body_total_t0) != 0) {
        local_ok = 0;
    }
    if (!FmAllRanksOk(local_ok)) goto cleanup;
    exit_code = 0;

cleanup:
    if (input) {
        if (sam_close(input) < 0) exit_code = 1;
        input = nullptr;
    }
    if (input_hfile) {
        hclose(input_hfile);
        input_hfile = nullptr;
    }
    if (header) sam_hdr_destroy(header);
    free(input_memory);
    free(middle.data);
    free(prefix.data);
    free(suffix.data);
    if (rank == 0) {
        printf("666fixmate total process cost %lf-----\n",
               GetTime() - total_t0);
    }
    return exit_code;
}
