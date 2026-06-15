#include "swbam_mpi.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <functional>
#include <limits>
#include <new>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

#include <mpi.h>

extern "C" {
    void slave_mpi_collate_extract();
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
int MpiCommonInitMemWriter(MemWriter &writer, size_t capacity);
int MpiCommonBuildBamHeaderMemory(
    sam_hdr_t *header, int compress_level,
    char **data, size_t *size);

namespace {

typedef MpiCollateRecordMetaShared CollateMeta;

const int kCollateCpes = 64;
const size_t kCollateControlReserve = 16ull * 1024 * 1024;
const size_t kCollateSimScratch = 8ull * 1024 * 1024;
const size_t kCollateExchangeRaw = 32ull * 1024 * 1024;
const size_t kCollateExchangeMeta = 4ull * 1024 * 1024;
const size_t kCollateCursorMeta = 64ull * 1024;
const size_t kCollateCursorRaw = 256ull * 1024;
const size_t kCollateConsolidateRaw = 8ull * 1024 * 1024;
const size_t kCollateConsolidateMeta = 1ull * 1024 * 1024;
const int kCollateTagHeader = 5600;
const int kCollateTagMeta = 5601;
const int kCollateTagRaw = 5602;

const unsigned char kCollateBgzfEofBlock[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
};

struct CollateBlockSet {
    bam_block *blocks;
    unsigned char *data;
    int count;

    CollateBlockSet() : blocks(nullptr), data(nullptr), count(0) {}
};

struct CollateExtractWorkspace {
    CollateBlockSet input;
    CollateBlockSet output;
    CollateMeta *records;

    CollateExtractWorkspace() : records(nullptr) {}
};

struct CollateMemorySegment {
    std::vector<CollateMeta> records;
    std::vector<unsigned char> raw;
};

static int CollateAllRanksOk(int local_ok) {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    return global_ok;
}

static double CollateReduceMax(double local) {
    double result = 0.0;
    MPI_Reduce(&local, &result, 1, MPI_DOUBLE,
               MPI_MAX, 0, MPI_COMM_WORLD);
    return result;
}

static uint32_t CollateReadLe32(const unsigned char *p) {
    return (uint32_t)p[0] |
           ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) |
           ((uint32_t)p[3] << 24);
}

static uint32_t CollateBgzfISize(const bam_block &block) {
    if (block.length < 4) return 0;
    return CollateReadLe32(block.data + block.length - 4);
}

static int CollateOwner(uint32_t bin, int bins, int ranks) {
    if (bins <= 0 || ranks <= 0) return 0;
    uint64_t owner = (uint64_t)bin * (uint64_t)ranks /
                     (uint64_t)bins;
    if (owner >= (uint64_t)ranks) owner = ranks - 1;
    return (int)owner;
}

static const unsigned char *CollateQname(
        const CollateMeta &meta,
        const unsigned char *raw, size_t raw_size) {
    if (!raw || meta.raw_offset > raw_size ||
        meta.raw_len > raw_size - (size_t)meta.raw_offset ||
        meta.raw_len < 36u + meta.qname_len + 1u) {
        return nullptr;
    }
    return raw + meta.raw_offset + 36;
}

static int CollateCompare(
        const CollateMeta &a, const unsigned char *a_raw,
        size_t a_raw_size,
        const CollateMeta &b, const unsigned char *b_raw,
        size_t b_raw_size) {
    if (a.bin != b.bin) return a.bin < b.bin ? -1 : 1;
    if (a.hash != b.hash) return a.hash < b.hash ? -1 : 1;
    const unsigned char *aq =
        CollateQname(a, a_raw, a_raw_size);
    const unsigned char *bq =
        CollateQname(b, b_raw, b_raw_size);
    if (!aq || !bq) {
        if (a.global_order == b.global_order) return 0;
        return a.global_order < b.global_order ? -1 : 1;
    }
    size_t common = std::min<size_t>(
        a.qname_len, b.qname_len);
    int cmp = memcmp(aq, bq, common);
    if (cmp != 0) return cmp < 0 ? -1 : 1;
    if (a.qname_len != b.qname_len) {
        return a.qname_len < b.qname_len ? -1 : 1;
    }
    if (a.flag_order != b.flag_order) {
        return a.flag_order < b.flag_order ? -1 : 1;
    }
    if (a.global_order != b.global_order) {
        return a.global_order < b.global_order ? -1 : 1;
    }
    return 0;
}

struct CollateLocalLess {
    const std::vector<unsigned char> *raw;

    bool operator()(const CollateMeta &a,
                    const CollateMeta &b) const {
        return CollateCompare(
                   a, raw->data(), raw->size(),
                   b, raw->data(), raw->size()) < 0;
    }
};

static int CollateAllocateBlockSet(
        CollateBlockSet *set, int count) {
    set->count = count;
    set->blocks = (bam_block *)aligned_alloc_custom(
        64, (size_t)count * sizeof(bam_block));
    set->data = aligned_alloc_custom(
        64, (size_t)count * BGZF_MAX_BLOCK_SIZE);
    if (!set->blocks || !set->data) return -1;
    memset(set->blocks, 0,
           (size_t)count * sizeof(bam_block));
    for (int i = 0; i < count; ++i) {
        set->blocks[i].data =
            set->data + (size_t)i * BGZF_MAX_BLOCK_SIZE;
        set->blocks[i].block_id = i;
    }
    return 0;
}

static void CollateFreeBlockSet(CollateBlockSet *set) {
    if (set->blocks) {
        aligned_free_custom(
            (unsigned char *)set->blocks);
    }
    if (set->data) aligned_free_custom(set->data);
    set->blocks = nullptr;
    set->data = nullptr;
    set->count = 0;
}

static int CollateAllocateExtract(
        CollateExtractWorkspace *ws) {
    if (CollateAllocateBlockSet(
            &ws->input, kCollateCpes) != 0 ||
        CollateAllocateBlockSet(
            &ws->output, kCollateCpes) != 0) {
        return -1;
    }
    ws->records = (CollateMeta *)aligned_alloc_custom(
        64, (size_t)kCollateCpes *
                MPI_RECORDS_PER_BLOCK *
                sizeof(CollateMeta));
    return ws->records ? 0 : -1;
}

static void CollateFreeExtract(
        CollateExtractWorkspace *ws) {
    CollateFreeBlockSet(&ws->input);
    CollateFreeBlockSet(&ws->output);
    if (ws->records) {
        aligned_free_custom(
            (unsigned char *)ws->records);
    }
    ws->records = nullptr;
}

static size_t CollateExtractWorkspaceBytes() {
    return 2ull * kCollateCpes *
               ((size_t)BGZF_MAX_BLOCK_SIZE +
                sizeof(bam_block)) +
           (size_t)kCollateCpes *
               MPI_RECORDS_PER_BLOCK *
               sizeof(CollateMeta);
}

static int CollateMemReadBlock(
        MemReader &reader, bam_block *block) {
    if (reader.pos >= reader.size ||
        reader.size - reader.pos <
            BLOCK_HEADER_LENGTH) {
        return -1;
    }
    const unsigned char *src =
        (const unsigned char *)reader.base + reader.pos;
    int length = (int)src[16] |
                 ((int)src[17] << 8);
    ++length;
    if (length <= 0 ||
        length > BGZF_MAX_BLOCK_SIZE ||
        (size_t)length > reader.size - reader.pos) {
        return -1;
    }
    memcpy(block->data, src, (size_t)length);
    block->length = (unsigned int)length;
    block->pos = 0;
    block->errcode = 0;
    reader.pos += (size_t)length;
    return 0;
}

static int CollateExtractAll(
        MemReader &reader, long long global_block_begin,
        int bins, std::vector<CollateMeta> *records,
        std::vector<unsigned char> *raw,
        MpiCollateStats *stats) {
    CollateExtractWorkspace ws;
    if (CollateAllocateExtract(&ws) != 0) {
        CollateFreeExtract(&ws);
        return -1;
    }
    MpiCollateExtractPara paras[kCollateCpes];
    long long local_block = 0;
    int ret = 0;
    while (reader.pos < reader.size) {
        int n_blocks = 0;
        size_t group_raw = 0;
        for (int b = 0; b < kCollateCpes; ++b) {
            if (CollateMemReadBlock(
                    reader, &ws.input.blocks[b]) != 0) {
                break;
            }
            uint32_t isize =
                CollateBgzfISize(ws.input.blocks[b]);
            if (isize > BGZF_MAX_BLOCK_SIZE) {
                ret = -1;
                break;
            }
            ws.output.blocks[b].data =
                ws.output.data +
                (size_t)b * BGZF_MAX_BLOCK_SIZE;
            ws.output.blocks[b].length = isize;
            group_raw += isize;
            ++n_blocks;
        }
        if (ret != 0 || n_blocks == 0) break;
        const size_t group_begin = raw->size();
        try {
            raw->resize(group_begin + group_raw);
        } catch (...) {
            ret = -1;
            break;
        }
        size_t raw_pos = group_begin;
        for (int b = 0; b < kCollateCpes; ++b) {
            memset(&paras[b], 0, sizeof(paras[b]));
            paras[b].block_id = b;
            paras[b].status = b < n_blocks ? 0 : -1;
            if (b >= n_blocks) continue;
            uint32_t isize =
                CollateBgzfISize(ws.input.blocks[b]);
            ws.output.blocks[b].data =
                raw->data() + raw_pos;
            ws.output.blocks[b].length = isize;
            paras[b].input_block =
                &ws.input.blocks[b];
            paras[b].un_comp_block =
                &ws.output.blocks[b];
            paras[b].raw_arena =
                raw->data() + raw_pos;
            paras[b].raw_capacity = isize;
            paras[b].raw_base_offset = raw_pos;
            paras[b].records =
                ws.records +
                (size_t)b * MPI_RECORDS_PER_BLOCK;
            paras[b].record_capacity =
                MPI_RECORDS_PER_BLOCK;
            paras[b].n_bins = bins;
            paras[b].global_block_index =
                global_block_begin + local_block + b;
            raw_pos += isize;
        }
        double t0 = GetTime();
        __real_athread_spawn(
            (void *)slave_mpi_collate_extract,
            paras, 1);
        athread_join();
        stats->t_extract += GetTime() - t0;
        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].status != 0 ||
                paras[b].raw_used !=
                    CollateBgzfISize(ws.input.blocks[b])) {
                fprintf(stderr,
                        "ERROR: collate extract failed at "
                        "global block %lld status=%d record=%d.\n",
                        global_block_begin + local_block + b,
                        paras[b].status,
                        paras[b].record_index);
                ret = -1;
                break;
            }
            try {
                records->insert(
                    records->end(), paras[b].records,
                    paras[b].records +
                        paras[b].n_records);
            } catch (...) {
                ret = -1;
                break;
            }
        }
        if (ret != 0) break;
        stats->input_blocks += n_blocks;
        local_block += n_blocks;
    }
    CollateFreeExtract(&ws);
    if (ret == 0) {
        stats->total_records =
            (long long)records->size();
    }
    return ret;
}

static size_t CollateChunkEnd(
        const std::vector<CollateMeta> &records,
        size_t begin, size_t end,
        size_t meta_capacity, size_t raw_capacity,
        size_t *raw_bytes) {
    size_t pos = begin;
    size_t bytes = 0;
    while (pos < end &&
           pos - begin < meta_capacity) {
        size_t len = records[pos].raw_len;
        if (len > raw_capacity - bytes) break;
        bytes += len;
        ++pos;
    }
    *raw_bytes = bytes;
    return pos;
}

static int CollatePackRange(
        const std::vector<CollateMeta> &records,
        const std::vector<unsigned char> &raw,
        size_t begin, size_t end,
        CollateMeta *out_meta,
        unsigned char *out_raw,
        size_t out_raw_capacity,
        size_t *out_raw_size) {
    size_t raw_pos = 0;
    for (size_t i = begin; i < end; ++i) {
        const CollateMeta &src = records[i];
        if (src.raw_offset > raw.size() ||
            src.raw_len >
                raw.size() - (size_t)src.raw_offset ||
            src.raw_len > out_raw_capacity - raw_pos) {
            return -1;
        }
        out_meta[i - begin] = src;
        out_meta[i - begin].raw_offset = raw_pos;
        memcpy(out_raw + raw_pos,
               raw.data() + src.raw_offset,
               src.raw_len);
        raw_pos += src.raw_len;
    }
    *out_raw_size = raw_pos;
    return 0;
}

static int CollateSendRecvBytes(
        int send_to, const char *send_data,
        size_t send_size, int recv_from,
        char *recv_data, size_t recv_size,
        int tag) {
    size_t sent = 0;
    size_t received = 0;
    while (sent < send_size ||
           received < recv_size) {
        int send_count = (int)std::min<size_t>(
            send_size - sent, INT_MAX);
        int recv_count = (int)std::min<size_t>(
            recv_size - received, INT_MAX);
        MPI_Status status;
        if (MPI_Sendrecv(
                send_count ? (void *)(send_data + sent) : nullptr,
                send_count, MPI_BYTE, send_to, tag,
                recv_count ? recv_data + received : nullptr,
                recv_count, MPI_BYTE, recv_from, tag,
                MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
            return -1;
        }
        sent += (size_t)send_count;
        received += (size_t)recv_count;
    }
    return 0;
}

static int CollateSendBytes(
        int dst, int tag, const char *data,
        long long size) {
    long long sent = 0;
    while (sent < size) {
        int count = (int)std::min<long long>(
            size - sent, INT_MAX);
        if (MPI_Send(
                (void *)(data + sent), count,
                MPI_BYTE, dst, tag,
                MPI_COMM_WORLD) != MPI_SUCCESS) {
            return -1;
        }
        sent += count;
    }
    return 0;
}

static int CollateRecvBytes(
        int src, int tag, char *data,
        long long size) {
    long long received = 0;
    while (received < size) {
        int count = (int)std::min<long long>(
            size - received, INT_MAX);
        MPI_Status status;
        if (MPI_Recv(
                data + received, count, MPI_BYTE,
                src, tag, MPI_COMM_WORLD,
                &status) != MPI_SUCCESS) {
            return -1;
        }
        received += count;
    }
    return 0;
}

static int CollateAppendMemoryChunk(
        const CollateMeta *records, size_t count,
        const unsigned char *raw, size_t raw_size,
        CollateMemorySegment *segment) {
    const uint64_t raw_base =
        (uint64_t)segment->raw.size();
    try {
        size_t old = segment->records.size();
        segment->records.resize(old + count);
        segment->raw.insert(
            segment->raw.end(), raw, raw + raw_size);
        for (size_t i = 0; i < count; ++i) {
            segment->records[old + i] = records[i];
            segment->records[old + i].raw_offset +=
                raw_base;
        }
    } catch (...) {
        return -1;
    }
    return 0;
}

static int CollateExchangeMemory(
        const std::vector<CollateMeta> &records,
        const std::vector<unsigned char> &raw,
        int rank, int comm_size, int bins,
        std::vector<CollateMemorySegment> *segments,
        MpiCollateStats *stats) {
    std::vector<size_t> bounds(
        (size_t)comm_size + 1, records.size());
    bounds[0] = 0;
    size_t pos = 0;
    for (int owner = 0; owner < comm_size; ++owner) {
        while (pos < records.size() &&
               CollateOwner(records[pos].bin,
                            bins, comm_size) == owner) {
            ++pos;
        }
        bounds[(size_t)owner + 1] = pos;
    }
    if (pos != records.size()) return -1;

    const size_t meta_capacity =
        kCollateExchangeMeta / sizeof(CollateMeta);
    std::vector<CollateMeta> send_meta(meta_capacity);
    std::vector<CollateMeta> recv_meta(meta_capacity);
    std::vector<unsigned char> send_raw(
        kCollateExchangeRaw);
    std::vector<unsigned char> recv_raw(
        kCollateExchangeRaw);
    segments->clear();
    segments->reserve((size_t)comm_size);

    double exchange_t0 = GetTime();
    for (int step = 0; step < comm_size; ++step) {
        int send_to = (rank + step) % comm_size;
        int recv_from =
            (rank - step + comm_size) % comm_size;
        size_t send_pos = bounds[(size_t)send_to];
        const size_t send_end =
            bounds[(size_t)send_to + 1];
        CollateMemorySegment segment;
        if (step == 0) {
            while (send_pos < send_end) {
                size_t raw_bytes = 0;
                size_t chunk_end = CollateChunkEnd(
                    records, send_pos, send_end,
                    meta_capacity, send_raw.size(),
                    &raw_bytes);
                if (chunk_end == send_pos ||
                    CollatePackRange(
                        records, raw, send_pos, chunk_end,
                        send_meta.data(), send_raw.data(),
                        send_raw.size(), &raw_bytes) != 0 ||
                    CollateAppendMemoryChunk(
                        send_meta.data(),
                        chunk_end - send_pos,
                        send_raw.data(), raw_bytes,
                        &segment) != 0) {
                    return -1;
                }
                send_pos = chunk_end;
            }
        } else {
            bool send_done = send_pos >= send_end;
            bool recv_done = false;
            while (!send_done || !recv_done) {
                size_t send_raw_bytes = 0;
                size_t chunk_end = send_pos;
                if (!send_done) {
                    chunk_end = CollateChunkEnd(
                        records, send_pos, send_end,
                        meta_capacity, send_raw.size(),
                        &send_raw_bytes);
                    if (chunk_end == send_pos ||
                        CollatePackRange(
                            records, raw, send_pos,
                            chunk_end, send_meta.data(),
                            send_raw.data(), send_raw.size(),
                            &send_raw_bytes) != 0) {
                        return -1;
                    }
                }
                unsigned long long send_header[3] = {
                    (unsigned long long)(chunk_end -
                                         send_pos),
                    (unsigned long long)send_raw_bytes,
                    chunk_end >= send_end ? 1ull : 0ull
                };
                unsigned long long recv_header[3] = {};
                MPI_Status status;
                double mpi_t0 = GetTime();
                if (MPI_Sendrecv(
                        send_header, 3,
                        MPI_UNSIGNED_LONG_LONG,
                        send_to, kCollateTagHeader,
                        recv_header, 3,
                        MPI_UNSIGNED_LONG_LONG,
                        recv_from, kCollateTagHeader,
                        MPI_COMM_WORLD, &status) !=
                    MPI_SUCCESS) {
                    return -1;
                }
                size_t recv_count =
                    (size_t)recv_header[0];
                size_t recv_raw_bytes =
                    (size_t)recv_header[1];
                if (recv_count > meta_capacity ||
                    recv_raw_bytes > recv_raw.size() ||
                    CollateSendRecvBytes(
                        send_to,
                        (const char *)send_meta.data(),
                        (chunk_end - send_pos) *
                            sizeof(CollateMeta),
                        recv_from,
                        (char *)recv_meta.data(),
                        recv_count *
                            sizeof(CollateMeta),
                        kCollateTagMeta) != 0 ||
                    CollateSendRecvBytes(
                        send_to,
                        (const char *)send_raw.data(),
                        send_raw_bytes, recv_from,
                        (char *)recv_raw.data(),
                        recv_raw_bytes,
                        kCollateTagRaw) != 0) {
                    return -1;
                }
                stats->t_mpi += GetTime() - mpi_t0;
                if (recv_count > 0 &&
                    CollateAppendMemoryChunk(
                        recv_meta.data(), recv_count,
                        recv_raw.data(), recv_raw_bytes,
                        &segment) != 0) {
                    return -1;
                }
                send_pos = chunk_end;
                send_done = send_header[2] != 0;
                recv_done = recv_header[2] != 0;
            }
        }
        if (!segment.records.empty()) {
            stats->received_records +=
                (long long)segment.records.size();
            segments->push_back(std::move(segment));
        }
    }
    stats->t_exchange += GetTime() - exchange_t0;
    return 0;
}

struct CollateMemoryCursor {
    const CollateMemorySegment *segment;
    size_t pos;
};

class CollateMemoryLoserTree {
public:
    explicit CollateMemoryLoserTree(
            const std::vector<CollateMemorySegment> *segments)
        : segments_(segments) {
        for (size_t i = 0; i < segments->size(); ++i) {
            if (!(*segments)[i].records.empty()) {
                CollateMemoryCursor cursor;
                cursor.segment = &(*segments)[i];
                cursor.pos = 0;
                cursors_.push_back(cursor);
            }
        }
        base_ = 1;
        while (base_ < cursors_.size()) base_ <<= 1;
        tree_.assign(base_ * 2, -1);
        for (size_t i = 0; i < cursors_.size(); ++i) {
            tree_[base_ + i] = (int)i;
        }
        for (size_t node = base_; node-- > 1;) {
            tree_[node] = Winner(
                tree_[node << 1],
                tree_[(node << 1) | 1]);
        }
    }

    int winner() const {
        return tree_.size() <= 1 ? -1 : tree_[1];
    }

    const CollateMeta &meta(int index) const {
        const CollateMemoryCursor &cursor =
            cursors_[(size_t)index];
        return cursor.segment->records[cursor.pos];
    }

    const unsigned char *raw(int index) const {
        const CollateMemoryCursor &cursor =
            cursors_[(size_t)index];
        return cursor.segment->raw.data() +
               meta(index).raw_offset;
    }

    void advance(int index) {
        ++cursors_[(size_t)index].pos;
        size_t node = base_ + (size_t)index;
        tree_[node] = Active(index) ? index : -1;
        while ((node >>= 1) > 0) {
            tree_[node] = Winner(
                tree_[node << 1],
                tree_[(node << 1) | 1]);
        }
    }

private:
    bool Active(int index) const {
        return index >= 0 &&
               cursors_[(size_t)index].pos <
                   cursors_[(size_t)index]
                       .segment->records.size();
    }

    bool Less(int lhs, int rhs) const {
        if (!Active(lhs)) return false;
        if (!Active(rhs)) return true;
        const CollateMemoryCursor &a =
            cursors_[(size_t)lhs];
        const CollateMemoryCursor &b =
            cursors_[(size_t)rhs];
        return CollateCompare(
                   a.segment->records[a.pos],
                   a.segment->raw.data(),
                   a.segment->raw.size(),
                   b.segment->records[b.pos],
                   b.segment->raw.data(),
                   b.segment->raw.size()) < 0;
    }

    int Winner(int lhs, int rhs) const {
        if (!Active(lhs)) return Active(rhs) ? rhs : -1;
        if (!Active(rhs)) return lhs;
        if (Less(lhs, rhs)) return lhs;
        if (Less(rhs, lhs)) return rhs;
        return lhs < rhs ? lhs : rhs;
    }

    const std::vector<CollateMemorySegment> *segments_;
    std::vector<CollateMemoryCursor> cursors_;
    std::vector<int> tree_;
    size_t base_;
};

static void CollateInitCompressPara(
        MpiSortRawCompressPara *para, int id) {
    memset(para, 0, sizeof(*para));
    para->block_id = id;
    para->status = -1;
    para->compress_level = 1;
}

static int CollateCompressStream(
        const std::function<int(
            const unsigned char **, uint32_t *)> &next,
        int compress_level, MemWriter *writer,
        MpiCollateStats *stats) {
    CollateBlockSet payload_a, payload_b;
    CollateBlockSet output_a, output_b;
    if (CollateAllocateBlockSet(
            &payload_a, kCollateCpes) != 0 ||
        CollateAllocateBlockSet(
            &payload_b, kCollateCpes) != 0 ||
        CollateAllocateBlockSet(
            &output_a, kCollateCpes) != 0 ||
        CollateAllocateBlockSet(
            &output_b, kCollateCpes) != 0) {
        CollateFreeBlockSet(&payload_a);
        CollateFreeBlockSet(&payload_b);
        CollateFreeBlockSet(&output_a);
        CollateFreeBlockSet(&output_b);
        return -1;
    }
    MpiSortRawCompressPara paras_a[kCollateCpes];
    MpiSortRawCompressPara paras_b[kCollateCpes];
    for (int i = 0; i < kCollateCpes; ++i) {
        CollateInitCompressPara(&paras_a[i], i);
        CollateInitCompressPara(&paras_b[i], i);
    }
    CollateBlockSet *payload_active = &payload_a;
    CollateBlockSet *payload_pending = &payload_b;
    CollateBlockSet *output_active = &output_a;
    CollateBlockSet *output_pending = &output_b;
    MpiSortRawCompressPara *active = paras_a;
    MpiSortRawCompressPara *pending = paras_b;
    int pending_count = 0;
    bool source_done = false;

    auto fill = [&](CollateBlockSet *payload,
                    CollateBlockSet *output,
                    MpiSortRawCompressPara *paras,
                    int *active_count) -> int {
        int count = 0;
        while (!source_done && count < kCollateCpes) {
            bam_block *block =
                &payload->blocks[count];
            block->pos = 0;
            block->length = 0;
            while (!source_done) {
                const unsigned char *record = nullptr;
                uint32_t length = 0;
                int ret = next(&record, &length);
                if (ret < 0) return -1;
                if (ret == 0) {
                    source_done = true;
                    break;
                }
                if (!record || length > BGZF_BLOCK_SIZE) {
                    return -1;
                }
                if (block->pos > 0 &&
                    block->pos + length >
                        BGZF_BLOCK_SIZE) {
                    return -2;
                }
                memcpy(block->data + block->pos,
                       record, length);
                block->pos += length;
                block->length = block->pos;
            }
            if (block->pos == 0) break;
            memset(&paras[count], 0,
                   sizeof(paras[count]));
            paras[count].block_id = count;
            paras[count].un_comp_block = block;
            paras[count].un_comp_size =
                (int)block->pos;
            paras[count].output_block =
                &output->blocks[count];
            paras[count].compress_level =
                compress_level;
            paras[count].status = 0;
            ++count;
        }
        for (int i = count; i < kCollateCpes; ++i) {
            CollateInitCompressPara(&paras[i], i);
        }
        *active_count = count;
        return 0;
    };

    /*
     * A one-record lookahead is needed when the next record does not fit
     * in the current BGZF payload. Wrap the source once so fill() never
     * consumes a record it cannot place.
     */
    const unsigned char *held_record = nullptr;
    uint32_t held_length = 0;
    std::function<int(const unsigned char **, uint32_t *)>
        original_next = next;
    auto buffered_next =
        [&](const unsigned char **record,
            uint32_t *length) -> int {
            if (held_record) {
                *record = held_record;
                *length = held_length;
                held_record = nullptr;
                held_length = 0;
                return 1;
            }
            return original_next(record, length);
        };

    auto fill_buffered =
        [&](CollateBlockSet *payload,
            CollateBlockSet *output,
            MpiSortRawCompressPara *paras,
            int *active_count) -> int {
            int count = 0;
            while (!source_done &&
                   count < kCollateCpes) {
                bam_block *block =
                    &payload->blocks[count];
                block->pos = 0;
                block->length = 0;
                while (!source_done) {
                    const unsigned char *record = nullptr;
                    uint32_t length = 0;
                    int ret =
                        buffered_next(&record, &length);
                    if (ret < 0) return -1;
                    if (ret == 0) {
                        source_done = true;
                        break;
                    }
                    if (!record ||
                        length > BGZF_BLOCK_SIZE) {
                        return -1;
                    }
                    if (block->pos > 0 &&
                        block->pos + length >
                            BGZF_BLOCK_SIZE) {
                        held_record = record;
                        held_length = length;
                        break;
                    }
                    memcpy(block->data + block->pos,
                           record, length);
                    block->pos += length;
                    block->length = block->pos;
                }
                if (block->pos == 0) break;
                memset(&paras[count], 0,
                       sizeof(paras[count]));
                paras[count].block_id = count;
                paras[count].un_comp_block = block;
                paras[count].un_comp_size =
                    (int)block->pos;
                paras[count].output_block =
                    &output->blocks[count];
                paras[count].compress_level =
                    compress_level;
                paras[count].status = 0;
                ++count;
            }
            for (int i = count;
                 i < kCollateCpes; ++i) {
                CollateInitCompressPara(
                    &paras[i], i);
            }
            *active_count = count;
            return 0;
        };

    auto flush_pending = [&]() -> int {
        double t0 = GetTime();
        for (int i = 0; i < pending_count; ++i) {
            if (pending[i].status != 0 ||
                !pending[i].output_block ||
                MpiWriteBlockToMem(
                    *writer,
                    pending[i].output_block) != 0) {
                return -1;
            }
            ++stats->bgzf_blocks;
        }
        stats->t_write += GetTime() - t0;
        pending_count = 0;
        return 0;
    };

    int active_count = 0;
    if (fill_buffered(
            payload_active, output_active,
            active, &active_count) != 0) {
        CollateFreeBlockSet(&payload_a);
        CollateFreeBlockSet(&payload_b);
        CollateFreeBlockSet(&output_a);
        CollateFreeBlockSet(&output_b);
        return -1;
    }
    while (active_count > 0) {
        double t0 = GetTime();
        __real_athread_spawn(
            (void *)slave_mpi_sort_compress_payload,
            active, 1);
        if (flush_pending() != 0) {
            athread_join();
            active_count = -1;
            break;
        }
        int next_count = 0;
        if (fill_buffered(
                payload_pending, output_pending,
                pending, &next_count) != 0) {
            athread_join();
            active_count = -1;
            break;
        }
        athread_join();
        stats->t_compress += GetTime() - t0;
        for (int i = 0; i < active_count; ++i) {
            if (active[i].status != 0) {
                active_count = -1;
                break;
            }
        }
        if (active_count < 0) break;
        pending_count = active_count;
        std::swap(payload_active, payload_pending);
        std::swap(output_active, output_pending);
        std::swap(active, pending);
        active_count = next_count;
    }
    int result = active_count < 0
        ? -1 : flush_pending();
    CollateFreeBlockSet(&payload_a);
    CollateFreeBlockSet(&payload_b);
    CollateFreeBlockSet(&output_a);
    CollateFreeBlockSet(&output_b);
    (void)fill;
    return result;
}

static void CollateCountGroup(
        const CollateMeta &meta,
        const unsigned char *raw,
        std::string *previous,
        uint32_t *previous_hash,
        MpiCollateStats *stats) {
    const char *qname =
        (const char *)(raw + 36);
    bool changed =
        previous->size() != meta.qname_len ||
        (meta.qname_len > 0 &&
         memcmp(previous->data(), qname,
                meta.qname_len) != 0);
    if (changed) {
        if (!previous->empty() &&
            meta.hash == *previous_hash) {
            ++stats->hash_collision_groups;
        }
        ++stats->qname_groups;
        previous->assign(qname, meta.qname_len);
        *previous_hash = meta.hash;
    }
}

static int CollateMergeMemory(
        const std::vector<CollateMemorySegment> &segments,
        int compress_level, MemWriter *writer,
        MpiCollateStats *stats) {
    CollateMemoryLoserTree tree(&segments);
    std::string previous;
    uint32_t previous_hash = 0;
    auto next = [&](const unsigned char **raw,
                    uint32_t *length) -> int {
        int winner = tree.winner();
        if (winner < 0) return 0;
        const CollateMeta meta = tree.meta(winner);
        const unsigned char *record =
            tree.raw(winner);
        CollateCountGroup(
            meta, record, &previous,
            &previous_hash, stats);
        *raw = record;
        *length = meta.raw_len;
        tree.advance(winner);
        return 1;
    };
    double t0 = GetTime();
    int ret = CollateCompressStream(
        next, compress_level, writer, stats);
    stats->t_merge += GetTime() - t0;
    return ret;
}

static int CollateParseMemory(
        const std::string &text, size_t *bytes) {
    *bytes = 0;
    if (text.empty()) return 0;
    char *end = nullptr;
    errno = 0;
    unsigned long long value =
        strtoull(text.c_str(), &end, 10);
    if (errno || end == text.c_str()) return -1;
    unsigned long long multiplier = 1;
    if (*end) {
        if (end[1] != '\0') return -1;
        if (*end == 'K' || *end == 'k') {
            multiplier = 1024ull;
        } else if (*end == 'M' || *end == 'm') {
            multiplier = 1024ull * 1024ull;
        } else if (*end == 'G' || *end == 'g') {
            multiplier =
                1024ull * 1024ull * 1024ull;
        } else {
            return -1;
        }
    }
    if (value > SIZE_MAX / multiplier) return -1;
    *bytes = (size_t)(value * multiplier);
    return *bytes ? 0 : -1;
}

static int CollateGatherOutput(
        const MemWriter &local_writer,
        sam_hdr_t *header, int compress_level,
        const std::string &output_path,
        int rank, int comm_size,
        MpiCollateStats *stats) {
    long long local_size =
        (long long)local_writer.size;
    std::vector<long long> sizes(
        (size_t)comm_size, 0);
    MPI_Allgather(
        &local_size, 1, MPI_LONG_LONG,
        sizes.data(), 1, MPI_LONG_LONG,
        MPI_COMM_WORLD);
    std::vector<long long> offsets(
        (size_t)comm_size, 0);
    long long body_size = 0;
    for (int i = 0; i < comm_size; ++i) {
        offsets[(size_t)i] = body_size;
        if (sizes[(size_t)i] < 0 ||
            body_size >
                LLONG_MAX - sizes[(size_t)i]) {
            return -1;
        }
        body_size += sizes[(size_t)i];
    }

    char *header_memory = nullptr;
    size_t header_size = 0;
    char *output_memory = nullptr;
    size_t output_size = 0;
    int local_ok = 1;
    if (rank == 0) {
        if (sam_hdr_update_hd(
                header, "SO", "unsorted",
                "GO", "query", NULL) < 0 ||
            sam_hdr_add_pg(
                header, "RabbitBAM-MPI",
                "PN", "RabbitBAM-MPI",
                "CL", "collate", NULL) != 0 ||
            MpiCommonBuildBamHeaderMemory(
                header, compress_level,
                &header_memory,
                &header_size) != 0) {
            local_ok = 0;
        }
        unsigned long long total =
            (unsigned long long)header_size +
            (unsigned long long)body_size +
            sizeof(kCollateBgzfEofBlock);
        if (local_ok && total <= SIZE_MAX) {
            output_size = (size_t)total;
            output_memory =
                (char *)malloc(output_size);
            if (!output_memory) local_ok = 0;
        } else if (local_ok) {
            local_ok = 0;
        }
    }
    if (!CollateAllRanksOk(local_ok)) {
        free(header_memory);
        free(output_memory);
        return -1;
    }

    double t0 = GetTime();
    if (rank == 0) {
        memcpy(output_memory, header_memory,
               header_size);
        if (local_writer.size) {
            memcpy(output_memory + header_size,
                   local_writer.data,
                   local_writer.size);
        }
        for (int src = 1; src < comm_size; ++src) {
            if (CollateRecvBytes(
                    src, 5700,
                    output_memory + header_size +
                        offsets[(size_t)src],
                    sizes[(size_t)src]) != 0) {
                local_ok = 0;
            }
        }
        memcpy(output_memory + header_size +
                   body_size,
               kCollateBgzfEofBlock,
               sizeof(kCollateBgzfEofBlock));
    } else if (CollateSendBytes(
                   0, 5700, local_writer.data,
                   local_size) != 0) {
        local_ok = 0;
    }
    stats->t_write += GetTime() - t0;
    if (rank == 0 && local_ok &&
        MpiCommonDumpMemoryToFile(
            output_path, output_memory,
            output_size) != 0) {
        local_ok = 0;
    }
    free(header_memory);
    free(output_memory);
    return CollateAllRanksOk(local_ok) ? 0 : -1;
}

/*
 * The strict external implementation follows below. It uses the same
 * collate key and exchange order as the memory path, while keeping runs
 * and received segments in uncompressed internal temporary files.
 */

struct CollateExtent {
    int fd;
    off_t meta_offset;
    off_t raw_offset;
    uint64_t record_count;
    uint64_t raw_size;

    CollateExtent()
        : fd(-1), meta_offset(0), raw_offset(0),
          record_count(0), raw_size(0) {}
};

struct CollateRun {
    CollateExtent extent;
    bool resident;

    CollateRun() : resident(false) {}
};

struct CollateSegment {
    std::vector<CollateExtent> extents;
    uint64_t record_count;
    uint64_t raw_size;

    CollateSegment() : record_count(0), raw_size(0) {}
};

struct CollateTempStore {
    int runs_fd;
    int segments_fd;
    off_t runs_end;
    off_t segments_end;

    CollateTempStore()
        : runs_fd(-1), segments_fd(-1),
          runs_end(0), segments_end(0) {}
};

struct CollateRunArena {
    unsigned char *data;
    size_t capacity;
    size_t raw_used;
    size_t meta_begin;
    size_t record_count;

    CollateRunArena()
        : data(nullptr), capacity(0), raw_used(0),
          meta_begin(0), record_count(0) {}

    CollateMeta *records() {
        return reinterpret_cast<CollateMeta *>(
            data + meta_begin);
    }

    const CollateMeta *records() const {
        return reinterpret_cast<const CollateMeta *>(
            data + meta_begin);
    }

    void reset() {
        raw_used = 0;
        meta_begin = capacity;
        record_count = 0;
    }
};

struct CollatePointerLess {
    const unsigned char *raw;
    size_t raw_size;

    bool operator()(const CollateMeta &a,
                    const CollateMeta &b) const {
        return CollateCompare(
                   a, raw, raw_size,
                   b, raw, raw_size) < 0;
    }
};

static size_t CollateCompressWorkspaceBytes() {
    return 4ull * kCollateCpes *
           ((size_t)BGZF_MAX_BLOCK_SIZE +
            sizeof(bam_block));
}

static int CollatePwriteAll(
        int fd, const void *data, size_t size,
        off_t offset, MpiCollateStats *stats) {
    const unsigned char *src =
        static_cast<const unsigned char *>(data);
    size_t done = 0;
    double t0 = GetTime();
    while (done < size) {
        ssize_t n = pwrite(
            fd, src + done, size - done,
            offset + (off_t)done);
        if (n < 0 && errno == EINTR) continue;
        if (n <= 0) {
            stats->t_temp_write_actual += GetTime() - t0;
            return -1;
        }
        done += (size_t)n;
    }
    stats->t_temp_write_actual += GetTime() - t0;
    stats->temp_write_bytes += (long long)size;
    return 0;
}

static int CollatePreadAll(
        int fd, void *data, size_t size,
        off_t offset, MpiCollateStats *stats) {
    unsigned char *dst =
        static_cast<unsigned char *>(data);
    size_t done = 0;
    double t0 = GetTime();
    while (done < size) {
        ssize_t n = pread(
            fd, dst + done, size - done,
            offset + (off_t)done);
        if (n < 0 && errno == EINTR) continue;
        if (n <= 0) {
            stats->t_temp_read_actual += GetTime() - t0;
            return -1;
        }
        done += (size_t)n;
    }
    stats->t_temp_read_actual += GetTime() - t0;
    stats->temp_read_bytes += (long long)size;
    return 0;
}

static void CollateSimulateCopy(
        const void *data, size_t size,
        unsigned char *scratch, size_t scratch_size,
        double *elapsed) {
    if (!data || !size || !scratch || !scratch_size) return;
    const unsigned char *src =
        static_cast<const unsigned char *>(data);
    volatile unsigned char guard = 0;
    double t0 = GetTime();
    size_t done = 0;
    while (done < size) {
        size_t n = std::min(scratch_size, size - done);
        memcpy(scratch, src + done, n);
        guard ^= scratch[n - 1];
        done += n;
    }
    *elapsed += GetTime() - t0;
    (void)guard;
}

static int CollateOpenOneTemp(
        const std::string &path, int *fd) {
    int opened = open(
        path.c_str(), O_RDWR | O_CREAT | O_EXCL, 0600);
    if (opened < 0) return -1;
    if (unlink(path.c_str()) != 0) {
        close(opened);
        return -1;
    }
    *fd = opened;
    return 0;
}

static int CollateOpenTemp(
        const std::string &prefix, int rank,
        CollateTempStore *store) {
    struct stat st;
    bool is_directory =
        stat(prefix.c_str(), &st) == 0 &&
        S_ISDIR(st.st_mode);
    for (int attempt = 0; attempt < 100; ++attempt) {
        char base[PATH_MAX];
        if (is_directory) {
            snprintf(
                base, sizeof(base),
                "%s%srabbitbam-collate.%ld.%d.%d",
                prefix.c_str(),
                !prefix.empty() && prefix.back() == '/'
                    ? "" : "/",
                (long)getpid(), rank, attempt);
        } else {
            snprintf(
                base, sizeof(base),
                "%s.%ld.%d.%d",
                prefix.c_str(), (long)getpid(),
                rank, attempt);
        }
        std::string runs =
            std::string(base) + ".runs.tmp";
        std::string segments =
            std::string(base) + ".segments.tmp";
        if (CollateOpenOneTemp(
                runs, &store->runs_fd) != 0) {
            if (errno == EEXIST) continue;
            return -1;
        }
        if (CollateOpenOneTemp(
                segments, &store->segments_fd) != 0) {
            close(store->runs_fd);
            store->runs_fd = -1;
            if (errno == EEXIST) continue;
            return -1;
        }
        return 0;
    }
    return -1;
}

static void CollateCloseTemp(CollateTempStore *store) {
    if (store->runs_fd >= 0) close(store->runs_fd);
    if (store->segments_fd >= 0) {
        close(store->segments_fd);
    }
    store->runs_fd = -1;
    store->segments_fd = -1;
}

static int CollateAppendExtent(
        int fd, off_t *file_end,
        const CollateMeta *records, size_t count,
        const unsigned char *raw, size_t raw_size,
        unsigned char *scratch, size_t scratch_size,
        MpiCollateStats *stats,
        CollateExtent *extent) {
    if (count > SIZE_MAX / sizeof(CollateMeta)) return -1;
    const size_t meta_bytes =
        count * sizeof(CollateMeta);
    const off_t off_max =
        std::numeric_limits<off_t>::max();
    if (*file_end < 0 ||
        meta_bytes > (size_t)off_max ||
        raw_size > (size_t)off_max ||
        *file_end > off_max - (off_t)meta_bytes ||
        *file_end + (off_t)meta_bytes >
            off_max - (off_t)raw_size) {
        return -1;
    }
    extent->fd = fd;
    extent->meta_offset = *file_end;
    extent->raw_offset =
        *file_end + (off_t)meta_bytes;
    extent->record_count = count;
    extent->raw_size = raw_size;
    if (CollatePwriteAll(
            fd, records, meta_bytes,
            extent->meta_offset, stats) != 0 ||
        CollatePwriteAll(
            fd, raw, raw_size,
            extent->raw_offset, stats) != 0) {
        return -1;
    }
    CollateSimulateCopy(
        records, meta_bytes, scratch, scratch_size,
        &stats->t_temp_write_sim);
    CollateSimulateCopy(
        raw, raw_size, scratch, scratch_size,
        &stats->t_temp_write_sim);
    *file_end =
        extent->raw_offset + (off_t)raw_size;
    return 0;
}

static int CollateArenaCanAppend(
        const CollateRunArena &arena,
        size_t raw_size, size_t records) {
    if (records >
        SIZE_MAX / sizeof(CollateMeta)) return 0;
    size_t meta_bytes =
        records * sizeof(CollateMeta);
    if (raw_size > arena.meta_begin - arena.raw_used ||
        meta_bytes >
            arena.meta_begin - arena.raw_used - raw_size) {
        return 0;
    }
    return 1;
}

static int CollateArenaAppendBlock(
        CollateRunArena *arena,
        const unsigned char *raw, size_t raw_size,
        const CollateMeta *records, size_t count) {
    if (!CollateArenaCanAppend(
            *arena, raw_size, count)) {
        return 1;
    }
    const size_t raw_base = arena->raw_used;
    memcpy(arena->data + raw_base, raw, raw_size);
    const size_t meta_bytes =
        count * sizeof(CollateMeta);
    arena->meta_begin -= meta_bytes;
    CollateMeta *dst = arena->records();
    for (size_t i = 0; i < count; ++i) {
        dst[i] = records[i];
        dst[i].raw_offset =
            raw_base + records[i].raw_offset;
    }
    arena->raw_used += raw_size;
    arena->record_count += count;
    return 0;
}

static void CollateSortArena(
        CollateRunArena *arena,
        MpiCollateStats *stats) {
    double t0 = GetTime();
    std::sort(
        arena->records(),
        arena->records() + arena->record_count,
        CollatePointerLess{
            arena->data, arena->raw_used});
    stats->t_local_sort += GetTime() - t0;
}

static int CollateSpillArena(
        CollateRunArena *arena,
        CollateTempStore *store,
        unsigned char *scratch, size_t scratch_size,
        MpiCollateStats *stats,
        std::vector<CollateRun> *runs) {
    if (arena->record_count == 0) return 0;
    CollateSortArena(arena, stats);
    CollateRun run;
    if (CollateAppendExtent(
            store->runs_fd, &store->runs_end,
            arena->records(), arena->record_count,
            arena->data, arena->raw_used,
            scratch, scratch_size, stats,
            &run.extent) != 0) {
        return -1;
    }
    runs->push_back(run);
    arena->reset();
    return 0;
}

static int CollateExtractRuns(
        MemReader &reader, long long global_block_begin,
        int bins, CollateRunArena *arena,
        CollateTempStore *store,
        unsigned char *scratch, size_t scratch_size,
        std::vector<CollateRun> *runs,
        MpiCollateStats *stats) {
    CollateExtractWorkspace ws;
    if (CollateAllocateExtract(&ws) != 0) {
        CollateFreeExtract(&ws);
        return -1;
    }
    MpiCollateExtractPara paras[kCollateCpes];
    long long local_block = 0;
    int result = 0;
    while (reader.pos < reader.size) {
        int n_blocks = 0;
        for (int b = 0; b < kCollateCpes; ++b) {
            if (CollateMemReadBlock(
                    reader, &ws.input.blocks[b]) != 0) {
                break;
            }
            uint32_t isize =
                CollateBgzfISize(ws.input.blocks[b]);
            if (isize > BGZF_MAX_BLOCK_SIZE) {
                result = -1;
                break;
            }
            ws.output.blocks[b].data =
                ws.output.data +
                (size_t)b * BGZF_MAX_BLOCK_SIZE;
            ws.output.blocks[b].length = isize;
            memset(&paras[b], 0, sizeof(paras[b]));
            paras[b].block_id = b;
            paras[b].input_block =
                &ws.input.blocks[b];
            paras[b].un_comp_block =
                &ws.output.blocks[b];
            paras[b].raw_arena =
                ws.output.blocks[b].data;
            paras[b].raw_capacity = isize;
            paras[b].raw_base_offset = 0;
            paras[b].records =
                ws.records +
                (size_t)b * MPI_RECORDS_PER_BLOCK;
            paras[b].record_capacity =
                MPI_RECORDS_PER_BLOCK;
            paras[b].n_bins = bins;
            paras[b].global_block_index =
                global_block_begin + local_block + b;
            ++n_blocks;
        }
        for (int b = n_blocks; b < kCollateCpes; ++b) {
            memset(&paras[b], 0, sizeof(paras[b]));
            paras[b].status = -1;
        }
        if (result != 0 || n_blocks == 0) break;
        double t0 = GetTime();
        __real_athread_spawn(
            (void *)slave_mpi_collate_extract,
            paras, 1);
        athread_join();
        stats->t_extract += GetTime() - t0;
        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].status != 0 ||
                paras[b].raw_used !=
                    CollateBgzfISize(ws.input.blocks[b])) {
                fprintf(
                    stderr,
                    "ERROR: collate external extract failed "
                    "at block %lld status=%d.\n",
                    global_block_begin + local_block + b,
                    paras[b].status);
                result = -1;
                break;
            }
            CollateMeta *block_records =
                ws.records +
                (size_t)b * MPI_RECORDS_PER_BLOCK;
            int append = CollateArenaAppendBlock(
                arena, ws.output.blocks[b].data,
                paras[b].raw_used, block_records,
                (size_t)paras[b].n_records);
            if (append == 1) {
                if (CollateSpillArena(
                        arena, store, scratch,
                        scratch_size, stats, runs) != 0) {
                    result = -1;
                    break;
                }
                append = CollateArenaAppendBlock(
                    arena, ws.output.blocks[b].data,
                    paras[b].raw_used, block_records,
                    (size_t)paras[b].n_records);
            }
            if (append != 0) {
                fprintf(
                    stderr,
                    "ERROR: collate -m is too small for one "
                    "BGZF block (raw=%zu records=%d arena=%zu).\n",
                    paras[b].raw_used,
                    paras[b].n_records,
                    arena->capacity);
                result = -1;
                break;
            }
            stats->total_records += paras[b].n_records;
        }
        if (result != 0) break;
        stats->input_blocks += n_blocks;
        local_block += n_blocks;
    }
    CollateFreeExtract(&ws);
    if (result != 0) return result;
    if (arena->record_count > 0) {
        CollateSortArena(arena, stats);
        CollateRun resident;
        resident.resident = true;
        resident.extent.record_count =
            arena->record_count;
        resident.extent.raw_size = arena->raw_used;
        runs->insert(runs->begin(), resident);
        stats->resident_run_records =
            (long long)arena->record_count;
        stats->resident_run_raw_bytes =
            (long long)arena->raw_used;
    }
    stats->runs = (long long)runs->size();
    return 0;
}

static int CollateLoadRun(
        const CollateRun &run,
        CollateRunArena *arena,
        unsigned char *scratch, size_t scratch_size,
        MpiCollateStats *stats) {
    if (run.resident) return 0;
    const uint64_t meta_bytes =
        run.extent.record_count * sizeof(CollateMeta);
    if (meta_bytes > arena->capacity ||
        run.extent.raw_size >
            arena->capacity - (size_t)meta_bytes) {
        return -1;
    }
    arena->raw_used = (size_t)run.extent.raw_size;
    arena->record_count =
        (size_t)run.extent.record_count;
    arena->meta_begin =
        arena->capacity - (size_t)meta_bytes;
    if (CollatePreadAll(
            run.extent.fd, arena->records(),
            (size_t)meta_bytes,
            run.extent.meta_offset, stats) != 0 ||
        CollatePreadAll(
            run.extent.fd, arena->data,
            arena->raw_used,
            run.extent.raw_offset, stats) != 0) {
        return -1;
    }
    CollateSimulateCopy(
        arena->records(), (size_t)meta_bytes,
        scratch, scratch_size,
        &stats->t_temp_read_sim);
    CollateSimulateCopy(
        arena->data, arena->raw_used,
        scratch, scratch_size,
        &stats->t_temp_read_sim);
    return 0;
}

static size_t CollateChunkEndPtr(
        const CollateMeta *records,
        size_t begin, size_t end,
        size_t meta_capacity, size_t raw_capacity,
        size_t *raw_bytes) {
    size_t pos = begin;
    size_t bytes = 0;
    while (pos < end &&
           pos - begin < meta_capacity) {
        size_t len = records[pos].raw_len;
        if (len > raw_capacity - bytes) break;
        bytes += len;
        ++pos;
    }
    *raw_bytes = bytes;
    return pos;
}

static int CollatePackRangePtr(
        const CollateMeta *records, size_t record_count,
        const unsigned char *raw, size_t raw_size,
        size_t begin, size_t end,
        CollateMeta *out_meta,
        unsigned char *out_raw,
        size_t out_raw_capacity,
        size_t *out_raw_size) {
    if (begin > end || end > record_count) return -1;
    size_t raw_pos = 0;
    for (size_t i = begin; i < end; ++i) {
        const CollateMeta &src = records[i];
        if (src.raw_offset > raw_size ||
            src.raw_len >
                raw_size - (size_t)src.raw_offset ||
            src.raw_len >
                out_raw_capacity - raw_pos) {
            return -1;
        }
        out_meta[i - begin] = src;
        out_meta[i - begin].raw_offset = raw_pos;
        memcpy(
            out_raw + raw_pos,
            raw + src.raw_offset, src.raw_len);
        raw_pos += src.raw_len;
    }
    *out_raw_size = raw_pos;
    return 0;
}

static int CollateAppendSegmentChunk(
        CollateTempStore *store,
        const CollateMeta *records, size_t count,
        const unsigned char *raw, size_t raw_size,
        unsigned char *scratch, size_t scratch_size,
        MpiCollateStats *stats,
        CollateSegment *segment) {
    CollateExtent extent;
    if (CollateAppendExtent(
            store->segments_fd,
            &store->segments_end,
            records, count, raw, raw_size,
            scratch, scratch_size, stats,
            &extent) != 0) {
        return -1;
    }
    segment->record_count += count;
    segment->raw_size += raw_size;
    segment->extents.push_back(extent);
    return 0;
}

static int CollateExchangeRuns(
        const std::vector<CollateRun> &runs,
        CollateRunArena *arena,
        CollateTempStore *store,
        int rank, int comm_size, int bins,
        unsigned char *scratch, size_t scratch_size,
        std::vector<CollateSegment> *segments,
        MpiCollateStats *stats) {
    unsigned long long local_runs = runs.size();
    unsigned long long max_runs = 0;
    MPI_Allreduce(
        &local_runs, &max_runs, 1,
        MPI_UNSIGNED_LONG_LONG, MPI_MAX,
        MPI_COMM_WORLD);
    const size_t raw_capacity =
        std::min<size_t>(
            kCollateExchangeRaw, 8ull * 1024 * 1024);
    const size_t meta_bytes =
        std::min<size_t>(
            kCollateExchangeMeta, 1ull * 1024 * 1024);
    const size_t meta_capacity =
        meta_bytes / sizeof(CollateMeta);
    std::vector<CollateMeta> send_meta(meta_capacity);
    std::vector<CollateMeta> recv_meta(meta_capacity);
    std::vector<unsigned char> send_raw(raw_capacity);
    std::vector<unsigned char> recv_raw(raw_capacity);

    for (unsigned long long round = 0;
         round < max_runs; ++round) {
        bool has_run = round < runs.size();
        if (has_run &&
            CollateLoadRun(
                runs[(size_t)round], arena,
                scratch, scratch_size, stats) != 0) {
            return -1;
        }
        const CollateMeta *records =
            has_run ? arena->records() : nullptr;
        size_t record_count =
            has_run ? arena->record_count : 0;
        std::vector<size_t> bounds(
            (size_t)comm_size + 1, record_count);
        bounds[0] = 0;
        size_t pos = 0;
        for (int owner = 0;
             owner < comm_size; ++owner) {
            while (pos < record_count &&
                   CollateOwner(
                       records[pos].bin,
                       bins, comm_size) == owner) {
                ++pos;
            }
            bounds[(size_t)owner + 1] = pos;
        }
        if (pos != record_count) return -1;

        double exchange_t0 = GetTime();
        for (int step = 0; step < comm_size; ++step) {
            const int send_to =
                (rank + step) % comm_size;
            const int recv_from =
                (rank - step + comm_size) %
                comm_size;
            size_t send_pos =
                bounds[(size_t)send_to];
            const size_t send_end =
                bounds[(size_t)send_to + 1];
            CollateSegment segment;
            bool send_done = send_pos >= send_end;
            bool recv_done = false;
            if (step == 0) {
                while (send_pos < send_end) {
                    size_t raw_bytes = 0;
                    size_t chunk_end =
                        CollateChunkEndPtr(
                            records, send_pos,
                            send_end, meta_capacity,
                            send_raw.size(),
                            &raw_bytes);
                    if (chunk_end == send_pos ||
                        CollatePackRangePtr(
                            records, record_count,
                            arena->data,
                            arena->raw_used,
                            send_pos, chunk_end,
                            send_meta.data(),
                            send_raw.data(),
                            send_raw.size(),
                            &raw_bytes) != 0 ||
                        CollateAppendSegmentChunk(
                            store, send_meta.data(),
                            chunk_end - send_pos,
                            send_raw.data(), raw_bytes,
                            scratch, scratch_size,
                            stats, &segment) != 0) {
                        return -1;
                    }
                    send_pos = chunk_end;
                }
                recv_done = true;
            } else {
                while (!send_done || !recv_done) {
                    size_t send_raw_bytes = 0;
                    size_t chunk_end = send_pos;
                    if (!send_done) {
                        chunk_end =
                            CollateChunkEndPtr(
                                records, send_pos,
                                send_end,
                                meta_capacity,
                                send_raw.size(),
                                &send_raw_bytes);
                        if (chunk_end == send_pos ||
                            CollatePackRangePtr(
                                records, record_count,
                                arena->data,
                                arena->raw_used,
                                send_pos, chunk_end,
                                send_meta.data(),
                                send_raw.data(),
                                send_raw.size(),
                                &send_raw_bytes) != 0) {
                            return -1;
                        }
                    }
                    unsigned long long send_header[3] = {
                        (unsigned long long)
                            (chunk_end - send_pos),
                        (unsigned long long)
                            send_raw_bytes,
                        chunk_end >= send_end
                            ? 1ull : 0ull
                    };
                    unsigned long long recv_header[3] = {};
                    MPI_Status status;
                    double mpi_t0 = GetTime();
                    if (MPI_Sendrecv(
                            send_header, 3,
                            MPI_UNSIGNED_LONG_LONG,
                            send_to, kCollateTagHeader,
                            recv_header, 3,
                            MPI_UNSIGNED_LONG_LONG,
                            recv_from, kCollateTagHeader,
                            MPI_COMM_WORLD, &status) !=
                        MPI_SUCCESS) {
                        return -1;
                    }
                    const size_t recv_count =
                        (size_t)recv_header[0];
                    const size_t recv_raw_bytes =
                        (size_t)recv_header[1];
                    if (recv_count > meta_capacity ||
                        recv_raw_bytes >
                            recv_raw.size() ||
                        CollateSendRecvBytes(
                            send_to,
                            reinterpret_cast<const char *>(
                                send_meta.data()),
                            (chunk_end - send_pos) *
                                sizeof(CollateMeta),
                            recv_from,
                            reinterpret_cast<char *>(
                                recv_meta.data()),
                            recv_count *
                                sizeof(CollateMeta),
                            kCollateTagMeta) != 0 ||
                        CollateSendRecvBytes(
                            send_to,
                            reinterpret_cast<const char *>(
                                send_raw.data()),
                            send_raw_bytes, recv_from,
                            reinterpret_cast<char *>(
                                recv_raw.data()),
                            recv_raw_bytes,
                            kCollateTagRaw) != 0) {
                        return -1;
                    }
                    stats->t_mpi += GetTime() - mpi_t0;
                    if (recv_count > 0 &&
                        CollateAppendSegmentChunk(
                            store, recv_meta.data(),
                            recv_count,
                            recv_raw.data(),
                            recv_raw_bytes,
                            scratch, scratch_size,
                            stats, &segment) != 0) {
                        return -1;
                    }
                    send_pos = chunk_end;
                    send_done = send_header[2] != 0;
                    recv_done = recv_header[2] != 0;
                }
            }
            if (segment.record_count > 0) {
                stats->received_records +=
                    (long long)segment.record_count;
                segments->push_back(
                    std::move(segment));
            }
        }
        stats->t_exchange += GetTime() - exchange_t0;
        if (has_run &&
            runs[(size_t)round].resident) {
            arena->reset();
        }
    }
    stats->segments = (long long)segments->size();
    return 0;
}

struct CollateFileCursor {
    const CollateSegment *segment;
    size_t extent_index;
    uint64_t extent_record_pos;
    CollateMeta *meta_buffer;
    size_t meta_capacity;
    size_t meta_count;
    size_t meta_pos;
    unsigned char *raw_buffer;
    size_t raw_capacity;
    uint64_t raw_begin;
    size_t raw_size;
    bool active;

    CollateFileCursor()
        : segment(nullptr), extent_index(0),
          extent_record_pos(0), meta_buffer(nullptr),
          meta_capacity(0), meta_count(0), meta_pos(0),
          raw_buffer(nullptr), raw_capacity(0),
          raw_begin(0), raw_size(0), active(false) {}
};

static int CollateCursorLoadMeta(
        CollateFileCursor *cursor,
        unsigned char *scratch, size_t scratch_size,
        MpiCollateStats *stats) {
    cursor->active = false;
    while (cursor->extent_index <
           cursor->segment->extents.size()) {
        const CollateExtent &extent =
            cursor->segment->extents[
                cursor->extent_index];
        if (cursor->extent_record_pos >=
            extent.record_count) {
            ++cursor->extent_index;
            cursor->extent_record_pos = 0;
            cursor->raw_size = 0;
            continue;
        }
        size_t count = (size_t)std::min<uint64_t>(
            extent.record_count -
                cursor->extent_record_pos,
            cursor->meta_capacity);
        size_t bytes = count * sizeof(CollateMeta);
        off_t offset = extent.meta_offset +
            (off_t)(cursor->extent_record_pos *
                    sizeof(CollateMeta));
        if (CollatePreadAll(
                extent.fd, cursor->meta_buffer,
                bytes, offset, stats) != 0) {
            return -1;
        }
        CollateSimulateCopy(
            cursor->meta_buffer, bytes,
            scratch, scratch_size,
            &stats->t_temp_read_sim);
        cursor->meta_count = count;
        cursor->meta_pos = 0;
        cursor->active = true;
        return 0;
    }
    return 0;
}

static const CollateMeta &CollateCursorMeta(
        const CollateFileCursor &cursor) {
    return cursor.meta_buffer[cursor.meta_pos];
}

static int CollateCursorRaw(
        CollateFileCursor *cursor,
        unsigned char *scratch, size_t scratch_size,
        MpiCollateStats *stats,
        const unsigned char **raw) {
    const CollateMeta &meta =
        CollateCursorMeta(*cursor);
    const CollateExtent &extent =
        cursor->segment->extents[
            cursor->extent_index];
    uint64_t raw_end =
        meta.raw_offset + meta.raw_len;
    if (raw_end > extent.raw_size ||
        meta.raw_len > cursor->raw_capacity) {
        return -1;
    }
    if (cursor->raw_size == 0 ||
        meta.raw_offset < cursor->raw_begin ||
        raw_end >
            cursor->raw_begin + cursor->raw_size) {
        cursor->raw_begin = meta.raw_offset;
        cursor->raw_size =
            (size_t)std::min<uint64_t>(
                cursor->raw_capacity,
                extent.raw_size -
                    meta.raw_offset);
        if (CollatePreadAll(
                extent.fd, cursor->raw_buffer,
                cursor->raw_size,
                extent.raw_offset +
                    (off_t)meta.raw_offset,
                stats) != 0) {
            return -1;
        }
        CollateSimulateCopy(
            cursor->raw_buffer,
            cursor->raw_size,
            scratch, scratch_size,
            &stats->t_temp_read_sim);
    }
    *raw = cursor->raw_buffer +
        (size_t)(meta.raw_offset -
                 cursor->raw_begin);
    if (meta.raw_len < 36 ||
        (uint64_t)CollateReadLe32(*raw) + 4u !=
            meta.raw_len) {
        return -1;
    }
    return 0;
}

static int CollateCursorAdvance(
        CollateFileCursor *cursor,
        unsigned char *scratch, size_t scratch_size,
        MpiCollateStats *stats) {
    ++cursor->meta_pos;
    ++cursor->extent_record_pos;
    if (cursor->meta_pos <
        cursor->meta_count) {
        return 0;
    }
    return CollateCursorLoadMeta(
        cursor, scratch, scratch_size, stats);
}

class CollateFileLoserTree {
public:
    explicit CollateFileLoserTree(
            std::vector<CollateFileCursor> *cursors)
        : cursors_(cursors), base_(1) {
        while (base_ < cursors_->size()) base_ <<= 1;
        tree_.assign(base_ * 2, -1);
        for (size_t i = 0; i < cursors_->size(); ++i) {
            if ((*cursors_)[i].active) {
                tree_[base_ + i] = (int)i;
            }
        }
        for (size_t node = base_; node-- > 1;) {
            tree_[node] = Winner(
                tree_[node << 1],
                tree_[(node << 1) | 1]);
        }
    }

    int winner() const {
        return tree_.size() > 1 ? tree_[1] : -1;
    }

    void replay(int player) {
        size_t node = base_ + (size_t)player;
        tree_[node] =
            Active(player) ? player : -1;
        while ((node >>= 1) > 0) {
            tree_[node] = Winner(
                tree_[node << 1],
                tree_[(node << 1) | 1]);
        }
    }

private:
    bool Active(int player) const {
        return player >= 0 &&
               (size_t)player < cursors_->size() &&
               (*cursors_)[(size_t)player].active;
    }

    bool Less(int lhs, int rhs) const {
        const CollateFileCursor &a =
            (*cursors_)[(size_t)lhs];
        const CollateFileCursor &b =
            (*cursors_)[(size_t)rhs];
        CollateMeta a_meta = CollateCursorMeta(a);
        CollateMeta b_meta = CollateCursorMeta(b);
        const unsigned char *a_raw =
            a.raw_buffer + (size_t)(
                a_meta.raw_offset - a.raw_begin);
        const unsigned char *b_raw =
            b.raw_buffer + (size_t)(
                b_meta.raw_offset - b.raw_begin);
        a_meta.raw_offset = 0;
        b_meta.raw_offset = 0;
        return CollateCompare(
                   a_meta, a_raw, a_meta.raw_len,
                   b_meta, b_raw, b_meta.raw_len) < 0;
    }

    int Winner(int lhs, int rhs) const {
        if (!Active(lhs)) return Active(rhs) ? rhs : -1;
        if (!Active(rhs)) return lhs;
        if (Less(lhs, rhs)) return lhs;
        if (Less(rhs, lhs)) return rhs;
        return lhs < rhs ? lhs : rhs;
    }

    std::vector<CollateFileCursor> *cursors_;
    size_t base_;
    std::vector<int> tree_;
};

static int CollateBuildCursors(
        const std::vector<CollateSegment> &segments,
        size_t begin, size_t end,
        unsigned char *cursor_memory,
        unsigned char *scratch, size_t scratch_size,
        MpiCollateStats *stats,
        std::vector<CollateFileCursor> *cursors) {
    const size_t meta_capacity =
        kCollateCursorMeta / sizeof(CollateMeta);
    const size_t stride =
        kCollateCursorMeta + kCollateCursorRaw;
    cursors->clear();
    cursors->reserve(end - begin);
    for (size_t i = begin; i < end; ++i) {
        if (!segments[i].record_count) continue;
        size_t slot = cursors->size();
        unsigned char *base =
            cursor_memory + slot * stride;
        CollateFileCursor cursor;
        cursor.segment = &segments[i];
        cursor.meta_buffer =
            reinterpret_cast<CollateMeta *>(base);
        cursor.meta_capacity = meta_capacity;
        cursor.raw_buffer =
            base + kCollateCursorMeta;
        cursor.raw_capacity = kCollateCursorRaw;
        if (CollateCursorLoadMeta(
                &cursor, scratch, scratch_size,
                stats) != 0) {
            return -1;
        }
        if (cursor.active) {
            const unsigned char *raw = nullptr;
            if (CollateCursorRaw(
                    &cursor, scratch, scratch_size,
                    stats, &raw) != 0) {
                return -1;
            }
            (void)raw;
        }
        cursors->push_back(cursor);
    }
    return 0;
}

static int CollateMergeGroupToTemp(
        const std::vector<CollateSegment> &segments,
        size_t begin, size_t end,
        unsigned char *cursor_memory,
        std::vector<CollateMeta> *out_meta,
        std::vector<unsigned char> *out_raw,
        CollateTempStore *store,
        unsigned char *scratch, size_t scratch_size,
        MpiCollateStats *stats,
        CollateSegment *output) {
    std::vector<CollateFileCursor> cursors;
    if (CollateBuildCursors(
            segments, begin, end, cursor_memory,
            scratch, scratch_size, stats,
            &cursors) != 0) {
        return -1;
    }
    CollateFileLoserTree tree(&cursors);
    size_t meta_count = 0;
    size_t raw_used = 0;
    auto flush = [&]() -> int {
        if (!meta_count) return 0;
        int ret = CollateAppendSegmentChunk(
            store, out_meta->data(), meta_count,
            out_raw->data(), raw_used,
            scratch, scratch_size, stats, output);
        meta_count = 0;
        raw_used = 0;
        return ret;
    };
    for (;;) {
        int winner = tree.winner();
        if (winner < 0) break;
        CollateFileCursor &cursor =
            cursors[(size_t)winner];
        const CollateMeta meta =
            CollateCursorMeta(cursor);
        const unsigned char *raw = nullptr;
        if (CollateCursorRaw(
                &cursor, scratch, scratch_size,
                stats, &raw) != 0) {
            return -1;
        }
        if (meta_count == out_meta->size() ||
            meta.raw_len >
                out_raw->size() - raw_used) {
            if (flush() != 0) return -1;
        }
        if (meta.raw_len > out_raw->size()) return -1;
        (*out_meta)[meta_count] = meta;
        (*out_meta)[meta_count].raw_offset = raw_used;
        memcpy(
            out_raw->data() + raw_used,
            raw, meta.raw_len);
        raw_used += meta.raw_len;
        ++meta_count;
        if (CollateCursorAdvance(
                &cursor, scratch, scratch_size,
                stats) != 0) {
            return -1;
        }
        if (cursor.active) {
            if (CollateCursorRaw(
                    &cursor, scratch,
                    scratch_size, stats,
                    &raw) != 0) {
                return -1;
            }
        }
        tree.replay(winner);
    }
    return flush();
}

static int CollateConsolidate(
        std::vector<CollateSegment> *segments,
        size_t fan_in,
        CollateTempStore *store,
        unsigned char *scratch, size_t scratch_size,
        MpiCollateStats *stats) {
    if (segments->size() <= fan_in) return 0;
    const size_t stride =
        kCollateCursorMeta + kCollateCursorRaw;
    std::vector<unsigned char> cursor_memory(
        fan_in * stride);
    std::vector<CollateMeta> out_meta(
        kCollateConsolidateMeta /
        sizeof(CollateMeta));
    std::vector<unsigned char> out_raw(
        kCollateConsolidateRaw);
    while (segments->size() > fan_in) {
        std::vector<CollateSegment> next;
        next.reserve(
            (segments->size() + fan_in - 1) /
            fan_in);
        double t0 = GetTime();
        for (size_t begin = 0;
             begin < segments->size();
             begin += fan_in) {
            size_t end = std::min(
                segments->size(),
                begin + fan_in);
            CollateSegment merged;
            if (CollateMergeGroupToTemp(
                    *segments, begin, end,
                    cursor_memory.data(),
                    &out_meta, &out_raw, store,
                    scratch, scratch_size,
                    stats, &merged) != 0) {
                return -1;
            }
            next.push_back(std::move(merged));
        }
        stats->t_merge += GetTime() - t0;
        ++stats->consolidation_passes;
        segments->swap(next);
    }
    return 0;
}

static int CollateMergeExternal(
        const std::vector<CollateSegment> &segments,
        int compress_level, size_t fan_in,
        unsigned char *scratch, size_t scratch_size,
        MemWriter *writer, MpiCollateStats *stats) {
    const size_t stride =
        kCollateCursorMeta + kCollateCursorRaw;
    std::vector<unsigned char> cursor_memory(
        fan_in * stride);
    std::vector<CollateFileCursor> cursors;
    if (CollateBuildCursors(
            segments, 0, segments.size(),
            cursor_memory.data(), scratch,
            scratch_size, stats, &cursors) != 0) {
        return -1;
    }
    CollateFileLoserTree tree(&cursors);
    int pending = -1;
    std::string previous;
    uint32_t previous_hash = 0;
    auto next = [&](const unsigned char **raw,
                    uint32_t *length) -> int {
        if (pending >= 0) {
            CollateFileCursor &cursor =
                cursors[(size_t)pending];
            if (CollateCursorAdvance(
                    &cursor, scratch, scratch_size,
                    stats) != 0) {
                return -1;
            }
            const unsigned char *unused = nullptr;
            if (cursor.active &&
                CollateCursorRaw(
                    &cursor, scratch,
                    scratch_size, stats,
                    &unused) != 0) {
                return -1;
            }
            tree.replay(pending);
            pending = -1;
        }
        int winner = tree.winner();
        if (winner < 0) return 0;
        CollateFileCursor &cursor =
            cursors[(size_t)winner];
        const CollateMeta &meta =
            CollateCursorMeta(cursor);
        if (CollateCursorRaw(
                &cursor, scratch, scratch_size,
                stats, raw) != 0) {
            return -1;
        }
        CollateCountGroup(
            meta, *raw, &previous,
            &previous_hash, stats);
        *length = meta.raw_len;
        pending = winner;
        return 1;
    };
    double t0 = GetTime();
    int ret = CollateCompressStream(
        next, compress_level, writer, stats);
    stats->t_merge += GetTime() - t0;
    return ret;
}

static int FusedBamExternalCollateMPI(
        MemReader &reader, MemWriter &writer,
        long long global_block_begin,
        int rank, int comm_size, int bins,
        int compress_level, size_t memory_limit,
        const std::string &temp_prefix,
        MpiCollateStats *stats) {
    const size_t exchange_workspace =
        2ull * (8ull * 1024 * 1024) +
        2ull * (1ull * 1024 * 1024);
    const size_t phase_workspace = std::max(
        CollateExtractWorkspaceBytes(),
        std::max(exchange_workspace,
                 CollateCompressWorkspaceBytes()));
    const size_t fixed =
        kCollateControlReserve +
        kCollateSimScratch +
        phase_workspace;
    const size_t minimum_arena =
        32ull * 1024 * 1024;
    if (memory_limit <= fixed + minimum_arena) {
        if (rank == 0) {
            fprintf(
                stderr,
                "ERROR: collate external mode needs at least "
                "%zu bytes per rank; requested %zu.\n",
                fixed + minimum_arena,
                memory_limit);
        }
        return -1;
    }

    CollateTempStore store;
    CollateRunArena arena;
    std::vector<CollateRun> runs;
    std::vector<CollateSegment> segments;
    unsigned char *scratch = nullptr;
    int result = -1;
    double wall_t0 = GetTime();
    arena.capacity = memory_limit - fixed;
    arena.capacity &= ~(size_t)63;
    arena.data = (unsigned char *)aligned_alloc_custom(
        64, arena.capacity);
    scratch = (unsigned char *)aligned_alloc_custom(
        64, kCollateSimScratch);
    int local_ready =
        arena.data && scratch &&
        CollateOpenTemp(
            temp_prefix, rank, &store) == 0;
    if (!CollateAllRanksOk(local_ready)) {
        goto cleanup;
    }
    arena.reset();
    stats->run_arena_bytes =
        (long long)arena.capacity;
    stats->tracked_peak_bytes =
        (long long)(arena.capacity + fixed);

    {
        int local_extract_ok =
            CollateExtractRuns(
                reader, global_block_begin, bins,
                &arena, &store, scratch,
                kCollateSimScratch, &runs,
                stats) == 0;
        if (!CollateAllRanksOk(local_extract_ok)) {
            goto cleanup;
        }
    }
    if (CollateExchangeRuns(
            runs, &arena, &store,
            rank, comm_size, bins,
            scratch, kCollateSimScratch,
            &segments, stats) != 0) {
        goto cleanup;
    }
    aligned_free_custom(arena.data);
    arena.data = nullptr;

    {
        const size_t stride =
            kCollateCursorMeta + kCollateCursorRaw;
        const size_t merge_fixed =
            kCollateControlReserve +
            kCollateSimScratch +
            CollateCompressWorkspaceBytes() +
            kCollateConsolidateMeta +
            kCollateConsolidateRaw;
        if (memory_limit <= merge_fixed + 2 * stride) {
            goto cleanup;
        }
        size_t fan_in =
            (memory_limit - merge_fixed) / stride;
        fan_in = std::max<size_t>(2, fan_in);
        fan_in = std::min<size_t>(
            fan_in,
            std::max<size_t>(2, segments.size()));
        stats->merge_fan_in = (long long)fan_in;
        stats->tracked_peak_bytes =
            std::max<long long>(
                stats->tracked_peak_bytes,
                (long long)(merge_fixed +
                            fan_in * stride));
        if (CollateConsolidate(
                &segments, fan_in, &store,
                scratch, kCollateSimScratch,
                stats) != 0 ||
            CollateMergeExternal(
                segments, compress_level, fan_in,
                scratch, kCollateSimScratch,
                &writer, stats) != 0) {
            goto cleanup;
        }
    }
    result = 0;

cleanup:
    stats->t_actual = GetTime() - wall_t0;
    stats->t_fused =
        stats->t_actual -
        stats->t_temp_read_actual -
        stats->t_temp_write_actual;
    if (stats->t_fused < 0.0) stats->t_fused = 0.0;
    CollateCloseTemp(&store);
    if (arena.data) aligned_free_custom(arena.data);
    if (scratch) aligned_free_custom(scratch);
    return result;
}

static int FusedBamMemoryCollateMPI(
        MemReader &reader, MemWriter &writer,
        long long global_block_begin,
        int rank, int comm_size, int bins,
        int compress_level, size_t memory_limit,
        MpiCollateStats *stats) {
    double wall_t0 = GetTime();
    std::vector<CollateMeta> records;
    std::vector<unsigned char> raw;
    std::vector<CollateMemorySegment> segments;
    int local_extract_ok = CollateExtractAll(
            reader, global_block_begin, bins,
            &records, &raw, stats) == 0;
    if (!CollateAllRanksOk(local_extract_ok)) {
        return -1;
    }
    double t0 = GetTime();
    std::sort(
        records.begin(), records.end(),
        CollateLocalLess{&raw});
    stats->t_local_sort += GetTime() - t0;
    if (CollateExchangeMemory(
            records, raw, rank, comm_size,
            bins, &segments, stats) != 0) {
        return -1;
    }
    size_t segment_bytes = 0;
    for (size_t i = 0; i < segments.size(); ++i) {
        if (segments[i].records.capacity() >
                (SIZE_MAX - segment_bytes) /
                    sizeof(CollateMeta)) {
            return -1;
        }
        segment_bytes +=
            segments[i].records.capacity() *
            sizeof(CollateMeta);
        if (segments[i].raw.capacity() >
            SIZE_MAX - segment_bytes) {
            return -1;
        }
        segment_bytes += segments[i].raw.capacity();
    }
    size_t tracked =
        records.capacity() * sizeof(CollateMeta);
    if (raw.capacity() > SIZE_MAX - tracked ||
        segment_bytes >
            SIZE_MAX - tracked - raw.capacity()) {
        return -1;
    }
    tracked += raw.capacity() + segment_bytes;
    if (tracked <=
        SIZE_MAX - kCollateControlReserve) {
        tracked += kCollateControlReserve;
    }
    stats->tracked_peak_bytes =
        tracked > (size_t)LLONG_MAX
            ? LLONG_MAX : (long long)tracked;
    if (memory_limit && tracked > memory_limit) {
        fprintf(
            stderr,
            "[rank %d] ERROR: collate memory mode "
            "tracked workset=%zu exceeds -m=%zu.\n",
            rank, tracked, memory_limit);
        return -1;
    }
    stats->segments = (long long)segments.size();
    if (CollateMergeMemory(
            segments, compress_level,
            &writer, stats) != 0) {
        return -1;
    }
    stats->t_actual = GetTime() - wall_t0;
    stats->t_fused = stats->t_actual;
    return 0;
}

static void CollatePrintStats(
        const MpiCollateStats &stats,
        int rank, int comm_size, int bins) {
    char local[4096] = {};
    snprintf(
        local, sizeof(local),
        "[rank %d] mode=%s blocks=%lld records=%lld "
        "received=%lld groups=%lld collisions=%lld bins=%d "
        "runs=%lld segments=%lld bgzf=%lld\n"
        "[rank %d] collate_time extract=%.3f sort=%.3f "
        "exchange=%.3f mpi=%.3f merge=%.3f compress=%.3f "
        "write=%.3f fused=%.3f actual=%.3f\n"
        "[rank %d] collate_temp write_actual=%.3f "
        "read_actual=%.3f write_sim=%.3f read_sim=%.3f "
        "write_bytes=%lld read_bytes=%lld\n"
        "[rank %d] collate_memory resident_records=%lld "
        "resident_raw=%lld peak=%lld arena=%lld fan_in=%lld "
        "consolidation_passes=%lld\n",
        rank, stats.mode ? "external" : "memory",
        stats.input_blocks, stats.total_records,
        stats.received_records, stats.qname_groups,
        stats.hash_collision_groups, bins,
        stats.runs, stats.segments,
        stats.bgzf_blocks,
        rank, stats.t_extract, stats.t_local_sort,
        stats.t_exchange, stats.t_mpi,
        stats.t_merge, stats.t_compress,
        stats.t_write, stats.t_fused,
        stats.t_actual,
        rank, stats.t_temp_write_actual,
        stats.t_temp_read_actual,
        stats.t_temp_write_sim,
        stats.t_temp_read_sim,
        stats.temp_write_bytes,
        stats.temp_read_bytes,
        rank, stats.resident_run_records,
        stats.resident_run_raw_bytes,
        stats.tracked_peak_bytes,
        stats.run_arena_bytes,
        stats.merge_fan_in,
        stats.consolidation_passes);
    const int width = 4096;
    std::vector<char> gathered;
    if (rank == 0) {
        gathered.resize((size_t)comm_size * width);
    }
    MPI_Gather(
        local, width, MPI_CHAR,
        rank == 0 ? gathered.data() : nullptr,
        width, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < comm_size; ++i) {
            fputs(
                gathered.data() +
                    (size_t)i * width,
                stdout);
        }
        fflush(stdout);
    }
}

} // namespace

int ProcessCollateMPI(CmdInfo *cmd_info) {
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
    MemWriter writer = {};
    MpiCollateStats stats = {};
    size_t memory_limit = 0;
    std::string temp_prefix;

    if (cmd_info->collate_bins_ <= 0 ||
        cmd_info->out_file_name_ == "-" ||
        CollateParseMemory(
            cmd_info->collate_memory_,
            &memory_limit) != 0) {
        if (rank == 0) {
            fprintf(
                stderr,
                "ERROR: invalid collate bins or memory limit.\n");
        }
        local_ok = 0;
    }
    if (!cmd_info->collate_temp_prefix_.empty()) {
        temp_prefix =
            cmd_info->collate_temp_prefix_;
    } else {
        size_t slash =
            cmd_info->out_file_name_.find_last_of('/');
        temp_prefix =
            slash == std::string::npos
                ? "./rabbitbam-collate"
                : cmd_info->out_file_name_.substr(
                      0, slash + 1) +
                      "rabbitbam-collate";
    }
    if (!CollateAllRanksOk(local_ok)) goto cleanup;
    {
        double cost =
            CollateReduceMax(GetTime() - total_t0);
        if (rank == 0) {
            printf(
                "111Complete the initialization cost %lf-----\n",
                cost);
        }
    }

    {
        double t0 = GetTime();
        if (MpiCommonLoadFileToMemory(
                cmd_info->in_file_name_,
                &input_memory, &input_size) != 0) {
            fprintf(
                stderr,
                "[rank %d] ERROR: cannot preload collate input %s.\n",
                rank,
                cmd_info->in_file_name_.c_str());
            local_ok = 0;
        }
        double cost =
            CollateReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf(
                "222Complete the memory cost %lf--\n",
                cost);
        }
    }
    if (!CollateAllRanksOk(local_ok)) goto cleanup;

    {
        double t0 = GetTime();
        input_hfile = hopen(
            "mem:", "rb:", input_memory, input_size);
        if (input_hfile) {
            input = (samFile *)hts_hopen(
                input_hfile, "data", "rb");
            if (input) input_hfile = nullptr;
        }
        if (!input) local_ok = 0;
        if (local_ok) header = sam_hdr_read(input);
        if (!header) local_ok = 0;
        if (local_ok &&
            input->format.format != bam &&
            input->format.format != binary_format) {
            if (rank == 0) {
                fprintf(
                    stderr,
                    "ERROR: RabbitBAM-MPI collate only supports BAM input.\n");
            }
            local_ok = 0;
        }
        if (local_ok) {
            body_start =
                (long long)input->fp.bgzf->block_address;
            if (body_start < 0 ||
                (unsigned long long)body_start >
                    input_size) {
                local_ok = 0;
            }
        }
        double cost =
            CollateReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf(
                "333Complete the head cost %lf---\n",
                cost);
        }
    }
    if (!CollateAllRanksOk(local_ok)) goto cleanup;

    {
        double t0 = GetTime();
        if (rank == 0) {
            printf(
                "Enable MPI BAM COLLATE mode "
                "(%d MPE + %d CPEs)!!!\n",
                comm_size, comm_size * 64);
            printf(
                "MPI BAM collate bins=%d compression=%d\n",
                cmd_info->collate_bins_,
                cmd_info->compress_level_);
            if (MpiCommonScanBgzfBlocksInMemory(
                    input_memory, input_size,
                    body_start, &block_offsets,
                    &block_lengths) != 0) {
                local_ok = 0;
            }
            n_blocks =
                (long long)block_offsets.size();
            if (n_blocks > INT_MAX) local_ok = 0;
        }
        MPI_Bcast(
            &local_ok, 1, MPI_INT, 0,
            MPI_COMM_WORLD);
        if (!local_ok) goto cleanup;
        MPI_Bcast(
            &n_blocks, 1, MPI_LONG_LONG, 0,
            MPI_COMM_WORLD);
        if (rank != 0) {
            block_offsets.resize((size_t)n_blocks);
            block_lengths.resize((size_t)n_blocks);
        }
        if (n_blocks > 0) {
            MPI_Bcast(
                block_offsets.data(), (int)n_blocks,
                MPI_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(
                block_lengths.data(), (int)n_blocks,
                MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        }
        local_begin =
            n_blocks * rank / comm_size;
        local_end =
            n_blocks * (rank + 1) / comm_size;
        if (MpiCommonSelectBlockRangeFromMemory(
                input_memory, input_size,
                block_offsets, block_lengths,
                local_begin, local_end,
                &rank_input,
                &rank_input_size) != 0) {
            local_ok = 0;
        }
        double cost =
            CollateReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf(
                "MPI BAM scan complete. data_blocks=%lld "
                "body_start=%lld\n",
                n_blocks, body_start);
            printf(
                "Complete the 4.1 scan/split/broadcast/select cost %lf\n",
                cost);
        }
    }
    if (!CollateAllRanksOk(local_ok)) goto cleanup;

    reader.base = rank_input;
    reader.size = rank_input_size;
    reader.pos = 0;
    if (MpiCommonInitMemWriter(
            writer, rank_input_size
                ? rank_input_size
                : 1024 * 1024) != 0) {
        local_ok = 0;
    }
    if (!CollateAllRanksOk(local_ok)) goto cleanup;

    {
        size_t estimate =
            rank_input_size >
                    (SIZE_MAX -
                     256ull * 1024 * 1024) / 6
                ? SIZE_MAX
                : rank_input_size * 6 +
                      256ull * 1024 * 1024;
        unsigned long long local_estimate =
            estimate == SIZE_MAX
                ? ULLONG_MAX
                : (unsigned long long)estimate;
        unsigned long long max_estimate = 0;
        MPI_Allreduce(
            &local_estimate, &max_estimate, 1,
            MPI_UNSIGNED_LONG_LONG, MPI_MAX,
            MPI_COMM_WORLD);
        int local_external =
            memory_limit != 0 &&
            estimate > memory_limit;
        int use_external = 0;
        MPI_Allreduce(
            &local_external, &use_external, 1,
            MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        stats.mode = use_external;
        if (rank == 0) {
            printf(
                "MPI BAM collate selected mode=%s "
                "memory_estimate_max=%llu memory_limit=%zu\n",
                use_external ? "external" : "memory",
                max_estimate, memory_limit);
            printf(
                "MPI BAM collate mode policy: no -m => memory; "
                "otherwise external when max(6 * rank_input_bytes "
                "+ 256MiB) exceeds -m\n");
        }
        double t0 = GetTime();
        int ret = use_external
            ? FusedBamExternalCollateMPI(
                  reader, writer, local_begin,
                  rank, comm_size,
                  cmd_info->collate_bins_,
                  cmd_info->compress_level_,
                  memory_limit, temp_prefix, &stats)
            : FusedBamMemoryCollateMPI(
                  reader, writer, local_begin,
                  rank, comm_size,
                  cmd_info->collate_bins_,
                  cmd_info->compress_level_,
                  memory_limit,
                  &stats);
        if (ret != 0) local_ok = 0;
        int global_ok =
            CollateAllRanksOk(local_ok);
        double modeled =
            use_external
                ? stats.t_fused
                : GetTime() - t0;
        double modeled_max =
            CollateReduceMax(modeled);
        double actual_max =
            CollateReduceMax(stats.t_actual);
        if (rank == 0 && global_ok) {
            printf(
                "Complete the 4.3 FusedBamCollateMPI cost %lf\n",
                modeled_max);
            printf(
                "FusedBamCollateMPI actual wall %lf\n",
                actual_max);
        }
        if (!global_ok) goto cleanup;
    }

    CollatePrintStats(
        stats, rank, comm_size,
        cmd_info->collate_bins_);
    if (CollateGatherOutput(
            writer, header,
            cmd_info->compress_level_,
            cmd_info->out_file_name_,
            rank, comm_size, &stats) != 0) {
        local_ok = 0;
    }
    if (!CollateAllRanksOk(local_ok)) goto cleanup;
    exit_code = 0;

cleanup:
    if (input) {
        if (sam_close(input) < 0) exit_code = 1;
        input = nullptr;
    }
    if (input_hfile &&
        hclose(input_hfile) != 0) {
        exit_code = 1;
    }
    if (header) sam_hdr_destroy(header);
    free(input_memory);
    free(writer.data);
    if (rank == 0) {
        printf(
            "666collate total process cost %lf-----\n",
            GetTime() - total_t0);
    }
    return exit_code;
}
