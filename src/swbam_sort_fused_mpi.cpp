#include "swbam_mpi.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <vector>

#include <mpi.h>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

extern "C" {
    void slave_mpi_sort_extract_raw();
    void slave_mpi_sort_compress_payload();
    void slave_mpi_sort_bucket_pack();
}

namespace {

const int kSortSampleFactor = 64;
const int kSortTagMeta = 7101;
const int kSortTagRaw = 7102;
const int kSortExchangeChunk = 64 * 1024 * 1024;
const int kSortCpeBlocks = 64;

//用于排序的结构体key
struct MpiSortKey {
    int32_t tid;
    int32_t pos;
    uint16_t flag;
    uint16_t pad;
    uint64_t global_order;
};

typedef MpiSortRecordMetaShared MpiSortRecordMeta;

struct MpiSortBlockSet {
    bam_block *blocks;
    unsigned char *data;
    int n;
};

struct MpiSortRawBuffer {
    unsigned char *data;
    size_t size;
    size_t capacity;

    MpiSortRawBuffer() : data(nullptr), size(0), capacity(0) {}
    ~MpiSortRawBuffer() {
        if (data) aligned_free_custom(data);
        data = nullptr;
        size = 0;
        capacity = 0;
    }

private:
    MpiSortRawBuffer(const MpiSortRawBuffer &);
    MpiSortRawBuffer &operator=(const MpiSortRawBuffer &);
};

struct MpiSortExtractWorkspace {
    MpiSortBlockSet input_blocks;
    MpiSortBlockSet un_blocks;
    unsigned char *raw_data;
    MpiSortRecordMeta *records;
    int n_blocks;
    int records_per_block;
    size_t raw_stride;
};

uint16_t MpiSortReadLe16(const unsigned char *p) {
    return (uint16_t)p[0] | ((uint16_t)p[1] << 8);
}

uint32_t MpiSortReadLe32(const unsigned char *p) {
    return (uint32_t)p[0] |
           ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) |
           ((uint32_t)p[3] << 24);
}

uint32_t MpiSortBgzfISize(const bam_block *block) {
    if (!block || block->length < BLOCK_FOOTER_LENGTH) return 0;
    return MpiSortReadLe32((const unsigned char *)block->data + block->length - 4);
}

int MpiSortRawReserve(MpiSortRawBuffer *buf, size_t needed) {
    if (!buf) return -1;
    if (needed <= buf->capacity) return 0;
    size_t new_capacity = buf->capacity ? buf->capacity + buf->capacity / 2 : (size_t)(64 * 1024 * 1024);
    if (new_capacity < needed) new_capacity = needed;
    unsigned char *new_data = aligned_alloc_custom(64, new_capacity);
    if (!new_data) return -1;
    if (buf->data && buf->size > 0) {
        memcpy(new_data, buf->data, buf->size);
    }
    if (buf->data) aligned_free_custom(buf->data);
    buf->data = new_data;
    buf->capacity = new_capacity;
    return 0;
}

void MpiSortRawRelease(MpiSortRawBuffer *buf) {
    if (!buf) return;
    if (buf->data) aligned_free_custom(buf->data);
    buf->data = nullptr;
    buf->size = 0;
    buf->capacity = 0;
}

uint64_t MpiSortTidKey(int32_t tid) {
    return tid < 0 ? UINT64_MAX : (uint64_t)(uint32_t)tid;
}

uint64_t MpiSortPosKey(int32_t pos) {
    return pos < 0 ? 0 : (uint64_t)pos + 1u;
}

/*用于排序的函数：
mapped records 在前，unmapped records 在后
mapped records 按 tid、pos 排
同 tid/pos 下，forward 在 reverse 前
再相同，就按原始输入顺序 global_order 排*/
int MpiSortCompareFields(int32_t atid, int32_t apos, uint16_t aflag, uint64_t aorder,
                         int32_t btid, int32_t bpos, uint16_t bflag, uint64_t border) {
    uint64_t av = MpiSortTidKey(atid);
    uint64_t bv = MpiSortTidKey(btid);
    if (av != bv) return av < bv ? -1 : 1;
    av = MpiSortPosKey(apos);
    bv = MpiSortPosKey(bpos);
    if (av != bv) return av < bv ? -1 : 1;
    av = (aflag & BAM_FREVERSE) ? 1u : 0u;
    bv = (bflag & BAM_FREVERSE) ? 1u : 0u;
    if (av != bv) return av < bv ? -1 : 1;
    if (aorder != border) return aorder < border ? -1 : 1;
    return 0;
}

int MpiSortCompareMetaKey(const MpiSortRecordMeta &a, const MpiSortKey &b) {
    return MpiSortCompareFields(a.tid, a.pos, a.flag, a.global_order,
                                b.tid, b.pos, b.flag, b.global_order);
}

struct MpiSortMetaLess {
    bool operator()(const MpiSortRecordMeta &a, const MpiSortRecordMeta &b) const {
        return MpiSortCompareFields(a.tid, a.pos, a.flag, a.global_order,
                                    b.tid, b.pos, b.flag, b.global_order) < 0;
    }
};

struct MpiSortKeyLess {
    bool operator()(const MpiSortKey &a, const MpiSortKey &b) const {
        return MpiSortCompareFields(a.tid, a.pos, a.flag, a.global_order,
                                    b.tid, b.pos, b.flag, b.global_order) < 0;
    }
};

MpiSortKey MpiSortKeyFromMeta(const MpiSortRecordMeta &m) {
    MpiSortKey key;
    key.tid = m.tid;
    key.pos = m.pos;
    key.flag = m.flag;
    key.pad = 0;
    key.global_order = m.global_order;
    return key;
}

//内存限制检查：
//用于检查当前阶段估算内存是否超过 -m 参数。
int MpiSortCheckMemoryLimit(size_t limit, size_t estimate, int rank, const char *stage) {
    if (limit == 0 || estimate <= limit) return 0;
    fprintf(stderr,
            "[rank %d] ERROR: sort v1 requires enough memory for %s. estimate=%zu limit=%zu; external merge is not implemented.\n",
            rank, stage, estimate, limit);
    return -1;
}

int MpiSortSendRecvBytes(int send_to, const char *send_data, long long send_len,
                         int recv_from, char *recv_data, long long recv_len,
                         int tag) {
    long long sent = 0;
    long long received = 0;
    MPI_Status status;
    while (sent < send_len || received < recv_len) {
        int send_chunk = 0;
        int recv_chunk = 0;
        if (sent < send_len) {
            send_chunk = (int)std::min<long long>(send_len - sent, kSortExchangeChunk);
        }
        if (received < recv_len) {
            recv_chunk = (int)std::min<long long>(recv_len - received, kSortExchangeChunk);
        }
        MPI_Sendrecv(send_chunk ? (void *)(send_data + sent) : nullptr, send_chunk, MPI_BYTE,
                     send_to, tag,
                     recv_chunk ? (void *)(recv_data + received) : nullptr, recv_chunk, MPI_BYTE,
                     recv_from, tag, MPI_COMM_WORLD, &status);
        sent += send_chunk;
        received += recv_chunk;
    }
    return 0;
}

int MpiSortAllocateBlockSet(MpiSortBlockSet *set, int n) {
    set->n = n;
    set->blocks = (bam_block *)aligned_alloc_custom(64, (size_t)n * sizeof(bam_block));
    set->data = aligned_alloc_custom(64, (size_t)n * BGZF_MAX_BLOCK_SIZE);
    if (!set->blocks || !set->data) {
        if (set->blocks) aligned_free_custom((unsigned char *)set->blocks);
        if (set->data) aligned_free_custom(set->data);
        set->blocks = nullptr;
        set->data = nullptr;
        set->n = 0;
        return -1;
    }
    memset(set->blocks, 0, (size_t)n * sizeof(bam_block));
    for (int i = 0; i < n; ++i) {
        set->blocks[i].data = set->data + (size_t)i * BGZF_MAX_BLOCK_SIZE;
        set->blocks[i].block_id = i;
    }
    return 0;
}

void MpiSortFreeBlockSet(MpiSortBlockSet *set) {
    if (set->blocks) aligned_free_custom((unsigned char *)set->blocks);
    if (set->data) aligned_free_custom(set->data);
    set->blocks = nullptr;
    set->data = nullptr;
    set->n = 0;
}

int MpiSortAllocateExtractWorkspace(MpiSortExtractWorkspace *ws,
                                    int n_blocks,
                                    int records_per_block,
                                    size_t raw_stride) {
    ws->n_blocks = n_blocks;
    ws->records_per_block = records_per_block;
    ws->raw_stride = raw_stride;
    ws->raw_data = nullptr;
    ws->records = nullptr;
    if (MpiSortAllocateBlockSet(&ws->input_blocks, n_blocks) != 0 ||
        MpiSortAllocateBlockSet(&ws->un_blocks, n_blocks) != 0) {
        return -1;
    }
    ws->records = (MpiSortRecordMeta *)aligned_alloc_custom(
        64, (size_t)n_blocks * records_per_block * sizeof(MpiSortRecordMeta));
    if (!ws->records) {
        MpiSortFreeBlockSet(&ws->input_blocks);
        MpiSortFreeBlockSet(&ws->un_blocks);
        if (ws->records) aligned_free_custom((unsigned char *)ws->records);
        ws->raw_data = nullptr;
        ws->records = nullptr;
        return -1;
    }
    return 0;
}

void MpiSortFreeExtractWorkspace(MpiSortExtractWorkspace *ws) {
    MpiSortFreeBlockSet(&ws->input_blocks);
    MpiSortFreeBlockSet(&ws->un_blocks);
    if (ws->raw_data) aligned_free_custom(ws->raw_data);
    if (ws->records) aligned_free_custom((unsigned char *)ws->records);
    ws->raw_data = nullptr;
    ws->records = nullptr;
    ws->n_blocks = 0;
    ws->records_per_block = 0;
    ws->raw_stride = 0;
}

//从内存里读取 BGZF block
int MpiSortMemReadBlock(char *base, size_t size, size_t &pos, bam_block *block) {
    if (pos >= size) return -1;
    if (pos + BLOCK_HEADER_LENGTH > size) return -1;
    uint16_t bsize = MpiSortReadLe16((const unsigned char *)base + pos + 16);
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

void MpiSortAccumulateExtractDetail(const MpiSortExtractPara *paras,
                                    int active_blocks,
                                    double wall,
                                    MpiSortStats *stats) {
    if (!stats || active_blocks <= 0 || wall <= 0.0) return;
    const MpiSortExtractPara *critical = nullptr;
    uint64_t critical_total = 0;
    for (int b = 0; b < active_blocks; ++b) {
        if (paras[b].decomp_total_cycles >= critical_total) {
            critical_total = paras[b].decomp_total_cycles;
            critical = &paras[b];
        }
    }
    if (!critical || critical_total == 0) {
        stats->t_extract_other += wall;
        return;
    }
    const double scale = wall / (double)critical_total;
    double alloc_time = scale * (double)critical->decomp_alloc_cycles;
    double inflate_time = scale * (double)critical->decomp_inflate_cycles;
    double crc_time = scale * (double)critical->decomp_crc_cycles;
    double parse_time = scale * (double)critical->decomp_parse_cycles;
    double other_time = wall - alloc_time - inflate_time - crc_time - parse_time;
    if (other_time < 0.0) other_time = 0.0;
    stats->t_extract_alloc += alloc_time;
    stats->t_extract_inflate += inflate_time;
    stats->t_extract_crc += crc_time;
    stats->t_extract_parse += parse_time;
    stats->t_extract_other += other_time;
}

void MpiSortAccumulateCompressDetail(const MpiSortRawCompressPara *paras,
                                     int active_blocks,
                                     double wall,
                                     MpiSortStats *stats) {
    if (!stats || active_blocks <= 0 || wall <= 0.0) return;
    const MpiSortRawCompressPara *critical = nullptr;
    uint64_t critical_total = 0;
    for (int b = 0; b < active_blocks; ++b) {
        if (paras[b].compress_total_cycles >= critical_total) {
            critical_total = paras[b].compress_total_cycles;
            critical = &paras[b];
        }
    }
    if (!critical || critical_total == 0) {
        stats->t_compress_other += wall;
        return;
    }
    const double scale = wall / (double)critical_total;
    double pack_time = scale * (double)critical->compress_pack_cycles;
    double alloc_time = scale * (double)critical->compress_alloc_cycles;
    double deflate_time = scale * (double)critical->compress_deflate_cycles;
    double footer_time = scale * (double)critical->compress_footer_cycles;
    double other_time = wall - pack_time - alloc_time - deflate_time - footer_time;
    if (other_time < 0.0) other_time = 0.0;
    stats->t_compress_pack += pack_time;
    stats->t_compress_alloc += alloc_time;
    stats->t_compress_deflate += deflate_time;
    stats->t_compress_footer += footer_time;
    stats->t_compress_other += other_time;
}

// CPE 批量提取 BAM records
int MpiSortExtractLocalRecordsCpe(MemReader &reader,
                                  long long global_block_begin,
                                  std::vector<MpiSortRecordMeta> *records,
                                  MpiSortRawBuffer *raw_data,
                                  MpiSortStats *stats) {
    // 一次处理 64 个 block。输入 block 使用双缓冲：CPE 处理当前批次时，
    // MPE 读取下一批；CPE 直接解压到最终 local_raw 位置，metadata scratch 复用一套。
    MpiSortExtractWorkspace ws = {};
    MpiSortBlockSet input_pending_set = {};
    MpiSortExtractPara paras[kSortCpeBlocks];
    const int records_per_block = (int)MPI_RECORDS_PER_BLOCK;
    const size_t raw_stride = 0;
    double setup_t0 = GetTime();
    if (MpiSortAllocateExtractWorkspace(&ws, kSortCpeBlocks, records_per_block, raw_stride) != 0 ||
        MpiSortAllocateBlockSet(&input_pending_set, kSortCpeBlocks) != 0) {
        fprintf(stderr, "ERROR: failed to allocate MPI sort extract workspace.\n");
        MpiSortFreeBlockSet(&input_pending_set);
        MpiSortFreeExtractWorkspace(&ws);
        return -1;
    }
    stats->t_setup += GetTime() - setup_t0;

    MpiSortBlockSet *input_active = &ws.input_blocks;
    MpiSortBlockSet *input_pending = &input_pending_set;

    auto do_read_group = [&](MpiSortBlockSet *input_blocks,
                             int *n_blocks,
                             bool unhidden) -> int {
        double read_t0 = GetTime();
        int count = 0;
        for (int b = 0; b < kSortCpeBlocks; ++b) {
            bam_block *blk = &input_blocks->blocks[b];
            int ret = MpiSortMemReadBlock(reader.base, reader.size, reader.pos, blk);
            if (ret < 0 || blk->length == 28) break;
            blk->block_id = b;
            blk->pos = 0;
            count++;
        }
        double read_wall = GetTime() - read_t0;
        *n_blocks = count;
        stats->t_extract_read += read_wall;
        if (unhidden) stats->t_extract_read_unhidden += read_wall;
        return 0;
    };

    size_t local_block = 0;
    int n_blocks = 0;
    if (do_read_group(input_active, &n_blocks, true) != 0) {
        MpiSortFreeBlockSet(&input_pending_set);
        MpiSortFreeExtractWorkspace(&ws);
        return -1;
    }
    while (n_blocks > 0) {
        stats->input_blocks += n_blocks;

        double prepare_t0 = GetTime();
        size_t block_raw_offsets[kSortCpeBlocks] = {};
        size_t block_raw_sizes[kSortCpeBlocks] = {};
        size_t group_raw_begin = raw_data->size;
        size_t group_raw_total = 0;
        for (int b = 0; b < n_blocks; ++b) {
            uint32_t isize = MpiSortBgzfISize(&input_active->blocks[b]);
            if (isize > BGZF_MAX_BLOCK_SIZE) {
                fprintf(stderr,
                        "ERROR: MPI sort extract saw invalid BGZF ISIZE=%u on input block %lld.\n",
                        isize, global_block_begin + (long long)local_block + b);
                double cleanup_t0 = GetTime();
                MpiSortFreeBlockSet(&input_pending_set);
                MpiSortFreeExtractWorkspace(&ws);
                stats->t_cleanup += GetTime() - cleanup_t0;
                return -1;
            }
            block_raw_offsets[b] = group_raw_begin + group_raw_total;
            block_raw_sizes[b] = (size_t)isize;
            group_raw_total += (size_t)isize;
        }
        if (MpiSortRawReserve(raw_data, group_raw_begin + group_raw_total) != 0) {
            fprintf(stderr, "ERROR: MPI sort failed to reserve final raw buffer for extract.\n");
            double cleanup_t0 = GetTime();
            MpiSortFreeBlockSet(&input_pending_set);
            MpiSortFreeExtractWorkspace(&ws);
            stats->t_cleanup += GetTime() - cleanup_t0;
            return -1;
        }

        for (int b = 0; b < kSortCpeBlocks; ++b) {
            if (b < n_blocks) {
                ws.un_blocks.blocks[b].data = raw_data->data + block_raw_offsets[b];
                ws.un_blocks.blocks[b].length = (int)block_raw_sizes[b];
                ws.un_blocks.blocks[b].pos = 0;
                ws.un_blocks.blocks[b].errcode = 0;
            }
            paras[b].block_id = b;
            paras[b].input_block = b < n_blocks ? &input_active->blocks[b] : nullptr;
            paras[b].un_comp_block = b < n_blocks ? &ws.un_blocks.blocks[b] : nullptr;
            paras[b].raw_arena = b < n_blocks ? raw_data->data + block_raw_offsets[b] : nullptr;
            paras[b].raw_capacity = b < n_blocks ? block_raw_sizes[b] : 0;
            paras[b].raw_used = 0;
            paras[b].raw_base_offset = b < n_blocks ? (uint64_t)block_raw_offsets[b] : 0;
            paras[b].records = ws.records + (size_t)b * records_per_block;
            paras[b].record_capacity = records_per_block;
            paras[b].n_records = 0;
            paras[b].global_block_index = global_block_begin + (long long)local_block + b;
            paras[b].status = b < n_blocks ? 0 : -1;
            paras[b].record_index = 0;
            paras[b].actual_value = 0;
            paras[b].limit_value = 0;
            paras[b].limit_id = BOUNDS_LIMIT_NONE;
            paras[b].decomp_alloc_cycles = 0;
            paras[b].decomp_inflate_cycles = 0;
            paras[b].decomp_crc_cycles = 0;
            paras[b].decomp_parse_cycles = 0;
            paras[b].decomp_total_cycles = 0;
        }
        stats->t_extract_prepare += GetTime() - prepare_t0;

        // 把解压、CRC、BAM record 解析交给 CPE；同时由 MPE 读取下一批输入。
        double t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_sort_extract_raw, paras, 1);
        int next_n_blocks = 0;
#ifdef ENABLE_MASKING
        if (do_read_group(input_pending, &next_n_blocks, false) != 0) {
            athread_join();
            MpiSortFreeBlockSet(&input_pending_set);
            MpiSortFreeExtractWorkspace(&ws);
            return -1;
        }
#endif
        athread_join();
        double wall = GetTime() - t0;
#ifndef ENABLE_MASKING
        if (do_read_group(input_pending, &next_n_blocks, true) != 0) {
            MpiSortFreeBlockSet(&input_pending_set);
            MpiSortFreeExtractWorkspace(&ws);
            return -1;
        }
#endif
        stats->t_extract += wall;
        MpiSortAccumulateExtractDetail(paras, n_blocks, wall, stats);

        double status_t0 = GetTime();
        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].status != 0) {
                if (paras[b].status == -3) {
                    fprintf(stderr,
                            "ERROR: MPI sort extract capacity exceeded on input block %lld. limit_id=%d limit=%lld actual=%lld record=%d.\n",
                            global_block_begin + (long long)local_block + b,
                            paras[b].limit_id, paras[b].limit_value,
                            paras[b].actual_value, paras[b].record_index);
                } else if (paras[b].status == -4) {
                    fprintf(stderr,
                            "ERROR: MPI sort v1 assumes BAM records do not cross BGZF blocks; boundary crossed near global block %lld.\n",
                            global_block_begin + (long long)local_block + b);
                } else {
                    fprintf(stderr, "ERROR: MPI sort extract failed on input block %lld with status %d.\n",
                            global_block_begin + (long long)local_block + b, paras[b].status);
                }
                stats->t_status_check += GetTime() - status_t0;
                double cleanup_t0 = GetTime();
                MpiSortFreeBlockSet(&input_pending_set);
                MpiSortFreeExtractWorkspace(&ws);
                stats->t_cleanup += GetTime() - cleanup_t0;
                return -1;
            }
        }
        stats->t_status_check += GetTime() - status_t0;

        // 从核每处理一个 BGZF block，会把结果先放在 workspace 的临时区域里；
        // raw bytes 已经直接写入最终 local_raw，这里只追加 metadata。
        double merge_t0 = GetTime();
        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].raw_used != block_raw_sizes[b]) {
                fprintf(stderr,
                        "ERROR: MPI sort extract raw size mismatch on input block %lld. expected=%zu actual=%zu.\n",
                        global_block_begin + (long long)local_block + b,
                        block_raw_sizes[b], paras[b].raw_used);
                double cleanup_t0 = GetTime();
                MpiSortFreeBlockSet(&input_pending_set);
                MpiSortFreeExtractWorkspace(&ws);
                stats->t_cleanup += GetTime() - cleanup_t0;
                return -1;
            }
            MpiSortRecordMeta *block_records = ws.records + (size_t)b * records_per_block;
            for (int r = 0; r < paras[b].n_records; ++r) {
                records->push_back(block_records[r]);
            }
        }
        raw_data->size = group_raw_begin + group_raw_total;
        stats->t_extract_merge += GetTime() - merge_t0;
        local_block += (size_t)n_blocks;
        std::swap(input_active, input_pending);
        n_blocks = next_n_blocks;
    }

    stats->local_records += (long long)records->size();
    double cleanup_t0 = GetTime();
    MpiSortFreeBlockSet(&input_pending_set);
    MpiSortFreeExtractWorkspace(&ws);
    stats->t_cleanup += GetTime() - cleanup_t0;
    return 0;
}

void MpiSortMakeLocalSamples(const std::vector<MpiSortRecordMeta> &records,
                             int comm_size,
                             std::vector<MpiSortKey> *samples) {
    samples->clear();
    if (records.empty()) return;
    size_t target = (size_t)comm_size * kSortSampleFactor;
    size_t n = std::min(target, records.size());
    samples->reserve(n);
    for (size_t i = 0; i < n; ++i) {
        size_t idx = ((i + 1) * records.size()) / (n + 1);
        if (idx >= records.size()) idx = records.size() - 1;
        samples->push_back(MpiSortKeyFromMeta(records[idx]));
    }
}

int MpiSortChooseSplitters(const std::vector<MpiSortKey> &local_samples,
                           int rank,
                           int comm_size,
                           std::vector<MpiSortKey> *splitters,
                           long long *global_sample_count) {
    int local_count = (int)local_samples.size();
    std::vector<int> counts;
    if (rank == 0) counts.resize((size_t)comm_size, 0);
    MPI_Gather(&local_count, 1, MPI_INT,
               rank == 0 ? counts.data() : nullptr, 1, MPI_INT,
               0, MPI_COMM_WORLD);

    std::vector<int> byte_counts;
    std::vector<int> byte_displs;
    std::vector<MpiSortKey> gathered;
    int total_count = 0;
    if (rank == 0) {
        byte_counts.resize((size_t)comm_size, 0);
        byte_displs.resize((size_t)comm_size, 0);
        for (int i = 0; i < comm_size; ++i) {
            byte_displs[(size_t)i] = total_count * (int)sizeof(MpiSortKey);
            byte_counts[(size_t)i] = counts[(size_t)i] * (int)sizeof(MpiSortKey);
            total_count += counts[(size_t)i];
        }
        gathered.resize((size_t)total_count);
    }

    MPI_Gatherv((void *)(local_samples.empty() ? nullptr : local_samples.data()),
                local_count * (int)sizeof(MpiSortKey), MPI_BYTE,
                rank == 0 ? (void *)gathered.data() : nullptr,
                rank == 0 ? byte_counts.data() : nullptr,
                rank == 0 ? byte_displs.data() : nullptr,
                MPI_BYTE, 0, MPI_COMM_WORLD);

    splitters->assign((size_t)std::max(0, comm_size - 1), MpiSortKey());
    if (rank == 0) {
        std::sort(gathered.begin(), gathered.end(), MpiSortKeyLess());
        for (int s = 1; s < comm_size; ++s) {
            if (!gathered.empty()) {
                size_t idx = (size_t)((long long)gathered.size() * s / comm_size);
                if (idx >= gathered.size()) idx = gathered.size() - 1;
                (*splitters)[(size_t)s - 1] = gathered[idx];
            }
        }
        *global_sample_count = (long long)gathered.size();
    }
    MPI_Bcast(global_sample_count, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    if (comm_size > 1) {
        MPI_Bcast(splitters->data(), (int)splitters->size() * (int)sizeof(MpiSortKey),
                  MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    return 0;
}

int MpiSortFindBucket(const MpiSortRecordMeta &record,
                      const std::vector<MpiSortKey> &splitters) {
    int bucket = 0;
    while (bucket < (int)splitters.size() &&
           MpiSortCompareMetaKey(record, splitters[(size_t)bucket]) > 0) {
        bucket++;
    }
    return bucket;
}

// MPI 分布式排序里的数据重分发阶段
int MpiSortExchangeBuckets(const std::vector<MpiSortRecordMeta> &local_records,
                           const MpiSortRawBuffer &local_raw,
                           const std::vector<MpiSortKey> &splitters,
                           int rank,
                           int comm_size,
                           std::vector<MpiSortRecordMeta> *recv_records,
                           std::vector<unsigned char> *recv_raw,
                           double *bucket_count_time,
                           double *bucket_pack_time,
                           double *exchange_time,
                           double *offset_fix_time,
                           double *cleanup_time,
                           long long *bucket_self_records,
                           long long *bucket_remote_records,
                           long long *bucket_self_raw_bytes,
                           long long *bucket_remote_raw_bytes) {
    // 第一遍计算每条 record 的 bucket，以及它在目标 bucket 内的精确偏移。
    // 这些 per-record 偏移让第二遍可以由 64 个 CPE 无锁并行写 send buffers。
    double count_t0 = GetTime();
    std::vector<int> bucket_ids(local_records.size(), 0);
    std::vector<uint64_t> record_meta_in_bucket(local_records.size(), 0);
    std::vector<uint64_t> record_raw_in_bucket(local_records.size(), 0);
    std::vector<long long> send_meta_counts((size_t)comm_size, 0);
    std::vector<long long> recv_meta_counts((size_t)comm_size, 0);
    std::vector<long long> send_raw_counts((size_t)comm_size, 0);
    std::vector<long long> recv_raw_counts((size_t)comm_size, 0);
    for (size_t i = 0; i < local_records.size(); ++i) {
        const MpiSortRecordMeta &src = local_records[i];
        if ((unsigned long long)src.raw_offset + src.raw_len > (unsigned long long)local_raw.size) return -1;
        int bucket = MpiSortFindBucket(src, splitters);
        bucket_ids[i] = bucket;
        record_meta_in_bucket[i] = (uint64_t)send_meta_counts[(size_t)bucket];
        record_raw_in_bucket[i] = (uint64_t)send_raw_counts[(size_t)bucket];
        send_meta_counts[(size_t)bucket]++;
        send_raw_counts[(size_t)bucket] += src.raw_len;
        if (bucket == rank) {
            (*bucket_self_records)++;
            (*bucket_self_raw_bytes) += src.raw_len;
        } else {
            (*bucket_remote_records)++;
            (*bucket_remote_raw_bytes) += src.raw_len;
        }
    }
    *bucket_count_time += GetTime() - count_t0;

    std::vector<long long> send_meta_displs((size_t)comm_size, 0);
    std::vector<long long> send_raw_displs((size_t)comm_size, 0);
    long long total_send_meta = 0;
    long long total_send_raw = 0;
    for (int i = 0; i < comm_size; ++i) {
        send_meta_displs[(size_t)i] = total_send_meta;
        send_raw_displs[(size_t)i] = total_send_raw;
        total_send_meta += send_meta_counts[(size_t)i];
        total_send_raw += send_raw_counts[(size_t)i];
    }
    if (total_send_meta < 0 || total_send_raw < 0 ||
        (unsigned long long)total_send_meta > (unsigned long long)(SIZE_MAX / sizeof(MpiSortRecordMeta)) ||
        (unsigned long long)total_send_raw > (unsigned long long)SIZE_MAX) {
        return -1;
    }

    // 第二遍写入一次性分配的连续发送缓冲。CPE 做变长 raw bytes 拷贝；
    // MPE 同时交换 counts 并分配接收缓冲，用 CPE pack 掩盖这部分主核工作。
    double pack_t0 = GetTime();
    std::vector<MpiSortRecordMeta> send_meta((size_t)total_send_meta);
    std::vector<unsigned char> send_raw((size_t)total_send_raw);
    MpiSortBucketPackPara pack_paras[kSortCpeBlocks];
    size_t n_records = local_records.size();
    size_t base_records = n_records / kSortCpeBlocks;
    size_t extra_records = n_records % kSortCpeBlocks;
    size_t begin = 0;
    for (int c = 0; c < kSortCpeBlocks; ++c) {
        size_t span = base_records + ((size_t)c < extra_records ? 1 : 0);
        MpiSortBucketPackPara &para = pack_paras[c];
        para.core_id = c;
        para.local_records = local_records.empty() ? nullptr : local_records.data();
        para.local_raw = local_raw.data;
        para.local_raw_size = local_raw.size;
        para.bucket_ids = bucket_ids.empty() ? nullptr : bucket_ids.data();
        para.record_meta_in_bucket = record_meta_in_bucket.empty() ? nullptr : record_meta_in_bucket.data();
        para.record_raw_in_bucket = record_raw_in_bucket.empty() ? nullptr : record_raw_in_bucket.data();
        para.send_meta_displs = send_meta_displs.data();
        para.send_raw_displs = send_raw_displs.data();
        para.send_meta = send_meta.empty() ? nullptr : send_meta.data();
        para.send_raw = send_raw.empty() ? nullptr : send_raw.data();
        para.send_meta_capacity = (uint64_t)send_meta.size();
        para.send_raw_capacity = (uint64_t)send_raw.size();
        para.record_begin = begin;
        para.record_end = begin + span;
        para.status = 0;
        para.record_index = -1;
        para.actual_value = 0;
        para.limit_value = 0;
        para.limit_id = BOUNDS_LIMIT_NONE;
        para.pack_cycles = 0;
        para.total_cycles = 0;
        begin += span;
    }
    __real_athread_spawn((void *)slave_mpi_sort_bucket_pack, pack_paras, 1);

    // counts 交换和 recv buffer 分配不依赖 send_meta/send_raw 内容，可与 CPE pack 并行。
    MPI_Alltoall(send_meta_counts.data(), 1, MPI_LONG_LONG,
                 recv_meta_counts.data(), 1, MPI_LONG_LONG, MPI_COMM_WORLD);
    MPI_Alltoall(send_raw_counts.data(), 1, MPI_LONG_LONG,
                 recv_raw_counts.data(), 1, MPI_LONG_LONG, MPI_COMM_WORLD);

    std::vector<long long> meta_displs((size_t)comm_size, 0);
    std::vector<long long> raw_displs((size_t)comm_size, 0);
    long long total_meta = 0;
    long long total_raw = 0;
    // 计算接收偏移并分配接收数组
    for (int i = 0; i < comm_size; ++i) {
        meta_displs[(size_t)i] = total_meta;
        raw_displs[(size_t)i] = total_raw;
        total_meta += recv_meta_counts[(size_t)i];
        total_raw += recv_raw_counts[(size_t)i];
    }
    int recv_layout_ok = 1;
    if (total_meta < 0 || total_raw < 0 ||
        (unsigned long long)total_meta > (unsigned long long)(SIZE_MAX / sizeof(MpiSortRecordMeta)) ||
        (unsigned long long)total_raw > (unsigned long long)SIZE_MAX) {
        recv_layout_ok = 0;
    }
    if (recv_layout_ok) {
        // 分配接收空间
        recv_records->assign((size_t)total_meta, MpiSortRecordMeta());
        recv_raw->assign((size_t)total_raw, 0);
    }
    athread_join();
    *bucket_pack_time += GetTime() - pack_t0;
    if (!recv_layout_ok) return -1;

    for (int c = 0; c < kSortCpeBlocks; ++c) {
        if (pack_paras[c].status != 0) {
            fprintf(stderr,
                    "ERROR: MPI sort bucket pack failed on CPE %d status=%d record=%d actual=%lld limit=%lld.\n",
                    c, pack_paras[c].status, pack_paras[c].record_index,
                    pack_paras[c].actual_value, pack_paras[c].limit_value);
            return -1;
        }
    }

    // MPI rank 间交换 metadata 和 raw bytes。counts 交换已在 CPE pack 期间完成，
    // 这里主要记录 self copy + 环形 Sendrecv 的非掩盖时间。
    double exch_t0 = GetTime();

    if (recv_meta_counts[(size_t)rank] > 0) {
        memcpy(recv_records->data() + meta_displs[(size_t)rank],
               send_meta.data() + send_meta_displs[(size_t)rank],
               (size_t)recv_meta_counts[(size_t)rank] * sizeof(MpiSortRecordMeta));
    }
    if (recv_raw_counts[(size_t)rank] > 0) {
        memcpy(recv_raw->data() + raw_displs[(size_t)rank],
               send_raw.data() + send_raw_displs[(size_t)rank],
               (size_t)recv_raw_counts[(size_t)rank]);
    }

    // 环形 Sendrecv 交换其他 rank 的数据
    for (int step = 1; step < comm_size; ++step) {
        int send_to = (rank + step) % comm_size;
        int recv_from = (rank - step + comm_size) % comm_size;
        long long send_meta_bytes = send_meta_counts[(size_t)send_to] * (long long)sizeof(MpiSortRecordMeta);
        long long recv_meta_bytes = recv_meta_counts[(size_t)recv_from] * (long long)sizeof(MpiSortRecordMeta);
        const char *send_meta_ptr = send_meta_bytes > 0
                                    ? (const char *)(send_meta.data() + send_meta_displs[(size_t)send_to])
                                    : nullptr;
        char *recv_meta_ptr = recv_meta_bytes > 0
                              ? (char *)(recv_records->data() + meta_displs[(size_t)recv_from])
                              : nullptr;
        if (MpiSortSendRecvBytes(send_to, send_meta_ptr, send_meta_bytes,
                                 recv_from, recv_meta_ptr, recv_meta_bytes,
                                 kSortTagMeta) != 0) {
            return -1;
        }

        long long send_raw_bytes = send_raw_counts[(size_t)send_to];
        long long recv_raw_bytes = recv_raw_counts[(size_t)recv_from];
        const char *send_raw_ptr = send_raw_bytes > 0
                                   ? (const char *)(send_raw.data() + send_raw_displs[(size_t)send_to])
                                   : nullptr;
        char *recv_raw_ptr = recv_raw_bytes > 0
                             ? (char *)(recv_raw->data() + raw_displs[(size_t)recv_from])
                             : nullptr;
        if (MpiSortSendRecvBytes(send_to, send_raw_ptr, send_raw_bytes,
                                 recv_from, recv_raw_ptr, recv_raw_bytes,
                                 kSortTagRaw) != 0) {
            return -1;
        }
    }
    *exchange_time += GetTime() - exch_t0;

    double offset_t0 = GetTime();
    for (int src = 0; src < comm_size; ++src) {
        long long begin = meta_displs[(size_t)src];
        long long end = begin + recv_meta_counts[(size_t)src];
        uint64_t raw_base = (uint64_t)raw_displs[(size_t)src];
        for (long long i = begin; i < end; ++i) {
            (*recv_records)[(size_t)i].raw_offset += raw_base;
        }
    }
    *offset_fix_time += GetTime() - offset_t0;

    // 显式释放大容量发送缓冲，使析构/释放成本进入日志，而不是落在未计时区间。
    double cleanup_t0 = GetTime();
    std::vector<int>().swap(bucket_ids);
    std::vector<MpiSortRecordMeta>().swap(send_meta);
    std::vector<unsigned char>().swap(send_raw);
    *cleanup_time += GetTime() - cleanup_t0;
    return 0;
}

void MpiSortInitEmptyRawCompressPara(MpiSortRawCompressPara *para, int block_id) {
    para->block_id = block_id;
    para->un_comp_block = nullptr;
    para->un_comp_size = 0;
    para->output_block = nullptr;
    para->output_size = 0;
    para->status = -1;
    para->compress_level = 1;
    para->compress_pack_cycles = 0;
    para->compress_alloc_cycles = 0;
    para->compress_deflate_cycles = 0;
    para->compress_footer_cycles = 0;
    para->compress_total_cycles = 0;
}

int MpiSortPreparePayloadBatch(const std::vector<MpiSortRecordMeta> &records,
                               const std::vector<unsigned char> &raw,
                               size_t *record_pos,
                               MpiSortBlockSet *payload_blocks,
                               MpiSortBlockSet *output_blocks,
                               MpiSortRawCompressPara *paras,
                               int compress_level,
                               int *active_blocks) {
    int active = 0;
    while (*record_pos < records.size() && active < kSortCpeBlocks) {
        bam_block *payload = &payload_blocks->blocks[active];
        payload->pos = 0;
        payload->length = 0;
        payload->errcode = 0;
        payload->block_id = active;
        while (*record_pos < records.size()) {
            const MpiSortRecordMeta &rec = records[*record_pos];
            if (rec.raw_len > BGZF_BLOCK_SIZE) return -1;
            if ((unsigned long long)rec.raw_offset + rec.raw_len > (unsigned long long)raw.size()) return -1;
            if (payload->pos > 0 && payload->pos + rec.raw_len > BGZF_BLOCK_SIZE) break;
            memcpy(payload->data + payload->pos, raw.data() + rec.raw_offset, rec.raw_len);
            payload->pos += rec.raw_len;
            payload->length = payload->pos;
            (*record_pos)++;
        }
        if (payload->pos == 0) return -1;

        paras[active].block_id = active;
        paras[active].un_comp_block = payload;
        paras[active].un_comp_size = (int)payload->pos;
        paras[active].output_block = &output_blocks->blocks[active];
        paras[active].output_size = 0;
        paras[active].status = 0;
        paras[active].compress_level = compress_level;
        paras[active].compress_pack_cycles = 0;
        paras[active].compress_alloc_cycles = 0;
        paras[active].compress_deflate_cycles = 0;
        paras[active].compress_footer_cycles = 0;
        paras[active].compress_total_cycles = 0;
        active++;
    }
    for (int i = active; i < kSortCpeBlocks; ++i) {
        MpiSortInitEmptyRawCompressPara(&paras[i], i);
    }
    *active_blocks = active;
    return 0;
}

int MpiSortCompressSortedRecordsCpe(const std::vector<MpiSortRecordMeta> &records,
                                    const std::vector<unsigned char> &raw,
                                    MemWriter &mem_writer,
                                    int compress_level,
                                    MpiSortStats *stats) {
    double setup_t0 = GetTime();
    MpiSortBlockSet payload_a = {}, payload_b = {}, out_a = {}, out_b = {};
    MpiSortRawCompressPara paras_a[kSortCpeBlocks], paras_b[kSortCpeBlocks];
    if (MpiSortAllocateBlockSet(&payload_a, kSortCpeBlocks) != 0 ||
        MpiSortAllocateBlockSet(&payload_b, kSortCpeBlocks) != 0 ||
        MpiSortAllocateBlockSet(&out_a, kSortCpeBlocks) != 0 ||
        MpiSortAllocateBlockSet(&out_b, kSortCpeBlocks) != 0) {
        fprintf(stderr, "ERROR: failed to allocate MPI sort compress workspace.\n");
        MpiSortFreeBlockSet(&payload_a);
        MpiSortFreeBlockSet(&payload_b);
        MpiSortFreeBlockSet(&out_a);
        MpiSortFreeBlockSet(&out_b);
        return -1;
    }
    for (int i = 0; i < kSortCpeBlocks; ++i) {
        MpiSortInitEmptyRawCompressPara(&paras_a[i], i);
        MpiSortInitEmptyRawCompressPara(&paras_b[i], i);
    }

    MpiSortBlockSet *payload_active = &payload_a;
    MpiSortBlockSet *payload_pending = &payload_b;
    MpiSortBlockSet *out_active = &out_a;
    MpiSortBlockSet *out_pending = &out_b;
    MpiSortRawCompressPara *paras_active = paras_a;
    MpiSortRawCompressPara *paras_pending = paras_b;
    bool has_pending = false;

    auto flush_pending = [&]() -> int {
        if (!has_pending) return 0;
        double t0 = GetTime();
        for (int i = 0; i < kSortCpeBlocks; ++i) {
            if (paras_pending[i].status == 0 && paras_pending[i].output_block) {
                if (MpiWriteBlockToMem(mem_writer, paras_pending[i].output_block) != 0) return -1;
                stats->bgzf_blocks++;
            }
            MpiSortInitEmptyRawCompressPara(&paras_pending[i], i);
        }
        has_pending = false;
        stats->t_write += GetTime() - t0;
        return 0;
    };

    size_t record_pos = 0;
    int active_blocks = 0;
    if (MpiSortPreparePayloadBatch(records, raw, &record_pos, payload_active, out_active,
                                   paras_active, compress_level, &active_blocks) != 0) {
        MpiSortFreeBlockSet(&payload_a);
        MpiSortFreeBlockSet(&payload_b);
        MpiSortFreeBlockSet(&out_a);
        MpiSortFreeBlockSet(&out_b);
        return -1;
    }
    stats->t_setup += GetTime() - setup_t0;

    int ret_code = 0;
    while (active_blocks > 0) {
        int next_active_blocks = 0;
        double t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_sort_compress_payload, paras_active, 1);
        if (flush_pending() != 0) {
            athread_join();
            ret_code = -1;
            break;
        }
        if (MpiSortPreparePayloadBatch(records, raw, &record_pos, payload_pending, out_pending,
                                       paras_pending, compress_level, &next_active_blocks) != 0) {
            athread_join();
            ret_code = -1;
            break;
        }
        athread_join();
        double wall = GetTime() - t0;
        stats->t_compress += wall;
        MpiSortAccumulateCompressDetail(paras_active, active_blocks, wall, stats);

        double status_t0 = GetTime();
        for (int i = 0; i < active_blocks; ++i) {
            if (paras_active[i].status != 0) {
                fprintf(stderr, "ERROR: MPI sort raw compress failed on block %d with status %d.\n",
                        i, paras_active[i].status);
                ret_code = -1;
                break;
            }
        }
        stats->t_status_check += GetTime() - status_t0;
        if (ret_code != 0) break;

        has_pending = active_blocks > 0;
        std::swap(payload_active, payload_pending);
        std::swap(out_active, out_pending);
        std::swap(paras_active, paras_pending);
        active_blocks = next_active_blocks;
    }

    int ret = ret_code == 0 ? flush_pending() : ret_code;
    double cleanup_t0 = GetTime();
    MpiSortFreeBlockSet(&payload_a);
    MpiSortFreeBlockSet(&payload_b);
    MpiSortFreeBlockSet(&out_a);
    MpiSortFreeBlockSet(&out_b);
    stats->t_cleanup += GetTime() - cleanup_t0;
    return ret;
}

} // namespace


int FusedBamSortMPI(MemReader &reader,
                    MemWriter &mem_writer,
                    long long global_block_begin,
                    int rank,
                    int comm_size,
                    int compress_level,
                    size_t memory_limit,
                    MpiSortStats *stats) {
    double fused_t0 = GetTime();
    double setup_t0 = GetTime();
    //1. 参数准备
    std::vector<MpiSortRecordMeta> local_records;
    MpiSortRawBuffer local_raw;
    std::vector<MpiSortKey> local_samples;
    std::vector<MpiSortKey> splitters;
    std::vector<MpiSortRecordMeta> received_records;
    std::vector<unsigned char> received_raw;

    if (reader.size > 0) {
        if (MpiSortRawReserve(&local_raw, reader.size * 2) != 0) {
            fprintf(stderr, "ERROR: MPI sort failed to reserve local raw buffer.\n");
            return -1;
        }
        local_records.reserve(reader.size / 64);
    }
    stats->t_setup += GetTime() - setup_t0;

    //2. 把本 rank 负责的 BGZF blocks 解析成：local_records 和 local_raw
    // local_records：record metadata，包括 tid/pos/flag/global_order/raw_offset/raw_len
    // local_raw：真实 BAM record bytes
    if (MpiSortExtractLocalRecordsCpe(reader, global_block_begin, &local_records, &local_raw, stats) != 0) {
        return -1;
    }
    //检查预估的大小是否超过输出的限制
    if (MpiSortCheckMemoryLimit(memory_limit,
                                local_raw.size + local_records.size() * sizeof(MpiSortRecordMeta),
                                rank, "extract") != 0) {
        return -1;
    }

    //3. 对本 rank 的 records 做本地排序。
    // 注意这里只排序 metadata，不移动真实 BAM record bytes，所以比移动完整 record 更轻。
    double t0 = GetTime();
    std::sort(local_records.begin(), local_records.end(), MpiSortMetaLess());
    stats->t_local_sort += GetTime() - t0;

    //4. 抽样划分数据
    t0 = GetTime();
    //每个 rank 从本地有序 records 中均匀抽样
    MpiSortMakeLocalSamples(local_records, comm_size, &local_samples);
    long long global_sample_count = 0;
    //rank 0 汇总所有 samples，排序后选出 comm_size - 1 个 splitter，再广播给所有 rank。这样全局排序空间被切成 comm_size 个区间。
    if (MpiSortChooseSplitters(local_samples, rank, comm_size, &splitters, &global_sample_count) != 0) {
        return -1;
    }
    stats->sample_records = (long long)local_samples.size();
    stats->t_sample += GetTime() - t0;


    //5. 按 splitter 分桶并交换 records，这一步完成后，当前 rank 拿到的是自己负责的全局坐标区间。
    double bucket_count_time = 0.0;
    double bucket_pack_time = 0.0;
    double exchange_time = 0.0;
    double offset_fix_time = 0.0;
    double exchange_cleanup_time = 0.0;
    long long bucket_self_records = 0;
    long long bucket_remote_records = 0;
    long long bucket_self_raw_bytes = 0;
    long long bucket_remote_raw_bytes = 0;
    if (MpiSortExchangeBuckets(local_records, local_raw, splitters, rank, comm_size,
                               &received_records, &received_raw,
                               &bucket_count_time, &bucket_pack_time,
                               &exchange_time, &offset_fix_time,
                               &exchange_cleanup_time,
                               &bucket_self_records, &bucket_remote_records,
                               &bucket_self_raw_bytes, &bucket_remote_raw_bytes) != 0) {
        return -1;
    }
    stats->bucket_self_records += bucket_self_records;
    stats->bucket_remote_records += bucket_remote_records;
    stats->bucket_self_raw_bytes += bucket_self_raw_bytes;
    stats->bucket_remote_raw_bytes += bucket_remote_raw_bytes;
    stats->t_bucket_count += bucket_count_time;
    stats->t_bucket_pack += bucket_pack_time;
    stats->t_partition += bucket_count_time + bucket_pack_time;
    stats->t_mpi_exchange += exchange_time;
    stats->t_exchange += exchange_time;
    stats->t_offset_fix += offset_fix_time;
    stats->t_cleanup += exchange_cleanup_time;
    stats->received_records = (long long)received_records.size();
    if (MpiSortCheckMemoryLimit(memory_limit,
                                received_raw.size() + received_records.size() * sizeof(MpiSortRecordMeta),
                                rank, "exchange") != 0) {
        return -1;
    }

    double cleanup_t0 = GetTime();
    std::vector<MpiSortRecordMeta>().swap(local_records);
    MpiSortRawRelease(&local_raw);
    stats->t_cleanup += GetTime() - cleanup_t0;

    //6. 对收到的全局数据再排序
    t0 = GetTime();
    std::sort(received_records.begin(), received_records.end(), MpiSortMetaLess());
    stats->t_final_sort += GetTime() - t0;

    //7. 把排序后的 records 重新压缩成 BGZF body
    if (MpiSortCompressSortedRecordsCpe(received_records, received_raw,
                                        mem_writer, compress_level, stats) != 0) {
        return -1;
    }
    stats->t_fused_total += GetTime() - fused_t0;
    (void)global_sample_count;
    return 0;
}
