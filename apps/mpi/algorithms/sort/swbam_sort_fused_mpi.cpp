#include "swbam_mpi.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <stdint.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

#include <mpi.h>
#include <libdeflate.h>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

extern "C" {
    void slave_mpi_sort_extract_raw();
    void slave_mpi_sort_compress_payload();
    void slave_mpi_sort_bucket_pack();
    void slave_mpi_sort_range_pack();
}

namespace {

const int kSortSampleFactor = 64;
const int kSortTagMeta = 7101;
const int kSortTagRaw = 7102;
const int kSortExchangeChunk = 64 * 1024 * 1024;
const int kSortCpeBlocks = 64;

#ifndef RABBITBAM_SORT_ENABLE_SANITY_CHECKS
#define RABBITBAM_SORT_ENABLE_SANITY_CHECKS 0
#endif

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

#if RABBITBAM_SORT_ENABLE_SANITY_CHECKS
int MpiSortValidatePayloadRecords(const bam_block *payload,
                                  const char *stage) {
    if (!payload || !payload->data) return -1;
    const size_t payload_len = (size_t)payload->pos;
    size_t pos = 0;
    int record = 0;
    while (pos < payload_len) {
        int rank = -1;
        if (payload_len - pos < 4) {
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            fprintf(stderr,
                    "[rank %d] ERROR: sort payload sanity failed "
                    "stage=%s block=%d record=%d offset=%zu "
                    "remaining=%zu reason=truncated_block_len\n",
                    rank, stage, payload->block_id, record, pos,
                    payload_len - pos);
            return -1;
        }
        const uint32_t block_len =
            MpiSortReadLe32((const unsigned char *)payload->data + pos);
        const uint64_t record_end = (uint64_t)pos + 4u + block_len;
        if (block_len < 32 || record_end > payload_len) {
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            fprintf(stderr,
                    "[rank %d] ERROR: sort payload sanity failed "
                    "stage=%s block=%d record=%d offset=%zu "
                    "payload_len=%zu block_len=%u record_end=%llu "
                    "reason=bad_record_extent\n",
                    rank, stage, payload->block_id, record, pos,
                    payload_len, block_len,
                    (unsigned long long)record_end);
            return -1;
        }
        const unsigned char *body =
            (const unsigned char *)payload->data + pos + 4;
        const uint32_t x2 = MpiSortReadLe32(body + 8);
        const uint32_t x3 = MpiSortReadLe32(body + 12);
        const uint32_t raw_l_qname = x2 & 0xffu;
        const uint32_t n_cigar = x3 & 0xffffu;
        const uint32_t l_qseq = MpiSortReadLe32(body + 16);
        const uint64_t data_len = (uint64_t)block_len - 32u;
        const uint64_t minimum_data =
            (uint64_t)raw_l_qname + ((uint64_t)n_cigar << 2) +
            (((uint64_t)l_qseq + 1u) >> 1) + (uint64_t)l_qseq;
        if (raw_l_qname == 0 || minimum_data > data_len) {
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            fprintf(stderr,
                    "[rank %d] ERROR: sort payload sanity failed "
                    "stage=%s block=%d record=%d offset=%zu "
                    "payload_len=%zu block_len=%u raw_l_qname=%u "
                    "n_cigar=%u l_qseq=%u min_data=%llu data_len=%llu "
                    "reason=bad_record_fields\n",
                    rank, stage, payload->block_id, record, pos,
                    payload_len, block_len, raw_l_qname, n_cigar,
                    l_qseq, (unsigned long long)minimum_data,
                    (unsigned long long)data_len);
            return -1;
        }
        pos = (size_t)record_end;
        ++record;
    }
    return 0;
}

int MpiSortValidateMetaRawBuffer(
        const std::vector<MpiSortRecordMeta> &records,
        const unsigned char *raw_data,
        size_t raw_size,
        const std::vector<long long> *source_counts,
        const char *stage,
        bool verbose = true) {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    size_t source = 0;
    long long source_end = source_counts && !source_counts->empty()
        ? (*source_counts)[0] : (long long)records.size();
    for (size_t i = 0; i < records.size(); ++i) {
        while (source_counts && source + 1 < source_counts->size() &&
               (long long)i >= source_end) {
            ++source;
            source_end += (*source_counts)[source];
        }
        const MpiSortRecordMeta &rec = records[i];
        const uint64_t end = rec.raw_offset + (uint64_t)rec.raw_len;
        uint32_t block_len = 0;
        if (rec.raw_len >= 4 && raw_data && end <= raw_size) {
            block_len = MpiSortReadLe32(raw_data + rec.raw_offset);
        }
        if (rec.raw_len < 36 || end > raw_size ||
            (uint64_t)block_len + 4u != rec.raw_len) {
            if (verbose) {
                fprintf(stderr,
                        "[rank %d] ERROR: sort meta/raw sanity failed "
                        "stage=%s source=%zu meta_index=%zu raw_offset=%llu "
                        "raw_len=%u raw_size=%zu block_len=%u expected_raw_len=%llu "
                        "tid=%d pos=%d flag=%u order=%llu\n",
                        rank, stage, source, i,
                        (unsigned long long)rec.raw_offset, rec.raw_len,
                        raw_size, block_len,
                        (unsigned long long)block_len + 4u,
                        rec.tid, rec.pos, rec.flag,
                        (unsigned long long)rec.global_order);
            }
            return -1;
        }
    }
    return 0;
}

int MpiSortValidateMetaRawRecords(
        const std::vector<MpiSortRecordMeta> &records,
        const std::vector<unsigned char> &raw,
        const std::vector<long long> *source_counts,
        const char *stage) {
    return MpiSortValidateMetaRawBuffer(records, raw.empty() ? nullptr : raw.data(),
                                        raw.size(), source_counts, stage);
}
#endif

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
            "[rank %d] ERROR: MPI sort memory mode exceeded -m during %s. estimate=%zu limit=%zu; use a smaller -m to force external mode or increase -m.\n",
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
    stats->cpe_calibration_cycles += (long long)critical_total;
    stats->t_cpe_calibration_wall += wall;
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

uint64_t MpiSortCompressCriticalCycles(const MpiSortRawCompressPara *paras,
                                       int active_blocks) {
    uint64_t critical_total = 0;
    for (int b = 0; b < active_blocks; ++b) {
        critical_total =
            std::max(critical_total, paras[b].compress_total_cycles);
    }
    return critical_total;
}

double MpiSortEstimateCompressSeconds(const MpiSortRawCompressPara *paras,
                                      int active_blocks,
                                      double fallback_wall,
                                      const MpiSortStats *stats) {
    const uint64_t critical_cycles =
        MpiSortCompressCriticalCycles(paras, active_blocks);
    if (critical_cycles > 0 && stats->cpe_calibration_cycles > 0 &&
        stats->t_cpe_calibration_wall > 0.0) {
        return (double)critical_cycles *
               stats->t_cpe_calibration_wall /
               (double)stats->cpe_calibration_cycles;
    }
    return fallback_wall;
}

class MpiSortBatchReader {
public:
    virtual ~MpiSortBatchReader() {}
    virtual int Read(MpiSortBlockSet *blocks, int *count) = 0;
};

class MpiSortMemoryBatchReader : public MpiSortBatchReader {
public:
    explicit MpiSortMemoryBatchReader(MemReader *reader)
        : reader_(reader) {}

    int Read(MpiSortBlockSet *blocks, int *count) {
        if (!reader_ || !blocks || !count) return -1;
        if (reader_->pos >= reader_->size) {
            *count = 0;
            return 0;
        }
        int n = 0;
        for (int b = 0; b < kSortCpeBlocks; ++b) {
            bam_block *block = &blocks->blocks[b];
            const int ret = MpiSortMemReadBlock(
                reader_->base, reader_->size, reader_->pos, block);
            if (ret < 0) return -1;
            if (block->length == 28) break;
            block->block_id = b;
            block->pos = 0;
            ++n;
            if (reader_->pos >= reader_->size) break;
        }
        *count = n;
        return 0;
    }

private:
    MemReader *reader_;
};

class MpiSortBackendBatchReader : public MpiSortBatchReader {
public:
    MpiSortBackendBatchReader(const swbam::BamInputBackend *input,
                              const swbam::BgzfBlockSpan *spans,
                              size_t count)
        : reader_(input, spans, count) {}

    int Read(MpiSortBlockSet *blocks, int *count) {
        if (!blocks || !count) return -1;
        swbam::BgzfBlockBatch batch;
        if (batch.Attach(blocks->blocks, kSortCpeBlocks) != 0) return -1;
        size_t n = 0;
        if (reader_.ReadNext(&batch, &n) != 0 || n > kSortCpeBlocks) {
            return -1;
        }
        *count = (int)n;
        return 0;
    }

private:
    swbam::BgzfSpanBatchReader reader_;
};

// CPE 批量提取 BAM records
int MpiSortExtractLocalRecordsCpeImpl(
        MpiSortBatchReader *block_reader,
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
        const int read_ret = block_reader->Read(input_blocks, &count);
        double read_wall = GetTime() - read_t0;
        *n_blocks = count;
        stats->t_extract_read += read_wall;
        if (unhidden) stats->t_extract_read_unhidden += read_wall;
        return read_ret;
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
                ws.un_blocks.blocks[b].data =
                    ws.un_blocks.data + (size_t)b * BGZF_MAX_BLOCK_SIZE;
                ws.un_blocks.blocks[b].length = (int)block_raw_sizes[b];
                ws.un_blocks.blocks[b].pos = 0;
                ws.un_blocks.blocks[b].errcode = 0;
            }
            paras[b].block_id = b;
            paras[b].input_block = b < n_blocks ? &input_active->blocks[b] : nullptr;
            paras[b].un_comp_block = b < n_blocks ? &ws.un_blocks.blocks[b] : nullptr;
            paras[b].raw_arena = b < n_blocks ? ws.un_blocks.blocks[b].data : nullptr;
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
            if (block_raw_sizes[b] > 0) {
                memcpy(raw_data->data + block_raw_offsets[b],
                       ws.un_blocks.blocks[b].data,
                       block_raw_sizes[b]);
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

int MpiSortExtractLocalRecordsCpe(
        MemReader &reader, long long global_block_begin,
        std::vector<MpiSortRecordMeta> *records,
        MpiSortRawBuffer *raw_data, MpiSortStats *stats) {
    MpiSortMemoryBatchReader block_reader(&reader);
    return MpiSortExtractLocalRecordsCpeImpl(
        &block_reader, global_block_begin, records, raw_data, stats);
}

int MpiSortExtractLocalRecordsCpe(
        const swbam::BamInputBackend &input,
        const swbam::BgzfBlockSpan *spans, size_t span_count,
        long long global_block_begin,
        std::vector<MpiSortRecordMeta> *records,
        MpiSortRawBuffer *raw_data, MpiSortStats *stats) {
    MpiSortBackendBatchReader block_reader(&input, spans, span_count);
    return MpiSortExtractLocalRecordsCpeImpl(
        &block_reader, global_block_begin, records, raw_data, stats);
}

#if RABBITBAM_SORT_ENABLE_SANITY_CHECKS
int MpiSortHostDecompressBlock(const bam_block *comp,
                               unsigned char *dst,
                               size_t capacity,
                               struct libdeflate_decompressor *decompressor,
                               size_t *out_len) {
    if (!comp || !dst || !decompressor || !out_len ||
        comp->length < BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH) {
        return -1;
    }
    const unsigned char *src = (const unsigned char *)comp->data;
    const uint32_t expected_crc =
        MpiSortReadLe32(src + comp->length - 8);
    size_t actual_len = capacity;
    int ret = libdeflate_deflate_decompress(
        decompressor,
        src + BLOCK_HEADER_LENGTH,
        (size_t)comp->length - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH,
        dst, capacity, &actual_len);
    if (ret != 0) return -1;
    const uint32_t crc = libdeflate_crc32(0, dst, actual_len);
    if (crc != expected_crc) return -1;
    *out_len = actual_len;
    return 0;
}

int MpiSortExtractLocalRecordsHost(MemReader &reader,
                                   long long global_block_begin,
                                   std::vector<MpiSortRecordMeta> *records,
                                   MpiSortRawBuffer *raw_data,
                                   MpiSortStats *stats) {
    MpiSortBlockSet input = {};
    if (MpiSortAllocateBlockSet(&input, 1) != 0) {
        fprintf(stderr, "ERROR: failed to allocate MPI sort host extract block.\n");
        return -1;
    }
    struct libdeflate_decompressor *decompressor =
        libdeflate_alloc_decompressor();
    if (!decompressor) {
        MpiSortFreeBlockSet(&input);
        fprintf(stderr, "ERROR: failed to allocate MPI sort host decompressor.\n");
        return -1;
    }

    int ret = 0;
    size_t local_block = 0;
    while (true) {
        bam_block *blk = &input.blocks[0];
        int read_ret = MpiSortMemReadBlock(reader.base, reader.size,
                                           reader.pos, blk);
        if (read_ret < 0 || blk->length == 28) break;
        const uint32_t isize = MpiSortBgzfISize(blk);
        if (isize > BGZF_MAX_BLOCK_SIZE) {
            fprintf(stderr,
                    "ERROR: MPI sort host extract saw invalid BGZF ISIZE=%u on input block %lld.\n",
                    isize, global_block_begin + (long long)local_block);
            ret = -1;
            break;
        }
        const size_t raw_begin = raw_data->size;
        if (MpiSortRawReserve(raw_data, raw_begin + (size_t)isize) != 0) {
            fprintf(stderr,
                    "ERROR: MPI sort host extract failed to reserve raw buffer.\n");
            ret = -1;
            break;
        }

        double t0 = GetTime();
        size_t out_len = (size_t)isize;
        if (MpiSortHostDecompressBlock(
                blk, raw_data->data + raw_begin, (size_t)isize,
                decompressor, &out_len) != 0 ||
            out_len != (size_t)isize) {
            fprintf(stderr,
                    "ERROR: MPI sort host extract failed to decompress input block %lld. expected=%u actual=%zu.\n",
                    global_block_begin + (long long)local_block,
                    isize, out_len);
            ret = -1;
            break;
        }
        if (stats) stats->t_extract += GetTime() - t0;

        size_t pos = 0;
        int count = 0;
        while (pos < out_len) {
            if (out_len - pos < 4) {
                fprintf(stderr,
                        "ERROR: MPI sort host extract saw truncated block_len at input block %lld offset=%zu.\n",
                        global_block_begin + (long long)local_block, pos);
                ret = -1;
                break;
            }
            const unsigned char *record = raw_data->data + raw_begin + pos;
            const uint32_t block_len = MpiSortReadLe32(record);
            const size_t raw_len = (size_t)block_len + 4u;
            if (block_len < 32 || raw_len > out_len - pos) {
                fprintf(stderr,
                        "ERROR: MPI sort host extract saw invalid record extent at input block %lld record=%d offset=%zu block_len=%u remaining=%zu.\n",
                        global_block_begin + (long long)local_block,
                        count, pos, block_len, out_len - pos);
                ret = -1;
                break;
            }
            const unsigned char *body = record + 4;
            const uint32_t x2 = MpiSortReadLe32(body + 8);
            const uint32_t x3 = MpiSortReadLe32(body + 12);
            const uint32_t raw_l_qname = x2 & 0xffu;
            const uint32_t n_cigar = x3 & 0xffffu;
            const uint32_t l_qseq = MpiSortReadLe32(body + 16);
            const uint64_t data_len = (uint64_t)block_len - 32u;
            const uint64_t minimum_data =
                (uint64_t)raw_l_qname + ((uint64_t)n_cigar << 2) +
                (((uint64_t)l_qseq + 1u) >> 1) + (uint64_t)l_qseq;
            if (raw_l_qname == 0 || minimum_data > data_len) {
                fprintf(stderr,
                        "ERROR: MPI sort host extract saw invalid record fields at input block %lld record=%d offset=%zu block_len=%u l_qname=%u n_cigar=%u l_qseq=%u.\n",
                        global_block_begin + (long long)local_block,
                        count, pos, block_len, raw_l_qname,
                        n_cigar, l_qseq);
                ret = -1;
                break;
            }

            MpiSortRecordMeta meta;
            meta.tid = (int32_t)MpiSortReadLe32(body);
            meta.pos = (int32_t)MpiSortReadLe32(body + 4);
            meta.flag = MpiSortReadLe16(body + 14);
            meta.pad = 0;
            meta.raw_len = (uint32_t)raw_len;
            meta.pad2 = 0;
            meta.raw_offset = (uint64_t)raw_begin + (uint64_t)pos;
            meta.global_order =
                ((uint64_t)(global_block_begin + (long long)local_block) << 32) |
                (uint32_t)count;
            records->push_back(meta);
            pos += raw_len;
            ++count;
        }
        if (ret != 0) break;
        raw_data->size = raw_begin + out_len;
        ++local_block;
    }

    if (ret == 0 && stats) {
        stats->local_records += (long long)records->size();
    }
    libdeflate_free_decompressor(decompressor);
    MpiSortFreeBlockSet(&input);
    return ret;
}
#endif

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

// MPI 分布式排序里的数据重分发阶段
int MpiSortExchangeBuckets(const std::vector<MpiSortRecordMeta> &local_records,
                           const MpiSortRawBuffer &local_raw,
                           const std::vector<MpiSortKey> &splitters,
                           int rank,
                           int comm_size,
                           std::vector<MpiSortRecordMeta> *recv_records,
                           std::vector<unsigned char> *recv_raw,
                           std::vector<long long> *recv_source_counts,
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
    int current_bucket = 0;
    const int n_splitters = (int)splitters.size();
    for (size_t i = 0; i < local_records.size(); ++i) {
        const MpiSortRecordMeta &src = local_records[i];
        if ((unsigned long long)src.raw_offset + src.raw_len > (unsigned long long)local_raw.size) return -1;
        while (current_bucket < n_splitters &&
               MpiSortCompareMetaKey(src, splitters[(size_t)current_bucket]) > 0) {
            current_bucket++;
        }
        int bucket = current_bucket;
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

    // 第二遍写入一次性分配的连续发送缓冲。CPE 负责变长 raw bytes pack；
    // 从核侧使用安全字节拷贝，允许 BAM record 起点和长度不对齐。
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
    if (recv_source_counts) {
        *recv_source_counts = recv_meta_counts;
    }
#if RABBITBAM_SORT_ENABLE_SANITY_CHECKS
    if (MpiSortValidateMetaRawRecords(*recv_records, *recv_raw,
                                      &recv_meta_counts,
                                      "memory-post-exchange") != 0) {
        return -1;
    }
#endif

    // 显式释放大容量发送缓冲，使析构/释放成本进入日志，而不是落在未计时区间。
    double cleanup_t0 = GetTime();
    std::vector<int>().swap(bucket_ids);
    std::vector<MpiSortRecordMeta>().swap(send_meta);
    std::vector<unsigned char>().swap(send_raw);
    *cleanup_time += GetTime() - cleanup_t0;
    return 0;
}

int MpiSortKWayMergeReceivedRecords(const std::vector<long long> &source_counts,
                                    std::vector<MpiSortRecordMeta> *records) {
    if (!records || records->empty()) return 0;
    if (source_counts.empty()) return -1;

    size_t total = 0;
    std::vector<size_t> begin(source_counts.size(), 0);
    std::vector<size_t> end(source_counts.size(), 0);
    std::vector<size_t> cursor(source_counts.size(), 0);
    for (size_t i = 0; i < source_counts.size(); ++i) {
        if (source_counts[i] < 0) return -1;
        begin[i] = total;
        cursor[i] = total;
        total += (size_t)source_counts[i];
        end[i] = total;
    }
    if (total != records->size()) return -1;

    std::vector<MpiSortRecordMeta> merged(records->size());
    MpiSortMetaLess less;
    for (size_t out = 0; out < merged.size(); ++out) {
        int best = -1;
        for (size_t src = 0; src < source_counts.size(); ++src) {
            if (cursor[src] >= end[src]) continue;
            if (best < 0 || less((*records)[cursor[src]], (*records)[cursor[(size_t)best]])) {
                best = (int)src;
            }
        }
        if (best < 0) return -1;
        merged[out] = (*records)[cursor[(size_t)best]++];
    }
    records->swap(merged);
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
#if RABBITBAM_SORT_ENABLE_SANITY_CHECKS
        if (MpiSortValidatePayloadRecords(payload, "memory-final") != 0) {
            return -1;
        }
#endif

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
                                    swbam::RankBodySink &body_sink,
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
                if (swbam::AppendBgzfBlock(
                        &body_sink, paras_pending[i].output_block) != 0) {
                    return -1;
                }
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

struct MpiSortExternalRunRange {
    size_t byte_begin;
    size_t byte_end;
    long long global_block_begin;
    int block_count;
};

struct MpiSortExternalRun {
    std::vector<MpiSortRecordMeta> records;
    std::vector<unsigned char> raw;
};

struct MpiSortExternalReceivedBatch {
    std::vector<MpiSortRecordMeta> records;
    std::vector<unsigned char> raw;
};

struct MpiSortExternalSegment {
    const MpiSortRecordMeta *records_begin;
    const MpiSortRecordMeta *records_end;
    const unsigned char *raw;
    size_t raw_size;
    size_t record_count;
};

struct MpiSortMergeCursor {
    const MpiSortRecordMeta *current;
    const MpiSortRecordMeta *end;
    const unsigned char *raw;
    size_t raw_size;
};

class MpiSortLoserTree {
public:
    explicit MpiSortLoserTree(std::vector<MpiSortMergeCursor> cursors)
        : cursors_(std::move(cursors)),
          losers_(cursors_.size(), (int)cursors_.size()),
          sentinel_((int)cursors_.size()) {
        for (int i = (int)cursors_.size() - 1; i >= 0; --i) {
            Adjust(i);
        }
    }

    bool empty() const {
        if (cursors_.empty()) return true;
        int winner = losers_[0];
        return winner == sentinel_ || !Active(winner);
    }

    const MpiSortMergeCursor &winner() const {
        return cursors_[(size_t)losers_[0]];
    }

    void advance() {
        int winner = losers_[0];
        cursors_[(size_t)winner].current++;
        Adjust(winner);
    }

private:
    bool Active(int player) const {
        return player >= 0 && player < sentinel_ &&
               cursors_[(size_t)player].current < cursors_[(size_t)player].end;
    }

    bool Greater(int lhs, int rhs) const {
        if (lhs == sentinel_) return false;
        if (rhs == sentinel_) return true;
        const bool lhs_active = Active(lhs);
        const bool rhs_active = Active(rhs);
        if (!lhs_active || !rhs_active) {
            if (lhs_active != rhs_active) return !lhs_active;
            return lhs > rhs;
        }
        MpiSortMetaLess less;
        const MpiSortRecordMeta &a = *cursors_[(size_t)lhs].current;
        const MpiSortRecordMeta &b = *cursors_[(size_t)rhs].current;
        if (less(b, a)) return true;
        if (less(a, b)) return false;
        return lhs > rhs;
    }

    void Adjust(int player) {
        int parent = (player + sentinel_) >> 1;
        while (parent > 0) {
            if (Greater(player, losers_[(size_t)parent])) {
                std::swap(player, losers_[(size_t)parent]);
            }
            parent >>= 1;
        }
        losers_[0] = player;
    }

    std::vector<MpiSortMergeCursor> cursors_;
    std::vector<int> losers_;
    int sentinel_;
};

size_t MpiSortEstimateBlockMemory(uint32_t isize) {
    size_t n = (size_t)(isize / 64u) + 1u;
    if (n > (size_t)MPI_RECORDS_PER_BLOCK) n = (size_t)MPI_RECORDS_PER_BLOCK;
    return (size_t)isize + n * sizeof(MpiSortRecordMeta);
}

size_t MpiSortExternalRunBudget(size_t memory_limit) {
    if (memory_limit == 0) return (size_t)-1 / 4;
    size_t budget = memory_limit / 3;
    if (budget < 256 * 1024 && memory_limit >= 256 * 1024) budget = 256 * 1024;
    if (budget == 0) budget = memory_limit;
    return budget;
}

int MpiSortBuildExternalRunRanges(MemReader &reader,
                                  long long global_block_begin,
                                  size_t memory_limit,
                                  std::vector<MpiSortExternalRunRange> *ranges) {
    ranges->clear();
    const size_t run_budget = MpiSortExternalRunBudget(memory_limit);
    size_t pos = 0;
    size_t run_begin = 0;
    size_t run_estimate = 0;
    int run_blocks = 0;
    long long block_index = 0;
    long long run_global_begin = global_block_begin;

    while (pos < reader.size) {
        if (pos + BLOCK_HEADER_LENGTH > reader.size) return -1;
        uint16_t bsize16 = MpiSortReadLe16((const unsigned char *)reader.base + pos + 16);
        size_t bsize = (size_t)bsize16 + 1u;
        if (bsize == 28) break;
        if (bsize < BLOCK_HEADER_LENGTH || pos + bsize > reader.size) return -1;
        uint32_t isize = MpiSortReadLe32((const unsigned char *)reader.base + pos + bsize - 4);
        if (isize > BGZF_MAX_BLOCK_SIZE) return -1;
        size_t block_estimate = MpiSortEstimateBlockMemory(isize);

        if (run_blocks > 0 && run_estimate + block_estimate > run_budget) {
            MpiSortExternalRunRange range;
            range.byte_begin = run_begin;
            range.byte_end = pos;
            range.global_block_begin = run_global_begin;
            range.block_count = run_blocks;
            ranges->push_back(range);
            run_begin = pos;
            run_global_begin = global_block_begin + block_index;
            run_estimate = 0;
            run_blocks = 0;
        }

        if (run_blocks == 0 && block_estimate > run_budget) {
            fprintf(stderr,
                    "ERROR: MPI sort external run needs at least one BGZF block in memory. block_estimate=%zu budget=%zu limit=%zu\n",
                    block_estimate, run_budget, memory_limit);
            return -2;
        }

        run_estimate += block_estimate;
        run_blocks++;
        pos += bsize;
        block_index++;
    }

    if (run_blocks > 0) {
        MpiSortExternalRunRange range;
        range.byte_begin = run_begin;
        range.byte_end = pos;
        range.global_block_begin = run_global_begin;
        range.block_count = run_blocks;
        ranges->push_back(range);
    }
    return 0;
}

int MpiSortStoreExternalRun(const std::vector<MpiSortRecordMeta> &records,
                            const MpiSortRawBuffer &raw,
                            MpiSortExternalRun *run,
                            MpiSortStats *stats) {
    double t0 = GetTime();
    run->records.assign(records.begin(), records.end());
    run->raw.assign(raw.data, raw.data + raw.size);
    stats->t_temp_write_sim += GetTime() - t0;
    return 0;
}

int MpiSortReadExternalRun(const MpiSortExternalRun &run,
                           std::vector<MpiSortRecordMeta> *records,
                           MpiSortRawBuffer *raw,
                           MpiSortStats *stats) {
    double t0 = GetTime();
    *records = run.records;
    if (MpiSortRawReserve(raw, run.raw.size()) != 0) return -1;
    if (!run.raw.empty()) memcpy(raw->data, run.raw.data(), run.raw.size());
    raw->size = run.raw.size();
    stats->t_temp_read_sim += GetTime() - t0;
    return 0;
}

int MpiSortAppendReceivedSegments(std::vector<MpiSortRecordMeta> *recv_records,
                                  std::vector<unsigned char> *recv_raw,
                                  const std::vector<long long> &source_counts,
                                  std::vector<MpiSortExternalReceivedBatch> *batches,
                                  std::vector<MpiSortExternalSegment> *segments,
                                  MpiSortStats *stats) {
    if (!recv_records || !recv_raw || !batches || !segments) return -1;

    size_t meta_begin = 0;
    for (size_t src = 0; src < source_counts.size(); ++src) {
        if (source_counts[src] < 0) return -1;
        size_t count = (size_t)source_counts[src];
        size_t meta_end = meta_begin + count;
        if (meta_end > recv_records->size()) return -1;
        if (count == 0) {
            meta_begin = meta_end;
            continue;
        }

        for (size_t i = meta_begin; i < meta_end; ++i) {
            const MpiSortRecordMeta &rec = (*recv_records)[i];
            uint64_t end = rec.raw_offset + (uint64_t)rec.raw_len;
            if (end > recv_raw->size()) return -1;
        }
        meta_begin = meta_end;
    }
    if (meta_begin != recv_records->size()) return -1;

    MpiSortExternalReceivedBatch batch;
    batch.records = std::move(*recv_records);
    batch.raw = std::move(*recv_raw);
    batches->push_back(std::move(batch));
    const MpiSortExternalReceivedBatch &stored = batches->back();

    meta_begin = 0;
    for (size_t src = 0; src < source_counts.size(); ++src) {
        size_t count = (size_t)source_counts[src];
        if (count > 0) {
            MpiSortExternalSegment seg;
            seg.records_begin = stored.records.data() + meta_begin;
            seg.records_end = seg.records_begin + count;
            seg.raw = stored.raw.data();
            seg.raw_size = stored.raw.size();
            seg.record_count = count;
            segments->push_back(seg);
        }
        meta_begin += count;
    }
    (void)stats;
    return 0;
}

int MpiSortPreparePayloadBatchFromSegments(MpiSortLoserTree *tree,
                                           MpiSortBlockSet *payload_blocks,
                                           MpiSortBlockSet *output_blocks,
                                           MpiSortRawCompressPara *paras,
                                           int compress_level,
                                           int *active_blocks) {
    if (!tree) return -1;

    int active = 0;
    while (!tree->empty() && active < kSortCpeBlocks) {
        bam_block *payload = &payload_blocks->blocks[active];
        payload->pos = 0;
        payload->length = 0;
        payload->errcode = 0;
        payload->block_id = active;

        while (!tree->empty()) {
            const MpiSortMergeCursor &cur = tree->winner();
            const MpiSortRecordMeta &rec = *cur.current;
            if (rec.raw_len > BGZF_BLOCK_SIZE) return -1;
            if ((unsigned long long)rec.raw_offset + rec.raw_len >
                (unsigned long long)cur.raw_size) return -1;
            if (payload->pos > 0 && payload->pos + rec.raw_len > BGZF_BLOCK_SIZE) {
                break;
            }
            memcpy(payload->data + payload->pos,
                   cur.raw + rec.raw_offset, rec.raw_len);
            payload->pos += rec.raw_len;
            payload->length = payload->pos;
            tree->advance();
        }
        if (payload->pos == 0) return -1;
#if RABBITBAM_SORT_ENABLE_SANITY_CHECKS
        if (MpiSortValidatePayloadRecords(payload, "external-memory") != 0) {
            return -1;
        }
#endif

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

int MpiSortCompressExternalSegmentsCpe(const std::vector<MpiSortExternalReceivedBatch> &batches,
                                       const std::vector<MpiSortExternalSegment> &segments,
                                       swbam::RankBodySink &body_sink,
                                       int compress_level,
                                       MpiSortStats *stats) {
    double setup_t0 = GetTime();
    MpiSortBlockSet payload_a = {}, payload_b = {}, out_a = {}, out_b = {};
    MpiSortRawCompressPara paras_a[kSortCpeBlocks], paras_b[kSortCpeBlocks];
    if (MpiSortAllocateBlockSet(&payload_a, kSortCpeBlocks) != 0 ||
        MpiSortAllocateBlockSet(&payload_b, kSortCpeBlocks) != 0 ||
        MpiSortAllocateBlockSet(&out_a, kSortCpeBlocks) != 0 ||
        MpiSortAllocateBlockSet(&out_b, kSortCpeBlocks) != 0) {
        fprintf(stderr, "ERROR: failed to allocate MPI sort external compress workspace.\n");
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

    std::vector<MpiSortMergeCursor> cursors;
    cursors.reserve(segments.size());
    for (size_t i = 0; i < segments.size(); ++i) {
        if (segments[i].record_count > 0) {
            MpiSortMergeCursor cur;
            cur.current = segments[i].records_begin;
            cur.end = segments[i].records_end;
            cur.raw = segments[i].raw;
            cur.raw_size = segments[i].raw_size;
            cursors.push_back(cur);
        }
    }
    MpiSortLoserTree tree(std::move(cursors));

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
                if (swbam::AppendBgzfBlock(
                        &body_sink, paras_pending[i].output_block) != 0) {
                    return -1;
                }
                stats->bgzf_blocks++;
            }
            MpiSortInitEmptyRawCompressPara(&paras_pending[i], i);
        }
        has_pending = false;
        stats->t_write += GetTime() - t0;
        return 0;
    };

    int active_blocks = 0;
    double merge_t0 = GetTime();
    stats->t_setup += GetTime() - setup_t0;

    if (MpiSortPreparePayloadBatchFromSegments(&tree,
                                               payload_active, out_active,
                                               paras_active, compress_level, &active_blocks) != 0) {
        MpiSortFreeBlockSet(&payload_a);
        MpiSortFreeBlockSet(&payload_b);
        MpiSortFreeBlockSet(&out_a);
        MpiSortFreeBlockSet(&out_b);
        return -1;
    }
    double merge_wall = GetTime() - merge_t0;
    stats->t_final_sort += merge_wall;
    stats->t_run_merge += merge_wall;
    stats->t_merge_unhidden += merge_wall;

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
        merge_t0 = GetTime();
        if (MpiSortPreparePayloadBatchFromSegments(&tree,
                                                   payload_pending, out_pending,
                                                   paras_pending, compress_level, &next_active_blocks) != 0) {
            athread_join();
            ret_code = -1;
            break;
        }
        merge_wall = GetTime() - merge_t0;
        stats->t_final_sort += merge_wall;
        stats->t_run_merge += merge_wall;
        athread_join();
        double wall = GetTime() - t0;
        stats->t_compress += wall;
        MpiSortAccumulateCompressDetail(paras_active, active_blocks, wall, stats);

        double status_t0 = GetTime();
        for (int i = 0; i < active_blocks; ++i) {
            if (paras_active[i].status != 0) {
                fprintf(stderr, "ERROR: MPI sort external raw compress failed on block %d with status %d.\n",
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

const size_t kSortStrictControlReserve = 16ull * 1024ull * 1024ull;
const size_t kSortStrictSimScratch = 8ull * 1024ull * 1024ull;
const size_t kSortStrictMetaBuffer = 64ull * 1024ull;
const size_t kSortStrictRawBuffer = 256ull * 1024ull;
const size_t kSortStrictConsolidateRaw = 4ull * 1024ull * 1024ull;
const size_t kSortStrictConsolidateMeta = 2ull * 1024ull * 1024ull;
const int kSortStrictTagHeader = 7200;
const int kSortStrictTagMeta = 7201;
const int kSortStrictTagRaw = 7202;
const int kSortStrictTagSimHeader = 7210;
const int kSortStrictTagSimMeta = 7211;
const int kSortStrictTagSimRaw = 7212;

struct MpiSortMemoryTracker {
    size_t limit;
    size_t current;
    size_t peak;

    explicit MpiSortMemoryTracker(size_t memory_limit)
        : limit(memory_limit), current(0), peak(0) {}

    bool acquire(size_t bytes) {
        if (bytes > SIZE_MAX - current) return false;
        if (limit != 0 && current + bytes > limit) return false;
        current += bytes;
        if (current > peak) peak = current;
        return true;
    }

    void release(size_t bytes) {
        current = bytes > current ? 0 : current - bytes;
    }
};

struct MpiSortStrictBuffer {
    unsigned char *data;
    size_t size;

    MpiSortStrictBuffer() : data(nullptr), size(0) {}
};

int MpiSortStrictAlloc(MpiSortMemoryTracker *tracker,
                       size_t bytes,
                       MpiSortStrictBuffer *buffer) {
    buffer->data = nullptr;
    buffer->size = 0;
    if (bytes == 0) return 0;
    if (!tracker->acquire(bytes)) return -1;
    buffer->data = aligned_alloc_custom(64, bytes);
    if (!buffer->data) {
        tracker->release(bytes);
        return -1;
    }
    buffer->size = bytes;
    return 0;
}

void MpiSortStrictFree(MpiSortMemoryTracker *tracker,
                       MpiSortStrictBuffer *buffer) {
    if (buffer->data) aligned_free_custom(buffer->data);
    tracker->release(buffer->size);
    buffer->data = nullptr;
    buffer->size = 0;
}

struct MpiSortStrictTempStore {
    int runs_fd;
    int segments_fd;
    off_t runs_end;
    off_t segments_end;

    MpiSortStrictTempStore()
        : runs_fd(-1), segments_fd(-1), runs_end(0), segments_end(0) {}
};

struct MpiSortStrictExtent {
    int fd;
    off_t meta_offset;
    off_t raw_offset;
    uint64_t record_count;
    uint64_t raw_size;
};

struct MpiSortStrictRun {
    MpiSortStrictExtent extent;
    bool memory_backed;

    MpiSortStrictRun()
        : extent{-1, 0, 0, 0, 0}, memory_backed(false) {}
};

struct MpiSortStrictSegment {
    std::vector<MpiSortStrictExtent> extents;
    uint64_t record_count;
    uint64_t raw_size;

    MpiSortStrictSegment() : record_count(0), raw_size(0) {}
};

std::string MpiSortStrictTempRoot(const char *temp_prefix) {
    std::string root = temp_prefix && temp_prefix[0] ? temp_prefix : "./rabbitbam-sort";
    struct stat st;
    if (stat(root.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
        if (!root.empty() && root[root.size() - 1] != '/') root += '/';
        root += "rabbitbam-sort";
    }
    return root;
}

int MpiSortStrictOpenTemp(const char *temp_prefix,
                          int rank,
                          MpiSortStrictTempStore *store) {
    const std::string root = MpiSortStrictTempRoot(temp_prefix);
    const long pid = (long)getpid();
    for (int attempt = 0; attempt < 1000; ++attempt) {
        char run_path[4096];
        char segment_path[4096];
        snprintf(run_path, sizeof(run_path), "%s.%ld.%d.%d.runs.tmp",
                 root.c_str(), pid, rank, attempt);
        snprintf(segment_path, sizeof(segment_path), "%s.%ld.%d.%d.segments.tmp",
                 root.c_str(), pid, rank, attempt);
        int run_fd = open(run_path, O_RDWR | O_CREAT | O_EXCL, 0600);
        if (run_fd < 0) {
            if (errno == EEXIST) continue;
            fprintf(stderr, "[rank %d] ERROR: cannot create sort temp file %s: %s\n",
                    rank, run_path, strerror(errno));
            return -1;
        }
        int segment_fd = open(segment_path, O_RDWR | O_CREAT | O_EXCL, 0600);
        if (segment_fd < 0) {
            int saved_errno = errno;
            close(run_fd);
            unlink(run_path);
            if (saved_errno == EEXIST) continue;
            fprintf(stderr, "[rank %d] ERROR: cannot create sort temp file %s: %s\n",
                    rank, segment_path, strerror(saved_errno));
            return -1;
        }
        int unlink_run_ret = unlink(run_path);
        int unlink_segment_ret = unlink(segment_path);
        if (unlink_run_ret != 0 || unlink_segment_ret != 0) {
            int saved_errno = errno;
            close(run_fd);
            close(segment_fd);
            if (unlink_run_ret != 0) unlink(run_path);
            if (unlink_segment_ret != 0) unlink(segment_path);
            fprintf(stderr, "[rank %d] ERROR: cannot unlink open sort temp files: %s\n",
                    rank, strerror(saved_errno));
            return -1;
        }
        store->runs_fd = run_fd;
        store->segments_fd = segment_fd;
        store->runs_end = 0;
        store->segments_end = 0;
        return 0;
    }
    fprintf(stderr, "[rank %d] ERROR: exhausted unique sort temp file attempts.\n", rank);
    return -1;
}

void MpiSortStrictCloseTemp(MpiSortStrictTempStore *store) {
    if (store->runs_fd >= 0) close(store->runs_fd);
    if (store->segments_fd >= 0) close(store->segments_fd);
    store->runs_fd = -1;
    store->segments_fd = -1;
}

void MpiSortStrictSimulateCopy(const void *src,
                               size_t bytes,
                               unsigned char *scratch,
                               size_t scratch_size,
                               double *counter) {
    if (!src || bytes == 0 || !scratch || scratch_size == 0) return;
    double t0 = GetTime();
    const unsigned char *p = (const unsigned char *)src;
    size_t pos = 0;
    volatile unsigned char guard = 0;
    while (pos < bytes) {
        size_t n = std::min(scratch_size, bytes - pos);
        memcpy(scratch, p + pos, n);
        guard ^= scratch[0];
        guard ^= scratch[n - 1];
        pos += n;
    }
    if (guard == 255 && bytes == 0) scratch[0] = guard;
    *counter += GetTime() - t0;
}

int MpiSortStrictPwriteAll(int fd,
                           const void *data,
                           size_t bytes,
                           off_t offset,
                           MpiSortStats *stats) {
    const unsigned char *p = (const unsigned char *)data;
    size_t done = 0;
    double t0 = GetTime();
    while (done < bytes) {
        size_t request = std::min(
            bytes - done, (size_t)std::numeric_limits<ssize_t>::max());
        ssize_t n = pwrite(fd, p + done, request, offset + (off_t)done);
        if (n < 0 && errno == EINTR) continue;
        if (n <= 0) {
            fprintf(stderr, "ERROR: sort temp pwrite failed: %s\n", strerror(errno));
            stats->t_temp_write_actual += GetTime() - t0;
            return -1;
        }
        done += (size_t)n;
    }
    stats->t_temp_write_actual += GetTime() - t0;
    stats->temp_write_bytes += (long long)bytes;
    return 0;
}

int MpiSortStrictPreadAll(int fd,
                          void *data,
                          size_t bytes,
                          off_t offset,
                          MpiSortStats *stats) {
    unsigned char *p = (unsigned char *)data;
    size_t done = 0;
    double t0 = GetTime();
    while (done < bytes) {
        size_t request = std::min(
            bytes - done, (size_t)std::numeric_limits<ssize_t>::max());
        ssize_t n = pread(fd, p + done, request, offset + (off_t)done);
        if (n < 0 && errno == EINTR) continue;
        if (n <= 0) {
            fprintf(stderr, "ERROR: sort temp pread failed or reached EOF: %s\n",
                    n < 0 ? strerror(errno) : "short read");
            stats->t_temp_read_actual += GetTime() - t0;
            return -1;
        }
        done += (size_t)n;
    }
    stats->t_temp_read_actual += GetTime() - t0;
    stats->temp_read_bytes += (long long)bytes;
    return 0;
}

int MpiSortStrictAppendExtent(int fd,
                              off_t *file_end,
                              const MpiSortRecordMeta *records,
                              size_t record_count,
                              const unsigned char *raw,
                              size_t raw_size,
                              unsigned char *sim_scratch,
                              size_t sim_scratch_size,
                              MpiSortStats *stats,
                              MpiSortStrictExtent *extent) {
    for (size_t i = 0; i < record_count; ++i) {
        if (records[i].raw_len < 36 ||
            records[i].raw_offset > raw_size ||
            records[i].raw_len >
                raw_size - (size_t)records[i].raw_offset ||
            (uint64_t)MpiSortReadLe32(raw + records[i].raw_offset) + 4u !=
                records[i].raw_len) {
            fprintf(stderr,
                    "ERROR: external sort attempted to spill an invalid "
                    "record. record=%zu raw_offset=%llu raw_len=%u "
                    "raw_size=%zu\n",
                    i, (unsigned long long)records[i].raw_offset,
                    records[i].raw_len, raw_size);
            return -1;
        }
    }
    const size_t meta_bytes = record_count * sizeof(MpiSortRecordMeta);
    extent->fd = fd;
    extent->meta_offset = *file_end;
    extent->record_count = (uint64_t)record_count;
    extent->raw_size = (uint64_t)raw_size;
    if (MpiSortStrictPwriteAll(fd, records, meta_bytes, *file_end, stats) != 0) return -1;
    MpiSortStrictSimulateCopy(records, meta_bytes, sim_scratch, sim_scratch_size,
                              &stats->t_temp_write_sim);
    *file_end += (off_t)meta_bytes;
    extent->raw_offset = *file_end;
    if (MpiSortStrictPwriteAll(fd, raw, raw_size, *file_end, stats) != 0) return -1;
    MpiSortStrictSimulateCopy(raw, raw_size, sim_scratch, sim_scratch_size,
                              &stats->t_temp_write_sim);
    *file_end += (off_t)raw_size;
    return 0;
}

struct MpiSortStrictRunArena {
    unsigned char *data;
    size_t capacity;
    size_t raw_used;
    size_t meta_begin;
    size_t record_count;

    MpiSortStrictRunArena()
        : data(nullptr), capacity(0), raw_used(0), meta_begin(0), record_count(0) {}

    void reset() {
        raw_used = 0;
        meta_begin = capacity - capacity % sizeof(MpiSortRecordMeta);
        record_count = 0;
    }

    MpiSortRecordMeta *records() {
        return (MpiSortRecordMeta *)(data + meta_begin);
    }

    bool can_append(size_t raw_bytes, size_t records_count) const {
        if (records_count > SIZE_MAX / sizeof(MpiSortRecordMeta)) return false;
        size_t meta_bytes = records_count * sizeof(MpiSortRecordMeta);
        if (meta_bytes > meta_begin) return false;
        size_t new_meta_begin = meta_begin - meta_bytes;
        return raw_bytes <= new_meta_begin - raw_used;
    }
};

int MpiSortStrictAppendBlock(MpiSortStrictRunArena *arena,
                             const unsigned char *raw,
                             size_t raw_size,
                             uint64_t source_raw_base,
                             const MpiSortRecordMeta *records,
                             size_t record_count) {
    if (!arena->can_append(raw_size, record_count)) return 1;
    size_t new_meta_begin = arena->meta_begin -
                            record_count * sizeof(MpiSortRecordMeta);
    MpiSortRecordMeta *dst =
        (MpiSortRecordMeta *)(arena->data + new_meta_begin);
    for (size_t i = 0; i < record_count; ++i) {
        if (records[i].raw_offset < source_raw_base ||
            records[i].raw_offset - source_raw_base + records[i].raw_len > raw_size) {
            return -1;
        }
        dst[i] = records[i];
        dst[i].raw_offset = (uint64_t)arena->raw_used +
                            records[i].raw_offset - source_raw_base;
    }
    memcpy(arena->data + arena->raw_used, raw, raw_size);
    arena->raw_used += raw_size;
    arena->meta_begin = new_meta_begin;
    arena->record_count += record_count;
    return 0;
}

void MpiSortStrictMakeSamples(const MpiSortRecordMeta *records,
                              size_t record_count,
                              int comm_size,
                              std::vector<MpiSortKey> *samples) {
    if (!records || record_count == 0) return;
    size_t target = (size_t)comm_size * kSortSampleFactor;
    size_t n = std::min(target, record_count);
    for (size_t i = 0; i < n; ++i) {
        size_t idx = ((i + 1) * record_count) / (n + 1);
        if (idx >= record_count) idx = record_count - 1;
        samples->push_back(MpiSortKeyFromMeta(records[idx]));
    }
}

int MpiSortStrictAppendSortedRun(MpiSortStrictRunArena *arena,
                                 MpiSortStrictTempStore *store,
                                 unsigned char *sim_scratch,
                                 size_t sim_scratch_size,
                                 MpiSortStats *stats,
                                 MpiSortStrictExtent *extent) {
    return MpiSortStrictAppendExtent(
        store->runs_fd, &store->runs_end,
        arena->records(), arena->record_count,
        arena->data, arena->raw_used,
        sim_scratch, sim_scratch_size, stats, extent);
}

void MpiSortStrictSortRun(MpiSortStrictRunArena *arena,
                          int comm_size,
                          std::vector<MpiSortKey> *samples,
                          MpiSortStats *stats) {
    double t0 = GetTime();
    std::sort(arena->records(), arena->records() + arena->record_count,
              MpiSortMetaLess());
    double wall = GetTime() - t0;
    stats->t_local_sort += wall;
    stats->t_run_sort += wall;
    MpiSortStrictMakeSamples(arena->records(), arena->record_count,
                             comm_size, samples);
}

int MpiSortStrictSpillRun(MpiSortStrictRunArena *arena,
                          MpiSortStrictTempStore *store,
                          unsigned char *sim_scratch,
                          size_t sim_scratch_size,
                          int comm_size,
                          std::vector<MpiSortStrictRun> *runs,
                          std::vector<MpiSortKey> *samples,
                          MpiSortStats *stats) {
    if (arena->record_count == 0) return 0;
    MpiSortStrictSortRun(arena, comm_size, samples, stats);
    MpiSortStrictRun run;
    if (MpiSortStrictAppendSortedRun(
            arena, store, sim_scratch, sim_scratch_size,
            stats, &run.extent) != 0) {
        return -1;
    }
    runs->push_back(run);
    arena->reset();
    return 0;
}

void MpiSortStrictRetainRun(MpiSortStrictRunArena *arena,
                            int comm_size,
                            std::vector<MpiSortStrictRun> *runs,
                            std::vector<MpiSortKey> *samples,
                            MpiSortStats *stats) {
    if (arena->record_count == 0) return;
    MpiSortStrictSortRun(arena, comm_size, samples, stats);
    MpiSortStrictRun run;
    run.memory_backed = true;
    run.extent.record_count = (uint64_t)arena->record_count;
    run.extent.raw_size = (uint64_t)arena->raw_used;
    stats->resident_run_records = (long long)arena->record_count;
    stats->resident_run_raw_bytes = (long long)arena->raw_used;
    runs->insert(runs->begin(), run);
}

size_t MpiSortStrictExtractWorkspaceBytes() {
    return 2ull * kSortCpeBlocks *
               (sizeof(bam_block) + (size_t)BGZF_MAX_BLOCK_SIZE) +
           (size_t)kSortCpeBlocks * MPI_RECORDS_PER_BLOCK *
               sizeof(MpiSortRecordMeta);
}

int MpiSortStrictGenerateRunsImpl(MpiSortBatchReader *block_reader,
                              long long global_block_begin,
                              int rank,
                              int comm_size,
                              MpiSortStrictRunArena *arena,
                              MpiSortStrictTempStore *store,
                              unsigned char *sim_scratch,
                              size_t sim_scratch_size,
                              MpiSortMemoryTracker *tracker,
                              std::vector<MpiSortStrictRun> *runs,
                              std::vector<MpiSortKey> *samples,
                              MpiSortStats *stats) {
    const size_t workspace_bytes = MpiSortStrictExtractWorkspaceBytes();
    if (!tracker->acquire(workspace_bytes)) {
        fprintf(stderr,
                "[rank %d] ERROR: -m cannot hold the minimum external extract workspace. required_at_least=%zu limit=%zu\n",
                rank, tracker->current + workspace_bytes, tracker->limit);
        return -1;
    }
    MpiSortExtractWorkspace ws = {};
    if (MpiSortAllocateExtractWorkspace(&ws, kSortCpeBlocks,
                                        (int)MPI_RECORDS_PER_BLOCK, 0) != 0) {
        tracker->release(workspace_bytes);
        return -1;
    }
    MpiSortExtractPara paras[kSortCpeBlocks];
    size_t local_block = 0;
    int ret = 0;
    while (true) {
        int n_blocks = 0;
        double read_t0 = GetTime();
        const int read_ret = block_reader->Read(
            &ws.input_blocks, &n_blocks);
        double read_wall = GetTime() - read_t0;
        stats->t_extract_read += read_wall;
        stats->t_extract_read_unhidden += read_wall;
        if (read_ret != 0) {
            ret = -1;
            break;
        }
        if (n_blocks == 0) break;
        stats->input_blocks += n_blocks;

        double prepare_t0 = GetTime();
        for (int b = 0; b < kSortCpeBlocks; ++b) {
            uint32_t isize = b < n_blocks ? MpiSortBgzfISize(&ws.input_blocks.blocks[b]) : 0;
            if (b < n_blocks && isize > BGZF_MAX_BLOCK_SIZE) {
                fprintf(stderr, "[rank %d] ERROR: invalid BGZF ISIZE=%u.\n", rank, isize);
                ret = -1;
                break;
            }
            if (b < n_blocks) {
                ws.un_blocks.blocks[b].data =
                    ws.un_blocks.data + (size_t)b * BGZF_MAX_BLOCK_SIZE;
                ws.un_blocks.blocks[b].length = (int)isize;
                ws.un_blocks.blocks[b].pos = 0;
                ws.un_blocks.blocks[b].errcode = 0;
            }
            paras[b].block_id = b;
            paras[b].input_block = b < n_blocks ? &ws.input_blocks.blocks[b] : nullptr;
            paras[b].un_comp_block = b < n_blocks ? &ws.un_blocks.blocks[b] : nullptr;
            paras[b].raw_arena = b < n_blocks ? ws.un_blocks.blocks[b].data : nullptr;
            paras[b].raw_capacity = b < n_blocks ? (size_t)isize : 0;
            paras[b].raw_used = 0;
            paras[b].raw_base_offset = (uint64_t)b * BGZF_MAX_BLOCK_SIZE;
            paras[b].records = ws.records + (size_t)b * MPI_RECORDS_PER_BLOCK;
            paras[b].record_capacity = (int)MPI_RECORDS_PER_BLOCK;
            paras[b].n_records = 0;
            paras[b].global_block_index =
                global_block_begin + (long long)local_block + b;
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
        if (ret != 0) break;

        double extract_t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_sort_extract_raw, paras, 1);
        athread_join();
        double extract_wall = GetTime() - extract_t0;
        stats->t_extract += extract_wall;
        MpiSortAccumulateExtractDetail(paras, n_blocks, extract_wall, stats);

        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].status != 0) {
                if (paras[b].status == -4) {
                    fprintf(stderr,
                            "[rank %d] ERROR: sort assumes BAM records do not cross BGZF blocks; global_block=%lld.\n",
                            rank, global_block_begin + (long long)local_block + b);
                } else {
                    fprintf(stderr,
                            "[rank %d] ERROR: external sort extract failed. global_block=%lld status=%d limit=%lld actual=%lld\n",
                            rank, global_block_begin + (long long)local_block + b,
                            paras[b].status, paras[b].limit_value, paras[b].actual_value);
                }
                ret = -1;
                break;
            }
            const size_t block_raw = paras[b].raw_used;
            const size_t block_records = (size_t)paras[b].n_records;
            const uint64_t block_base = (uint64_t)b * BGZF_MAX_BLOCK_SIZE;
            double append_t0 = GetTime();
            int append_ret = MpiSortStrictAppendBlock(
                arena, ws.un_blocks.blocks[b].data, block_raw, block_base,
                ws.records + (size_t)b * MPI_RECORDS_PER_BLOCK,
                block_records);
            stats->t_extract_merge += GetTime() - append_t0;
            if (append_ret == 1 && arena->record_count > 0) {
                if (MpiSortStrictSpillRun(arena, store, sim_scratch,
                                          sim_scratch_size, comm_size,
                                          runs, samples, stats) != 0) {
                    ret = -1;
                    break;
                }
                append_t0 = GetTime();
                append_ret = MpiSortStrictAppendBlock(
                    arena, ws.un_blocks.blocks[b].data, block_raw, block_base,
                    ws.records + (size_t)b * MPI_RECORDS_PER_BLOCK,
                    block_records);
                stats->t_extract_merge += GetTime() - append_t0;
            }
            if (append_ret != 0) {
                size_t minimum = kSortStrictControlReserve + kSortStrictSimScratch +
                                 workspace_bytes + block_raw +
                                 block_records * sizeof(MpiSortRecordMeta);
                fprintf(stderr,
                        "[rank %d] ERROR: one BGZF block cannot fit in the external run arena. minimum_memory=%zu limit=%zu block_raw=%zu records=%zu\n",
                        rank, minimum, tracker->limit, block_raw, block_records);
                ret = -1;
                break;
            }
            stats->local_records += (long long)block_records;
        }
        if (ret != 0) break;
        local_block += (size_t)n_blocks;
    }
    if (ret == 0) {
        MpiSortStrictRetainRun(arena, comm_size, runs, samples, stats);
    }
    MpiSortFreeExtractWorkspace(&ws);
    tracker->release(workspace_bytes);
    return ret;
}

int MpiSortStrictGenerateRuns(
        MemReader &reader, long long global_block_begin,
        int rank, int comm_size,
        MpiSortStrictRunArena *arena,
        MpiSortStrictTempStore *store,
        unsigned char *sim_scratch, size_t sim_scratch_size,
        MpiSortMemoryTracker *tracker,
        std::vector<MpiSortStrictRun> *runs,
        std::vector<MpiSortKey> *samples,
        MpiSortStats *stats) {
    MpiSortMemoryBatchReader block_reader(&reader);
    return MpiSortStrictGenerateRunsImpl(
        &block_reader, global_block_begin, rank, comm_size,
        arena, store, sim_scratch, sim_scratch_size, tracker,
        runs, samples, stats);
}

int MpiSortStrictGenerateRuns(
        const swbam::BamInputBackend &input,
        const swbam::BgzfBlockSpan *spans, size_t span_count,
        long long global_block_begin,
        int rank, int comm_size,
        MpiSortStrictRunArena *arena,
        MpiSortStrictTempStore *store,
        unsigned char *sim_scratch, size_t sim_scratch_size,
        MpiSortMemoryTracker *tracker,
        std::vector<MpiSortStrictRun> *runs,
        std::vector<MpiSortKey> *samples,
        MpiSortStats *stats) {
    MpiSortBackendBatchReader block_reader(&input, spans, span_count);
    return MpiSortStrictGenerateRunsImpl(
        &block_reader, global_block_begin, rank, comm_size,
        arena, store, sim_scratch, sim_scratch_size, tracker,
        runs, samples, stats);
}

int MpiSortStrictLoadRun(const MpiSortStrictRun &run,
                         MpiSortStrictRunArena *arena,
                         MpiSortStrictTempStore *store,
                         unsigned char *sim_scratch,
                         size_t sim_scratch_size,
                         MpiSortStats *stats) {
    if (run.memory_backed) return -1;
    arena->reset();
    size_t count = (size_t)run.extent.record_count;
    size_t raw_size = (size_t)run.extent.raw_size;
    size_t meta_bytes = count * sizeof(MpiSortRecordMeta);
    if (!arena->can_append(raw_size, count)) return -1;
    arena->meta_begin -= meta_bytes;
    arena->record_count = count;
    arena->raw_used = raw_size;
    if (MpiSortStrictPreadAll(run.extent.fd, arena->records(), meta_bytes,
                              run.extent.meta_offset, stats) != 0 ||
        MpiSortStrictPreadAll(run.extent.fd, arena->data, raw_size,
                              run.extent.raw_offset, stats) != 0) {
        return -1;
    }
    MpiSortStrictSimulateCopy(arena->records(), meta_bytes, sim_scratch,
                              sim_scratch_size, &stats->t_temp_read_sim);
    MpiSortStrictSimulateCopy(arena->data, raw_size, sim_scratch,
                              sim_scratch_size, &stats->t_temp_read_sim);
    (void)store;
    return 0;
}

struct MpiSortStrictExchangeWorkspace {
    MpiSortStrictBuffer send_meta_buffer;
    MpiSortStrictBuffer recv_meta_buffer;
    MpiSortStrictBuffer send_raw_buffer;
    MpiSortStrictBuffer recv_raw_buffer;
    MpiSortStrictBuffer prefix_buffer;
    size_t meta_capacity;
    size_t raw_capacity;
    size_t prefix_capacity;

    MpiSortStrictExchangeWorkspace()
        : meta_capacity(0), raw_capacity(0), prefix_capacity(0) {}
};

int MpiSortStrictAllocateExchange(MpiSortMemoryTracker *tracker,
                                  size_t budget,
                                  MpiSortStrictExchangeWorkspace *ws) {
    size_t raw_each = budget * 3 / 8;
    size_t meta_each = budget / 16;
    size_t prefix_bytes = budget - 2 * raw_each - 2 * meta_each;
    if (raw_each < BGZF_MAX_BLOCK_SIZE ||
        meta_each < 64 * sizeof(MpiSortRecordMeta) ||
        prefix_bytes < 64 * sizeof(uint64_t)) {
        return -1;
    }
    if (MpiSortStrictAlloc(tracker, raw_each, &ws->send_raw_buffer) != 0 ||
        MpiSortStrictAlloc(tracker, raw_each, &ws->recv_raw_buffer) != 0 ||
        MpiSortStrictAlloc(tracker, meta_each, &ws->send_meta_buffer) != 0 ||
        MpiSortStrictAlloc(tracker, meta_each, &ws->recv_meta_buffer) != 0 ||
        MpiSortStrictAlloc(tracker, prefix_bytes, &ws->prefix_buffer) != 0) {
        MpiSortStrictFree(tracker, &ws->send_meta_buffer);
        MpiSortStrictFree(tracker, &ws->recv_meta_buffer);
        MpiSortStrictFree(tracker, &ws->send_raw_buffer);
        MpiSortStrictFree(tracker, &ws->recv_raw_buffer);
        MpiSortStrictFree(tracker, &ws->prefix_buffer);
        return -1;
    }
    ws->meta_capacity = meta_each / sizeof(MpiSortRecordMeta);
    ws->raw_capacity = raw_each;
    ws->prefix_capacity = prefix_bytes / sizeof(uint64_t);
    return 0;
}

void MpiSortStrictFreeExchange(MpiSortMemoryTracker *tracker,
                               MpiSortStrictExchangeWorkspace *ws) {
    MpiSortStrictFree(tracker, &ws->send_meta_buffer);
    MpiSortStrictFree(tracker, &ws->recv_meta_buffer);
    MpiSortStrictFree(tracker, &ws->send_raw_buffer);
    MpiSortStrictFree(tracker, &ws->recv_raw_buffer);
    MpiSortStrictFree(tracker, &ws->prefix_buffer);
}

size_t MpiSortStrictChunkEnd(const MpiSortRecordMeta *records,
                             size_t begin,
                             size_t end,
                             size_t meta_capacity,
                             size_t prefix_capacity,
                             size_t raw_capacity,
                             size_t *raw_bytes) {
    const size_t record_limit = std::min(meta_capacity, prefix_capacity);
    size_t pos = begin;
    size_t raw = 0;
    while (pos < end && pos - begin < record_limit) {
        if (records[pos].raw_len > raw_capacity - raw) break;
        raw += records[pos].raw_len;
        pos++;
    }
    *raw_bytes = raw;
    return pos;
}

int MpiSortStrictPackRange(const MpiSortRecordMeta *records,
                           size_t record_count,
                           const unsigned char *raw,
                           size_t raw_size,
                           size_t begin,
                           size_t end,
                           MpiSortStrictExchangeWorkspace *ws,
                           MpiSortStats *stats) {
    const size_t count = end - begin;
    uint64_t *prefix = (uint64_t *)ws->prefix_buffer.data;
    size_t raw_pos = 0;
    for (size_t i = 0; i < count; ++i) {
        prefix[i] = (uint64_t)raw_pos;
        raw_pos += records[begin + i].raw_len;
    }
    if (count > ws->meta_capacity || count > ws->prefix_capacity ||
        raw_pos > ws->raw_capacity) {
        return -1;
    }
    MpiSortRangePackPara paras[kSortCpeBlocks];
    for (int c = 0; c < kSortCpeBlocks; ++c) {
        size_t c_begin = begin + count * (size_t)c / kSortCpeBlocks;
        size_t c_end = begin + count * (size_t)(c + 1) / kSortCpeBlocks;
        paras[c].core_id = c;
        paras[c].local_records = records;
        paras[c].local_raw = raw;
        paras[c].local_raw_size = raw_size;
        paras[c].record_raw_offsets = prefix;
        paras[c].send_meta = (MpiSortRecordMeta *)ws->send_meta_buffer.data;
        paras[c].send_raw = ws->send_raw_buffer.data;
        paras[c].send_meta_capacity = ws->meta_capacity;
        paras[c].send_raw_capacity = ws->raw_capacity;
        paras[c].record_begin = c_begin;
        paras[c].record_end = c_end;
        paras[c].output_record_begin = c_begin - begin;
        paras[c].status = 0;
        paras[c].record_index = 0;
        paras[c].actual_value = 0;
        paras[c].limit_value = 0;
        paras[c].limit_id = BOUNDS_LIMIT_NONE;
        paras[c].pack_cycles = 0;
        paras[c].total_cycles = 0;
    }
    double t0 = GetTime();
    __real_athread_spawn((void *)slave_mpi_sort_range_pack, paras, 1);
    athread_join();
    stats->t_bucket_pack += GetTime() - t0;
    bool cpe_ok = true;
    for (int c = 0; c < kSortCpeBlocks; ++c) {
        if (paras[c].status != 0) cpe_ok = false;
    }
    MpiSortRecordMeta *packed =
        (MpiSortRecordMeta *)ws->send_meta_buffer.data;
    for (size_t i = 0; cpe_ok && i < count; ++i) {
        const MpiSortRecordMeta &src = records[begin + i];
        const MpiSortRecordMeta &dst = packed[i];
        const uint64_t offset = prefix[i];
        if (dst.tid != src.tid || dst.pos != src.pos ||
            dst.flag != src.flag || dst.raw_len != src.raw_len ||
            dst.global_order != src.global_order ||
            dst.raw_offset != offset ||
            offset > raw_pos || dst.raw_len > raw_pos - offset ||
            (uint64_t)MpiSortReadLe32(
                ws->send_raw_buffer.data + offset) + 4u != dst.raw_len) {
            cpe_ok = false;
        }
    }
    if (!cpe_ok) {
        for (size_t i = 0; i < count; ++i) {
            const MpiSortRecordMeta &src = records[begin + i];
            const uint64_t offset = prefix[i];
            if (src.raw_offset > raw_size ||
                src.raw_len > raw_size - (size_t)src.raw_offset ||
                offset > ws->raw_capacity ||
                src.raw_len > ws->raw_capacity - (size_t)offset) {
                return -1;
            }
            packed[i] = src;
            packed[i].raw_offset = offset;
            memcpy(ws->send_raw_buffer.data + offset,
                   raw + src.raw_offset, src.raw_len);
        }
    }
    (void)record_count;
    return 0;
}

int MpiSortStrictAppendReceivedChunk(MpiSortStrictTempStore *store,
                                     const MpiSortRecordMeta *records,
                                     size_t record_count,
                                     const unsigned char *raw,
                                     size_t raw_size,
                                     unsigned char *sim_scratch,
                                     size_t sim_scratch_size,
                                     MpiSortStats *stats,
                                     MpiSortStrictSegment *segment) {
    MpiSortStrictExtent extent;
    if (MpiSortStrictAppendExtent(store->segments_fd, &store->segments_end,
                                  records, record_count, raw, raw_size,
                                  sim_scratch, sim_scratch_size,
                                  stats, &extent) != 0) {
        return -1;
    }
    segment->extents.push_back(extent);
    segment->record_count += record_count;
    segment->raw_size += raw_size;
    return 0;
}

int MpiSortStrictExchangeRuns(const std::vector<MpiSortStrictRun> &runs,
                              const std::vector<MpiSortKey> &splitters,
                              MpiSortStrictRunArena *arena,
                              MpiSortStrictExchangeWorkspace *ws,
                              MpiSortStrictTempStore *store,
                              unsigned char *sim_scratch,
                              size_t sim_scratch_size,
                              int rank,
                              int comm_size,
                              std::vector<MpiSortStrictSegment> *segments,
                              MpiSortStats *stats) {
    long long local_runs = (long long)runs.size();
    long long max_runs = 0;
    MPI_Allreduce(&local_runs, &max_runs, 1, MPI_LONG_LONG, MPI_MAX,
                  MPI_COMM_WORLD);
    for (long long round = 0; round < max_runs; ++round) {
        int load_ok = 1;
        if (round < local_runs) {
            const MpiSortStrictRun &run = runs[(size_t)round];
            if (run.memory_backed) {
                if (round != 0 ||
                    arena->record_count != (size_t)run.extent.record_count ||
                    arena->raw_used != (size_t)run.extent.raw_size) {
                    load_ok = 0;
                }
            } else if (MpiSortStrictLoadRun(
                           run, arena, store, sim_scratch,
                           sim_scratch_size, stats) != 0) {
                load_ok = 0;
            }
        } else {
            arena->reset();
        }
        int global_load_ok = 0;
        MPI_Allreduce(&load_ok, &global_load_ok, 1, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD);
        if (!global_load_ok) {
            return -1;
        }
        MpiSortRecordMeta *records = arena->records();
        size_t record_count = arena->record_count;
        std::vector<size_t> bounds((size_t)comm_size + 1, record_count);
        std::vector<size_t> bucket_raw_bytes((size_t)comm_size, 0);
        bounds[0] = 0;
        double count_t0 = GetTime();
        size_t pos = 0;
        for (int bucket = 0; bucket < comm_size - 1; ++bucket) {
            while (pos < record_count &&
                   MpiSortCompareMetaKey(records[pos], splitters[(size_t)bucket]) <= 0) {
                pos++;
            }
            bounds[(size_t)bucket + 1] = pos;
        }
        bounds[(size_t)comm_size] = record_count;
        double count_wall = GetTime() - count_t0;
        stats->t_bucket_count += count_wall;

        for (int bucket = 0; bucket < comm_size; ++bucket) {
            size_t braw = 0;
            for (size_t i = bounds[(size_t)bucket];
                 i < bounds[(size_t)bucket + 1]; ++i) {
                braw += records[i].raw_len;
            }
            bucket_raw_bytes[(size_t)bucket] = braw;
            long long n = (long long)(bounds[(size_t)bucket + 1] -
                                      bounds[(size_t)bucket]);
            if (bucket == rank) {
                stats->bucket_self_records += n;
                stats->bucket_self_raw_bytes += (long long)braw;
            } else {
                stats->bucket_remote_records += n;
                stats->bucket_remote_raw_bytes += (long long)braw;
            }
        }

        double exchange_round_t0 = GetTime();
        int round_io_ok = 1;
        for (int step = 0; step < comm_size; ++step) {
            int send_to = (rank + step) % comm_size;
            int recv_from = (rank - step + comm_size) % comm_size;
            size_t send_pos = bounds[(size_t)send_to];
            size_t send_end = bounds[(size_t)send_to + 1];
            MpiSortStrictSegment received_segment;

            if (step == 0) {
                while (send_pos < send_end) {
                    size_t send_raw_bytes = 0;
                    size_t send_chunk_end = MpiSortStrictChunkEnd(
                        records, send_pos, send_end, ws->meta_capacity,
                        ws->prefix_capacity, ws->raw_capacity,
                        &send_raw_bytes);
                    if (send_chunk_end == send_pos ||
                        MpiSortStrictPackRange(
                            records, record_count, arena->data,
                            arena->raw_used, send_pos, send_chunk_end,
                            ws, stats) != 0 ||
                        MpiSortStrictAppendReceivedChunk(
                            store,
                            (MpiSortRecordMeta *)ws->send_meta_buffer.data,
                            send_chunk_end - send_pos,
                            ws->send_raw_buffer.data, send_raw_bytes,
                            sim_scratch, sim_scratch_size,
                            stats, &received_segment) != 0) {
                        round_io_ok = 0;
                        break;
                    }
                    send_pos = send_chunk_end;
                }
            } else {
                bool send_done = send_pos >= send_end;
                bool recv_done = false;
                while (!send_done || !recv_done) {
                    size_t send_raw_bytes = 0;
                    size_t send_chunk_end = send_pos;
                    if (!send_done) {
                        send_chunk_end = MpiSortStrictChunkEnd(
                            records, send_pos, send_end, ws->meta_capacity,
                            ws->prefix_capacity, ws->raw_capacity,
                            &send_raw_bytes);
                        if (send_chunk_end == send_pos ||
                            MpiSortStrictPackRange(records, record_count,
                                                   arena->data, arena->raw_used,
                                                   send_pos, send_chunk_end,
                                                   ws, stats) != 0) {
                            return -1;
                        }
                    }
                    unsigned long long send_header[3] = {
                        (unsigned long long)(send_chunk_end - send_pos),
                        (unsigned long long)send_raw_bytes,
                        send_chunk_end >= send_end ? 1ull : 0ull
                    };
                    unsigned long long recv_header[3] = {};
                    MPI_Status mpi_status;
                    double mpi_actual_t0 = GetTime();
                    MPI_Sendrecv(send_header, 3, MPI_UNSIGNED_LONG_LONG,
                                 send_to, kSortStrictTagHeader,
                                 recv_header, 3, MPI_UNSIGNED_LONG_LONG,
                                 recv_from, kSortStrictTagHeader,
                                 MPI_COMM_WORLD, &mpi_status);
                    size_t recv_count = (size_t)recv_header[0];
                    size_t recv_raw_bytes = (size_t)recv_header[1];
                    if (recv_count > ws->meta_capacity ||
                        recv_raw_bytes > ws->raw_capacity) {
                        return -1;
                    }
                    if (MpiSortSendRecvBytes(
                            send_to, (char *)ws->send_meta_buffer.data,
                            (long long)(send_header[0] * sizeof(MpiSortRecordMeta)),
                            recv_from, (char *)ws->recv_meta_buffer.data,
                            (long long)(recv_count * sizeof(MpiSortRecordMeta)),
                            kSortStrictTagMeta) != 0 ||
                        MpiSortSendRecvBytes(
                            send_to, (char *)ws->send_raw_buffer.data,
                            (long long)send_raw_bytes,
                            recv_from, (char *)ws->recv_raw_buffer.data,
                            (long long)recv_raw_bytes,
                            kSortStrictTagRaw) != 0) {
                        return -1;
                    }
                    stats->t_mpi_exchange += GetTime() - mpi_actual_t0;

                    // Repeat the same transfer before either side writes its
                    // received chunk. This isolates MPI cost from temp-file
                    // stalls that otherwise propagate to the peer.
                    unsigned long long sim_recv_header[3] = {};
                    double mpi_sim_t0 = GetTime();
                    MPI_Sendrecv(send_header, 3, MPI_UNSIGNED_LONG_LONG,
                                 send_to, kSortStrictTagSimHeader,
                                 sim_recv_header, 3, MPI_UNSIGNED_LONG_LONG,
                                 recv_from, kSortStrictTagSimHeader,
                                 MPI_COMM_WORLD, &mpi_status);
                    if (sim_recv_header[0] != recv_header[0] ||
                        sim_recv_header[1] != recv_header[1] ||
                        sim_recv_header[2] != recv_header[2]) {
                        return -1;
                    }
                    if (MpiSortSendRecvBytes(
                            send_to, (char *)ws->send_meta_buffer.data,
                            (long long)(send_header[0] *
                                        sizeof(MpiSortRecordMeta)),
                            recv_from, (char *)ws->recv_meta_buffer.data,
                            (long long)(recv_count *
                                        sizeof(MpiSortRecordMeta)),
                            kSortStrictTagSimMeta) != 0 ||
                        MpiSortSendRecvBytes(
                            send_to, (char *)ws->send_raw_buffer.data,
                            (long long)send_raw_bytes,
                            recv_from, (char *)ws->recv_raw_buffer.data,
                            (long long)recv_raw_bytes,
                            kSortStrictTagSimRaw) != 0) {
                        return -1;
                    }
                    stats->t_mpi_simulated += GetTime() - mpi_sim_t0;
                    if (round_io_ok && recv_count > 0 &&
                        MpiSortStrictAppendReceivedChunk(
                            store, (MpiSortRecordMeta *)ws->recv_meta_buffer.data,
                            recv_count, ws->recv_raw_buffer.data,
                            recv_raw_bytes, sim_scratch, sim_scratch_size,
                            stats, &received_segment) != 0) {
                        round_io_ok = 0;
                    }
                    send_pos = send_chunk_end;
                    send_done = send_header[2] != 0;
                    recv_done = recv_header[2] != 0;
                }
            }
            if (round_io_ok && received_segment.record_count > 0) {
                stats->received_records +=
                    (long long)received_segment.record_count;
                segments->push_back(std::move(received_segment));
            }
        }
        double exchange_wall = GetTime() - exchange_round_t0;
        stats->t_exchange += exchange_wall;
        int global_round_io_ok = 0;
        MPI_Allreduce(&round_io_ok, &global_round_io_ok, 1, MPI_INT,
                      MPI_MIN, MPI_COMM_WORLD);
        if (!global_round_io_ok) return -1;
    }
    stats->t_partition = stats->t_bucket_count + stats->t_bucket_pack;
    stats->t_run_bucket = stats->t_partition;
    stats->t_run_exchange = stats->t_exchange;
    return 0;
}

struct MpiSortStrictCursor {
    const MpiSortStrictSegment *segment;
    size_t extent_index;
    uint64_t extent_record_pos;
    MpiSortRecordMeta *meta_buffer;
    size_t meta_capacity;
    size_t meta_count;
    size_t meta_pos;
    unsigned char *raw_buffer;
    size_t raw_capacity;
    uint64_t raw_buffer_begin;
    size_t raw_buffer_size;
    bool active;

    MpiSortStrictCursor()
        : segment(nullptr), extent_index(0), extent_record_pos(0),
          meta_buffer(nullptr), meta_capacity(0), meta_count(0), meta_pos(0),
          raw_buffer(nullptr), raw_capacity(0), raw_buffer_begin(0),
          raw_buffer_size(0), active(false) {}
};

int MpiSortStrictCursorLoadMeta(MpiSortStrictCursor *cursor,
                                MpiSortStrictTempStore *store,
                                unsigned char *sim_scratch,
                                size_t sim_scratch_size,
                                MpiSortStats *stats) {
    cursor->active = false;
    while (cursor->extent_index < cursor->segment->extents.size()) {
        const MpiSortStrictExtent &extent =
            cursor->segment->extents[cursor->extent_index];
        if (cursor->extent_record_pos >= extent.record_count) {
            cursor->extent_index++;
            cursor->extent_record_pos = 0;
            cursor->raw_buffer_size = 0;
            continue;
        }
        size_t count = (size_t)std::min<uint64_t>(
            extent.record_count - cursor->extent_record_pos,
            cursor->meta_capacity);
        size_t bytes = count * sizeof(MpiSortRecordMeta);
        off_t offset = extent.meta_offset +
            (off_t)(cursor->extent_record_pos * sizeof(MpiSortRecordMeta));
        if (MpiSortStrictPreadAll(extent.fd, cursor->meta_buffer,
                                  bytes, offset, stats) != 0) {
            return -1;
        }
        MpiSortStrictSimulateCopy(cursor->meta_buffer, bytes, sim_scratch,
                                  sim_scratch_size, &stats->t_temp_read_sim);
        cursor->meta_count = count;
        cursor->meta_pos = 0;
        cursor->active = true;
        (void)store;
        return 0;
    }
    (void)store;
    return 0;
}

int MpiSortStrictCursorInit(MpiSortStrictCursor *cursor,
                            const MpiSortStrictSegment *segment,
                            MpiSortRecordMeta *meta_buffer,
                            size_t meta_capacity,
                            unsigned char *raw_buffer,
                            size_t raw_capacity,
                            MpiSortStrictTempStore *store,
                            unsigned char *sim_scratch,
                            size_t sim_scratch_size,
                            MpiSortStats *stats) {
    cursor->segment = segment;
    cursor->meta_buffer = meta_buffer;
    cursor->meta_capacity = meta_capacity;
    cursor->raw_buffer = raw_buffer;
    cursor->raw_capacity = raw_capacity;
    return MpiSortStrictCursorLoadMeta(cursor, store, sim_scratch,
                                       sim_scratch_size, stats);
}

const MpiSortRecordMeta &MpiSortStrictCursorMeta(
        const MpiSortStrictCursor &cursor) {
    return cursor.meta_buffer[cursor.meta_pos];
}

int MpiSortStrictCursorRaw(MpiSortStrictCursor *cursor,
                           MpiSortStrictTempStore *store,
                           unsigned char *sim_scratch,
                           size_t sim_scratch_size,
                           MpiSortStats *stats,
                           const unsigned char **raw_ptr) {
    const MpiSortRecordMeta &meta = MpiSortStrictCursorMeta(*cursor);
    const MpiSortStrictExtent &extent =
        cursor->segment->extents[cursor->extent_index];
    uint64_t raw_end = meta.raw_offset + meta.raw_len;
    if (raw_end > extent.raw_size || meta.raw_len > cursor->raw_capacity) {
        return -1;
    }
    if (cursor->raw_buffer_size == 0 ||
        meta.raw_offset < cursor->raw_buffer_begin ||
        raw_end > cursor->raw_buffer_begin + cursor->raw_buffer_size) {
        cursor->raw_buffer_begin = meta.raw_offset;
        cursor->raw_buffer_size = (size_t)std::min<uint64_t>(
            cursor->raw_capacity, extent.raw_size - meta.raw_offset);
        if (MpiSortStrictPreadAll(
                extent.fd, cursor->raw_buffer,
                cursor->raw_buffer_size,
                extent.raw_offset + (off_t)meta.raw_offset, stats) != 0) {
            return -1;
        }
        MpiSortStrictSimulateCopy(cursor->raw_buffer,
                                  cursor->raw_buffer_size,
                                  sim_scratch, sim_scratch_size,
                                  &stats->t_temp_read_sim);
    }
    *raw_ptr = cursor->raw_buffer +
        (size_t)(meta.raw_offset - cursor->raw_buffer_begin);
    if (meta.raw_len < 36 ||
        (uint64_t)MpiSortReadLe32(*raw_ptr) + 4u != meta.raw_len) {
        fprintf(stderr,
                "ERROR: external sort temp record length mismatch. "
                "extent=%zu record=%llu raw_offset=%llu raw_len=%u\n",
                cursor->extent_index,
                (unsigned long long)cursor->extent_record_pos,
                (unsigned long long)meta.raw_offset, meta.raw_len);
        return -1;
    }
    (void)store;
    return 0;
}

int MpiSortStrictCursorAdvance(MpiSortStrictCursor *cursor,
                               MpiSortStrictTempStore *store,
                               unsigned char *sim_scratch,
                               size_t sim_scratch_size,
                               MpiSortStats *stats) {
    cursor->meta_pos++;
    cursor->extent_record_pos++;
    if (cursor->meta_pos < cursor->meta_count) return 0;
    return MpiSortStrictCursorLoadMeta(cursor, store, sim_scratch,
                                       sim_scratch_size, stats);
}

class MpiSortStrictLoserTree {
public:
    explicit MpiSortStrictLoserTree(std::vector<MpiSortStrictCursor> *cursors)
        : cursors_(cursors), losers_(cursors->size(), (int)cursors->size()),
          sentinel_((int)cursors->size()) {
        for (int i = sentinel_ - 1; i >= 0; --i) Adjust(i);
    }

    bool empty() const {
        int winner = losers_.empty() ? sentinel_ : losers_[0];
        return winner == sentinel_ || !Active(winner);
    }

    int winner_index() const {
        return losers_[0];
    }

    void replay(int player) {
        Adjust(player);
    }

private:
    bool Active(int player) const {
        return player >= 0 && player < sentinel_ &&
               (*cursors_)[(size_t)player].active;
    }

    bool Greater(int lhs, int rhs) const {
        if (lhs == sentinel_) return false;
        if (rhs == sentinel_) return true;
        bool lhs_active = Active(lhs);
        bool rhs_active = Active(rhs);
        if (lhs_active != rhs_active) return !lhs_active;
        if (!lhs_active) return lhs > rhs;
        MpiSortMetaLess less;
        const MpiSortRecordMeta &a =
            MpiSortStrictCursorMeta((*cursors_)[(size_t)lhs]);
        const MpiSortRecordMeta &b =
            MpiSortStrictCursorMeta((*cursors_)[(size_t)rhs]);
        if (less(b, a)) return true;
        if (less(a, b)) return false;
        return lhs > rhs;
    }

    void Adjust(int player) {
        int parent = (player + sentinel_) >> 1;
        while (parent > 0) {
            if (Greater(player, losers_[(size_t)parent])) {
                std::swap(player, losers_[(size_t)parent]);
            }
            parent >>= 1;
        }
        if (!losers_.empty()) losers_[0] = player;
    }

    std::vector<MpiSortStrictCursor> *cursors_;
    std::vector<int> losers_;
    int sentinel_;
};

size_t MpiSortStrictCompressWorkspaceBytes() {
    return 4ull * kSortCpeBlocks *
           (sizeof(bam_block) + (size_t)BGZF_MAX_BLOCK_SIZE);
}

int MpiSortStrictBuildCursors(
        const std::vector<MpiSortStrictSegment> &segments,
        size_t begin,
        size_t end,
        MpiSortStrictBuffer *cursor_memory,
        MpiSortStrictTempStore *store,
        unsigned char *sim_scratch,
        size_t sim_scratch_size,
        MpiSortStats *stats,
        std::vector<MpiSortStrictCursor> *cursors) {
    const size_t meta_capacity =
        kSortStrictMetaBuffer / sizeof(MpiSortRecordMeta);
    const size_t stride = kSortStrictMetaBuffer + kSortStrictRawBuffer;
    cursors->clear();
    cursors->reserve(end - begin);
    for (size_t i = begin; i < end; ++i) {
        if (segments[i].record_count == 0) continue;
        size_t slot = cursors->size();
        unsigned char *base = cursor_memory->data + slot * stride;
        MpiSortStrictCursor cursor;
        if (MpiSortStrictCursorInit(
                &cursor, &segments[i],
                (MpiSortRecordMeta *)base, meta_capacity,
                base + kSortStrictMetaBuffer, kSortStrictRawBuffer,
                store, sim_scratch, sim_scratch_size, stats) != 0) {
            return -1;
        }
        cursors->push_back(cursor);
    }
    return 0;
}

int MpiSortStrictMergeGroupToTemp(
        const std::vector<MpiSortStrictSegment> &segments,
        size_t begin,
        size_t end,
        MpiSortStrictBuffer *cursor_memory,
        MpiSortStrictBuffer *output_meta_memory,
        MpiSortStrictBuffer *output_raw_memory,
        MpiSortStrictTempStore *store,
        unsigned char *sim_scratch,
        size_t sim_scratch_size,
        MpiSortStats *stats,
        MpiSortStrictSegment *output) {
    std::vector<MpiSortStrictCursor> cursors;
    if (MpiSortStrictBuildCursors(segments, begin, end, cursor_memory,
                                  store, sim_scratch, sim_scratch_size,
                                  stats, &cursors) != 0) {
        return -1;
    }
    MpiSortStrictLoserTree tree(&cursors);
    MpiSortRecordMeta *out_meta =
        (MpiSortRecordMeta *)output_meta_memory->data;
    size_t meta_capacity =
        output_meta_memory->size / sizeof(MpiSortRecordMeta);
    size_t meta_count = 0;
    size_t raw_used = 0;
    auto flush = [&]() -> int {
        if (meta_count == 0) return 0;
        int ret = MpiSortStrictAppendReceivedChunk(
            store, out_meta, meta_count, output_raw_memory->data, raw_used,
            sim_scratch, sim_scratch_size, stats, output);
        meta_count = 0;
        raw_used = 0;
        return ret;
    };
    while (!tree.empty()) {
        int winner = tree.winner_index();
        MpiSortStrictCursor &cursor = cursors[(size_t)winner];
        const MpiSortRecordMeta meta = MpiSortStrictCursorMeta(cursor);
        const unsigned char *raw_ptr = nullptr;
        if (MpiSortStrictCursorRaw(&cursor, store, sim_scratch,
                                   sim_scratch_size, stats, &raw_ptr) != 0) {
            return -1;
        }
        if (meta_count == meta_capacity ||
            meta.raw_len > output_raw_memory->size - raw_used) {
            if (flush() != 0) return -1;
        }
        if (meta.raw_len > output_raw_memory->size) return -1;
        out_meta[meta_count] = meta;
        out_meta[meta_count].raw_offset = raw_used;
        memcpy(output_raw_memory->data + raw_used, raw_ptr, meta.raw_len);
        raw_used += meta.raw_len;
        meta_count++;
        if (MpiSortStrictCursorAdvance(&cursor, store, sim_scratch,
                                       sim_scratch_size, stats) != 0) {
            return -1;
        }
        tree.replay(winner);
    }
    return flush();
}

int MpiSortStrictConsolidate(
        std::vector<MpiSortStrictSegment> *segments,
        size_t fan_in,
        MpiSortMemoryTracker *tracker,
        MpiSortStrictTempStore *store,
        unsigned char *sim_scratch,
        size_t sim_scratch_size,
        MpiSortStats *stats) {
    if (segments->size() <= fan_in) return 0;
    const size_t cursor_bytes =
        fan_in * (kSortStrictMetaBuffer + kSortStrictRawBuffer);
    MpiSortStrictBuffer cursor_memory;
    MpiSortStrictBuffer output_meta;
    MpiSortStrictBuffer output_raw;
    if (MpiSortStrictAlloc(tracker, cursor_bytes, &cursor_memory) != 0 ||
        MpiSortStrictAlloc(tracker, kSortStrictConsolidateMeta,
                           &output_meta) != 0 ||
        MpiSortStrictAlloc(tracker, kSortStrictConsolidateRaw,
                           &output_raw) != 0) {
        MpiSortStrictFree(tracker, &cursor_memory);
        MpiSortStrictFree(tracker, &output_meta);
        MpiSortStrictFree(tracker, &output_raw);
        return -1;
    }
    while (segments->size() > fan_in) {
        std::vector<MpiSortStrictSegment> next;
        next.reserve((segments->size() + fan_in - 1) / fan_in);
        double read_actual_before = stats->t_temp_read_actual;
        double write_actual_before = stats->t_temp_write_actual;
        double read_sim_before = stats->t_temp_read_sim;
        double write_sim_before = stats->t_temp_write_sim;
        double t0 = GetTime();
        for (size_t begin = 0; begin < segments->size(); begin += fan_in) {
            size_t end = std::min(segments->size(), begin + fan_in);
            MpiSortStrictSegment merged;
            if (MpiSortStrictMergeGroupToTemp(
                    *segments, begin, end, &cursor_memory,
                    &output_meta, &output_raw, store,
                    sim_scratch, sim_scratch_size, stats, &merged) != 0) {
                MpiSortStrictFree(tracker, &cursor_memory);
                MpiSortStrictFree(tracker, &output_meta);
                MpiSortStrictFree(tracker, &output_raw);
                return -1;
            }
            next.push_back(std::move(merged));
        }
        double wall = GetTime() - t0;
        double actual_io =
            (stats->t_temp_read_actual - read_actual_before) +
            (stats->t_temp_write_actual - write_actual_before);
        double simulated_wall = wall - actual_io;
        if (simulated_wall < 0.0) simulated_wall = 0.0;
        stats->t_run_merge += wall;
        stats->t_consolidation_simulated += simulated_wall;
        stats->t_consolidation_temp_read_sim +=
            stats->t_temp_read_sim - read_sim_before;
        stats->t_consolidation_temp_write_sim +=
            stats->t_temp_write_sim - write_sim_before;
        stats->consolidation_passes++;
        segments->swap(next);
    }
    MpiSortStrictFree(tracker, &cursor_memory);
    MpiSortStrictFree(tracker, &output_meta);
    MpiSortStrictFree(tracker, &output_raw);
    return 0;
}

int MpiSortStrictPreparePayloadBatch(
        MpiSortStrictLoserTree *tree,
        std::vector<MpiSortStrictCursor> *cursors,
        MpiSortStrictTempStore *store,
        unsigned char *sim_scratch,
        size_t sim_scratch_size,
        MpiSortBlockSet *payload_blocks,
        MpiSortBlockSet *output_blocks,
        MpiSortRawCompressPara *paras,
        int compress_level,
        MpiSortStats *stats,
        int *active_blocks) {
    int active = 0;
    while (!tree->empty() && active < kSortCpeBlocks) {
        bam_block *payload = &payload_blocks->blocks[active];
        payload->pos = 0;
        payload->length = 0;
        payload->errcode = 0;
        payload->block_id = active;
        while (!tree->empty()) {
            int winner = tree->winner_index();
            MpiSortStrictCursor &cursor = (*cursors)[(size_t)winner];
            const MpiSortRecordMeta meta = MpiSortStrictCursorMeta(cursor);
            if (meta.raw_len > BGZF_BLOCK_SIZE) return -1;
            if (payload->pos > 0 &&
                payload->pos + meta.raw_len > BGZF_BLOCK_SIZE) {
                break;
            }
            const unsigned char *raw_ptr = nullptr;
            if (MpiSortStrictCursorRaw(&cursor, store, sim_scratch,
                                       sim_scratch_size, stats,
                                       &raw_ptr) != 0) {
                return -1;
            }
            memcpy(payload->data + payload->pos, raw_ptr, meta.raw_len);
            payload->pos += meta.raw_len;
            payload->length = payload->pos;
            if (MpiSortStrictCursorAdvance(&cursor, store, sim_scratch,
                                           sim_scratch_size, stats) != 0) {
                return -1;
            }
            tree->replay(winner);
        }
        if (payload->pos == 0) return -1;
#if RABBITBAM_SORT_ENABLE_SANITY_CHECKS
        if (MpiSortValidatePayloadRecords(payload, "external-strict") != 0) {
            return -1;
        }
#endif
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

int MpiSortStrictMergeCompress(
        const std::vector<MpiSortStrictSegment> &segments,
        MpiSortMemoryTracker *tracker,
        MpiSortStrictTempStore *store,
        unsigned char *sim_scratch,
        size_t sim_scratch_size,
        swbam::RankBodySink &body_sink,
        int compress_level,
        MpiSortStats *stats) {
    const size_t cursor_stride =
        kSortStrictMetaBuffer + kSortStrictRawBuffer;
    const size_t cursor_bytes = segments.size() * cursor_stride;
    const size_t compress_bytes = MpiSortStrictCompressWorkspaceBytes();
    MpiSortStrictBuffer cursor_memory;
    bool compress_acquired = tracker->acquire(compress_bytes);
    if (!compress_acquired ||
        MpiSortStrictAlloc(tracker, cursor_bytes, &cursor_memory) != 0) {
        if (compress_acquired) tracker->release(compress_bytes);
        return -1;
    }
    MpiSortBlockSet payload_a = {}, payload_b = {}, out_a = {}, out_b = {};
    MpiSortRawCompressPara paras_a[kSortCpeBlocks], paras_b[kSortCpeBlocks];
    if (MpiSortAllocateBlockSet(&payload_a, kSortCpeBlocks) != 0 ||
        MpiSortAllocateBlockSet(&payload_b, kSortCpeBlocks) != 0 ||
        MpiSortAllocateBlockSet(&out_a, kSortCpeBlocks) != 0 ||
        MpiSortAllocateBlockSet(&out_b, kSortCpeBlocks) != 0) {
        MpiSortStrictFree(tracker, &cursor_memory);
        tracker->release(compress_bytes);
        return -1;
    }
    for (int i = 0; i < kSortCpeBlocks; ++i) {
        MpiSortInitEmptyRawCompressPara(&paras_a[i], i);
        MpiSortInitEmptyRawCompressPara(&paras_b[i], i);
    }
    std::vector<MpiSortStrictCursor> cursors;
    double merge_read_actual_before = stats->t_temp_read_actual;
    double merge_read_sim_before = stats->t_temp_read_sim;
    double merge_t0 = GetTime();
    if (MpiSortStrictBuildCursors(segments, 0, segments.size(),
                                  &cursor_memory, store,
                                  sim_scratch, sim_scratch_size,
                                  stats, &cursors) != 0) {
        MpiSortFreeBlockSet(&payload_a);
        MpiSortFreeBlockSet(&payload_b);
        MpiSortFreeBlockSet(&out_a);
        MpiSortFreeBlockSet(&out_b);
        MpiSortStrictFree(tracker, &cursor_memory);
        tracker->release(compress_bytes);
        return -1;
    }
    double merge_wall = GetTime() - merge_t0;
    double merge_actual_io =
        stats->t_temp_read_actual - merge_read_actual_before;
    double merge_simulated = merge_wall - merge_actual_io;
    if (merge_simulated < 0.0) merge_simulated = 0.0;
    stats->t_final_sort += merge_wall;
    stats->t_run_merge += merge_wall;
    stats->t_merge_simulated += merge_simulated;
    stats->t_merge_temp_read_sim +=
        stats->t_temp_read_sim - merge_read_sim_before;
    stats->t_merge_unhidden += merge_simulated;
    MpiSortStrictLoserTree tree(&cursors);
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
            if (paras_pending[i].status == 0 &&
                paras_pending[i].output_block) {
                if (swbam::AppendBgzfBlock(
                        &body_sink,
                        paras_pending[i].output_block) != 0) {
                    return -1;
                }
                stats->bgzf_blocks++;
            }
            MpiSortInitEmptyRawCompressPara(&paras_pending[i], i);
        }
        has_pending = false;
        stats->t_write += GetTime() - t0;
        return 0;
    };

    int active_blocks = 0;
    merge_read_actual_before = stats->t_temp_read_actual;
    merge_read_sim_before = stats->t_temp_read_sim;
    merge_t0 = GetTime();
    if (MpiSortStrictPreparePayloadBatch(
            &tree, &cursors, store, sim_scratch, sim_scratch_size,
            payload_active, out_active, paras_active,
            compress_level, stats, &active_blocks) != 0) {
        MpiSortFreeBlockSet(&payload_a);
        MpiSortFreeBlockSet(&payload_b);
        MpiSortFreeBlockSet(&out_a);
        MpiSortFreeBlockSet(&out_b);
        MpiSortStrictFree(tracker, &cursor_memory);
        tracker->release(compress_bytes);
        return -1;
    }
    merge_wall = GetTime() - merge_t0;
    merge_actual_io = stats->t_temp_read_actual - merge_read_actual_before;
    merge_simulated = merge_wall - merge_actual_io;
    if (merge_simulated < 0.0) merge_simulated = 0.0;
    stats->t_final_sort += merge_wall;
    stats->t_run_merge += merge_wall;
    stats->t_merge_simulated += merge_simulated;
    stats->t_merge_temp_read_sim +=
        stats->t_temp_read_sim - merge_read_sim_before;
    stats->t_merge_unhidden += merge_simulated;

    int ret = 0;
    while (active_blocks > 0) {
        int next_active = 0;
        double t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_sort_compress_payload,
                             paras_active, 1);
        if (flush_pending() != 0) {
            athread_join();
            ret = -1;
            break;
        }
        merge_read_actual_before = stats->t_temp_read_actual;
        merge_read_sim_before = stats->t_temp_read_sim;
        merge_t0 = GetTime();
        if (MpiSortStrictPreparePayloadBatch(
                &tree, &cursors, store, sim_scratch, sim_scratch_size,
                payload_pending, out_pending, paras_pending,
                compress_level, stats, &next_active) != 0) {
            athread_join();
            ret = -1;
            break;
        }
        merge_wall = GetTime() - merge_t0;
        merge_actual_io =
            stats->t_temp_read_actual - merge_read_actual_before;
        merge_simulated = merge_wall - merge_actual_io;
        if (merge_simulated < 0.0) merge_simulated = 0.0;
        stats->t_final_sort += merge_wall;
        stats->t_run_merge += merge_wall;
        stats->t_merge_simulated += merge_simulated;
        stats->t_merge_temp_read_sim +=
            stats->t_temp_read_sim - merge_read_sim_before;
        athread_join();
        double wall = GetTime() - t0;
        stats->t_compress += wall;
        double compress_simulated = MpiSortEstimateCompressSeconds(
            paras_active, active_blocks, wall, stats);
        stats->t_compress_simulated += compress_simulated;
        MpiSortAccumulateCompressDetail(paras_active, active_blocks,
                                        compress_simulated, stats);
        for (int i = 0; i < active_blocks; ++i) {
            if (paras_active[i].status != 0) {
                ret = -1;
                break;
            }
        }
        if (ret != 0) break;
        has_pending = active_blocks > 0;
        std::swap(payload_active, payload_pending);
        std::swap(out_active, out_pending);
        std::swap(paras_active, paras_pending);
        active_blocks = next_active;
    }
    if (ret == 0) ret = flush_pending();
    MpiSortFreeBlockSet(&payload_a);
    MpiSortFreeBlockSet(&payload_b);
    MpiSortFreeBlockSet(&out_a);
    MpiSortFreeBlockSet(&out_b);
    MpiSortStrictFree(tracker, &cursor_memory);
    tracker->release(compress_bytes);
    return ret;
}

} // namespace


static int FusedBamSortMPIImpl(
        MemReader *reader,
        const swbam::BamInputBackend *input,
        const swbam::BgzfBlockSpan *spans,
        size_t span_count,
        size_t compressed_input_size,
        swbam::RankBodySink &body_sink,
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
    std::vector<long long> received_source_counts;

    if (compressed_input_size > 0) {
        const size_t raw_reserve = compressed_input_size > SIZE_MAX / 2
            ? SIZE_MAX : compressed_input_size * 2;
        if (raw_reserve == SIZE_MAX ||
            MpiSortRawReserve(&local_raw, raw_reserve) != 0) {
            fprintf(stderr, "ERROR: MPI sort failed to reserve local raw buffer.\n");
            return -1;
        }
        local_records.reserve(compressed_input_size / 64);
    }
    stats->t_setup += GetTime() - setup_t0;

    //2. 把本 rank 负责的 BGZF blocks 解析成：local_records 和 local_raw
    // local_records：record metadata，包括 tid/pos/flag/global_order/raw_offset/raw_len
    // local_raw：真实 BAM record bytes
#if RABBITBAM_SORT_ENABLE_SANITY_CHECKS
    const size_t extract_reader_pos = reader ? reader->pos : 0;
#endif
    const int extract_ret = reader
        ? MpiSortExtractLocalRecordsCpe(
              *reader, global_block_begin,
              &local_records, &local_raw, stats)
        : MpiSortExtractLocalRecordsCpe(
              *input, spans, span_count, global_block_begin,
              &local_records, &local_raw, stats);
    if (extract_ret != 0) {
        return -1;
    }
#if RABBITBAM_SORT_ENABLE_SANITY_CHECKS
    if (MpiSortValidateMetaRawBuffer(local_records, local_raw.data,
                                     local_raw.size, nullptr,
                                     "memory-post-extract", false) != 0) {
        fprintf(stderr,
                "[rank %d] WARNING: MPI sort CPE extract produced inconsistent raw metadata; retrying extract on MPE host path.\n",
                rank);
        stats->local_records -= (long long)local_records.size();
        local_records.clear();
        MpiSortRawRelease(&local_raw);
        if (!reader) return -1;
        reader->pos = extract_reader_pos;
        if (MpiSortExtractLocalRecordsHost(*reader, global_block_begin,
                                           &local_records, &local_raw,
                                           stats) != 0) {
            return -1;
        }
        if (MpiSortValidateMetaRawBuffer(local_records, local_raw.data,
                                         local_raw.size, nullptr,
                                         "memory-post-extract-host") != 0) {
            return -1;
        }
    }
#endif
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
                               &received_source_counts,
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
    if (MpiSortKWayMergeReceivedRecords(received_source_counts, &received_records) != 0) {
        fprintf(stderr, "[rank %d] ERROR: MPI sort k-way merge failed.\n", rank);
        return -1;
    }
    stats->t_final_sort += GetTime() - t0;

    //7. 把排序后的 records 重新压缩成 BGZF body
    if (MpiSortCompressSortedRecordsCpe(received_records, received_raw,
                                        body_sink, compress_level, stats) != 0) {
        return -1;
    }
    stats->t_fused_total += GetTime() - fused_t0;
    (void)global_sample_count;
    return 0;
}

int FusedBamSortMPI(MemReader &reader,
                    swbam::RankBodySink &body_sink,
                    long long global_block_begin,
                    int rank,
                    int comm_size,
                    int compress_level,
                    size_t memory_limit,
                    MpiSortStats *stats) {
    return FusedBamSortMPIImpl(
        &reader, nullptr, nullptr, 0, reader.size,
        body_sink, global_block_begin, rank, comm_size,
        compress_level, memory_limit, stats);
}

int FusedBamSortMPI(const swbam::BamInputBackend &input,
                    const swbam::BgzfBlockSpan *spans,
                    size_t span_count,
                    swbam::RankBodySink &body_sink,
                    long long global_block_begin,
                    int rank,
                    int comm_size,
                    int compress_level,
                    size_t memory_limit,
                    MpiSortStats *stats) {
    size_t compressed_size = 0;
    for (size_t i = 0; i < span_count; ++i) {
        if (spans[i].compressed_size > SIZE_MAX - compressed_size) return -1;
        compressed_size += spans[i].compressed_size;
    }
    return FusedBamSortMPIImpl(
        nullptr, &input, spans, span_count, compressed_size,
        body_sink, global_block_begin, rank, comm_size,
        compress_level, memory_limit, stats);
}

static int FusedBamExternalSortMPIImpl(
        MemReader *reader,
        const swbam::BamInputBackend *input,
        const swbam::BgzfBlockSpan *spans,
        size_t span_count,
        swbam::RankBodySink &body_sink,
        long long global_block_begin,
        int rank,
        int comm_size,
        int compress_level,
        size_t memory_limit,
        const char *temp_prefix,
        MpiSortStats *stats) {
    double fused_t0 = GetTime();
    stats->sort_mode = 1;
    MpiSortMemoryTracker tracker(memory_limit);
    MpiSortStrictTempStore temp_store;
    MpiSortStrictBuffer sim_scratch;
    MpiSortStrictBuffer run_arena_memory;
    MpiSortStrictRunArena run_arena;
    MpiSortStrictExchangeWorkspace exchange_ws;
    std::vector<MpiSortStrictRun> runs;
    std::vector<MpiSortStrictSegment> segments;
    std::vector<MpiSortKey> local_samples;
    std::vector<MpiSortKey> splitters;
    int ret = -1;
    bool exchange_allocated = false;
    bool control_acquired = false;
    bool run_arena_allocated = false;
    auto all_ranks_ok = [](int local_ok) -> bool {
        int global_ok = 0;
        MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD);
        return global_ok != 0;
    };

    do {
        int setup_ok = 1;
        if (memory_limit == 0) {
            fprintf(stderr,
                    "[rank %d] ERROR: strict external sort requires a non-zero -m limit.\n",
                    rank);
            setup_ok = 0;
        }
        const size_t extract_workspace = MpiSortStrictExtractWorkspaceBytes();
        const size_t extract_minimum =
            kSortStrictControlReserve + kSortStrictSimScratch +
            extract_workspace + 2 * BGZF_MAX_BLOCK_SIZE +
            128 * sizeof(MpiSortRecordMeta);
        const size_t merge_minimum =
            kSortStrictControlReserve + kSortStrictSimScratch +
            MpiSortStrictCompressWorkspaceBytes() +
            kSortStrictConsolidateMeta + kSortStrictConsolidateRaw +
            2 * (kSortStrictMetaBuffer + kSortStrictRawBuffer);
        const size_t minimum_memory = std::max(extract_minimum, merge_minimum);
        if (setup_ok && memory_limit < minimum_memory) {
            fprintf(stderr,
                    "[rank %d] ERROR: -m is too small for strict external sort. minimum_memory=%zu limit=%zu\n",
                    rank, minimum_memory, memory_limit);
            setup_ok = 0;
        }
        double setup_t0 = GetTime();
        if (setup_ok) {
            if (!tracker.acquire(kSortStrictControlReserve)) {
                setup_ok = 0;
            } else {
                control_acquired = true;
            }
        }
        if (setup_ok &&
            MpiSortStrictAlloc(&tracker, kSortStrictSimScratch,
                               &sim_scratch) != 0) {
            setup_ok = 0;
        }
        if (setup_ok) {
            double temp_open_t0 = GetTime();
            if (MpiSortStrictOpenTemp(temp_prefix, rank, &temp_store) != 0) {
                setup_ok = 0;
            }
            stats->t_temp_open_actual += GetTime() - temp_open_t0;
        }

        const size_t stage_budget = setup_ok
            ? memory_limit - kSortStrictControlReserve - kSortStrictSimScratch
            : 0;
        size_t run_arena_bytes = stage_budget * 3 / 5;
        const size_t exchange_min = 8ull * 1024ull * 1024ull;
        if (setup_ok && run_arena_bytes + exchange_min > stage_budget) {
            run_arena_bytes = stage_budget - exchange_min;
        }
        if (setup_ok && run_arena_bytes + extract_workspace > stage_budget) {
            run_arena_bytes = stage_budget - extract_workspace;
        }
        run_arena_bytes -= run_arena_bytes % 64;
        if (setup_ok && run_arena_bytes < 2 * BGZF_MAX_BLOCK_SIZE) {
            fprintf(stderr,
                    "[rank %d] ERROR: -m leaves no usable run arena. arena=%zu limit=%zu\n",
                    rank, run_arena_bytes, memory_limit);
            setup_ok = 0;
        }
        if (setup_ok &&
            MpiSortStrictAlloc(&tracker, run_arena_bytes,
                               &run_arena_memory) != 0) {
            setup_ok = 0;
        }
        if (setup_ok) {
            run_arena_allocated = true;
            run_arena.data = run_arena_memory.data;
            run_arena.capacity = run_arena_memory.size;
            run_arena.reset();
            stats->run_arena_bytes = (long long)run_arena_bytes;
        }
        stats->t_setup +=
            std::max(0.0,
                     GetTime() - setup_t0 - stats->t_temp_open_actual);
        if (!all_ranks_ok(setup_ok)) break;

        int generate_ret = reader
            ? MpiSortStrictGenerateRuns(
                  *reader, global_block_begin, rank, comm_size,
                  &run_arena, &temp_store, sim_scratch.data,
                  sim_scratch.size, &tracker, &runs, &local_samples,
                  stats)
            : MpiSortStrictGenerateRuns(
                  *input, spans, span_count, global_block_begin,
                  rank, comm_size, &run_arena, &temp_store,
                  sim_scratch.data, sim_scratch.size, &tracker,
                  &runs, &local_samples, stats);
        int generate_ok = generate_ret == 0 ? 1 : 0;
        if (!all_ranks_ok(generate_ok)) break;
        stats->external_runs = (long long)runs.size();

        double sample_t0 = GetTime();
        long long global_sample_count = 0;
        if (MpiSortChooseSplitters(local_samples, rank, comm_size,
                                   &splitters,
                                   &global_sample_count) != 0) {
            break;
        }
        stats->sample_records = (long long)local_samples.size();
        stats->t_sample += GetTime() - sample_t0;

        const size_t exchange_budget = stage_budget - run_arena_bytes;
        int exchange_setup_ok =
            MpiSortStrictAllocateExchange(&tracker, exchange_budget,
                                          &exchange_ws) == 0 ? 1 : 0;
        if (!exchange_setup_ok) {
            fprintf(stderr,
                    "[rank %d] ERROR: -m cannot allocate bounded exchange buffers. exchange_budget=%zu limit=%zu\n",
                    rank, exchange_budget, memory_limit);
        }
        exchange_allocated = exchange_setup_ok != 0;
        if (!all_ranks_ok(exchange_setup_ok)) break;
        if (MpiSortStrictExchangeRuns(
                runs, splitters, &run_arena, &exchange_ws, &temp_store,
                sim_scratch.data, sim_scratch.size, rank, comm_size,
                &segments, stats) != 0) {
            break;
        }
        stats->external_segments = (long long)segments.size();
        MpiSortStrictFreeExchange(&tracker, &exchange_ws);
        exchange_allocated = false;
        MpiSortStrictFree(&tracker, &run_arena_memory);
        run_arena_allocated = false;

        const size_t cursor_stride =
            kSortStrictMetaBuffer + kSortStrictRawBuffer;
        const size_t merge_available =
            memory_limit - tracker.current;
        const size_t merge_fixed =
            MpiSortStrictCompressWorkspaceBytes() +
            kSortStrictConsolidateMeta + kSortStrictConsolidateRaw;
        if (merge_available <= merge_fixed + 2 * cursor_stride) {
            fprintf(stderr,
                    "[rank %d] ERROR: -m cannot support a two-way external merge. minimum_additional=%zu available=%zu\n",
                    rank, merge_fixed + 2 * cursor_stride,
                    merge_available);
            break;
        }
        size_t fan_in = (merge_available - merge_fixed) / cursor_stride;
        if (fan_in < 2) fan_in = 2;
        stats->merge_fan_in = (long long)fan_in;
        if (MpiSortStrictConsolidate(
                &segments, fan_in, &tracker, &temp_store,
                sim_scratch.data, sim_scratch.size, stats) != 0) {
            break;
        }
        stats->external_segments = (long long)segments.size();
        if (MpiSortStrictMergeCompress(
                segments, &tracker, &temp_store, sim_scratch.data,
                sim_scratch.size, body_sink, compress_level,
                stats) != 0) {
            break;
        }
        (void)global_sample_count;
        ret = 0;
    } while (false);

    double cleanup_t0 = GetTime();
    if (exchange_allocated) {
        MpiSortStrictFreeExchange(&tracker, &exchange_ws);
    }
    if (run_arena_allocated) {
        MpiSortStrictFree(&tracker, &run_arena_memory);
    }
    MpiSortStrictCloseTemp(&temp_store);
    MpiSortStrictFree(&tracker, &sim_scratch);
    if (control_acquired) tracker.release(kSortStrictControlReserve);
    stats->t_cleanup += GetTime() - cleanup_t0;
    stats->tracked_peak_bytes = (long long)tracker.peak;
    stats->t_fused_actual = GetTime() - fused_t0;
    double non_merge_read_sim =
        stats->t_temp_read_sim -
        stats->t_merge_temp_read_sim -
        stats->t_consolidation_temp_read_sim;
    double non_consolidation_write_sim =
        stats->t_temp_write_sim -
        stats->t_consolidation_temp_write_sim;
    if (non_merge_read_sim < 0.0) non_merge_read_sim = 0.0;
    if (non_consolidation_write_sim < 0.0) {
        non_consolidation_write_sim = 0.0;
    }
    const double final_pipeline_sim =
        std::max(stats->t_merge_simulated,
                 stats->t_compress_simulated);
    stats->t_fused_total =
        stats->t_setup +
        stats->t_extract_read_unhidden +
        stats->t_extract +
        stats->t_extract_prepare +
        stats->t_extract_merge +
        stats->t_local_sort +
        stats->t_sample +
        stats->t_bucket_count +
        stats->t_bucket_pack +
        stats->t_mpi_simulated +
        stats->t_offset_fix +
        non_merge_read_sim +
        non_consolidation_write_sim +
        stats->t_consolidation_simulated +
        final_pipeline_sim +
        stats->t_status_check;
    return ret;
}

int FusedBamExternalSortMPI(MemReader &reader,
                            swbam::RankBodySink &body_sink,
                            long long global_block_begin,
                            int rank,
                            int comm_size,
                            int compress_level,
                            size_t memory_limit,
                            const char *temp_prefix,
                            MpiSortStats *stats) {
    return FusedBamExternalSortMPIImpl(
        &reader, nullptr, nullptr, 0, body_sink,
        global_block_begin, rank, comm_size, compress_level,
        memory_limit, temp_prefix, stats);
}

int FusedBamExternalSortMPI(const swbam::BamInputBackend &input,
                            const swbam::BgzfBlockSpan *spans,
                            size_t span_count,
                            swbam::RankBodySink &body_sink,
                            long long global_block_begin,
                            int rank,
                            int comm_size,
                            int compress_level,
                            size_t memory_limit,
                            const char *temp_prefix,
                            MpiSortStats *stats) {
    return FusedBamExternalSortMPIImpl(
        nullptr, &input, spans, span_count, body_sink,
        global_block_begin, rank, comm_size, compress_level,
        memory_limit, temp_prefix, stats);
}

int FusedBamSortMPI(MemReader &reader,
                    MemWriter &mem_writer,
                    long long global_block_begin,
                    int rank,
                    int comm_size,
                    int compress_level,
                    size_t memory_limit,
                    MpiSortStats *stats) {
    swbam::MemoryRankBodySink body_sink(&mem_writer);
    return FusedBamSortMPI(
        reader, body_sink, global_block_begin, rank, comm_size,
        compress_level, memory_limit, stats);
}

int FusedBamExternalSortMPI(MemReader &reader,
                            MemWriter &mem_writer,
                            long long global_block_begin,
                            int rank,
                            int comm_size,
                            int compress_level,
                            size_t memory_limit,
                            const char *temp_prefix,
                            MpiSortStats *stats) {
    swbam::MemoryRankBodySink body_sink(&mem_writer);
    return FusedBamExternalSortMPI(
        reader, body_sink, global_block_begin, rank, comm_size,
        compress_level, memory_limit, temp_prefix, stats);
}
