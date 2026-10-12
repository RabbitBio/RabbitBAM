#include "swbam_mpi.h"
#include "swbam/io.h"
#include "swbam/mpi_runtime.h"

#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <libdeflate.h>
#include <mpi.h>
#include <sys/stat.h>
#include <vector>

namespace {

const unsigned char kMpiBgzfEofBlock[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
};

bool MpiIsBamLikeFormat(int format) {
    return format == bam || format == binary_format;
}

bool MpiIsSamLikeFormat(int format) {
    return format == sam || format == text_format;
}

int MpiNormalizeFormat(int format) {
    if (MpiIsBamLikeFormat(format)) return bam;
    if (MpiIsSamLikeFormat(format)) return sam;
    return format;
}

bool MpiHasSuffix(const std::string &name, const char *suffix) {
    size_t suffix_len = strlen(suffix);
    return name.size() >= suffix_len &&
           name.compare(name.size() - suffix_len, suffix_len, suffix) == 0;
}

int MpiOutputFormatFromName(const std::string &out_name) {
    if (MpiHasSuffix(out_name, ".sam")) return sam;
    return bam;
}

bool MpiHasBamFilterRequest(const CmdInfo *cmd_info) {
    return cmd_info->min_mapq_ >= 0 ||
           cmd_info->max_mapq_ >= 0 ||
           cmd_info->require_flag_ != 0 ||
           cmd_info->exclude_flag_ != 0 ||
           !cmd_info->ref_name_.empty() ||
           cmd_info->min_read_len_ >= 0 ||
           cmd_info->max_read_len_ >= 0;
}

BamFilterOptions MpiBuildBamFilterOptions(const CmdInfo *cmd_info) {
    BamFilterOptions filter;
    filter.min_mapq = cmd_info->min_mapq_;
    filter.max_mapq = cmd_info->max_mapq_;
    filter.require_flag = cmd_info->require_flag_;
    filter.exclude_flag = cmd_info->exclude_flag_;
    filter.ref_tid = -2;
    filter.min_read_len = cmd_info->min_read_len_;
    filter.max_read_len = cmd_info->max_read_len_;
    return filter;
}

int MpiAllRanksOk(int local_ok) {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return global_ok;
}

double MpiReduceMaxCost(double local_cost) {
    double max_cost = 0.0;
    MPI_Reduce(&local_cost, &max_cost, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return max_cost;
}

int MpiLoadFileToMemory(const std::string &path, char **data, size_t *size) {
    return swbam::LoadFileToMemory(path, data, size);
}

int MpiBgzfBlockLengthAt(const char *base, size_t size,
                         size_t pos, size_t *block_len);

int MpiScanBgzfBlocksInMemory(const char *base, size_t size, long long body_start,
                              std::vector<long long> *offsets,
                              std::vector<long long> *lengths) {
    return swbam::ScanBgzfBlocksInMemory(
        base, size, body_start, offsets, lengths);
}

int MpiSelectBlockRangeFromMemory(char *base, size_t input_size,
                                  const std::vector<long long> &offsets,
                                  const std::vector<long long> &lengths,
                                  long long begin,
                                  long long end,
                                  char **data,
                                  size_t *size) {
    return swbam::SelectBlockRangeFromMemory(
        base, input_size, offsets, lengths,
        begin, end, data, size);
}

int MpiDumpMemoryToFile(const std::string &path, const char *data, size_t size) {
    swbam::PosixBamOutput output;
    if (output.Open(path) != 0 ||
        output.Write(data, size) != 0 ||
        output.Flush() != 0) {
        return -1;
    }
    return output.Close();
}

int MpiSendBytes(int dst, int tag, const char *data, long long len) {
    long long sent = 0;
    while (sent < len) {
        int chunk = (int)std::min<long long>(len - sent, INT_MAX);
        if (MPI_Send((void *)(data + sent), chunk, MPI_BYTE, dst, tag, MPI_COMM_WORLD) != MPI_SUCCESS) {
            return -1;
        }
        sent += chunk;
    }
    return 0;
}

int MpiRecvBytes(int src, int tag, char *data, long long len) {
    long long received = 0;
    while (received < len) {
        int chunk = (int)std::min<long long>(len - received, INT_MAX);
        MPI_Status status;
        if (MPI_Recv(data + received, chunk, MPI_BYTE, src, tag, MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
            return -1;
        }
        received += chunk;
    }
    return 0;
}

int MpiEnsureMemWriterCapacity(MemWriter &w, size_t add) {
    if (w.size + add <= w.capacity) return 0;
    size_t next = w.capacity ? w.capacity : 64 * 1024 * 1024;
    while (w.size + add > next) {
        size_t prev = next;
        next *= 2;
        if (next < prev) return -1;
    }
    char *new_data = (char *)realloc(w.data, next);
    if (!new_data) return -1;
    w.data = new_data;
    w.capacity = next;
    return 0;
}

int MpiInitMemWriter(MemWriter &w, size_t cap) {
    if (cap == 0) cap = 64 * 1024 * 1024;
    w.data = (char *)malloc(cap);
    w.capacity = w.data ? cap : 0;
    w.size = 0;
    return w.data ? 0 : -1;
}

void MpiPackLe16(unsigned char *buffer, uint16_t value) {
    buffer[0] = (unsigned char)(value & 0xff);
    buffer[1] = (unsigned char)((value >> 8) & 0xff);
}

void MpiPackLe32(unsigned char *buffer, uint32_t value) {
    buffer[0] = (unsigned char)(value & 0xff);
    buffer[1] = (unsigned char)((value >> 8) & 0xff);
    buffer[2] = (unsigned char)((value >> 16) & 0xff);
    buffer[3] = (unsigned char)((value >> 24) & 0xff);
}

uint32_t MpiReadLe32(const unsigned char *buffer) {
    return (uint32_t)buffer[0] |
           ((uint32_t)buffer[1] << 8) |
           ((uint32_t)buffer[2] << 16) |
           ((uint32_t)buffer[3] << 24);
}

int MpiBgzfBlockLengthAt(const char *base, size_t size,
                         size_t pos, size_t *block_len) {
    if (!base || !block_len || pos > size ||
        size - pos < BLOCK_HEADER_LENGTH) {
        return -1;
    }
    const unsigned char *block =
        (const unsigned char *)(base + pos);
    if (block[0] != 0x1f || block[1] != 0x8b ||
        block[2] != 0x08 || block[3] != 0x04 ||
        block[10] != 0x06 || block[12] != 'B' ||
        block[13] != 'C' || block[14] != 0x02 ||
        block[15] != 0x00) {
        return -1;
    }
    int len = (int)block[16] | ((int)block[17] << 8);
    len++;
    if (len <= 0 || len > BGZF_MAX_BLOCK_SIZE ||
        (size_t)len > size - pos) {
        return -1;
    }
    *block_len = (size_t)len;
    return 0;
}

void MpiVectorPutLe32(std::vector<unsigned char> *out, uint32_t value) {
    unsigned char buf[4];
    MpiPackLe32(buf, value);
    out->insert(out->end(), buf, buf + 4);
}

int MpiInflateBgzfBlockForHeader(
        const char *base, size_t size, size_t pos,
        struct libdeflate_decompressor *decompressor,
        std::vector<unsigned char> *raw,
        size_t *block_len) {
    if (!base || !decompressor || !raw || !block_len) {
        return -1;
    }
    size_t len = 0;
    if (MpiBgzfBlockLengthAt(base, size, pos, &len) != 0 ||
        len <= BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH) {
        return -1;
    }
    const unsigned char *block =
        (const unsigned char *)(base + pos);
    uint32_t isize = MpiReadLe32(block + len - 4);
    if (isize > BGZF_MAX_BLOCK_SIZE) return -1;
    size_t old_size = raw->size();
    raw->resize(old_size + isize);
    size_t actual = 0;
    enum libdeflate_result ret =
        libdeflate_deflate_decompress(
            decompressor,
            block + BLOCK_HEADER_LENGTH,
            len - BLOCK_HEADER_LENGTH -
                BLOCK_FOOTER_LENGTH,
            raw->data() + old_size, isize, &actual);
    if (ret != LIBDEFLATE_SUCCESS || actual != isize) {
        return -1;
    }
    uint32_t crc = libdeflate_crc32(
        0, raw->data() + old_size, isize);
    if (crc != MpiReadLe32(block + len - 8)) {
        return -1;
    }
    *block_len = len;
    return 0;
}

int MpiBgzfBlockStartsWithBamRecord(
        const char *base, size_t size, size_t pos,
        struct libdeflate_decompressor *decompressor) {
    std::vector<unsigned char> raw;
    raw.reserve(BGZF_BLOCK_SIZE);
    size_t block_len = 0;
    if (MpiInflateBgzfBlockForHeader(
            base, size, pos, decompressor,
            &raw, &block_len) != 0) {
        return 0;
    }
    if (raw.empty()) return 0;
    if (raw.size() < 36) return 0;
    int32_t record_len = (int32_t)MpiReadLe32(raw.data());
    if (record_len < 32 ||
        (uint64_t)record_len + 4u > raw.size()) {
        return 0;
    }
    uint32_t x2 = MpiReadLe32(raw.data() + 12);
    uint32_t x3 = MpiReadLe32(raw.data() + 16);
    uint32_t l_qname = x2 & 0xffu;
    uint32_t n_cigar = x3 & 0xffffu;
    uint32_t l_qseq = MpiReadLe32(raw.data() + 20);
    uint64_t min_payload =
        (uint64_t)l_qname + ((uint64_t)n_cigar << 2) +
        ((uint64_t)l_qseq + 1u) / 2u + (uint64_t)l_qseq;
    if (l_qname == 0 || min_payload >
        (uint64_t)record_len - 32u) {
        return 0;
    }
    return 1;
}

int MpiFindNextBamBodyBgzfBlock(
        const char *base, size_t size, size_t start,
        struct libdeflate_decompressor *decompressor,
        size_t *body_start) {
    if (!base || !decompressor || !body_start || start > size) {
        return -1;
    }
    for (size_t pos = start; pos + BLOCK_HEADER_LENGTH <= size;
         ++pos) {
        size_t len = 0;
        if (MpiBgzfBlockLengthAt(base, size, pos, &len) != 0) {
            continue;
        }
        if (MpiBgzfBlockStartsWithBamRecord(
                base, size, pos, decompressor)) {
            *body_start = pos;
            return 0;
        }
    }
    return -1;
}

int MpiParseBamHeaderRaw(const std::vector<unsigned char> &raw,
                         size_t *header_raw_size,
                         uint32_t *text_len) {
    if (raw.size() < 8) return 1;
    if (memcmp(raw.data(), "BAM\1", 4) != 0) return -1;
    uint32_t l_text = MpiReadLe32(raw.data() + 4);
    size_t pos = 8u + (size_t)l_text;
    if (pos > raw.size()) return 1;
    if (raw.size() - pos < 4) return 1;
    uint32_t n_ref = MpiReadLe32(raw.data() + pos);
    pos += 4;
    for (uint32_t i = 0; i < n_ref; ++i) {
        if (raw.size() - pos < 4) return 1;
        uint32_t l_name = MpiReadLe32(raw.data() + pos);
        pos += 4;
        if ((size_t)l_name > raw.size() - pos) return 1;
        pos += (size_t)l_name;
        if (raw.size() - pos < 4) return 1;
        pos += 4;
    }
    *header_raw_size = pos;
    *text_len = l_text;
    return 0;
}

int MpiFindSamBodyStartInMemory(const char *base, size_t size, long long *body_start) {
    return swbam::FindSamBodyStartInMemory(base, size, body_start);
}

int MpiSplitSamRangesInMemory(const char *base, size_t size, long long body_start,
                              int comm_size,
                              std::vector<long long> *offsets,
                              std::vector<long long> *lengths) {
    return swbam::SplitSamRangesInMemory(
        base, size, body_start, comm_size, offsets, lengths);
}

int MpiAppendBgzfPayloadToMem(MemWriter &w, const unsigned char *src, size_t src_len,
                              int compress_level) {
    if (src_len == 0) return 0;
    if (src_len > BGZF_BLOCK_SIZE) return -1;
    if (compress_level != 0 && compress_level != 1 && compress_level != 6) return -1;

    unsigned char out[BGZF_MAX_BLOCK_SIZE];
    size_t block_len = 0;
    if (compress_level == 0) {
        if (src_len + 5 + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH > BGZF_MAX_BLOCK_SIZE) return -1;
        out[BLOCK_HEADER_LENGTH] = 1;
        MpiPackLe16(out + BLOCK_HEADER_LENGTH + 1, (uint16_t)src_len);
        MpiPackLe16(out + BLOCK_HEADER_LENGTH + 3, (uint16_t)~(uint16_t)src_len);
        memcpy(out + BLOCK_HEADER_LENGTH + 5, src, src_len);
        block_len = src_len + 5 + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    } else {
        struct libdeflate_compressor *compressor = libdeflate_alloc_compressor(compress_level);
        if (!compressor) return -1;
        size_t clen = libdeflate_deflate_compress(
            compressor, src, src_len,
            out + BLOCK_HEADER_LENGTH,
            BGZF_MAX_BLOCK_SIZE - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH);
        libdeflate_free_compressor(compressor);
        if (clen == 0) return -1;
        block_len = clen + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    }
    if (block_len > BGZF_MAX_BLOCK_SIZE) return -1;

    memcpy(out, g_magic, BLOCK_HEADER_LENGTH);
    MpiPackLe16(out + 16, (uint16_t)(block_len - 1));
    uint32_t crc = libdeflate_crc32(0, src, src_len);
    MpiPackLe32(out + block_len - 8, crc);
    MpiPackLe32(out + block_len - 4, (uint32_t)src_len);
    return MpiWriteBytesToMem(w, (const char *)out, block_len);
}

int MpiBuildBamHeaderMemory(sam_hdr_t *hdr, int compress_level, char **data, size_t *size) {
    *data = nullptr;
    *size = 0;
    if (!hdr) return -1;

    const char *text = sam_hdr_str(hdr);
    size_t text_len = sam_hdr_length(hdr);
    if (text_len > 0 && !text) return -1;
    if ((unsigned long long)text_len > (unsigned long long)UINT32_MAX) return -1;
    if (hdr->n_targets < 0) return -1;

    std::vector<unsigned char> raw;
    size_t reserve_len = 12 + text_len + (size_t)hdr->n_targets * 32;
    raw.reserve(reserve_len);
    const unsigned char bam_magic[4] = {'B', 'A', 'M', 1};
    raw.insert(raw.end(), bam_magic, bam_magic + 4);
    MpiVectorPutLe32(&raw, (uint32_t)text_len);
    if (text_len > 0) raw.insert(raw.end(), text, text + text_len);
    MpiVectorPutLe32(&raw, (uint32_t)hdr->n_targets);
    for (int32_t i = 0; i < hdr->n_targets; ++i) {
        if (!hdr->target_name || !hdr->target_name[i] || !hdr->target_len) return -1;
        size_t name_len = strlen(hdr->target_name[i]) + 1;
        if ((unsigned long long)name_len > (unsigned long long)UINT32_MAX) return -1;
        MpiVectorPutLe32(&raw, (uint32_t)name_len);
        raw.insert(raw.end(), hdr->target_name[i], hdr->target_name[i] + name_len);
        MpiVectorPutLe32(&raw, hdr->target_len[i]);
    }

    MemWriter header_writer = {};
    if (MpiInitMemWriter(header_writer, raw.size() + 64 * 1024) != 0) return -1;
    size_t pos = 0;
    while (pos < raw.size()) {
        size_t n = std::min((size_t)BGZF_BLOCK_SIZE, raw.size() - pos);
        if (MpiAppendBgzfPayloadToMem(header_writer, raw.data() + pos, n, compress_level) != 0) {
            free(header_writer.data);
            return -1;
        }
        pos += n;
    }

    *data = header_writer.data;
    *size = header_writer.size;
    return 0;
}

void MpiPrintRankBufferedLines(int rank, int comm_size, const char *local_lines) {
    const int kLineBytes = 4096;
    std::vector<char> gathered;
    if (rank == 0) gathered.resize((size_t)comm_size * kLineBytes);
    MPI_Gather((void *)local_lines, kLineBytes, MPI_CHAR,
               rank == 0 ? gathered.data() : nullptr, kLineBytes, MPI_CHAR,
               0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int r = 0; r < comm_size; ++r) {
            fputs(gathered.data() + (size_t)r * kLineBytes, stdout);
        }
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void MpiPrintRankStats(int rank, int comm_size, const MpiBamToBamStats &stats, size_t body_size) {
    char local_lines[4096] = {};
    double keep_ratio = stats.total_records > 0
                        ? (double)stats.kept_records / (double)stats.total_records
                        : 1.0;
    double core_stage = stats.t_decomp_filter + stats.t_pack + stats.t_compress;
    snprintf(local_lines, sizeof(local_lines),
             "[rank %d] blocks=%lld groups=%lld records=%lld kept=%lld dropped=%lld bgzf=%lld body=%zu keep_ratio=%.6f\n"
             "[rank %d] decomp_filter_slave=%.3f  pack=%.3f  compress_stage=%.3f  read=%.3f  write=%.3f  mpi_write=%.3f\n"
             "[rank %d] decomp_detail alloc=%.3f  inflate=%.3f  crc=%.3f  parse_filter=%.3f  other=%.3f\n"
             "[rank %d] compress_detail serialize=%.3f  alloc=%.3f  deflate=%.3f  footer=%.3f  other=%.3f\n"
             "[rank %d] optimized_total=%.6f  core_stage=%.6f\n",
             rank, stats.input_blocks, stats.group_count, stats.total_records,
             stats.kept_records, stats.dropped_records, stats.bgzf_blocks, body_size, keep_ratio,
             rank, stats.t_decomp_filter, stats.t_pack, stats.t_compress,
             stats.t_read, stats.t_write, stats.t_mpi_write,
             rank, stats.t_decomp_alloc, stats.t_decomp_inflate,
             stats.t_decomp_crc, stats.t_decomp_parse, stats.t_decomp_other,
             rank, stats.t_compress_serialize, stats.t_compress_alloc,
             stats.t_compress_deflate, stats.t_compress_footer, stats.t_compress_other,
             rank, stats.t_optimized_total, core_stage);
    MpiPrintRankBufferedLines(rank, comm_size, local_lines);
}

void MpiPrintRankBamToSamStats(int rank, int comm_size,
                               const MpiBamToSamStats &stats,
                               size_t body_size) {
    char local_lines[4096] = {};
    double core_stage = stats.t_decomp + stats.t_collect + stats.t_format;
    snprintf(local_lines, sizeof(local_lines),
             "[rank %d] blocks=%lld groups=%lld records=%lld format_tiles=%lld body=%zu\n"
             "[rank %d] decomp_slave=%.3f  format_slave=%.3f  collect=%.3f  read=%.3f  write=%.3f  gather=%.3f\n"
             "[rank %d] decomp_detail alloc=%.3f  inflate=%.3f  crc=%.3f  parse=%.3f  other=%.3f\n"
             "[rank %d] optimized_total=%.6f  core_stage=%.6f  alloc_init=%.6f  free_workspace=%.6f\n",
             rank, stats.input_blocks, stats.group_count, stats.total_records,
             stats.format_tiles, body_size,
             rank, stats.t_decomp, stats.t_format, stats.t_collect,
             stats.t_read, stats.t_write, stats.t_gather,
             rank, stats.t_decomp_alloc, stats.t_decomp_inflate,
             stats.t_decomp_crc, stats.t_decomp_parse, stats.t_decomp_other,
             rank, stats.t_optimized_total, core_stage,
             stats.t_alloc_init, stats.t_free_workspace);
    MpiPrintRankBufferedLines(rank, comm_size, local_lines);
}

void MpiPrintRankSamToBamStats(int rank, int comm_size,
                               const MpiSamToBamStats &stats,
                               size_t body_size) {
    char local_lines[4096] = {};
    double core_stage = stats.t_copy_count + stats.t_parse + stats.t_pack + stats.t_compress;
    // double measured_extra = stats.t_split + stats.t_alloc_init + stats.t_setup_reset +
    //                         stats.t_compress_setup + stats.t_status_check +
    //                         stats.t_free_workspace;
    // double residual = stats.t_optimized_total - core_stage - measured_extra;
    snprintf(local_lines, sizeof(local_lines),
             "[rank %d] chunks=%lld chunk_groups=%lld records=%lld compress_groups=%lld bgzf=%lld body=%zu\n"
             "[rank %d] parse_fast=%lld  parse_fallback=%lld\n"
             "[rank %d] split=%.3f  copy_count_slave=%.3f  parse_slave=%.3f  pack=%.3f  compress_slave=%.3f  write=%.3f  gather=%.3f\n"
             "[rank %d] parse_detail core=%.3f  aux=%.3f  cg=%.3f  fallback=%.3f  other=%.3f\n"
             "[rank %d] compress_detail serialize=%.3f  alloc=%.3f  deflate=%.3f  footer=%.3f  other=%.3f\n"
             "[rank %d] optimized_total=%.6f  core_stage=%.6f\n",
             rank, stats.input_chunks, stats.chunk_groups, stats.total_records,
             stats.compress_groups, stats.bgzf_blocks, body_size,
             rank, stats.parse_fast_records, stats.parse_fallback_records,
             rank, stats.t_split, stats.t_copy_count, stats.t_parse,
             stats.t_pack, stats.t_compress, stats.t_write, stats.t_gather,
             rank, stats.t_parse_core, stats.t_parse_aux,
             stats.t_parse_cg, stats.t_parse_fallback, stats.t_parse_other,
             rank, stats.t_compress_serialize, stats.t_compress_alloc,
             stats.t_compress_deflate, stats.t_compress_footer, stats.t_compress_other,
             rank, stats.t_optimized_total, core_stage);
    MpiPrintRankBufferedLines(rank, comm_size, local_lines);
}

} // namespace

int MpiCommonLoadFileToMemory(const std::string &path,
                              char **data, size_t *size) {
    return MpiLoadFileToMemory(path, data, size);
}

int MpiCommonScanBgzfBlocksInMemory(
        const char *base, size_t size, long long body_start,
        std::vector<long long> *offsets,
        std::vector<long long> *lengths) {
    return MpiScanBgzfBlocksInMemory(
        base, size, body_start, offsets, lengths);
}

int MpiCommonSelectBlockRangeFromMemory(
        char *base, size_t input_size,
        const std::vector<long long> &offsets,
        const std::vector<long long> &lengths,
        long long begin, long long end,
        char **data, size_t *size) {
    return MpiSelectBlockRangeFromMemory(
        base, input_size, offsets, lengths,
        begin, end, data, size);
}

int MpiCommonDumpMemoryToFile(const std::string &path,
                              const char *data, size_t size) {
    return MpiDumpMemoryToFile(path, data, size);
}

int MpiCommonInitMemWriter(MemWriter &writer, size_t capacity) {
    return MpiInitMemWriter(writer, capacity);
}

int MpiCommonBuildBamHeaderMemory(sam_hdr_t *header,
                                  int compress_level,
                                  char **data, size_t *size) {
    return MpiBuildBamHeaderMemory(
        header, compress_level, data, size);
}

int MpiCommonReadBamHeaderFromMemory(const char *data, size_t size,
                                     sam_hdr_t **header,
                                     long long *body_start) {
    if (header) *header = nullptr;
    if (body_start) *body_start = 0;
    if (!data || !header || !body_start || size == 0) {
        return -1;
    }

    std::vector<unsigned char> raw;
    raw.reserve(64 * 1024);
    struct libdeflate_decompressor *decompressor =
        libdeflate_alloc_decompressor();
    if (!decompressor) return -1;

    int ret = -1;
    size_t pos = 0;
    while (pos < size) {
        size_t block_len = 0;
        if (MpiInflateBgzfBlockForHeader(
                data, size, pos, decompressor,
                &raw, &block_len) != 0) {
            break;
        }
        pos += block_len;

        size_t header_raw_size = 0;
        uint32_t text_len = 0;
        int parse_state = MpiParseBamHeaderRaw(
            raw, &header_raw_size, &text_len);
        if (parse_state < 0) break;
        if (parse_state > 0) continue;

        if (raw.size() != header_raw_size) break;
        sam_hdr_t *parsed = sam_hdr_parse(
            (size_t)text_len, (const char *)raw.data() + 8);
        if (!parsed) break;
        size_t aligned_body_start = 0;
        if (MpiFindNextBamBodyBgzfBlock(
                data, size, pos, decompressor,
                &aligned_body_start) != 0) {
            sam_hdr_destroy(parsed);
            break;
        }
        *header = parsed;
        *body_start = (long long)aligned_body_start;
        ret = 0;
        break;
    }

    libdeflate_free_decompressor(decompressor);
    return ret;
}

void MpiMemoryBamFree(MpiMemoryBam *bam) {
    if (!bam) return;
    free(bam->data);
    bam->data = nullptr;
    bam->size = 0;
}

int MpiBroadcastMemoryBam(MpiMemoryBam *bam, int root) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int local_ok = bam != nullptr ? 1 : 0;
    unsigned long long wire_size = 0;
    if (rank == root && local_ok) {
        wire_size = (unsigned long long)bam->size;
        if ((size_t)wire_size != bam->size ||
            (bam->size > 0 && bam->data == nullptr)) {
            local_ok = 0;
        }
    }

    MPI_Bcast(&wire_size, 1, MPI_UNSIGNED_LONG_LONG,
              root, MPI_COMM_WORLD);
    if (rank != root && local_ok) {
        if (wire_size > (unsigned long long)SIZE_MAX) {
            local_ok = 0;
        } else {
            MpiMemoryBamFree(bam);
            bam->size = (size_t)wire_size;
            bam->data = bam->size ? (char *)malloc(bam->size) : nullptr;
            if (bam->size > 0 && !bam->data) {
                local_ok = 0;
            }
        }
    }

    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    if (!global_ok) return -1;

    const unsigned long long kBroadcastChunk = 16ull * 1024ull * 1024ull;
    unsigned long long done = 0;
    local_ok = 1;
    while (done < wire_size) {
        int chunk = (int)std::min<unsigned long long>(
            wire_size - done, kBroadcastChunk);
        if (MPI_Bcast(bam->data + done, chunk, MPI_BYTE,
                      root, MPI_COMM_WORLD) != MPI_SUCCESS) {
            local_ok = 0;
        }
        done += (unsigned long long)chunk;
    }
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    return global_ok ? 0 : -1;
}

int MpiWriteBlockToMem(MemWriter &w, bam_block *block) {
    if (MpiEnsureMemWriterCapacity(w, block->length) != 0) return -1;
    memcpy(w.data + w.size, block->data, block->length);
    w.size += block->length;
    return 0;
}

int MpiWriteBytesToMem(MemWriter &w, const char *data, size_t len) {
    if (len == 0) return 0;
    if (MpiEnsureMemWriterCapacity(w, len) != 0) return -1;
    memcpy(w.data + w.size, data, len);
    w.size += len;
    return 0;
}

int ProcessSwBamMPI(CmdInfo *cmd_info) {
    const double body_start_time = GetTime();
    int rank = 0, ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    swbam::mpi::MpiBamInput input;
    swbam::mpi::MpiBamInputPlan plan;
    swbam::PosixSamInput native_sam;
    swbam::mpi::MpiIoSamInput mpiio_sam;
    swbam::SamInputBackend *sam_input = nullptr;
    swbam_usage::BodySink body;
    sam_hdr_t *header = nullptr;
    char *prefix = nullptr;
    size_t prefix_size = 0;
    const int input_format = swbam_usage::InputFormat(cmd_info->in_file_name_);
    const bool from_bam = MpiIsBamLikeFormat(input_format);
    const bool to_sam = MpiHasSuffix(cmd_info->out_file_name_, ".sam");
    int status = 1;
    const auto run = [&]() -> int {
        int ok = 1;
        if (!MpiAllRanksOk(from_bam || MpiIsSamLikeFormat(input_format))) {
            if (rank == 0) fprintf(stderr, "ERROR: conversion requires BAM or uncompressed SAM input\n");
            return 1;
        }
        if (from_bam) {
            ok = swbam_usage::OpenInput(input, cmd_info->in_file_name_,
                cmd_info->io_backend_, cmd_info->io_memory_limit_) == 0;
            if (ok) header = sam_hdr_dup(input.backend()->header());
        } else {
            sam_input = swbam_usage::IsRam()
                ? swbam_usage::CurrentSession()->sam_input()
                : cmd_info->io_backend_ == "mpiio"
                    ? static_cast<swbam::SamInputBackend *>(&mpiio_sam)
                    : static_cast<swbam::SamInputBackend *>(&native_sam);
            ok = sam_input && sam_input->Open(cmd_info->in_file_name_) == 0;
            if (ok) header = sam_hdr_dup(sam_input->header());
        }
        if (!MpiAllRanksOk(ok && header)) return 1;
        if ((!from_bam && to_sam) ||
            (MpiHasBamFilterRequest(cmd_info) && (!from_bam || to_sam))) {
            if (rank == 0) fprintf(stderr, "ERROR: filtering requires BAM -> BAM; SAM -> SAM is unsupported\n");
            return 1;
        }
        size_t input_bytes = 0;
        MemReader reader = {};
        std::vector<char> sam_storage;
        if (from_bam) {
            ok = input.PreparePlan(&plan) == 0;
            if (ok)
                for (size_t i = plan.rank_begin; i < plan.rank_end; ++i)
                    input_bytes += plan.blocks[i].compressed_size;
        } else {
            std::vector<long long> offsets(ranks), lengths(ranks);
            if (rank == 0) ok = sam_input->SplitRanges(ranks, &offsets, &lengths) == 0;
            if (!MpiAllRanksOk(ok)) return 1;
            MPI_Bcast(offsets.data(), ranks, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(lengths.data(), ranks, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
            input_bytes = lengths[rank];
            reader.base = sam_input->ResidentRange(offsets[rank], input_bytes);
            if (!reader.base && input_bytes) {
                ok = sam_input->ReadRange(offsets[rank], input_bytes, &sam_storage) == 0;
                reader.base = sam_storage.empty() ? nullptr : sam_storage.data();
            }
            reader.size = input_bytes;
        }
        uint64_t limit = 0;
        ok = ok && swbam::ParseByteSize(cmd_info->rank_body_memory_limit_, &limit) == 0 &&
            body.Open(cmd_info->rank_body_backend_, limit,
                      std::min<size_t>(input_bytes, 64u * 1024u * 1024u),
                      "swbam-conversion.body", cmd_info->rank_body_temp_dir_) == 0;
        if (!MpiAllRanksOk(ok)) return 1;
        BamFilterOptions filter = MpiBuildBamFilterOptions(cmd_info);
        if (!cmd_info->ref_name_.empty()) {
            filter.ref_tid = sam_hdr_name2tid(header, cmd_info->ref_name_.c_str());
            ok = filter.ref_tid >= 0;
        }
        if (cmd_info->max_mapq_ >= 0 && cmd_info->min_mapq_ > cmd_info->max_mapq_) ok = 0;
        if (cmd_info->max_read_len_ >= 0 && cmd_info->min_read_len_ > cmd_info->max_read_len_) ok = 0;
        if (cmd_info->validate_bounds_) ok = 0;
        if (!MpiAllRanksOk(ok)) {
            if (rank == 0) fprintf(stderr, "ERROR: invalid filter bounds/reference or unsupported --validate-bounds (use io-check)\n");
            return 1;
        }
        MpiBamToBamStats bm = {};
        MpiBamToSamStats bs = {};
        MpiSamToBamStats sb = {};
        MPI_Barrier(MPI_COMM_WORLD);
        const double core_start = GetTime();
        if (!from_bam)
            ok = OptimizedSamToBamMPI(reader, body, header, cmd_info->compress_level_, &sb) == 0;
        else if (to_sam)
            ok = OptimizedBamToSamMPI(*input.backend(), plan.rank_spans(),
                plan.rank_block_count(), body, header, &bs) == 0;
        else
            ok = OptimizedBamToBamMPI(*input.backend(), plan.rank_spans(),
                plan.rank_block_count(), body, filter, cmd_info->compress_level_, &bm) == 0;
        const double core = MpiReduceMaxCost(GetTime() - core_start);
        if (!MpiAllRanksOk(ok && body.Flush() == 0)) return 1;
        long long local[3] = {
            from_bam ? to_sam ? bs.total_records : bm.total_records : sb.total_records,
            from_bam && !to_sam ? bm.kept_records :
                from_bam ? bs.total_records : sb.total_records,
            (long long)body.size()
        }, totals[3] = {};
        MPI_Reduce(local, totals, 3, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        const char *name = !from_bam ? "OptimizedSamToBamMPI" :
            to_sam ? "OptimizedBamToSamMPI" : "OptimizedBamToBamMPI";
        if (rank == 0) {
            printf("Complete the 4.3 %s cost %lf\n", name, core);
            printf("%s finished. ranks=%d total_records=%lld kept_records=%lld body_bytes=%lld\n",
                   name, ranks, totals[0], totals[1], totals[2]);
            if (!to_sam) ok = MpiBuildBamHeaderMemory(header, cmd_info->compress_level_, &prefix, &prefix_size) == 0;
        }
        if (!MpiAllRanksOk(ok)) return 1;
        const void *header_bytes = to_sam ? sam_hdr_str(header) : prefix;
        const size_t header_size = to_sam ? sam_hdr_length(header) : prefix_size;
        const double output_start = GetTime();
        ok = swbam_usage::WriteOutput(cmd_info->out_file_name_, body,
            header_size, header_bytes, to_sam ? 0 : sizeof(kMpiBgzfEofBlock),
            to_sam ? nullptr : kMpiBgzfEofBlock) == 0;
        const double output_time = MpiReduceMaxCost(GetTime() - output_start);
        if (rank == 0) printf("Output commit cost %lf\n", output_time);
        if (!from_bam) MpiPrintRankSamToBamStats(rank, ranks, sb, body.size());
        else if (to_sam) MpiPrintRankBamToSamStats(rank, ranks, bs, body.size());
        else MpiPrintRankStats(rank, ranks, bm, body.size());
        return MpiAllRanksOk(ok) ? 0 : 1;
    };
    status = run();
    free(prefix);
    if (header) sam_hdr_destroy(header);
    if (sam_input) sam_input->Close();
    input.Close();
    if (body.Close()) status = 1;
    const double elapsed = MpiReduceMaxCost(GetTime() - body_start_time);
    if (rank == 0) printf("444Complete the total body cost %lf\n", elapsed);
    return status;
}

int ProcessBamToBamMPI(CmdInfo *cmd_info) {
    return ProcessSwBamMPI(cmd_info);
}
