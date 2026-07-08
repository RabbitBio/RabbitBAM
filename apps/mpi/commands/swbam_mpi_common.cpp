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
             "[rank %d] fused_total=%.6f  core_stage=%.6f\n",
             rank, stats.input_blocks, stats.group_count, stats.total_records,
             stats.kept_records, stats.dropped_records, stats.bgzf_blocks, body_size, keep_ratio,
             rank, stats.t_decomp_filter, stats.t_pack, stats.t_compress,
             stats.t_read, stats.t_write, stats.t_mpi_write,
             rank, stats.t_decomp_alloc, stats.t_decomp_inflate,
             stats.t_decomp_crc, stats.t_decomp_parse, stats.t_decomp_other,
             rank, stats.t_compress_serialize, stats.t_compress_alloc,
             stats.t_compress_deflate, stats.t_compress_footer, stats.t_compress_other,
             rank, stats.t_fused_total, core_stage);
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
             "[rank %d] fused_total=%.6f  core_stage=%.6f  alloc_init=%.6f  free_workspace=%.6f\n",
             rank, stats.input_blocks, stats.group_count, stats.total_records,
             stats.format_tiles, body_size,
             rank, stats.t_decomp, stats.t_format, stats.t_collect,
             stats.t_read, stats.t_write, stats.t_gather,
             rank, stats.t_decomp_alloc, stats.t_decomp_inflate,
             stats.t_decomp_crc, stats.t_decomp_parse, stats.t_decomp_other,
             rank, stats.t_fused_total, core_stage,
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
    // double residual = stats.t_fused_total - core_stage - measured_extra;
    snprintf(local_lines, sizeof(local_lines),
             "[rank %d] chunks=%lld chunk_groups=%lld records=%lld compress_groups=%lld bgzf=%lld body=%zu\n"
             "[rank %d] parse_fast=%lld  parse_fallback=%lld\n"
             "[rank %d] split=%.3f  copy_count_slave=%.3f  parse_slave=%.3f  pack=%.3f  compress_slave=%.3f  write=%.3f  gather=%.3f\n"
             "[rank %d] parse_detail core=%.3f  aux=%.3f  cg=%.3f  fallback=%.3f  other=%.3f\n"
             "[rank %d] compress_detail serialize=%.3f  alloc=%.3f  deflate=%.3f  footer=%.3f  other=%.3f\n"
             "[rank %d] fused_total=%.6f  core_stage=%.6f\n",
             rank, stats.input_chunks, stats.chunk_groups, stats.total_records,
             stats.compress_groups, stats.bgzf_blocks, body_size,
             rank, stats.parse_fast_records, stats.parse_fallback_records,
             rank, stats.t_split, stats.t_copy_count, stats.t_parse,
             stats.t_pack, stats.t_compress, stats.t_write, stats.t_gather,
             rank, stats.t_parse_core, stats.t_parse_aux,
             stats.t_parse_cg, stats.t_parse_fallback, stats.t_parse_other,
             rank, stats.t_compress_serialize, stats.t_compress_alloc,
             stats.t_compress_deflate, stats.t_compress_footer, stats.t_compress_other,
             rank, stats.t_fused_total, core_stage);
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

    //1.相关数据的初始化
    double t_init = GetTime();

    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    double init_cost, init_cost_max;
    double t_total = 0.0;
    int exit_code = 1;
    int local_ok = 1;
    samFile *sin = nullptr;
    sam_hdr_t *hdr = nullptr;
    swbam::mpi::MpiBamInput backend_input_handle;
    swbam::BamInputBackend *backend_input = nullptr;
    swbam::mpi::MpiBamInputPlan backend_plan;
    swbam::PosixSamInput posix_sam_input;
    swbam::mpi::MpiIoSamInput mpiio_sam_input;
    swbam::SamInputBackend *sam_input = nullptr;
    std::vector<char> rank_sam_storage;
    hFILE *input_mem_hfile = nullptr;
    char *input_file_mem = nullptr;
    size_t input_file_size = 0;
    char *rank_input_mem = nullptr;
    char *output_file_mem = nullptr;
    char *bam_header_mem = nullptr;
    char *simulated_write_mem = nullptr;
    size_t bam_header_size = 0;
    size_t output_body_start = 0;
    size_t output_file_size = 0;
    size_t simulated_write_size = 0;
    const char *sam_header_text = nullptr;
    size_t sam_header_len = 0;
    MemReader reader = {};
    swbam::AdaptiveRankBodySink body_sink;
    uint64_t rank_body_memory_limit = 0;
    MpiBamToBamStats stats = {};
    MpiBamToSamStats sam_stats = {};
    MpiSamToBamStats sam2bam_stats = {};
    BamFilterOptions filter = MpiBuildBamFilterOptions(cmd_info);
    bool filter_requested = MpiHasBamFilterRequest(cmd_info);
    std::vector<long long> block_offsets;
    std::vector<long long> block_lengths;
    std::vector<long long> body_sizes;
    std::vector<long long> body_prefixes;
    long long body_start = 0;
    long long n_blocks = 0;
    long long local_body_size = 0;
    long long local_prefix = 0;
    long long total_body_size = 0;
    volatile unsigned long long simulated_write_guard = 0;
    int input_format = -1;
    int output_format = MpiOutputFormatFromName(cmd_info->out_file_name_);
    bool bam_to_bam = false;
    bool bam_to_sam = false;
    bool sam_to_bam = false;
    bool backend_bam_input =
        (output_format == bam || output_format == sam) &&
        MpiHasSuffix(cmd_info->in_file_name_, ".bam");
    bool streaming_sam_input =
        output_format == bam &&
        MpiHasSuffix(cmd_info->in_file_name_, ".sam") &&
        (cmd_info->io_backend_ == "posix" ||
         cmd_info->io_backend_ == "mpiio" ||
         cmd_info->io_backend_ == "auto");
    bool mpiio_output = cmd_info->io_output_backend_ == "mpiio";

    if (swbam::ParseByteSize(
            cmd_info->rank_body_memory_limit_,
            &rank_body_memory_limit) != 0 ||
        rank_body_memory_limit == 0) {
        if (rank == 0) {
            fprintf(stderr,
                    "ERROR: invalid rank body memory limit '%s'.\n",
                    cmd_info->rank_body_memory_limit_.c_str());
        }
        local_ok = 0;
    }

    if (cmd_info->validate_bounds_) {
        if (rank == 0) {
            fprintf(stderr, "ERROR: --validate-bounds is not supported by RabbitBAM-MPI. Use RabbitBAM-X for checked paths.\n");
        }
        local_ok = 0;
    }
    if (cmd_info->min_mapq_ > cmd_info->max_mapq_ && cmd_info->max_mapq_ >= 0) {
        if (rank == 0) fprintf(stderr, "ERROR: --min-mapq cannot be greater than --max-mapq.\n");
        local_ok = 0;
    }
    if (cmd_info->min_read_len_ > cmd_info->max_read_len_ && cmd_info->max_read_len_ >= 0) {
        if (rank == 0) fprintf(stderr, "ERROR: --min-read-len cannot be greater than --max-read-len.\n");
        local_ok = 0;
    }
    if (!MpiAllRanksOk(local_ok)) goto cleanup;

    init_cost = GetTime() - t_init;
    init_cost_max = MpiReduceMaxCost(init_cost);
    if (rank == 0) {
        // t_total += init_cost_max;
        printf("111Complete the initialization cost %lf-----\n", init_cost_max);
    }


    //2.把文件加载到所有rank的内存中，一共6份，每个rank一份，后续每个rank从内存中读取自己的部分进行处理；不计入处理时间
    {
        double preload_t0 = GetTime();
        if (backend_bam_input) {
            if (backend_input_handle.Open(
                    cmd_info->in_file_name_, cmd_info->io_backend_,
                    cmd_info->io_memory_limit_) != 0) {
                fprintf(stderr, "[rank %d] ERROR: cannot open BAM input backend for %s\n",
                        rank, cmd_info->in_file_name_.c_str());
                local_ok = 0;
            } else {
                backend_input = backend_input_handle.backend();
                input_file_size = backend_input->size();
                swbam::MemoryBamInput *memory_input =
                    dynamic_cast<swbam::MemoryBamInput *>(backend_input);
                if (memory_input) input_file_mem = memory_input->data();
            }
        } else if (streaming_sam_input) {
            sam_input = cmd_info->io_backend_ == "mpiio"
                ? static_cast<swbam::SamInputBackend *>(&mpiio_sam_input)
                : static_cast<swbam::SamInputBackend *>(&posix_sam_input);
            if (sam_input->Open(cmd_info->in_file_name_) != 0) {
                fprintf(stderr, "[rank %d] ERROR: cannot open streaming SAM input %s\n",
                        rank, cmd_info->in_file_name_.c_str());
                local_ok = 0;
            } else {
                input_file_size = sam_input->size();
            }
        } else if (MpiLoadFileToMemory(
                       cmd_info->in_file_name_, &input_file_mem,
                       &input_file_size) != 0) {
            fprintf(stderr, "[rank %d] ERROR: cannot preload input %s into memory\n",
                    rank, cmd_info->in_file_name_.c_str());
            local_ok = 0;
        }
        double preload_cost = backend_bam_input
            ? backend_input_handle.data_open_cost()
            : GetTime() - preload_t0;
        double preload_cost_max = MpiReduceMaxCost(preload_cost);
        if (rank == 0 && local_ok) {
            if (backend_bam_input) {
                printf("MPI BAM input backend=%s auto_memory_budget=%llu\n",
                       backend_input_handle.selected_backend().c_str(),
                       (unsigned long long)backend_input_handle.auto_memory_budget());
            } else if (streaming_sam_input) {
                printf("MPI SAM input backend=%s\n",
                       cmd_info->io_backend_ == "mpiio"
                           ? "mpiio" : "posix");
            }
            printf("222Complete the memory/input open cost %lf--\n", preload_cost_max);
        }
    }
    if (!MpiAllRanksOk(local_ok)) goto cleanup;


    //3.从内存创建 HTSlib 输入句柄并读取 BAM header，顺便判断输入输出格式是否合法，记录BAM/SAM body的起始位置
    {
        double header_t0 = GetTime();

        if (backend_bam_input) {
            if (!backend_input || backend_input->format() != bam ||
                !backend_input->header()) {
                local_ok = 0;
            } else {
                hdr = sam_hdr_dup(backend_input->header());
                if (!hdr) local_ok = 0;
                input_format = bam;
                bam_to_bam = output_format == bam;
                bam_to_sam = output_format == sam;
                if (!bam_to_bam && !bam_to_sam) local_ok = 0;
                if (filter_requested && !bam_to_bam) {
                    if (rank == 0) {
                        fprintf(stderr, "ERROR: BAM filtering options are only supported for MPI BAM -> BAM.\n");
                    }
                    local_ok = 0;
                }
                if (bam_to_sam && cmd_info->compress_level_ != 1 &&
                    rank == 0) {
                    printf("NOTE: --compress-level=%d is ignored for MPI BAM2SAM because SAM output is plain text.\n",
                           cmd_info->compress_level_);
                }
                body_start = (long long)backend_input->body_offset();
            }
        } else if (streaming_sam_input) {
            if (!sam_input || !sam_input->header()) {
                local_ok = 0;
            } else {
                hdr = sam_hdr_dup(sam_input->header());
                if (!hdr) local_ok = 0;
                input_format = sam;
                sam_to_bam = true;
                body_start = (long long)sam_input->body_offset();
            }
        } else {
            input_mem_hfile = hopen("mem:", "rb:", input_file_mem, input_file_size);
        }
        if (!input_mem_hfile) {
            if (!backend_bam_input && !streaming_sam_input) {
                fprintf(stderr, "[rank %d] ERROR: cannot open preloaded BAM memory for %s\n",
                        rank, cmd_info->in_file_name_.c_str());
                local_ok = 0;
            }
        }
        if (local_ok && !backend_bam_input && !streaming_sam_input) {
            sin = (samFile *)hts_hopen(input_mem_hfile, "data", "rb");
            if (!sin) {
                fprintf(stderr, "[rank %d] ERROR: cannot create HTS input handle from memory\n", rank);
                if (hclose(input_mem_hfile) != 0) {
                    fprintf(stderr, "[rank %d] ERROR: closing failed HTS memory handle failed.\n", rank);
                }
                input_mem_hfile = nullptr;
                input_file_mem = nullptr;
                input_file_size = 0;
                local_ok = 0;
            } else {
                input_mem_hfile = nullptr;
            }
        }

        if (local_ok && !backend_bam_input && !streaming_sam_input) {
            hdr = sam_hdr_read(sin);
            if (!hdr) {
                fprintf(stderr, "[rank %d] ERROR: cannot read header from %s\n",
                        rank, cmd_info->in_file_name_.c_str());
                local_ok = 0;
            }
        }
        if (local_ok && !backend_bam_input && !streaming_sam_input) {
            input_format = MpiNormalizeFormat(sin->format.format);
            bam_to_bam = input_format == bam && output_format == bam;
            bam_to_sam = input_format == bam && output_format == sam;
            sam_to_bam = input_format == sam && output_format == bam;

            if (!bam_to_bam && !bam_to_sam && !sam_to_bam) {
                if (rank == 0) {
                    fprintf(stderr, "Unsupported MPI conversion: input format %d -> output format %d\n",
                            input_format, output_format);
                }
                local_ok = 0;
            }
            if (filter_requested && !bam_to_bam) {
                if (rank == 0) {
                    fprintf(stderr, "ERROR: BAM filtering options are only supported for MPI BAM -> BAM.\n");
                }
                local_ok = 0;
            }
            if (local_ok && bam_to_sam && cmd_info->compress_level_ != 1 && rank == 0) {
                printf("NOTE: --compress-level=%d is ignored for MPI BAM2SAM because SAM output is plain text.\n",
                       cmd_info->compress_level_);
            }
            if (local_ok && bam_to_bam && !cmd_info->ref_name_.empty()) {
                int ref_tid = sam_hdr_name2tid(hdr, cmd_info->ref_name_.c_str());
                if (ref_tid < 0) {
                    fprintf(stderr, "[rank %d] ERROR: reference name '%s' does not exist in the BAM header.\n",
                            rank, cmd_info->ref_name_.c_str());
                    local_ok = 0;
                } else {
                    filter.ref_tid = ref_tid;
                }
            }
            if (local_ok) {
                if (input_format == bam) {
                    body_start = (long long)sin->fp.bgzf->block_address;
                    if (body_start < 0 ||
                        (unsigned long long)body_start > (unsigned long long)input_file_size) {
                        fprintf(stderr, "[rank %d] ERROR: invalid BAM body start offset %lld.\n", rank, body_start);
                        local_ok = 0;
                    }
                } else if (input_format == sam) {
                    if (MpiFindSamBodyStartInMemory(input_file_mem, input_file_size, &body_start) != 0) {
                        fprintf(stderr, "[rank %d] ERROR: failed to find SAM body start in memory.\n", rank);
                        local_ok = 0;
                    }
                }
            }
        }
        if (local_ok && backend_bam_input && bam_to_bam &&
            !cmd_info->ref_name_.empty()) {
            int ref_tid = sam_hdr_name2tid(hdr, cmd_info->ref_name_.c_str());
            if (ref_tid < 0) {
                fprintf(stderr, "[rank %d] ERROR: reference name '%s' does not exist in the BAM header.\n",
                        rank, cmd_info->ref_name_.c_str());
                local_ok = 0;
            } else {
                filter.ref_tid = ref_tid;
            }
        }

        double header_cost = backend_bam_input
            ? backend_input_handle.header_open_cost()
            : GetTime() - header_t0;
        double header_cost_max = MpiReduceMaxCost(header_cost);
        if (rank == 0 && local_ok) {
            // t_total += header_cost_max;
            printf("333Complete the head cost %lf---\n", header_cost_max);
        }
    }
    if (!MpiAllRanksOk(local_ok)) goto cleanup;


    //4.核心处理阶段
    {
        double body_total_t0 = GetTime();
        double body_t0 = GetTime();
        double stage41_t0 = GetTime();
        //4.1 rank 0 扫描 BGZF blocks
        if (rank == 0) {
            const char *mode_name = sam_to_bam ? "SAM2BAM" : (bam_to_sam ? "BAM2SAM" : "BAM2BAM");
            printf("Enable MPI %s mode (%d MPE + %d CPEs)!!!\n",
                   mode_name, comm_size, comm_size * 64);
            if (bam_to_bam && filter_requested) printf("Enable MPI BAM filtering options.\n");
            if (bam_to_bam || sam_to_bam) {
                printf("MPI BAM output compression level=%d\n", cmd_info->compress_level_);
            }
            printf("MPI output backend=%s\n",
                   mpiio_output ? "mpiio" : "memory");

            if (input_format == bam && !backend_bam_input) {
                if (MpiScanBgzfBlocksInMemory(input_file_mem, input_file_size, body_start,
                                              &block_offsets, &block_lengths) != 0) {
                    fprintf(stderr, "ERROR: failed to scan input BGZF blocks.\n");
                    local_ok = 0;
                }
                n_blocks = (long long)block_offsets.size();
                if (local_ok && n_blocks > (long long)INT_MAX) {
                    fprintf(stderr, "ERROR: too many BGZF blocks for MPI_Bcast in RabbitBAM-MPI v1.\n");
                    local_ok = 0;
                }
                if (local_ok) {
                    printf("MPI BAM scan complete. data_blocks=%lld body_start=%lld header_end=%lld\n",
                           n_blocks, body_start, (long long)sin->fp.bgzf->block_address);
                }
            } else if (input_format == sam) {
                int split_ret = streaming_sam_input
                    ? sam_input->SplitRanges(
                          comm_size, &block_offsets, &block_lengths)
                    : MpiSplitSamRangesInMemory(
                          input_file_mem, input_file_size, body_start,
                          comm_size, &block_offsets, &block_lengths);
                if (split_ret != 0) {
                    fprintf(stderr, "ERROR: failed to split SAM body ranges for MPI ranks.\n");
                    local_ok = 0;
                }
                n_blocks = (long long)block_offsets.size();
                if (local_ok) {
                    printf("MPI SAM split complete. rank_ranges=%lld body_start=%lld file_size=%zu\n",
                           n_blocks, body_start, input_file_size);
                }
            }
        }

        if (backend_bam_input) {
            if (swbam::mpi::PrepareMpiBamInputPlan(
                    *backend_input, &backend_plan) != 0) {
                fprintf(stderr, "[rank %d] ERROR: failed to prepare MPI BAM input plan.\n", rank);
                local_ok = 0;
            } else {
                body_start = (long long)backend_plan.body_offset;
                n_blocks = (long long)backend_plan.blocks.size();
                if (rank == 0) {
                    printf("MPI BAM scan complete. data_blocks=%lld body_start=%lld header_end=%lld backend=%s\n",
                           n_blocks, body_start, body_start,
                           backend_input_handle.selected_backend().c_str());
                }
            }
        } else {
            // 把 block 信息广播给所有 rank。
            MPI_Bcast(&local_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (!local_ok) goto cleanup;
            MPI_Bcast(&body_start, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(&n_blocks, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
            if (rank != 0) {
                block_offsets.resize((size_t)n_blocks);
                block_lengths.resize((size_t)n_blocks);
            }
            if (n_blocks > 0) {
                MPI_Bcast(block_offsets.data(), (int)n_blocks, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
                MPI_Bcast(block_lengths.data(), (int)n_blocks, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
            }
        }
        if (!MpiAllRanksOk(local_ok)) goto cleanup;

        // 按 rank 切分 block 范围。memory 保留连续窗口快速路径，
        // POSIX 只计算当前 rank 的压缩字节数，4.3 内按 batch 读取。
        size_t rank_input_size = 0;
        if (backend_bam_input) {
            for (size_t i = backend_plan.rank_begin;
                 i < backend_plan.rank_end; ++i) {
                const size_t block_size =
                    (size_t)backend_plan.blocks[i].compressed_size;
                if (block_size > SIZE_MAX - rank_input_size) {
                    fprintf(stderr, "[rank %d] ERROR: assigned BAM byte count overflows size_t.\n", rank);
                    local_ok = 0;
                    break;
                }
                rank_input_size += block_size;
            }
            if (local_ok && input_file_mem &&
                backend_plan.rank_begin < backend_plan.rank_end) {
                const swbam::BgzfBlockSpan &first =
                    backend_plan.blocks[backend_plan.rank_begin];
                const swbam::BgzfBlockSpan &last =
                    backend_plan.blocks[backend_plan.rank_end - 1];
                const uint64_t range_end =
                    last.offset + (uint64_t)last.compressed_size;
                if (range_end < first.offset ||
                    range_end > (uint64_t)input_file_size ||
                    range_end - first.offset != (uint64_t)rank_input_size) {
                    fprintf(stderr, "[rank %d] ERROR: memory backend assigned BAM blocks are not contiguous.\n", rank);
                    local_ok = 0;
                } else {
                    rank_input_mem = input_file_mem + (size_t)first.offset;
                }
            }
        } else if (input_format == bam) {
            long long begin = n_blocks * rank / comm_size;
            long long end = n_blocks * (rank + 1) / comm_size;
            if (MpiSelectBlockRangeFromMemory(input_file_mem, input_file_size,
                                              block_offsets, block_lengths,
                                              begin, end, &rank_input_mem, &rank_input_size) != 0) {
                fprintf(stderr, "[rank %d] ERROR: failed to select assigned BGZF block range [%lld, %lld) from memory.\n",
                        rank, begin, end);
                local_ok = 0;
            }
        } else {
            if (rank >= (int)block_offsets.size()) {
                fprintf(stderr, "[rank %d] ERROR: missing assigned SAM range.\n", rank);
                local_ok = 0;
            } else {
                long long start = block_offsets[(size_t)rank];
                long long len = block_lengths[(size_t)rank];
                if (start < 0 || len < 0 ||
                    (unsigned long long)start + (unsigned long long)len > (unsigned long long)input_file_size) {
                    fprintf(stderr, "[rank %d] ERROR: invalid assigned SAM range start=%lld len=%lld.\n",
                            rank, start, len);
                    local_ok = 0;
                } else {
                    if (streaming_sam_input) {
                        if (sam_input->ReadRange(
                                (uint64_t)start, (size_t)len,
                                &rank_sam_storage) != 0) {
                            fprintf(stderr, "[rank %d] ERROR: failed to read assigned POSIX SAM range.\n", rank);
                            local_ok = 0;
                        } else {
                            rank_input_mem = rank_sam_storage.empty()
                                ? nullptr : rank_sam_storage.data();
                            rank_input_size = rank_sam_storage.size();
                        }
                    } else {
                        rank_input_mem = input_file_mem + start;
                        rank_input_size = (size_t)len;
                    }
                }
            }
        }
        double stage41_cost = GetTime() - stage41_t0;
        double stage41_cost_max = MpiReduceMaxCost(stage41_cost);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.1 scan/split/broadcast/select cost %lf\n", stage41_cost_max);
        }

        //4.2 初始化内存 reader 和 writer
        double stage42_t0 = GetTime();
        if (local_ok) {
            if (!backend_bam_input || input_file_mem) {
                reader.base = rank_input_mem;
                reader.size = rank_input_size;
                reader.pos = 0;
            }
            size_t writer_capacity = rank_input_size ? rank_input_size : 64 * 1024 * 1024;
            if (bam_to_sam && rank_input_size > 0) {
                if (rank_input_size > SIZE_MAX / 5) {
                    fprintf(stderr, "[rank %d] ERROR: MPI BAM2SAM memory writer estimate overflow.\n", rank);
                    local_ok = 0;
                } else {
                    writer_capacity = rank_input_size * 5;
                }
            } else if (sam_to_bam && rank_input_size > 0) {
                writer_capacity = rank_input_size / 2;
                if (writer_capacity == 0) writer_capacity = 64 * 1024 * 1024;
            }
            char spool_prefix[96];
            snprintf(spool_prefix, sizeof(spool_prefix),
                     "rabbitbam-convert-rank-%d.body", rank);
            if (local_ok && body_sink.Open(
                    cmd_info->rank_body_backend_,
                    rank_body_memory_limit, writer_capacity,
                    spool_prefix,
                    cmd_info->rank_body_temp_dir_) != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: failed to open rank body sink.\n",
                        rank);
                local_ok = 0;
            }
        }
        int stage42_ok = MpiAllRanksOk(local_ok);
        double stage42_cost = GetTime() - stage42_t0;
        double stage42_cost_max = MpiReduceMaxCost(stage42_cost);
        if (rank == 0 && stage42_ok) {
            printf("Complete the 4.2 init reader/writer cost %lf\n", stage42_cost_max);
        }
        if (!stage42_ok) goto cleanup;

        //4.3 核心处理：FusedBamToBamMPI / FusedBamToSamMPI / FusedSamToBamMPI
        const char *fused_name = sam_to_bam ? "FusedSamToBamMPI" :
                                 (bam_to_sam ? "FusedBamToSamMPI" : "FusedBamToBamMPI");
        double FusedBamToBamMPI_t0 = GetTime();
        if (bam_to_sam) {
            int ret = backend_bam_input && !input_file_mem
                ? FusedBamToSamMPI(
                      *backend_input, backend_plan.rank_spans(),
                      backend_plan.rank_block_count(), body_sink,
                      hdr, &sam_stats)
                : FusedBamToSamMPI(
                      reader, body_sink, hdr, &sam_stats);
            if (ret != 0) {
                fprintf(stderr, "[rank %d] ERROR: MPI bam2sam fused 1CG body failed.\n", rank);
                local_ok = 0;
            }
        } else if (sam_to_bam) {
            if (FusedSamToBamMPI(
                    reader, body_sink, hdr,
                    cmd_info->compress_level_,
                    &sam2bam_stats) != 0) {
                fprintf(stderr, "[rank %d] ERROR: MPI sam2bam fused 1CG body failed.\n", rank);
                local_ok = 0;
            }
        } else {
            int ret = backend_bam_input && !input_file_mem
                ? FusedBamToBamMPI(
                      *backend_input, backend_plan.rank_spans(),
                      backend_plan.rank_block_count(), body_sink,
                      filter, cmd_info->compress_level_, &stats)
                : FusedBamToBamMPI(
                      reader, body_sink, filter,
                      cmd_info->compress_level_, &stats);
            if (ret != 0) {
                fprintf(stderr, "[rank %d] ERROR: MPI bam2bam fused 1CG body failed.\n", rank);
                local_ok = 0;
            }
        }

        int global_ok = MpiAllRanksOk(local_ok);
        double FusedBamToBamMPI_cost = GetTime() - FusedBamToBamMPI_t0;
        double FusedBamToBamMPI_cost_max = MpiReduceMaxCost(FusedBamToBamMPI_cost);
        if (rank == 0 && global_ok) {
            printf("Complete the 4.3 %s cost %lf\n", fused_name, FusedBamToBamMPI_cost_max);
        }

        if (!global_ok) goto cleanup;

        {
            const int local_spool = body_sink.is_memory() ? 0 : 1;
            int spool_ranks = 0;
            MPI_Reduce(&local_spool, &spool_ranks, 1, MPI_INT,
                       MPI_SUM, 0, MPI_COMM_WORLD);
            if (rank == 0) {
                printf("MPI rank body backend requested=%s "
                       "memory_ranks=%d spool_ranks=%d\n",
                       cmd_info->rank_body_backend_.c_str(),
                       comm_size - spool_ranks, spool_ranks);
            }
        }

        //4.4 收集每个 rank 的输出大小
        double stage44_t0 = GetTime();
        if (body_sink.Flush() != 0 ||
            body_sink.size() > (uint64_t)LLONG_MAX) {
            local_ok = 0;
        }
        if (!MpiAllRanksOk(local_ok)) goto cleanup;
        local_body_size = (long long)body_sink.size();
        body_sizes.assign((size_t)comm_size, 0);
        body_prefixes.assign((size_t)comm_size, 0);
        MPI_Allgather(&local_body_size, 1, MPI_LONG_LONG,
                      body_sizes.data(), 1, MPI_LONG_LONG, MPI_COMM_WORLD);
        local_prefix = 0;
        total_body_size = 0;
        for (int i = 0; i < comm_size; ++i) {
            body_prefixes[(size_t)i] = total_body_size;
            if (i < rank) local_prefix += body_sizes[(size_t)i];
            total_body_size += body_sizes[(size_t)i];
        }
        double stage44_cost = GetTime() - stage44_t0;
        double stage44_cost_max = MpiReduceMaxCost(stage44_cost);
        if (rank == 0) {
            printf("Complete the 4.4 gather body sizes cost %lf\n", stage44_cost_max);
        }

        //4.5a 处理输出 header 和整体布局，这部分计入处理时间
        double stage45_header_t0 = GetTime();
        long long output_body_start_ll = 0;
        long long output_file_size_ll = 0;
        if (rank == 0) {
            unsigned long long final_size = 0;
            sam_header_text = nullptr;
            sam_header_len = 0;
            if (bam_to_sam) {
                sam_header_text = sam_hdr_str(hdr);
                sam_header_len = sam_hdr_length(hdr);
                if (sam_header_len > 0 && !sam_header_text) {
                    fprintf(stderr, "ERROR: failed to get SAM header text for MPI BAM2SAM output.\n");
                    local_ok = 0;
                }
                output_body_start = sam_header_len;
                final_size = (unsigned long long)output_body_start +
                             (unsigned long long)total_body_size;
            } else if (sam_to_bam) {
                if (MpiBuildBamHeaderMemory(hdr, cmd_info->compress_level_, &bam_header_mem, &bam_header_size) != 0) {
                    fprintf(stderr, "ERROR: failed to build MPI SAM2BAM output BAM header in memory.\n");
                    local_ok = 0;
                }
                output_body_start = bam_header_size;
                final_size = (unsigned long long)output_body_start +
                             (unsigned long long)total_body_size +
                             (unsigned long long)sizeof(kMpiBgzfEofBlock);
            } else {
                if (backend_bam_input && !input_file_mem) {
                    if (MpiBuildBamHeaderMemory(
                            hdr, cmd_info->compress_level_,
                            &bam_header_mem, &bam_header_size) != 0) {
                        fprintf(stderr, "ERROR: failed to build MPI BAM2BAM output header in memory.\n");
                        local_ok = 0;
                    }
                    output_body_start = bam_header_size;
                } else {
                    output_body_start = (size_t)body_start;
                }
                final_size = (unsigned long long)output_body_start +
                             (unsigned long long)total_body_size +
                             (unsigned long long)sizeof(kMpiBgzfEofBlock);
            }
            if (final_size > (unsigned long long)SIZE_MAX ||
                final_size > (unsigned long long)LLONG_MAX ||
                (unsigned long long)output_body_start > (unsigned long long)LLONG_MAX) {
                fprintf(stderr, "ERROR: MPI output buffer would exceed addressable memory.\n");
                local_ok = 0;
            } else {
                output_file_size = (size_t)final_size;
                output_body_start_ll = (long long)output_body_start;
                output_file_size_ll = (long long)output_file_size;
            }
        }

        int stage45_header_ok = MpiAllRanksOk(local_ok);
        MPI_Bcast(&output_body_start_ll, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&output_file_size_ll, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        if (stage45_header_ok) {
            output_body_start = (size_t)output_body_start_ll;
            output_file_size = (size_t)output_file_size_ll;
        }
        double stage45_header_cost = GetTime() - stage45_header_t0;
        double stage45_header_cost_max = MpiReduceMaxCost(stage45_header_cost);
        if (rank == 0 && stage45_header_ok) {
            printf("Complete the 4.5a header/layout cost %lf\n", stage45_header_cost_max);
        }
        if (!stage45_header_ok) goto cleanup;

        //4.5b 分配每个 rank 的本地模拟写内存，这部分不计入处理时间
        double stage45_malloc_t0 = GetTime();
        simulated_write_size = body_sink.is_memory() &&
            local_body_size > 0 ? (size_t)local_body_size : 0;
        if (body_sink.is_memory() && rank == 0) {
            if (output_body_start > SIZE_MAX - simulated_write_size) {
                fprintf(stderr, "[rank %d] ERROR: simulated output memory estimate overflow.\n", rank);
                local_ok = 0;
            } else {
                simulated_write_size += output_body_start;
            }
        }
        if (local_ok && simulated_write_size > 0) {
            simulated_write_mem = (char *)malloc(simulated_write_size);
            if (!simulated_write_mem) {
                fprintf(stderr, "[rank %d] ERROR: failed to allocate simulated output memory. size=%zu\n",
                        rank, simulated_write_size);
                local_ok = 0;
            }
        }
        int stage45_malloc_ok = MpiAllRanksOk(local_ok);
        double stage45_malloc_cost = GetTime() - stage45_malloc_t0;
        double stage45_malloc_cost_max = MpiReduceMaxCost(stage45_malloc_cost);
        if (rank == 0 && stage45_malloc_ok) {
            printf("Complete the 4.5b malloc simulated write memory cost %lf\n", stage45_malloc_cost_max);
        }
        if (!stage45_malloc_ok) goto cleanup;

        //4.6 各 rank 根据逻辑偏移模拟写入高速磁盘内存，只统计 memcpy 时间
        double sim_write_t0 = GetTime();
        long long simulated_write_offset = (long long)output_body_start + local_prefix;
        (void)simulated_write_offset;
        size_t simulated_pos = 0;
        if (body_sink.is_memory() &&
            rank == 0 && output_body_start > 0) {
            if (bam_to_sam) {
                if (sam_header_len > 0) {
                    memcpy(simulated_write_mem + simulated_pos, sam_header_text, sam_header_len);
                    simulated_pos += sam_header_len;
                }
            } else if (sam_to_bam) {
                if (bam_header_size > 0) {
                    memcpy(simulated_write_mem + simulated_pos, bam_header_mem, bam_header_size);
                    simulated_pos += bam_header_size;
                }
            } else {
                const char *header_data =
                    backend_bam_input && !input_file_mem
                        ? bam_header_mem : input_file_mem;
                if (output_body_start > 0) {
                    memcpy(simulated_write_mem + simulated_pos,
                           header_data, output_body_start);
                    simulated_pos += output_body_start;
                }
            }
        }
        if (body_sink.is_memory() && local_body_size > 0) {
            if (body_sink.ReadAt(
                    0, simulated_write_mem + simulated_pos,
                    (size_t)local_body_size) != 0) {
                local_ok = 0;
            } else {
                simulated_pos += (size_t)local_body_size;
            }
        }
        if (simulated_pos > 0) {
            unsigned char *guard_ptr = (unsigned char *)simulated_write_mem;
            simulated_write_guard += guard_ptr[0];
            simulated_write_guard += guard_ptr[simulated_pos - 1];
        }
        double sim_write_cost = GetTime() - sim_write_t0;
        if (bam_to_sam) {
            sam_stats.t_gather += sim_write_cost;
        } else if (sam_to_bam) {
            sam2bam_stats.t_gather += sim_write_cost;
        } else {
            stats.t_mpi_write += sim_write_cost;
        }

        global_ok = MpiAllRanksOk(local_ok);
        double stage46_cost_max = MpiReduceMaxCost(sim_write_cost);
        if (rank == 0 && global_ok) {
            printf("Complete the 4.6 simulated header/body distributed write cost %lf\n", stage46_cost_max);
        }
        if (simulated_write_guard == (unsigned long long)-1 && rank < 0) {
            fprintf(stderr, "unused simulated write guard %llu\n", simulated_write_guard);
        }
        double body_cost = GetTime() - body_t0;
        double body_cost_max = MpiReduceMaxCost(body_cost);
        if (rank == 0 && global_ok) {
            double body_counted_cost = stage41_cost_max + stage42_cost_max +
                                       FusedBamToBamMPI_cost_max + stage44_cost_max +
                                       stage45_header_cost_max + stage46_cost_max;
            t_total += body_counted_cost;
            // printf("Complete the counted body cost %lf\n", body_counted_cost);
            printf("Complete the total (4.1~4.6) cost %lf-----\n", t_total);
        }
        if (!global_ok) goto cleanup;


        //4.7 全局同步，统计每个 rank 的处理时间和统计数据，rank 0 汇总并打印最终统计结果
        double stage47_t0 = GetTime();
        if (bam_to_sam) {
            long long local_long_stats[4] = {
                sam_stats.input_blocks,
                sam_stats.group_count,
                sam_stats.total_records,
                sam_stats.format_tiles
            };
            long long global_long_stats[4] = {};
            MPI_Reduce(local_long_stats, global_long_stats, 4, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            double local_double_stats[15] = {
                sam_stats.t_decomp,
                sam_stats.t_decomp_alloc,
                sam_stats.t_decomp_inflate,
                sam_stats.t_decomp_crc,
                sam_stats.t_decomp_parse,
                sam_stats.t_decomp_other,
                sam_stats.t_format,
                sam_stats.t_collect,
                sam_stats.t_read,
                sam_stats.t_write,
                sam_stats.t_gather,
                sam_stats.t_fused_total,
                sam_stats.t_decomp + sam_stats.t_collect + sam_stats.t_format,
                sam_stats.t_alloc_init,
                sam_stats.t_free_workspace
            };
            double global_double_stats[15] = {};
            MPI_Reduce(local_double_stats, global_double_stats, 15, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            MpiPrintRankBamToSamStats(
                rank, comm_size, sam_stats, (size_t)body_sink.size());
            if (rank == 0) {
                printf("FusedBamToSamMPI finished. ranks=%d in_blocks=%lld groups=%lld total_records=%lld format_tiles=%lld body_bytes=%lld\n",
                       comm_size, global_long_stats[0], global_long_stats[1],
                       global_long_stats[2], global_long_stats[3], total_body_size);
                printf("  decomp_slave_sum=%.3f  format_slave_sum=%.3f  collect_sum=%.3f  read_sum=%.3f  write_sum=%.3f  gather_sum=%.3f\n",
                       global_double_stats[0], global_double_stats[6], global_double_stats[7],
                       global_double_stats[8], global_double_stats[9], global_double_stats[10]);
                printf("  decomp_detail_sum alloc=%.3f  inflate=%.3f  crc=%.3f  parse=%.3f  other=%.3f\n",
                       global_double_stats[1], global_double_stats[2], global_double_stats[3],
                       global_double_stats[4], global_double_stats[5]);
                printf("  fused_total_sum=%.3f  core_stage_sum=%.3f  alloc_init_sum=%.3f  free_workspace_sum=%.3f\n",
                       global_double_stats[11], global_double_stats[12],
                       global_double_stats[13], global_double_stats[14]);
            }
        } else if (sam_to_bam) {
            long long local_long_stats[7] = {
                sam2bam_stats.input_chunks,
                sam2bam_stats.chunk_groups,
                sam2bam_stats.total_records,
                sam2bam_stats.compress_groups,
                sam2bam_stats.bgzf_blocks,
                sam2bam_stats.parse_fast_records,
                sam2bam_stats.parse_fallback_records
            };
            long long global_long_stats[7] = {};
            MPI_Reduce(local_long_stats, global_long_stats, 7, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            double local_double_stats[22] = {
                sam2bam_stats.t_split,
                sam2bam_stats.t_copy_count,
                sam2bam_stats.t_parse,
                sam2bam_stats.t_parse_core,
                sam2bam_stats.t_parse_aux,
                sam2bam_stats.t_parse_cg,
                sam2bam_stats.t_parse_fallback,
                sam2bam_stats.t_parse_other,
                sam2bam_stats.t_pack,
                sam2bam_stats.t_compress,
                sam2bam_stats.t_compress_serialize,
                sam2bam_stats.t_compress_alloc,
                sam2bam_stats.t_compress_deflate,
                sam2bam_stats.t_compress_footer,
                sam2bam_stats.t_compress_other,
                sam2bam_stats.t_write,
                sam2bam_stats.t_gather,
                sam2bam_stats.t_fused_total,
                sam2bam_stats.t_setup_reset,
                sam2bam_stats.t_compress_setup,
                sam2bam_stats.t_status_check,
                sam2bam_stats.t_alloc_init + sam2bam_stats.t_free_workspace
            };
            double global_double_stats[22] = {};
            MPI_Reduce(local_double_stats, global_double_stats, 22, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


            MpiPrintRankSamToBamStats(
                rank, comm_size, sam2bam_stats, (size_t)body_sink.size());
            if (rank == 0) {
                printf("FusedSamToBamMPI finished. ranks=%d chunks=%lld chunk_groups=%lld total_records=%lld compress_groups=%lld bgzf_blocks=%lld body_bytes=%lld\n",
                       comm_size, global_long_stats[0], global_long_stats[1],
                       global_long_stats[2], global_long_stats[3], global_long_stats[4],
                       total_body_size);
                printf("  parse_fast_records=%lld  parse_fallback_records=%lld\n",
                       global_long_stats[5], global_long_stats[6]);
                printf("  split_sum=%.3f  copy_count_slave_sum=%.3f  parse_slave_sum=%.3f  pack_sum=%.3f  compress_slave_sum=%.3f  write_sum=%.3f  gather_sum=%.3f\n",
                       global_double_stats[0], global_double_stats[1], global_double_stats[2],
                       global_double_stats[8], global_double_stats[9], global_double_stats[15],
                       global_double_stats[16]);
                printf("  parse_detail_sum core=%.3f  aux=%.3f  cg=%.3f  fallback=%.3f  other=%.3f\n",
                       global_double_stats[3], global_double_stats[4], global_double_stats[5],
                       global_double_stats[6], global_double_stats[7]);
                printf("  compress_detail_sum serialize=%.3f  alloc=%.3f  deflate=%.3f  footer=%.3f  other=%.3f\n",
                       global_double_stats[10], global_double_stats[11], global_double_stats[12],
                       global_double_stats[13], global_double_stats[14]);
                printf("  fused_total_sum=%.3f  core_stage_sum=%.3f  setup_reset_sum=%.3f  compress_setup_sum=%.3f  status_check_sum=%.3f  alloc_free_sum=%.3f\n",
                       global_double_stats[17],
                       global_double_stats[1] + global_double_stats[2] +
                       global_double_stats[8] + global_double_stats[9],
                       global_double_stats[18], global_double_stats[19],
                       global_double_stats[20], global_double_stats[21]);
            }
        } else {
            long long local_long_stats[7] = {
                stats.input_blocks,
                stats.group_count,
                stats.total_records,
                stats.kept_records,
                stats.dropped_records,
                stats.bgzf_blocks,
                stats.pack_records
            };
            long long global_long_stats[7] = {};
            MPI_Reduce(local_long_stats, global_long_stats, 7, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            double local_double_stats[16] = {
                stats.t_decomp_filter,
                stats.t_decomp_alloc,
                stats.t_decomp_inflate,
                stats.t_decomp_crc,
                stats.t_decomp_parse,
                stats.t_decomp_other,
                stats.t_pack,
                stats.t_compress,
                stats.t_compress_serialize,
                stats.t_compress_alloc,
                stats.t_compress_deflate,
                stats.t_compress_footer,
                stats.t_compress_other,
                stats.t_read,
                stats.t_write,
                stats.t_mpi_write
            };
            double global_double_stats[16] = {};
            MPI_Reduce(local_double_stats, global_double_stats, 16, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            double local_detail_stats[3] = {
                stats.t_fused_total,
                stats.t_decomp_filter + stats.t_pack + stats.t_compress,
                stats.t_prepare_decomp,
            };
            double global_detail_stats[3] = {};
            MPI_Reduce(local_detail_stats, global_detail_stats, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            MpiPrintRankStats(
                rank, comm_size, stats, (size_t)body_sink.size());
            if (rank == 0) {
                double keep_ratio = global_long_stats[2] > 0
                                    ? (double)global_long_stats[3] / (double)global_long_stats[2]
                                    : 1.0;
                printf("FusedBamToBamMPI finished. ranks=%d in_blocks=%lld groups=%lld total_records=%lld kept_records=%lld dropped_records=%lld bgzf_blocks=%lld body_bytes=%lld\n",
                       comm_size, global_long_stats[0], global_long_stats[1], global_long_stats[2],
                       global_long_stats[3], global_long_stats[4], global_long_stats[5],
                       total_body_size);
                printf("  decomp_filter_slave_sum=%.3f  pack_sum=%.3f  compress_slave_sum=%.3f  read_sum=%.3f  write_sum=%.3f  mpi_write_sum=%.3f\n",
                       global_double_stats[0], global_double_stats[6], global_double_stats[7],
                       global_double_stats[13], global_double_stats[14], global_double_stats[15]);
                printf("  decomp_detail_sum alloc=%.3f  inflate=%.3f  crc=%.3f  parse_filter=%.3f  other=%.3f\n",
                       global_double_stats[1], global_double_stats[2], global_double_stats[3],
                       global_double_stats[4], global_double_stats[5]);
                printf("  compress_detail_sum serialize=%.3f  alloc=%.3f  deflate=%.3f  footer=%.3f  other=%.3f\n",
                       global_double_stats[8], global_double_stats[9], global_double_stats[10],
                       global_double_stats[11], global_double_stats[12]);
                printf("  fused_total_sum=%.3f  core_stage_sum=%.3f  prepare_decomp_sum=%.3f\n",
                       global_detail_stats[0], global_detail_stats[1], global_detail_stats[2]);
                printf("  pack_records=%lld  keep_ratio=%.6f  filter_mode=%s\n",
                       global_long_stats[6], keep_ratio, bam_filter_is_noop(filter) ? "passthrough" : "filtered");
            }
        }
        double stage47_cost = GetTime() - stage47_t0;
        double stage47_cost_max = MpiReduceMaxCost(stage47_cost);
        if (rank == 0 && global_ok) {
            printf("Complete the 4.7 rank stats reduce/print cost %lf\n", stage47_cost_max);
        }

        double body_total_cost = GetTime() - body_total_t0;
        double body_total_cost_max = MpiReduceMaxCost(body_total_cost);
        if (rank == 0 && global_ok) {
            printf("444Complete the total body cost %lf\n", body_total_cost_max);
        }
    }

    //5. 输出不计入核心处理时间。MPI-IO 路径跳过 rank0 全量 gather。
    if (mpiio_output) {
        double mpiio_t0 = GetTime();
        swbam::mpi::MpiFileOutput output;
        if (output.Open(cmd_info->out_file_name_,
                        (uint64_t)output_file_size) != 0) {
            fprintf(stderr, "[rank %d] ERROR: failed to open MPI-IO output %s.\n",
                    rank, cmd_info->out_file_name_.c_str());
            local_ok = 0;
        }
        int mpiio_open_ok = MpiAllRanksOk(local_ok);
        if (mpiio_open_ok) {
            if (rank == 0 && output_body_start > 0) {
                const char *header_data = nullptr;
                if (bam_to_sam) {
                    header_data = sam_header_text;
                } else if (sam_to_bam ||
                           (backend_bam_input && !input_file_mem)) {
                    header_data = bam_header_mem;
                } else {
                    header_data = input_file_mem;
                }
                if (output.WriteAt(0, header_data,
                                   (uint64_t)output_body_start) != 0) {
                    fprintf(stderr, "ERROR: failed to write MPI-IO output header.\n");
                    local_ok = 0;
                }
            }
            const size_t chunk_capacity = 8u * 1024u * 1024u;
            std::vector<unsigned char> chunk;
            if (local_body_size > 0) {
                try {
                    chunk.resize((size_t)std::min<long long>(
                        local_body_size, (long long)chunk_capacity));
                } catch (...) {
                    local_ok = 0;
                }
            }
            uint64_t copied = 0;
            while (local_ok && copied < (uint64_t)local_body_size) {
                const size_t count = (size_t)std::min<uint64_t>(
                    (uint64_t)local_body_size - copied,
                    chunk.size());
                if (body_sink.ReadAt(
                        copied, chunk.data(), count) != 0 ||
                    output.WriteAt(
                        (uint64_t)output_body_start +
                            (uint64_t)local_prefix + copied,
                        chunk.data(), count) != 0) {
                    fprintf(stderr,
                            "[rank %d] ERROR: failed to write MPI-IO "
                            "output body.\n", rank);
                    local_ok = 0;
                }
                copied += count;
            }
            if (rank == 0 && !bam_to_sam &&
                output.WriteAt(
                    (uint64_t)output_body_start + (uint64_t)total_body_size,
                    kMpiBgzfEofBlock,
                    (uint64_t)sizeof(kMpiBgzfEofBlock)) != 0) {
                fprintf(stderr, "ERROR: failed to write MPI-IO BGZF EOF.\n");
                local_ok = 0;
            }
        }
        int mpiio_write_ok = MpiAllRanksOk(local_ok);
        if (mpiio_write_ok && output.Sync() != 0) {
            fprintf(stderr, "[rank %d] ERROR: failed to sync MPI-IO output.\n", rank);
            local_ok = 0;
        }
        MpiAllRanksOk(local_ok);
        if (output.Close() != 0) {
            fprintf(stderr, "[rank %d] ERROR: failed to close MPI-IO output.\n", rank);
            local_ok = 0;
        }
        int mpiio_close_ok = MpiAllRanksOk(local_ok);
        double mpiio_cost = GetTime() - mpiio_t0;
        double mpiio_cost_max = MpiReduceMaxCost(mpiio_cost);
        if (rank == 0 && mpiio_close_ok) {
            printf("555Distributed MPI-IO output cost %lf--\n", mpiio_cost_max);
        }
        if (!mpiio_open_ok || !mpiio_write_ok || !mpiio_close_ok) goto cleanup;
    } else {
        double verify_alloc_t0 = GetTime();
        if (rank == 0) {
            output_file_mem = output_file_size ? (char *)malloc(output_file_size) : nullptr;
            if (output_file_size > 0 && !output_file_mem) {
                fprintf(stderr, "ERROR: failed to allocate final MPI output buffer.\n");
                local_ok = 0;
            }
            if (local_ok) {
                if (bam_to_sam) {
                    if (sam_header_len > 0) memcpy(output_file_mem, sam_header_text, sam_header_len);
                } else if (sam_to_bam) {
                    if (bam_header_size > 0) memcpy(output_file_mem, bam_header_mem, bam_header_size);
                } else {
                    const char *header_data =
                        backend_bam_input && !input_file_mem
                            ? bam_header_mem : input_file_mem;
                    if (output_body_start > 0) {
                        memcpy(output_file_mem, header_data,
                               output_body_start);
                    }
                }
            }
        }
        double verify_alloc_cost = GetTime() - verify_alloc_t0;
        double verify_alloc_cost_max = MpiReduceMaxCost(verify_alloc_cost);
        if (rank == 0 && local_ok) {
            printf("555Prepare verification output memory cost %lf--\n", verify_alloc_cost_max);
        }
        if (!MpiAllRanksOk(local_ok)) goto cleanup;

        double verify_gather_t0 = GetTime();
        if (rank == 0) {
            if (local_body_size > 0) {
                if (body_sink.ReadAt(
                        0,
                        output_file_mem + output_body_start + local_prefix,
                        (size_t)local_body_size) != 0) {
                    local_ok = 0;
                }
            }
            for (int src = 1; src < comm_size; ++src) {
                long long recv_size = body_sizes[(size_t)src];
                if (recv_size <= 0) continue;
                long long received = 0;
                while (received < recv_size) {
                    const long long count = std::min<long long>(
                        recv_size - received, 8ll * 1024ll * 1024ll);
                    if (MpiRecvBytes(
                            src, 0,
                            output_file_mem + output_body_start +
                                body_prefixes[(size_t)src] + received,
                            count) != 0) {
                        fprintf(stderr,
                                "ERROR: failed to gather MPI rank %d "
                                "body into output memory.\n", src);
                        local_ok = 0;
                        break;
                    }
                    received += count;
                }
                if (!local_ok) break;
            }
            if (local_ok && !bam_to_sam) {
                memcpy(output_file_mem + output_body_start + total_body_size,
                       kMpiBgzfEofBlock, sizeof(kMpiBgzfEofBlock));
            }
        } else if (local_body_size > 0) {
            std::vector<unsigned char> chunk;
            try {
                chunk.resize((size_t)std::min<long long>(
                    local_body_size, 8ll * 1024ll * 1024ll));
            } catch (...) {
                local_ok = 0;
            }
            uint64_t sent = 0;
            while (local_ok && sent < (uint64_t)local_body_size) {
                const size_t count = (size_t)std::min<uint64_t>(
                    (uint64_t)local_body_size - sent,
                    chunk.size());
                if (body_sink.ReadAt(
                        sent, chunk.data(), count) != 0 ||
                    MpiSendBytes(
                        0, 0,
                        reinterpret_cast<const char *>(chunk.data()),
                        count) != 0) {
                    fprintf(stderr,
                            "[rank %d] ERROR: failed to send output "
                            "body to rank 0 for verification dump.\n",
                            rank);
                    local_ok = 0;
                }
                sent += count;
            }
        }
        double verify_gather_cost = GetTime() - verify_gather_t0;
        double verify_gather_cost_max = MpiReduceMaxCost(verify_gather_cost);
        if (rank == 0 && local_ok) {
            printf("555Gather verification output memory cost %lf--\n", verify_gather_cost_max);
        }
        if (!MpiAllRanksOk(local_ok)) goto cleanup;

        double dump_t0 = 0.0;
        double dump_cost = 0.0;
        if (rank == 0) {
            dump_t0 = GetTime();
            if (MpiDumpMemoryToFile(cmd_info->out_file_name_, output_file_mem, output_file_size) != 0) {
                fprintf(stderr, "ERROR: failed to dump MPI output memory to %s\n",
                        cmd_info->out_file_name_.c_str());
                local_ok = 0;
            }
            dump_cost = GetTime() - dump_t0;
        }
        double dump_cost_max = MpiReduceMaxCost(dump_cost);
        if (rank == 0 && local_ok) printf("555Dump memory to output file cost %lf--\n", dump_cost_max);
        if (!MpiAllRanksOk(local_ok)) goto cleanup;
    }

    exit_code = 0;

    //6. 清理资源
cleanup:
    {
        double close_t0 = GetTime();
        body_sink.Close();
        if (output_file_mem) free(output_file_mem);
        if (bam_header_mem) free(bam_header_mem);
        if (simulated_write_mem) free(simulated_write_mem);
        if (hdr) sam_hdr_destroy(hdr);
        if (backend_bam_input) {
            input_file_mem = nullptr;
            backend_input = nullptr;
            backend_input_handle.Close();
        } else if (streaming_sam_input) {
            if (sam_input) sam_input->Close();
            sam_input = nullptr;
        } else if (sin) {
            int ret = hts_close(sin);
            if (ret < 0) fprintf(stderr, "[rank %d] ERROR: closing input failed.\n", rank);
            input_file_mem = nullptr;
        } else if (input_mem_hfile) {
            if (hclose(input_mem_hfile) != 0) {
                fprintf(stderr, "[rank %d] ERROR: closing memory hFILE failed.\n", rank);
            }
            input_file_mem = nullptr;
        } else if (input_file_mem) {
            free(input_file_mem);
            input_file_mem = nullptr;
        }
        double close_cost = GetTime() - close_t0;
        double close_cost_max = MpiReduceMaxCost(close_cost);
        if (rank == 0) printf("666close the files cost %lf-----\n", close_cost_max);
    }
    return exit_code;
}

int ProcessBamToBamMPI(CmdInfo *cmd_info) {
    return ProcessSwBamMPI(cmd_info);
}
