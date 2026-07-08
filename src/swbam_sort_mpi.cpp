#include "swbam_mpi.h"
#include "swbam/mpi_runtime.h"

#include <algorithm>
#include <cctype>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <string>
#include <sys/types.h>
#include <vector>

#include <libdeflate.h>
#include <mpi.h>

namespace {

const int kSortExchangeChunk = 64 * 1024 * 1024;

const unsigned char kSortBgzfEofBlock[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
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

void MpiSortPackLe16(unsigned char *p, uint16_t value) {
    p[0] = (unsigned char)(value & 0xff);
    p[1] = (unsigned char)((value >> 8) & 0xff);
}

void MpiSortPackLe32(unsigned char *p, uint32_t value) {
    p[0] = (unsigned char)(value & 0xff);
    p[1] = (unsigned char)((value >> 8) & 0xff);
    p[2] = (unsigned char)((value >> 16) & 0xff);
    p[3] = (unsigned char)((value >> 24) & 0xff);
}

void MpiSortVectorPutLe32(std::vector<unsigned char> *out, uint32_t value) {
    unsigned char buf[4];
    MpiSortPackLe32(buf, value);
    out->insert(out->end(), buf, buf + 4);
}

bool MpiSortIsBamLikeFormat(int format) {
    return format == bam || format == binary_format;
}

int MpiSortNormalizeFormat(int format) {
    return MpiSortIsBamLikeFormat(format) ? bam : format;
}

int MpiSortAllRanksOk(int local_ok) {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return global_ok;
}

double MpiSortReduceMaxCost(double local_cost) {
    double max_cost = 0.0;
    MPI_Reduce(&local_cost, &max_cost, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return max_cost;
}

int MpiSortLoadFileToMemory(const std::string &path, char **data, size_t *size) {
    *data = nullptr;
    *size = 0;
    FILE *fp = fopen(path.c_str(), "rb");
    if (!fp) return -1;
    if (fseeko(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        return -1;
    }
    off_t end = ftello(fp);
    if (end < 0 || fseeko(fp, 0, SEEK_SET) != 0) {
        fclose(fp);
        return -1;
    }
    if ((unsigned long long)end > (unsigned long long)SIZE_MAX) {
        fclose(fp);
        return -1;
    }
    *size = (size_t)end;
    *data = *size ? (char *)malloc(*size) : nullptr;
    if (*size > 0 && !*data) {
        fclose(fp);
        return -1;
    }
    if (*size > 0 && fread(*data, 1, *size, fp) != *size) {
        fclose(fp);
        free(*data);
        *data = nullptr;
        *size = 0;
        return -1;
    }
    fclose(fp);
    return 0;
}

int MpiSortDumpMemoryToFile(const std::string &path, const char *data, size_t size) {
    FILE *fp = fopen(path.c_str(), "wb");
    if (!fp) return -1;
    size_t written = 0;
    while (written < size) {
        size_t n = fwrite(data + written, 1, size - written, fp);
        if (n == 0) {
            fclose(fp);
            return -1;
        }
        written += n;
    }
    return fclose(fp) == 0 ? 0 : -1;
}

int MpiSortSendBytes(int dst, int tag, const char *data, long long len) {
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

int MpiSortRecvBytes(int src, int tag, char *data, long long len) {
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

int MpiSortInitMemWriter(MemWriter &w, size_t cap) {
    if (cap == 0) cap = 64 * 1024 * 1024;
    w.data = (char *)malloc(cap);
    w.capacity = w.data ? cap : 0;
    w.size = 0;
    return w.data ? 0 : -1;
}

int MpiSortSendRecvBytes(int send_to, const char *send_data, long long send_len,
                         int recv_from, char *recv_data, long long recv_len,
                         int tag) {
    long long sent = 0;
    long long received = 0;
    char dummy = 0;
    while (sent < send_len || received < recv_len) {
        int send_chunk = 0;
        int recv_chunk = 0;
        if (sent < send_len) {
            send_chunk = (int)std::min<long long>(send_len - sent, kSortExchangeChunk);
        }
        if (received < recv_len) {
            recv_chunk = (int)std::min<long long>(recv_len - received, kSortExchangeChunk);
        }
        MPI_Status status;
        const char *send_ptr = send_chunk > 0 ? send_data + sent : &dummy;
        char *recv_ptr = recv_chunk > 0 ? recv_data + received : &dummy;
        if (MPI_Sendrecv((void *)send_ptr, send_chunk, MPI_BYTE, send_to, tag,
                         recv_ptr, recv_chunk, MPI_BYTE, recv_from, tag,
                         MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
            return -1;
        }
        sent += send_chunk;
        received += recv_chunk;
    }
    return 0;
}

int MpiSortScanBgzfBlocksInMemory(const char *base, size_t size, long long body_start,
                                  std::vector<long long> *offsets,
                                  std::vector<long long> *lengths) {
    if (!base || body_start < 0 || (unsigned long long)body_start > (unsigned long long)size) return -1;
    long long pos = body_start;
    while ((unsigned long long)pos < (unsigned long long)size) {
        if ((unsigned long long)pos + BLOCK_HEADER_LENGTH > (unsigned long long)size) return -1;
        const unsigned char *header = (const unsigned char *)(base + pos);
        int block_len = (int)MpiSortReadLe16(header + 16) + 1;
        if (block_len <= 0 || (unsigned long long)pos + (unsigned long long)block_len > (unsigned long long)size) {
            return -1;
        }
        bool is_eof = block_len == (int)sizeof(kSortBgzfEofBlock) &&
                      memcmp(base + pos, kSortBgzfEofBlock, sizeof(kSortBgzfEofBlock)) == 0;
        if (is_eof) break;
        offsets->push_back(pos);
        lengths->push_back(block_len);
        pos += block_len;
    }
    return 0;
}

int MpiSortSelectBlockRangeFromMemory(char *base, size_t input_size,
                                      const std::vector<long long> &offsets,
                                      const std::vector<long long> &lengths,
                                      long long begin,
                                      long long end,
                                      char **data,
                                      size_t *size) {
    *data = nullptr;
    *size = 0;
    if (begin >= end) return 0;
    long long start = offsets[(size_t)begin];
    long long stop = offsets[(size_t)(end - 1)] + lengths[(size_t)(end - 1)];
    if (start < 0 || stop < start || (unsigned long long)stop > (unsigned long long)input_size) return -1;
    *data = base + start;
    *size = (size_t)(stop - start);
    return 0;
}

int MpiSortAppendBgzfPayloadToMem(MemWriter &w,
                                  const unsigned char *src,
                                  size_t src_len,
                                  int compress_level,
                                  struct libdeflate_compressor *compressor) {
    if (src_len == 0) return 0;
    if (src_len > BGZF_BLOCK_SIZE) return -1;
    if (compress_level != 0 && compress_level != 1 && compress_level != 6) return -1;
    if (compress_level != 0 && !compressor) return -1;

    unsigned char out[BGZF_MAX_BLOCK_SIZE];
    size_t block_len = 0;
    if (compress_level == 0) {
        if (src_len + 5 + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH > BGZF_MAX_BLOCK_SIZE) return -1;
        out[BLOCK_HEADER_LENGTH] = 1;
        MpiSortPackLe16(out + BLOCK_HEADER_LENGTH + 1, (uint16_t)src_len);
        MpiSortPackLe16(out + BLOCK_HEADER_LENGTH + 3, (uint16_t)~(uint16_t)src_len);
        memcpy(out + BLOCK_HEADER_LENGTH + 5, src, src_len);
        block_len = src_len + 5 + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    } else {
        size_t clen = libdeflate_deflate_compress(
            compressor, src, src_len,
            out + BLOCK_HEADER_LENGTH,
            BGZF_MAX_BLOCK_SIZE - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH);
        if (clen == 0) return -1;
        block_len = clen + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    }
    if (block_len > BGZF_MAX_BLOCK_SIZE) return -1;

    memcpy(out, g_magic, BLOCK_HEADER_LENGTH);
    MpiSortPackLe16(out + 16, (uint16_t)(block_len - 1));
    uint32_t crc = libdeflate_crc32(0, src, src_len);
    MpiSortPackLe32(out + block_len - 8, crc);
    MpiSortPackLe32(out + block_len - 4, (uint32_t)src_len);
    return MpiWriteBytesToMem(w, (const char *)out, block_len);
}

int MpiSortBuildBamHeaderMemory(sam_hdr_t *hdr, int compress_level, char **data, size_t *size) {
    *data = nullptr;
    *size = 0;
    if (!hdr) return -1;

    const char *text = sam_hdr_str(hdr);
    size_t text_len = sam_hdr_length(hdr);
    if (text_len > 0 && !text) return -1;
    if ((unsigned long long)text_len > (unsigned long long)UINT32_MAX) return -1;
    if (hdr->n_targets < 0) return -1;

    std::vector<unsigned char> raw;
    raw.reserve(12 + text_len + (size_t)hdr->n_targets * 32);
    const unsigned char magic[4] = {'B', 'A', 'M', 1};
    raw.insert(raw.end(), magic, magic + 4);
    MpiSortVectorPutLe32(&raw, (uint32_t)text_len);
    if (text_len > 0) raw.insert(raw.end(), text, text + text_len);
    MpiSortVectorPutLe32(&raw, (uint32_t)hdr->n_targets);
    for (int32_t i = 0; i < hdr->n_targets; ++i) {
        if (!hdr->target_name || !hdr->target_name[i] || !hdr->target_len) return -1;
        size_t name_len = strlen(hdr->target_name[i]) + 1;
        if ((unsigned long long)name_len > (unsigned long long)UINT32_MAX) return -1;
        MpiSortVectorPutLe32(&raw, (uint32_t)name_len);
        raw.insert(raw.end(), hdr->target_name[i], hdr->target_name[i] + name_len);
        MpiSortVectorPutLe32(&raw, hdr->target_len[i]);
    }

    MemWriter header_writer = {};
    if (MpiSortInitMemWriter(header_writer, raw.size() + 64 * 1024) != 0) return -1;
    struct libdeflate_compressor *compressor =
        compress_level == 0 ? nullptr : libdeflate_alloc_compressor(compress_level);
    if (compress_level != 0 && !compressor) {
        free(header_writer.data);
        return -1;
    }
    size_t pos = 0;
    while (pos < raw.size()) {
        size_t n = std::min((size_t)BGZF_BLOCK_SIZE, raw.size() - pos);
        if (MpiSortAppendBgzfPayloadToMem(header_writer, raw.data() + pos, n,
                                          compress_level, compressor) != 0) {
            if (compressor) libdeflate_free_compressor(compressor);
            free(header_writer.data);
            return -1;
        }
        pos += n;
    }
    if (compressor) libdeflate_free_compressor(compressor);
    *data = header_writer.data;
    *size = header_writer.size;
    return 0;
}

int MpiSortParseMemoryLimit(const std::string &text, size_t *bytes) {
    *bytes = 0;
    if (text.empty()) return 0;
    const char *s = text.c_str();
    char *end = nullptr;
    unsigned long long value = strtoull(s, &end, 10);
    if (end == s) return -1;
    while (*end && isspace((unsigned char)*end)) end++;
    unsigned long long mult = 1;
    if (*end) {
        char c = (char)tolower((unsigned char)*end);
        if (c == 'k') mult = 1024ull;
        else if (c == 'm') mult = 1024ull * 1024ull;
        else if (c == 'g') mult = 1024ull * 1024ull * 1024ull;
        else if (c == 't') mult = 1024ull * 1024ull * 1024ull * 1024ull;
        else return -1;
        end++;
        if (*end == 'b' || *end == 'B') end++;
        while (*end && isspace((unsigned char)*end)) end++;
        if (*end) return -1;
    }
    if (value != 0 && mult > ULLONG_MAX / value) return -1;
    unsigned long long result = value * mult;
    if (result > (unsigned long long)SIZE_MAX) return -1;
    *bytes = (size_t)result;
    return 0;
}

int MpiSortCheckMemoryLimit(size_t limit, size_t estimate, int rank, const char *stage) {
    if (limit == 0 || estimate <= limit) return 0;
    fprintf(stderr,
            "[rank %d] ERROR: MPI sort memory check failed. "
            "stage=%s estimate=%zu limit=%zu\n",
            rank, stage, estimate, limit);
    return -1;
}

size_t MpiSortEstimateInMemoryNeed(size_t rank_input_size) {
    const size_t margin = 256ull * 1024ull * 1024ull;
    if (rank_input_size > (SIZE_MAX - margin) / 3u) return SIZE_MAX;
    return rank_input_size * 3u + margin;
}

int MpiSortUpdateHeaderCoordinate(sam_hdr_t *hdr) {
    if (!hdr) return -1;
    if (sam_hdr_update_hd(hdr, "SO", "coordinate") < 0) {
        if (sam_hdr_add_line(hdr, "HD", "VN", "1.6", "SO", "coordinate", NULL) < 0) {
            return -1;
        }
    }
    (void)sam_hdr_remove_tag_hd(hdr, "GO");
    (void)sam_hdr_remove_tag_hd(hdr, "SS");
    return 0;
}

void MpiSortPrintRankStats(int rank, int comm_size,
                           const MpiSortStats &stats,
                           size_t body_size) {
    const double core_stage = stats.sort_mode == 1
        ? stats.t_fused_total
        : stats.t_setup + stats.t_extract_read_unhidden +
          stats.t_extract + stats.t_extract_prepare +
          stats.t_extract_merge + stats.t_local_sort + stats.t_sample +
          stats.t_partition + stats.t_exchange + stats.t_offset_fix +
          stats.t_final_sort + stats.t_compress +
          stats.t_status_check + stats.t_cleanup;
    const double unaccounted = stats.t_fused_total - core_stage;
    const long long bucket_total_raw =
        stats.bucket_self_raw_bytes + stats.bucket_remote_raw_bytes;
    const double bucket_self_raw_ratio =
        bucket_total_raw > 0 ? 100.0 * (double)stats.bucket_self_raw_bytes / (double)bucket_total_raw : 0.0;
    const char *mode = stats.sort_mode == 1 ? "external" : "memory";
    char local_lines[12288] = {};
    snprintf(local_lines, sizeof(local_lines),
             "[rank %d] mode=%s blocks=%lld local_records=%lld received_records=%lld samples=%lld bgzf=%lld body=%zu\n"
             "[rank %d] extract=%.3f local_sort=%.3f sample=%.3f partition=%.3f exchange=%.3f final_sort=%.3f compress=%.3f write=%.3f\n"
             "[rank %d] pipeline_detail setup=%.3f read=%.3f read_unhidden=%.3f extract_prepare=%.3f extract_merge=%.3f bucket_count=%.3f bucket_pack=%.3f mpi_exchange=%.3f offset_fix=%.3f status=%.3f cleanup=%.3f unaccounted=%.3f\n"
             "[rank %d] bucket_dist self_records=%lld remote_records=%lld self_raw=%lld remote_raw=%lld self_raw_ratio=%.2f%%\n"
             "[rank %d] external_sim mpi=%.3f merge=%.3f compress=%.3f consolidation=%.3f temp_write=%.3f temp_read=%.3f\n"
             "[rank %d] external_actual runs=%lld segments=%lld temp_open=%.3f temp_write=%.3f temp_read=%.3f exchange_wall=%.3f merge_wall=%.3f compress_pipeline_wall=%.3f\n"
             "[rank %d] external_memory temp_write_bytes=%lld temp_read_bytes=%lld resident_records=%lld resident_raw=%lld tracked_peak=%lld limit_arena=%lld merge_fan_in=%lld consolidation_passes=%lld\n"
             "[rank %d] extract_detail alloc=%.3f inflate=%.3f crc=%.3f parse=%.3f other=%.3f\n"
             "[rank %d] compress_detail pack=%.3f alloc=%.3f deflate=%.3f footer=%.3f other=%.3f\n"
             "[rank %d] fused_total=%.6f actual_wall_including_probe=%.6f core_stage=%.6f\n",
             rank, mode, stats.input_blocks, stats.local_records, stats.received_records,
             stats.sample_records, stats.bgzf_blocks, body_size,
             rank, stats.t_extract, stats.t_local_sort, stats.t_sample,
             stats.t_partition, stats.t_exchange, stats.t_final_sort,
             stats.t_compress, stats.t_write,
             rank, stats.t_setup, stats.t_extract_read, stats.t_extract_read_unhidden,
             stats.t_extract_prepare, stats.t_extract_merge,
             stats.t_bucket_count, stats.t_bucket_pack, stats.t_mpi_exchange,
             stats.t_offset_fix, stats.t_status_check, stats.t_cleanup, unaccounted,
             rank, stats.bucket_self_records, stats.bucket_remote_records,
             stats.bucket_self_raw_bytes, stats.bucket_remote_raw_bytes,
             bucket_self_raw_ratio,
             rank, stats.t_mpi_simulated, stats.t_merge_simulated,
             stats.t_compress_simulated,
             stats.t_consolidation_simulated,
             stats.t_temp_write_sim, stats.t_temp_read_sim,
             rank, stats.external_runs, stats.external_segments,
             stats.t_temp_open_actual,
             stats.t_temp_write_actual, stats.t_temp_read_actual,
             stats.t_exchange, stats.t_final_sort, stats.t_compress,
             rank, stats.temp_write_bytes, stats.temp_read_bytes,
             stats.resident_run_records, stats.resident_run_raw_bytes,
             stats.tracked_peak_bytes, stats.run_arena_bytes,
             stats.merge_fan_in, stats.consolidation_passes,
             rank, stats.t_extract_alloc, stats.t_extract_inflate, stats.t_extract_crc,
             stats.t_extract_parse, stats.t_extract_other,
             rank, stats.t_compress_pack, stats.t_compress_alloc, stats.t_compress_deflate,
             stats.t_compress_footer, stats.t_compress_other,
             rank, stats.t_fused_total, stats.t_fused_actual, core_stage);

    const int kLineBytes = 12288;
    std::vector<char> gathered;
    if (rank == 0) gathered.resize((size_t)comm_size * kLineBytes);
    MPI_Gather(local_lines, kLineBytes, MPI_CHAR,
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

} // namespace

int MpiSortMemoryToMemory(CmdInfo *cmd_info,
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
    char *input_file_mem = const_cast<char *>(input_memory_const);
    sam_hdr_t *hdr = nullptr;
    char *rank_input_mem = nullptr;
    char *bam_header_mem = nullptr;
    char *simulated_write_mem = nullptr;
    size_t bam_header_size = 0;
    size_t output_file_size = 0;
    size_t output_body_start = 0;
    MemReader reader = {};
    MemWriter mem_writer = {};
    MpiSortStats stats = {};
    std::vector<long long> block_offsets;
    std::vector<long long> block_lengths;
    std::vector<long long> body_sizes;
    std::vector<long long> body_prefixes;
    long long body_start = 0;
    long long n_blocks = 0;
    long long local_block_begin = 0;
    long long local_block_end = 0;
    long long local_body_size = 0;
    long long total_body_size = 0;
    size_t memory_limit = 0;
    volatile unsigned long long simulated_write_guard = 0;
    double body_total_t0 = 0.0;
    double stage41_cost_max = 0.0;
    double stage42_cost_max = 0.0;
    double fused_cost_max = 0.0;
    double stage44_cost_max = 0.0;
    double stage45_header_cost_max = 0.0;
    double stage46_cost_max = 0.0;

    if (output_bam) {
        output_bam->data = nullptr;
        output_bam->size = 0;
    }
    if (core_cost) *core_cost = 0.0;

    if (!input_file_mem || input_size == 0 || !output_bam ||
        MpiSortParseMemoryLimit(cmd_info->sort_memory_,
                                &memory_limit) != 0) {
        if (rank == 0) {
            fprintf(stderr,
                    "ERROR: invalid dedup-pipeline sort input or memory limit.\n");
        }
        local_ok = 0;
    }
    if (!MpiSortAllRanksOk(local_ok)) goto cleanup;
    {
        double init_cost_max =
            MpiSortReduceMaxCost(GetTime() - total_t0);
        if (rank == 0) {
            printf("111Complete the initialization cost %lf-----\n",
                   init_cost_max);
        }
    }

    {
        double header_t0 = GetTime();
        if (rank == 0 &&
            MpiCommonReadBamHeaderFromMemory(
                input_file_mem, input_size,
                &hdr, &body_start) != 0) {
            fprintf(stderr,
                    "ERROR: failed to parse sort pipeline BAM header.\n");
            local_ok = 0;
        }
        MPI_Bcast(&local_ok, 1, MPI_INT, 0,
                  MPI_COMM_WORLD);
        MPI_Bcast(&body_start, 1, MPI_LONG_LONG, 0,
                  MPI_COMM_WORLD);
        double header_cost_max =
            MpiSortReduceMaxCost(GetTime() - header_t0);
        if (rank == 0 && local_ok) {
            printf("333Complete the head cost %lf---\n",
                   header_cost_max);
        }
    }
    if (!MpiSortAllRanksOk(local_ok)) goto cleanup;

    body_total_t0 = GetTime();
    {
        double stage41_t0 = GetTime();
        if (rank == 0) {
            printf("Enable MPI BAM SORT mode (%d MPE + %d CPEs)!!!\n",
                   comm_size, comm_size * 64);
            printf("MPI BAM output compression level=%d\n",
                   cmd_info->compress_level_);
            if (MpiCommonScanBgzfBlocksInMemory(
                    input_file_mem, input_size, body_start,
                    &block_offsets, &block_lengths) != 0) {
                local_ok = 0;
            }
            n_blocks = (long long)block_offsets.size();
            if (n_blocks > (long long)INT_MAX) local_ok = 0;
            if (local_ok) {
                printf("MPI BAM scan complete. data_blocks=%lld "
                       "body_start=%lld header_end=%lld\n",
                       n_blocks, body_start,
                       body_start);
            }
        }
        MPI_Bcast(&local_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (!local_ok) goto cleanup;
        MPI_Bcast(&body_start, 1, MPI_LONG_LONG, 0,
                  MPI_COMM_WORLD);
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
        local_block_end = n_blocks * (rank + 1) / comm_size;
        size_t rank_input_size = 0;
        if (MpiSortSelectBlockRangeFromMemory(
                input_file_mem, input_size,
                block_offsets, block_lengths,
                local_block_begin, local_block_end,
                &rank_input_mem, &rank_input_size) != 0) {
            local_ok = 0;
        }
        stage41_cost_max =
            MpiSortReduceMaxCost(GetTime() - stage41_t0);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.1 scan/split/broadcast/select cost %lf\n",
                   stage41_cost_max);
        }
        if (!MpiSortAllRanksOk(local_ok)) goto cleanup;

        double stage42_t0 = GetTime();
        reader.base = rank_input_mem;
        reader.size = rank_input_size;
        reader.pos = 0;
        if (MpiSortInitMemWriter(
                mem_writer,
                rank_input_size ? rank_input_size
                                : 64 * 1024 * 1024) != 0) {
            local_ok = 0;
        }
        int stage42_ok = MpiSortAllRanksOk(local_ok);
        stage42_cost_max =
            MpiSortReduceMaxCost(GetTime() - stage42_t0);
        if (rank == 0 && stage42_ok) {
            printf("Complete the 4.2 init reader/writer cost %lf\n",
                   stage42_cost_max);
        }
        if (!stage42_ok) goto cleanup;

        size_t estimated =
            MpiSortEstimateInMemoryNeed(rank_input_size);
        unsigned long long local_estimate_u64 =
            estimated == SIZE_MAX ? ULLONG_MAX
                                  : (unsigned long long)estimated;
        unsigned long long max_estimate_u64 = 0;
        MPI_Allreduce(&local_estimate_u64, &max_estimate_u64,
                      1, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                      MPI_COMM_WORLD);
        int local_external =
            (memory_limit != 0 && estimated > memory_limit)
                ? 1 : 0;
        int use_external = 0;
        MPI_Allreduce(&local_external, &use_external, 1,
                      MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        stats.sort_mode = 0;
        if (rank == 0) {
            printf("MPI BAM sort selected mode=memory "
                   "memory_estimate_max=%llu memory_limit=%zu\n",
                   max_estimate_u64, memory_limit);
            printf("MPI BAM sort pipeline policy: memory only; "
                   "external mode is disabled.\n");
        }
        if (use_external) {
            if (rank == 0) {
                fprintf(stderr,
                        "ERROR: dedup-pipeline sort memory estimate exceeds -m; external intermediate storage is disabled in v1.\n");
            }
            local_ok = 0;
        }
        if (!MpiSortAllRanksOk(local_ok)) goto cleanup;

        double fused_t0 = GetTime();
        int fused_ret =
            FusedBamSortMPI(reader, mem_writer,
                            local_block_begin, rank,
                            comm_size, cmd_info->compress_level_,
                            memory_limit, &stats);
        if (fused_ret != 0) local_ok = 0;
        int global_ok = MpiSortAllRanksOk(local_ok);
        stats.t_fused_actual = GetTime() - fused_t0;
        fused_cost_max =
            MpiSortReduceMaxCost(stats.t_fused_actual);
        if (rank == 0 && global_ok) {
            printf("Complete the 4.3 FusedBamSortMPI cost %lf\n",
                   fused_cost_max);
            printf("FusedBamSortMPI reported cost model=measured-memory-mode\n");
            printf("FusedBamSortMPI actual wall %lf\n",
                   fused_cost_max);
        }
        if (!global_ok) goto cleanup;

        double stage44_t0 = GetTime();
        local_body_size = (long long)mem_writer.size;
        body_sizes.assign((size_t)comm_size, 0);
        body_prefixes.assign((size_t)comm_size, 0);
        MPI_Allgather(&local_body_size, 1, MPI_LONG_LONG,
                      body_sizes.data(), 1, MPI_LONG_LONG,
                      MPI_COMM_WORLD);
        total_body_size = 0;
        for (int i = 0; i < comm_size; ++i) {
            body_prefixes[(size_t)i] = total_body_size;
            if (body_sizes[(size_t)i] < 0 ||
                total_body_size >
                    LLONG_MAX - body_sizes[(size_t)i]) {
                local_ok = 0;
            } else {
                total_body_size += body_sizes[(size_t)i];
            }
        }
        stage44_cost_max =
            MpiSortReduceMaxCost(GetTime() - stage44_t0);
        if (rank == 0) {
            printf("Complete the 4.4 gather body sizes cost %lf\n",
                   stage44_cost_max);
        }
        if (!MpiSortAllRanksOk(local_ok)) goto cleanup;

        double stage45_header_t0 = GetTime();
        long long output_body_start_ll = 0;
        long long output_file_size_ll = 0;
        if (rank == 0) {
            if (MpiSortUpdateHeaderCoordinate(hdr) != 0 ||
                MpiSortBuildBamHeaderMemory(
                    hdr, cmd_info->compress_level_,
                    &bam_header_mem, &bam_header_size) != 0) {
                local_ok = 0;
            }
            if (local_ok) {
                unsigned long long final_size =
                    (unsigned long long)bam_header_size +
                    (unsigned long long)total_body_size +
                    (unsigned long long)sizeof(kSortBgzfEofBlock);
                if (final_size > (unsigned long long)SIZE_MAX ||
                    final_size > (unsigned long long)LLONG_MAX ||
                    (unsigned long long)bam_header_size >
                        (unsigned long long)LLONG_MAX) {
                    local_ok = 0;
                } else {
                    output_body_start = bam_header_size;
                    output_file_size = (size_t)final_size;
                    output_body_start_ll =
                        (long long)output_body_start;
                    output_file_size_ll =
                        (long long)output_file_size;
                }
            }
        }
        int stage45_header_ok = MpiSortAllRanksOk(local_ok);
        MPI_Bcast(&output_body_start_ll, 1, MPI_LONG_LONG,
                  0, MPI_COMM_WORLD);
        MPI_Bcast(&output_file_size_ll, 1, MPI_LONG_LONG,
                  0, MPI_COMM_WORLD);
        if (stage45_header_ok) {
            output_body_start = (size_t)output_body_start_ll;
            output_file_size = (size_t)output_file_size_ll;
        }
        stage45_header_cost_max =
            MpiSortReduceMaxCost(GetTime() - stage45_header_t0);
        if (rank == 0 && stage45_header_ok) {
            printf("Complete the 4.5a header/layout cost %lf\n",
                   stage45_header_cost_max);
        }
        if (!stage45_header_ok) goto cleanup;

        double stage45_malloc_t0 = GetTime();
        size_t simulated_write_size =
            local_body_size > 0 ? (size_t)local_body_size : 0;
        if (rank == 0) {
            if (output_body_start >
                SIZE_MAX - simulated_write_size) {
                local_ok = 0;
            } else {
                simulated_write_size += output_body_start;
            }
        }
        if (local_ok && simulated_write_size > 0) {
            simulated_write_mem =
                (char *)malloc(simulated_write_size);
            if (!simulated_write_mem) local_ok = 0;
        }
        int stage45_malloc_ok = MpiSortAllRanksOk(local_ok);
        double stage45_malloc_cost_max =
            MpiSortReduceMaxCost(GetTime() - stage45_malloc_t0);
        if (rank == 0 && stage45_malloc_ok) {
            printf("Complete the 4.5b malloc simulated write memory cost %lf\n",
                   stage45_malloc_cost_max);
        }
        if (!stage45_malloc_ok) goto cleanup;

        double sim_write_t0 = GetTime();
        size_t simulated_pos = 0;
        if (rank == 0 && output_body_start > 0) {
            memcpy(simulated_write_mem + simulated_pos,
                   bam_header_mem, bam_header_size);
            simulated_pos += bam_header_size;
        }
        if (local_body_size > 0) {
            memcpy(simulated_write_mem + simulated_pos,
                   mem_writer.data, (size_t)local_body_size);
            simulated_pos += (size_t)local_body_size;
        }
        if (simulated_pos > 0) {
            unsigned char *guard_ptr =
                (unsigned char *)simulated_write_mem;
            simulated_write_guard += guard_ptr[0];
            simulated_write_guard += guard_ptr[simulated_pos - 1];
        }
        stats.t_write += GetTime() - sim_write_t0;
        stage46_cost_max =
            MpiSortReduceMaxCost(stats.t_write);
        int stage46_ok = MpiSortAllRanksOk(local_ok);
        if (rank == 0 && stage46_ok) {
            printf("Complete the 4.6 simulated header/body distributed write cost %lf\n",
                   stage46_cost_max);
            printf("Complete the total (4.1~4.6) cost %lf-----\n",
                   stage41_cost_max + stage42_cost_max +
                   fused_cost_max + stage44_cost_max +
                   stage45_header_cost_max + stage46_cost_max);
        }
        if (simulated_write_guard == (unsigned long long)-1 &&
            rank < 0) {
            fprintf(stderr, "unused simulated write guard %llu\n",
                    simulated_write_guard);
        }
        if (!stage46_ok) goto cleanup;

        double stage47_t0 = GetTime();
        MpiSortPrintRankStats(rank, comm_size,
                              stats, mem_writer.size);
        double stage47_cost_max =
            MpiSortReduceMaxCost(GetTime() - stage47_t0);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.7 rank stats reduce/print cost %lf\n",
                   stage47_cost_max);
        }

        double body_total_cost_max =
            MpiSortReduceMaxCost(GetTime() - body_total_t0);
        if (rank == 0 && local_ok) {
            printf("444Complete the total body cost %lf\n",
                   body_total_cost_max);
        }
    }

    {
        double verify_alloc_t0 = GetTime();
        if (rank == 0) {
            output_bam->data =
                output_file_size ? (char *)malloc(output_file_size)
                                 : nullptr;
            output_bam->size = output_file_size;
            if (output_file_size > 0 && !output_bam->data) {
                local_ok = 0;
            }
            if (local_ok && bam_header_size > 0) {
                memcpy(output_bam->data, bam_header_mem,
                       bam_header_size);
            }
        }
        double verify_alloc_cost_max =
            MpiSortReduceMaxCost(GetTime() - verify_alloc_t0);
        if (rank == 0 && local_ok) {
            printf("555Prepare verification output memory cost %lf--\n",
                   verify_alloc_cost_max);
        }
        if (!MpiSortAllRanksOk(local_ok)) goto cleanup;

        double verify_gather_t0 = GetTime();
        if (rank == 0) {
            if (local_body_size > 0) {
                memcpy(output_bam->data + output_body_start +
                           body_prefixes[(size_t)rank],
                       mem_writer.data, (size_t)local_body_size);
            }
            for (int src = 1; src < comm_size; ++src) {
                long long recv_size = body_sizes[(size_t)src];
                if (recv_size <= 0) continue;
                if (MpiSortRecvBytes(
                        src, 0,
                        output_bam->data + output_body_start +
                            body_prefixes[(size_t)src],
                        recv_size) != 0) {
                    local_ok = 0;
                    break;
                }
            }
            if (local_ok) {
                memcpy(output_bam->data + output_body_start +
                           total_body_size,
                       kSortBgzfEofBlock,
                       sizeof(kSortBgzfEofBlock));
            }
        } else if (local_body_size > 0) {
            if (MpiSortSendBytes(0, 0, mem_writer.data,
                                 local_body_size) != 0) {
                local_ok = 0;
            }
        }
        double verify_gather_cost_max =
            MpiSortReduceMaxCost(GetTime() - verify_gather_t0);
        if (rank == 0 && local_ok) {
            printf("555Gather verification output memory cost %lf--\n",
                   verify_gather_cost_max);
            printf("555Keep output memory cost 0.000000--\n");
        }
        if (!MpiSortAllRanksOk(local_ok)) goto cleanup;
    }

    if (core_cost) {
        *core_cost = stage41_cost_max + stage42_cost_max +
                     fused_cost_max + stage44_cost_max +
                     stage45_header_cost_max + stage46_cost_max;
    }
    exit_code = 0;

cleanup:
    if (mem_writer.data) free(mem_writer.data);
    if (bam_header_mem) free(bam_header_mem);
    if (simulated_write_mem) free(simulated_write_mem);
    if (hdr) sam_hdr_destroy(hdr);
    if (rank == 0) {
        printf("666sort total process cost %lf-----\n",
               GetTime() - total_t0);
    }
    return exit_code;
}

int ProcessSortMPI(CmdInfo *cmd_info) {

    //1.相关数据的初始化
    double t_init = GetTime();

    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int exit_code = 1;
    int local_ok = 1;
    swbam::mpi::MpiBamInput input_handle;
    swbam::BamInputBackend *input_backend = nullptr;
    swbam::MemoryBamInput *memory_input = nullptr;
    swbam::mpi::MpiBamInputPlan input_plan;
    sam_hdr_t *hdr = nullptr;
    char *rank_input_mem = nullptr;
    char *bam_header_mem = nullptr;
    size_t bam_header_size = 0;
    char *simulated_write_mem = nullptr;
    char *output_file_mem = nullptr;
    size_t output_file_size = 0;
    size_t output_body_start = 0;
    size_t simulated_write_size = 0;
    MemReader reader = {};
    swbam::AdaptiveRankBodySink rank_body_sink;
    MpiSortStats stats = {};
    std::vector<long long> body_sizes;
    std::vector<long long> body_prefixes;
    std::vector<unsigned char> body_chunk;
    long long body_start = 0;
    long long n_blocks = 0;
    long long local_block_begin = 0;
    long long local_block_end = 0;
    long long local_body_size = 0;
    long long local_prefix = 0;
    long long total_body_size = 0;
    size_t memory_limit = 0;
    uint64_t rank_body_memory_limit = 0;
    std::string sort_temp_prefix;
    volatile unsigned long long simulated_write_guard = 0;

    if (MpiSortParseMemoryLimit(cmd_info->sort_memory_, &memory_limit) != 0) {
        if (rank == 0) fprintf(stderr, "ERROR: invalid sort memory limit '%s'.\n", cmd_info->sort_memory_.c_str());
        local_ok = 0;
    }
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
    if (!cmd_info->sort_temp_prefix_.empty()) {
        sort_temp_prefix = cmd_info->sort_temp_prefix_;
    } else {
        size_t slash = cmd_info->out_file_name_.find_last_of('/');
        sort_temp_prefix = slash == std::string::npos
            ? "./rabbitbam-sort"
            : cmd_info->out_file_name_.substr(0, slash + 1) + "rabbitbam-sort";
    }
    if (!MpiSortAllRanksOk(local_ok)) goto cleanup;

    {
        double init_cost = GetTime() - t_init;
        double init_cost_max = MpiSortReduceMaxCost(init_cost);
        if (rank == 0) printf("111Complete the initialization cost %lf-----\n", init_cost_max);
    }


    //2. 通过公共 backend 打开 BAM；memory 仍是默认快速路径。
    {
        if (input_handle.Open(
                cmd_info->in_file_name_, cmd_info->io_backend_,
                cmd_info->io_memory_limit_) != 0) {
            fprintf(stderr, "[rank %d] ERROR: cannot open sort input %s\n",
                    rank, cmd_info->in_file_name_.c_str());
            local_ok = 0;
        }
        input_backend = input_handle.backend();
        if (local_ok && (!input_backend ||
                         MpiSortNormalizeFormat(input_backend->format()) != bam)) {
            local_ok = 0;
        }
        memory_input = dynamic_cast<swbam::MemoryBamInput *>(input_backend);
        double open_cost_max = MpiSortReduceMaxCost(
            input_handle.data_open_cost());
        if (rank == 0 && local_ok) {
            printf("MPI BAM input backend=%s auto_memory_budget=%llu\n",
                   input_handle.selected_backend().c_str(),
                   (unsigned long long)input_handle.auto_memory_budget());
            printf("222Complete the input open cost %lf--\n", open_cost_max);
        }
    }
    if (!MpiSortAllRanksOk(local_ok)) goto cleanup;


    //3. 从 backend 复制 header，并记录 BAM body 起始位置。
    {
        double header_t0 = GetTime();
        if (local_ok && input_backend->header()) {
            hdr = sam_hdr_dup(input_backend->header());
        }
        if (!hdr) local_ok = 0;
        if (local_ok) {
            body_start = (long long)input_backend->body_offset();
            if (body_start < 0 || (unsigned long long)body_start >
                    (unsigned long long)input_backend->size()) {
                fprintf(stderr, "[rank %d] ERROR: invalid BAM body start offset %lld.\n", rank, body_start);
                local_ok = 0;
            }
        }
        double header_cost = input_handle.header_open_cost() +
                             GetTime() - header_t0;
        double header_cost_max = MpiSortReduceMaxCost(header_cost);
        if (rank == 0 && local_ok) printf("333Complete the head cost %lf---\n", header_cost_max);
    }
    if (!MpiSortAllRanksOk(local_ok)) goto cleanup;


    //4.核心处理阶段
    {
        double body_total_t0 = GetTime();
        double body_t0 = GetTime();

        //4.1 公共 runtime 扫描 BGZF blocks 并按 rank 分区。
        double stage41_t0 = GetTime();
        if (rank == 0) {
            printf("Enable MPI BAM SORT mode (%d MPE + %d CPEs)!!!\n",
                   comm_size, comm_size * 64);
            printf("MPI BAM output compression level=%d\n", cmd_info->compress_level_);
        }
        if (swbam::mpi::PrepareMpiBamInputPlan(
                *input_backend, &input_plan) != 0) {
            local_ok = 0;
        }
        n_blocks = (long long)input_plan.blocks.size();
        local_block_begin = (long long)input_plan.rank_begin;
        local_block_end = (long long)input_plan.rank_end;
        size_t rank_input_size = 0;
        for (size_t i = input_plan.rank_begin;
             local_ok && i < input_plan.rank_end; ++i) {
            const size_t block_size = input_plan.blocks[i].compressed_size;
            if (block_size > SIZE_MAX - rank_input_size) {
                local_ok = 0;
            } else {
                rank_input_size += block_size;
            }
        }
        if (local_ok && memory_input && input_plan.rank_block_count() > 0) {
            const swbam::BgzfBlockSpan *spans = input_plan.rank_spans();
            const uint64_t start = spans[0].offset;
            uint64_t expected = start;
            for (size_t i = 0;
                 local_ok && i < input_plan.rank_block_count(); ++i) {
                if (spans[i].offset != expected) {
                    local_ok = 0;
                } else {
                    expected += spans[i].compressed_size;
                }
            }
            if (local_ok &&
                (start > memory_input->size() ||
                 rank_input_size > memory_input->size() - (size_t)start)) {
                local_ok = 0;
            }
            if (local_ok) {
                rank_input_mem = memory_input->data() + (size_t)start;
                reader.base = rank_input_mem;
                reader.size = rank_input_size;
                reader.pos = 0;
            }
        }
        double stage41_cost = GetTime() - stage41_t0;
        double stage41_cost_max = MpiSortReduceMaxCost(stage41_cost);
        if (rank == 0 && local_ok) {
            printf("MPI BAM scan complete. data_blocks=%lld body_start=%lld\n",
                   n_blocks, body_start);
            printf("Complete the 4.1 scan/split/broadcast/select cost %lf\n", stage41_cost_max);
        }
        if (!MpiSortAllRanksOk(local_ok)) goto cleanup;

        //4.2 初始化每个 rank 的自适应 body sink。
        double stage42_t0 = GetTime();
        char spool_prefix[96];
        snprintf(spool_prefix, sizeof(spool_prefix),
                 "rabbitbam-sort-rank-%d.body", rank);
        if (rank_body_sink.Open(
                cmd_info->rank_body_backend_, rank_body_memory_limit,
                rank_input_size ? rank_input_size : 64 * 1024 * 1024,
                spool_prefix, cmd_info->rank_body_temp_dir_) != 0) {
            fprintf(stderr,
                    "[rank %d] ERROR: failed to open sort rank body sink "
                    "backend=%s temp_dir=%s.\n",
                    rank, cmd_info->rank_body_backend_.c_str(),
                    cmd_info->rank_body_temp_dir_.c_str());
            local_ok = 0;
        }
        int stage42_ok = MpiSortAllRanksOk(local_ok);
        double stage42_cost = GetTime() - stage42_t0;
        double stage42_cost_max = MpiSortReduceMaxCost(stage42_cost);
        if (rank == 0 && stage42_ok) {
            printf("Complete the 4.2 init reader/writer cost %lf\n", stage42_cost_max);
        }
        if (!stage42_ok) goto cleanup;

        //4.3 核心处理：内存足够时走内排序；-m 不足时自动切到外排序。
        size_t estimated_in_memory_need = MpiSortEstimateInMemoryNeed(rank_input_size);
        unsigned long long local_estimate_u64 =
            estimated_in_memory_need == SIZE_MAX ? ULLONG_MAX : (unsigned long long)estimated_in_memory_need;
        unsigned long long max_estimate_u64 = 0;
        MPI_Allreduce(&local_estimate_u64, &max_estimate_u64, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        int local_external = (memory_limit != 0 && estimated_in_memory_need > memory_limit) ? 1 : 0;
        int use_external = 0;
        MPI_Allreduce(&local_external, &use_external, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        stats.sort_mode = use_external ? 1 : 0;
        if (rank == 0) {
            printf("MPI BAM sort selected mode=%s memory_estimate_max=%llu memory_limit=%zu\n",
                   use_external ? "external" : "memory",
                   max_estimate_u64,
                   memory_limit);
            printf("MPI BAM sort mode policy: no -m => memory; otherwise all ranks use external when max(3 * rank_input_bytes + 256MiB) exceeds -m\n");
        }

        double fused_t0 = GetTime();
        int fused_ret = 0;
        if (memory_input) {
            fused_ret = use_external
                ? FusedBamExternalSortMPI(
                      reader, rank_body_sink, local_block_begin,
                      rank, comm_size, cmd_info->compress_level_,
                      memory_limit, sort_temp_prefix.c_str(), &stats)
                : FusedBamSortMPI(
                      reader, rank_body_sink, local_block_begin,
                      rank, comm_size, cmd_info->compress_level_,
                      memory_limit, &stats);
        } else {
            fused_ret = use_external
                ? FusedBamExternalSortMPI(
                      *input_backend, input_plan.rank_spans(),
                      input_plan.rank_block_count(), rank_body_sink,
                      local_block_begin, rank, comm_size,
                      cmd_info->compress_level_, memory_limit,
                      sort_temp_prefix.c_str(), &stats)
                : FusedBamSortMPI(
                      *input_backend, input_plan.rank_spans(),
                      input_plan.rank_block_count(), rank_body_sink,
                      local_block_begin, rank, comm_size,
                      cmd_info->compress_level_, memory_limit, &stats);
        }
        if (fused_ret != 0) {
            fprintf(stderr, "[rank %d] ERROR: MPI BAM sort fused body failed.\n", rank);
            local_ok = 0;
        }
        int global_ok = MpiSortAllRanksOk(local_ok);
        double fused_actual_cost_with_probe = GetTime() - fused_t0;
        if (!use_external) {
            stats.t_fused_actual = fused_actual_cost_with_probe;
        }
        double fused_actual_cost = fused_actual_cost_with_probe -
            (use_external ? stats.t_mpi_simulated : 0.0);
        if (fused_actual_cost < 0.0) fused_actual_cost = 0.0;
        double fused_cost = use_external
            ? stats.t_fused_total
            : fused_actual_cost;
        double fused_cost_max = MpiSortReduceMaxCost(fused_cost);
        double fused_actual_cost_max = MpiSortReduceMaxCost(fused_actual_cost);
        double fused_probe_cost_max = MpiSortReduceMaxCost(
            use_external ? stats.t_mpi_simulated : 0.0);
        if (rank == 0 && global_ok) {
            printf("Complete the 4.3 FusedBamSortMPI cost %lf\n", fused_cost_max);
            printf("FusedBamSortMPI reported cost model=%s\n",
                   use_external ? "simulated-high-speed-temp-io"
                                : "measured-memory-mode");
            printf("FusedBamSortMPI actual wall %lf\n", fused_actual_cost_max);
            if (use_external) {
                printf("FusedBamSortMPI timing probe overhead %lf\n",
                       fused_probe_cost_max);
            }
        }
        if (!global_ok) goto cleanup;

        //4.4 收集每个 rank 的输出大小
        double stage44_t0 = GetTime();
        if (rank_body_sink.Flush() != 0 ||
            rank_body_sink.size() > (uint64_t)LLONG_MAX) {
            local_ok = 0;
        }
        if (!MpiSortAllRanksOk(local_ok)) goto cleanup;
        local_body_size = (long long)rank_body_sink.size();
        {
            const int local_spool = rank_body_sink.is_memory() ? 0 : 1;
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
        double stage44_cost_max = MpiSortReduceMaxCost(stage44_cost);
        if (rank == 0) printf("Complete the 4.4 gather body sizes cost %lf\n", stage44_cost_max);

        //4.5a 处理输出 header 和整体布局，这部分计入处理时间
        double stage45_header_t0 = GetTime();
        long long output_body_start_ll = 0;
        long long output_file_size_ll = 0;
        if (rank == 0) {
            if (MpiSortUpdateHeaderCoordinate(hdr) != 0) {
                fprintf(stderr, "ERROR: failed to update BAM header SO:coordinate for sort output.\n");
                local_ok = 0;
            }
            if (local_ok &&
                MpiSortBuildBamHeaderMemory(hdr, cmd_info->compress_level_,
                                            &bam_header_mem, &bam_header_size) != 0) {
                fprintf(stderr, "ERROR: failed to build MPI SORT output BAM header in memory.\n");
                local_ok = 0;
            }
            if (local_ok) {
                unsigned long long final_size = (unsigned long long)bam_header_size +
                                               (unsigned long long)total_body_size +
                                               (unsigned long long)sizeof(kSortBgzfEofBlock);
                if (final_size > (unsigned long long)SIZE_MAX ||
                    final_size > (unsigned long long)LLONG_MAX ||
                    (unsigned long long)bam_header_size > (unsigned long long)LLONG_MAX) {
                    fprintf(stderr, "ERROR: MPI sort output buffer would exceed addressable memory.\n");
                    local_ok = 0;
                } else {
                    output_body_start = bam_header_size;
                    output_file_size = (size_t)final_size;
                    output_body_start_ll = (long long)output_body_start;
                    output_file_size_ll = (long long)output_file_size;
                }
            }
        }
        int stage45_header_ok = MpiSortAllRanksOk(local_ok);
        MPI_Bcast(&output_body_start_ll, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&output_file_size_ll, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        if (stage45_header_ok) {
            output_body_start = (size_t)output_body_start_ll;
            output_file_size = (size_t)output_file_size_ll;
        }
        double stage45_header_cost = GetTime() - stage45_header_t0;
        double stage45_header_cost_max = MpiSortReduceMaxCost(stage45_header_cost);
        if (rank == 0 && stage45_header_ok) {
            printf("Complete the 4.5a header/layout cost %lf\n", stage45_header_cost_max);
        }
        if (!stage45_header_ok) goto cleanup;

        //4.5b 分配每个 rank 的本地模拟写内存，这部分不计入处理时间
        double stage45_malloc_t0 = GetTime();
        const size_t simulation_chunk = 8u * 1024u * 1024u;
        simulated_write_size = rank_body_sink.is_memory()
            ? ((local_body_size > 0) ? (size_t)local_body_size : 0)
            : (size_t)std::min<long long>(local_body_size,
                                          (long long)simulation_chunk);
        if (!rank_body_sink.is_memory() && local_body_size > 0) {
            try {
                body_chunk.resize(simulated_write_size);
            } catch (...) {
                local_ok = 0;
            }
        }
        if (rank == 0) {
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
        int stage45_malloc_ok = MpiSortAllRanksOk(local_ok);
        double stage45_malloc_cost = GetTime() - stage45_malloc_t0;
        double stage45_malloc_cost_max = MpiSortReduceMaxCost(stage45_malloc_cost);
        if (rank == 0 && stage45_malloc_ok) {
            printf("Complete the 4.5b malloc simulated write memory cost %lf\n", stage45_malloc_cost_max);
        }
        if (!stage45_malloc_ok) goto cleanup;

        //4.6 各 rank 根据逻辑偏移模拟写入高速磁盘内存，只统计 memcpy 时间
        double sim_write_cost = 0.0;
        double body_source_read_cost = 0.0;
        size_t simulated_pos = 0;
        if (rank == 0 && output_body_start > 0) {
            double copy_t0 = GetTime();
            memcpy(simulated_write_mem + simulated_pos, bam_header_mem, bam_header_size);
            sim_write_cost += GetTime() - copy_t0;
            simulated_pos += bam_header_size;
        }
        if (local_body_size > 0) {
            const MemWriter *memory_writer = rank_body_sink.memory_writer();
            if (memory_writer) {
                double copy_t0 = GetTime();
                memcpy(simulated_write_mem + simulated_pos,
                       memory_writer->data, (size_t)local_body_size);
                sim_write_cost += GetTime() - copy_t0;
                simulated_pos += (size_t)local_body_size;
            } else {
                uint64_t copied = 0;
                const size_t destination_offset = simulated_pos;
                while (local_ok && copied < rank_body_sink.size()) {
                    const size_t count = (size_t)std::min<uint64_t>(
                        rank_body_sink.size() - copied, body_chunk.size());
                    double read_t0 = GetTime();
                    if (rank_body_sink.ReadAt(
                            copied, body_chunk.data(), count) != 0) {
                        local_ok = 0;
                        break;
                    }
                    body_source_read_cost += GetTime() - read_t0;
                    double copy_t0 = GetTime();
                    memcpy(simulated_write_mem + destination_offset,
                           body_chunk.data(), count);
                    sim_write_cost += GetTime() - copy_t0;
                    copied += count;
                }
                simulated_pos += simulated_write_size;
            }
        }
        if (simulated_pos > 0) {
            unsigned char *guard_ptr = (unsigned char *)simulated_write_mem;
            simulated_write_guard += guard_ptr[0];
            simulated_write_guard += guard_ptr[simulated_pos - 1];
        }
        stats.t_write += sim_write_cost;
        global_ok = MpiSortAllRanksOk(local_ok);
        double stage46_cost_max = MpiSortReduceMaxCost(sim_write_cost);
        double body_source_read_cost_max =
            MpiSortReduceMaxCost(body_source_read_cost);
        if (rank == 0 && global_ok) {
            printf("Complete the 4.6 simulated header/body distributed write cost %lf\n", stage46_cost_max);
            if (body_source_read_cost_max > 0.0) {
                printf("Rank body source read diagnostic cost %lf\n",
                       body_source_read_cost_max);
            }
        }
        if (simulated_write_guard == (unsigned long long)-1 && rank < 0) {
            fprintf(stderr, "unused simulated write guard %llu\n", simulated_write_guard);
        }
        if (!global_ok) goto cleanup;

        double body_cost = GetTime() - body_t0;
        (void)body_cost;
        if (rank == 0 && global_ok) {
            double body_counted_cost = stage41_cost_max + stage42_cost_max +
                                       fused_cost_max + stage44_cost_max +
                                       stage45_header_cost_max + stage46_cost_max;
            printf("Complete the total (4.1~4.6) cost %lf-----\n", body_counted_cost);
            double body_actual_counted_cost =
                stage41_cost_max + stage42_cost_max +
                fused_actual_cost_max + stage44_cost_max +
                stage45_header_cost_max + stage46_cost_max;
            printf("total actual wall (4.1~4.6) %lf-----\n",
                   body_actual_counted_cost);
        }

        //4.7 全局同步，统计每个 rank 的处理时间和统计数据，rank 0 汇总并打印最终统计结果
        double stage47_t0 = GetTime();
        long long local_long_stats[19] = {
            stats.input_blocks,
            stats.local_records,
            stats.received_records,
            stats.sample_records,
            stats.bgzf_blocks,
            stats.bucket_self_records,
            stats.bucket_remote_records,
            stats.bucket_self_raw_bytes,
            stats.bucket_remote_raw_bytes,
            stats.external_runs,
            stats.external_segments,
            stats.temp_read_bytes,
            stats.temp_write_bytes,
            stats.tracked_peak_bytes,
            stats.run_arena_bytes,
            stats.merge_fan_in,
            stats.consolidation_passes,
            stats.resident_run_records,
            stats.resident_run_raw_bytes
        };
        long long global_long_stats[19] = {};
        long long max_long_stats[19] = {};
        MPI_Reduce(local_long_stats, global_long_stats, 19, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(local_long_stats, max_long_stats, 19, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
        double local_double_stats[39] = {
            stats.t_setup,
            stats.t_extract_read,
            stats.t_extract_read_unhidden,
            stats.t_extract,
            stats.t_extract_prepare,
            stats.t_extract_merge,
            stats.t_local_sort,
            stats.t_sample,
            stats.t_partition,
            stats.t_bucket_count,
            stats.t_bucket_pack,
            stats.t_exchange,
            stats.t_mpi_exchange,
            stats.t_offset_fix,
            stats.t_final_sort,
            stats.t_compress,
            stats.t_write,
            stats.t_status_check,
            stats.t_cleanup,
            stats.t_fused_total,
            stats.t_temp_write_sim,
            stats.t_temp_read_sim,
            stats.t_run_sort,
            stats.t_run_bucket,
            stats.t_run_exchange,
            stats.t_run_merge,
            stats.t_merge_unhidden,
            stats.t_temp_write_actual,
            stats.t_temp_read_actual,
            stats.t_fused_actual,
            stats.t_mpi_simulated,
            stats.t_merge_simulated,
            stats.t_compress_simulated,
            stats.t_consolidation_simulated,
            stats.t_temp_open_actual,
            stats.t_merge_temp_read_sim,
            stats.t_consolidation_temp_read_sim,
            stats.t_consolidation_temp_write_sim,
            stats.t_cpe_calibration_wall
        };
        double global_double_stats[39] = {};
        MPI_Reduce(local_double_stats, global_double_stats, 39, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        MpiSortPrintRankStats(rank, comm_size, stats,
                              (size_t)rank_body_sink.size());
        if (rank == 0) {
            printf("FusedBamSortMPI finished. mode=%s ranks=%d in_blocks=%lld local_records=%lld received_records=%lld samples=%lld bgzf_blocks=%lld body_bytes=%lld\n",
                   stats.sort_mode == 1 ? "external" : "memory",
                   comm_size, global_long_stats[0], global_long_stats[1],
                   global_long_stats[2], global_long_stats[3], global_long_stats[4],
                   total_body_size);
            printf("  extract_sum=%.3f  local_sort_sum=%.3f  sample_sum=%.3f  partition_sum=%.3f  exchange_sum=%.3f  final_sort_sum=%.3f  compress_sum=%.3f  write_sum=%.3f\n",
                   global_double_stats[3], global_double_stats[6],
                   global_double_stats[7], global_double_stats[8],
                   global_double_stats[11], global_double_stats[14],
                   global_double_stats[15], global_double_stats[16]);
            printf("  pipeline_detail_sum setup=%.3f read=%.3f read_unhidden=%.3f extract_prepare=%.3f extract_merge=%.3f bucket_count=%.3f bucket_pack=%.3f mpi_exchange=%.3f offset_fix=%.3f status=%.3f cleanup=%.3f\n",
                   global_double_stats[0], global_double_stats[1], global_double_stats[2],
                   global_double_stats[4], global_double_stats[5],
                   global_double_stats[9], global_double_stats[10],
                   global_double_stats[12], global_double_stats[13],
                   global_double_stats[17], global_double_stats[18]);
            const long long global_bucket_raw = global_long_stats[7] + global_long_stats[8];
            const double global_self_raw_ratio =
                global_bucket_raw > 0 ? 100.0 * (double)global_long_stats[7] / (double)global_bucket_raw : 0.0;
            printf("  bucket_dist_sum self_records=%lld remote_records=%lld self_raw=%lld remote_raw=%lld self_raw_ratio=%.2f%%\n",
                   global_long_stats[5], global_long_stats[6],
                   global_long_stats[7], global_long_stats[8],
                   global_self_raw_ratio);
            printf("  external_sim_sum mpi=%.3f merge=%.3f compress=%.3f consolidation=%.3f temp_write=%.3f temp_read=%.3f merge_temp_read=%.3f consolidation_temp_read=%.3f consolidation_temp_write=%.3f\n",
                   global_double_stats[30], global_double_stats[31],
                   global_double_stats[32], global_double_stats[33],
                   global_double_stats[20], global_double_stats[21],
                   global_double_stats[35], global_double_stats[36],
                   global_double_stats[37]);
            printf("  external_actual_sum runs=%lld segments=%lld temp_open=%.3f temp_write=%.3f temp_read=%.3f exchange_wall=%.3f mpi_calls=%.3f merge_wall=%.3f compress_pipeline_wall=%.3f\n",
                   global_long_stats[9], global_long_stats[10],
                   global_double_stats[34],
                   global_double_stats[27], global_double_stats[28],
                   global_double_stats[11], global_double_stats[12],
                   global_double_stats[14], global_double_stats[15]);
            printf("  external_memory_sum temp_read_bytes=%lld temp_write_bytes=%lld resident_records=%lld resident_raw=%lld tracked_peak_max=%lld run_arena_max=%lld merge_fan_in_max=%lld consolidation_passes_max=%lld actual_wall_sum=%.3f\n",
                   global_long_stats[11], global_long_stats[12],
                   global_long_stats[17], global_long_stats[18],
                   max_long_stats[13], max_long_stats[14],
                   stats.sort_mode == 1 ? max_long_stats[15] : 0,
                   max_long_stats[16], global_double_stats[29]);
            const double core_stage_sum = stats.sort_mode == 1
                ? global_double_stats[19]
                : global_double_stats[0] + global_double_stats[2] +
                  global_double_stats[3] + global_double_stats[4] +
                  global_double_stats[5] + global_double_stats[6] +
                  global_double_stats[7] + global_double_stats[8] +
                  global_double_stats[11] + global_double_stats[13] +
                  global_double_stats[14] + global_double_stats[15] +
                  global_double_stats[17] + global_double_stats[18];
            printf("  fused_total_sum=%.3f  core_stage_sum=%.3f  unaccounted_sum=%.3f\n",
                   global_double_stats[19], core_stage_sum,
                   global_double_stats[19] - core_stage_sum);
        }
        double stage47_cost = GetTime() - stage47_t0;
        double stage47_cost_max = MpiSortReduceMaxCost(stage47_cost);
        if (rank == 0 && global_ok) {
            printf("Complete the 4.7 rank stats reduce/print cost %lf\n", stage47_cost_max);
        }

        double body_total_cost = GetTime() - body_total_t0;
        double body_total_cost_max = MpiSortReduceMaxCost(body_total_cost);
        if (rank == 0 && global_ok) {
            printf("444Complete the total body cost %lf\n", body_total_cost_max);
        }
    }


    //5. 最终输出不计入核心时间：memory 汇总到 rank 0，mpiio 分布式直写。
    {
        if (cmd_info->io_output_backend_ == "mpiio") {
            double write_t0 = GetTime();
            swbam::mpi::MpiFileOutput output;
            if (output.Open(cmd_info->out_file_name_,
                            (uint64_t)output_file_size) != 0) {
                local_ok = 0;
            }
            if (MpiSortAllRanksOk(local_ok)) {
                if (rank == 0 && bam_header_size > 0 &&
                    output.WriteAt(0, bam_header_mem,
                                   bam_header_size) != 0) {
                    local_ok = 0;
                }
                const size_t chunk_capacity = (size_t)std::min<uint64_t>(
                    rank_body_sink.size(), 8ull * 1024ull * 1024ull);
                try {
                    body_chunk.resize(chunk_capacity);
                } catch (...) {
                    local_ok = 0;
                }
                uint64_t copied = 0;
                while (local_ok && copied < rank_body_sink.size()) {
                    const size_t count = (size_t)std::min<uint64_t>(
                        rank_body_sink.size() - copied, body_chunk.size());
                    if (rank_body_sink.ReadAt(
                            copied, body_chunk.data(), count) != 0 ||
                        output.WriteAt(
                            (uint64_t)output_body_start +
                            (uint64_t)body_prefixes[(size_t)rank] + copied,
                            body_chunk.data(), count) != 0) {
                        local_ok = 0;
                    }
                    copied += count;
                }
                if (rank == 0 &&
                    output.WriteAt(
                        (uint64_t)output_body_start +
                        (uint64_t)total_body_size,
                        kSortBgzfEofBlock,
                        sizeof(kSortBgzfEofBlock)) != 0) {
                    local_ok = 0;
                }
            }
            const int writes_ok = MpiSortAllRanksOk(local_ok);
            if (writes_ok && output.Sync() != 0) local_ok = 0;
            if (output.Close() != 0) local_ok = 0;
            const double write_cost_max = MpiSortReduceMaxCost(
                GetTime() - write_t0);
            if (rank == 0 && local_ok) {
                printf("555MPI-IO distributed output cost %lf--\n",
                       write_cost_max);
            }
            if (!MpiSortAllRanksOk(local_ok)) goto cleanup;
        } else {
        double verify_alloc_t0 = GetTime();
        if (rank == 0) {
            output_file_mem = output_file_size ? (char *)malloc(output_file_size) : nullptr;
            if (output_file_size > 0 && !output_file_mem) {
                fprintf(stderr, "ERROR: failed to allocate final MPI sort output buffer.\n");
                local_ok = 0;
            }
            if (local_ok && bam_header_size > 0) {
                memcpy(output_file_mem, bam_header_mem, bam_header_size);
            }
        }
        double verify_alloc_cost = GetTime() - verify_alloc_t0;
        double verify_alloc_cost_max = MpiSortReduceMaxCost(verify_alloc_cost);
        if (rank == 0 && local_ok) {
            printf("555Prepare verification output memory cost %lf--\n", verify_alloc_cost_max);
        }
        if (!MpiSortAllRanksOk(local_ok)) goto cleanup;

        double verify_gather_t0 = GetTime();
        if (rank == 0) {
            if (local_body_size > 0) {
                if (rank_body_sink.ReadAt(
                        0, output_file_mem + output_body_start + local_prefix,
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
                    if (MpiSortRecvBytes(
                            src, 0,
                            output_file_mem + output_body_start +
                                body_prefixes[(size_t)src] + received,
                            count) != 0) {
                        fprintf(stderr,
                                "ERROR: failed to gather MPI sort rank %d "
                                "body into output memory.\n", src);
                        local_ok = 0;
                        break;
                    }
                    received += count;
                }
                if (!local_ok) break;
            }
            if (local_ok) {
                memcpy(output_file_mem + output_body_start + total_body_size,
                       kSortBgzfEofBlock, sizeof(kSortBgzfEofBlock));
            }
        } else if (local_body_size > 0) {
            const size_t chunk_capacity = (size_t)std::min<long long>(
                local_body_size, 8ll * 1024ll * 1024ll);
            try {
                body_chunk.resize(chunk_capacity);
            } catch (...) {
                local_ok = 0;
            }
            uint64_t sent = 0;
            while (local_ok && sent < rank_body_sink.size()) {
                const size_t count = (size_t)std::min<uint64_t>(
                    rank_body_sink.size() - sent, body_chunk.size());
                if (rank_body_sink.ReadAt(
                        sent, body_chunk.data(), count) != 0 ||
                    MpiSortSendBytes(
                        0, 0, (const char *)body_chunk.data(),
                        (long long)count) != 0) {
                    fprintf(stderr,
                            "[rank %d] ERROR: failed to send sort output "
                            "body to rank 0.\n", rank);
                    local_ok = 0;
                }
                sent += count;
            }
        }
        double verify_gather_cost = GetTime() - verify_gather_t0;
        double verify_gather_cost_max = MpiSortReduceMaxCost(verify_gather_cost);
        if (rank == 0 && local_ok) {
            printf("555Gather verification output memory cost %lf--\n", verify_gather_cost_max);
        }
        if (!MpiSortAllRanksOk(local_ok)) goto cleanup;

        double dump_t0 = 0.0;
        double dump_cost = 0.0;
        if (rank == 0) {
            dump_t0 = GetTime();
            if (MpiSortDumpMemoryToFile(cmd_info->out_file_name_, output_file_mem, output_file_size) != 0) {
                fprintf(stderr, "ERROR: failed to dump MPI sort output memory to %s\n",
                        cmd_info->out_file_name_.c_str());
                local_ok = 0;
            }
            dump_cost = GetTime() - dump_t0;
        }
        double dump_cost_max = MpiSortReduceMaxCost(dump_cost);
        if (rank == 0 && local_ok) printf("555Dump memory to output file cost %lf--\n", dump_cost_max);
        if (!MpiSortAllRanksOk(local_ok)) goto cleanup;
        }
    }

    exit_code = 0;


    //6. 清理资源
cleanup:
    {
        double close_t0 = GetTime();
        rank_body_sink.Close();
        if (output_file_mem) free(output_file_mem);
        if (bam_header_mem) free(bam_header_mem);
        if (simulated_write_mem) free(simulated_write_mem);
        if (hdr) sam_hdr_destroy(hdr);
        input_handle.Close();
        double close_cost = GetTime() - close_t0;
        double close_cost_max = MpiSortReduceMaxCost(close_cost);
        if (rank == 0) printf("666close the files cost %lf-----\n", close_cost_max);
    }
    return exit_code;
}
