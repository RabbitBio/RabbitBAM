#include "codec_ratio_probe.h"
#include "swbam/swbam.h"
#include "../swbam-sdk/examples/storage-native.h"
#include <libdeflate.h>
#include <athread.h>

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

namespace {

uint32_t Read32(const unsigned char *p) {
    return uint32_t(p[0]) | (uint32_t(p[1]) << 8) |
           (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24);
}

struct Body {
    MemWriter writer = {nullptr, 0, 0};
    ~Body() { std::free(writer.data); }
};

bool Verify(struct libdeflate_decompressor *decoder, const bam_block &compressed,
            const bam_block &decoded, std::vector<unsigned char> *scratch) {
    const size_t size = compressed.length;
    if (size < 26 || (size_t(compressed.data[16]) |
                     (size_t(compressed.data[17]) << 8)) + 1 != size ||
        Read32(compressed.data + size - 4) != decoded.length) return false;
    size_t actual = 0;
    if (libdeflate_deflate_decompress(decoder, compressed.data + 18, size - 26,
            scratch->data(), scratch->size(), &actual) != LIBDEFLATE_SUCCESS ||
        actual != decoded.length || std::memcmp(scratch->data(), decoded.data, actual)) return false;
    return libdeflate_crc32(0, scratch->data(), actual) == Read32(compressed.data + size - 8);
}

int Run(const char *path, unsigned long long expected, int rank, int ranks) {
    swbam::mpi::MpiBamInput input;
    swbam::mpi::MpiBamInputPlan plan;
    int ok = input.Open(path, "posix", "") == 0;
    if (!swbam::mpi::AllRanksOk(ok) || input.PreparePlan(&plan)) return 1;
    unsigned long long local_input = 0;
    for (size_t i = 0; i < plan.rank_block_count(); ++i)
        local_input += plan.rank_spans()[i].compressed_size;
    Body body;
    const size_t reserve = static_cast<size_t>(local_input + local_input / 2 + (64u << 20));
    body.writer.data = static_cast<char *>(std::malloc(reserve));
    body.writer.capacity = body.writer.data ? reserve : 0;
    if (!swbam::mpi::AllRanksOk(body.writer.data != nullptr)) return 1;
    swbam::MemoryRankBodySink sink(&body.writer);
    swbam::operators::BamTransformMetrics metrics = {};
    ok = swbam::operators::RunBamTransformPipeline(*input.backend(),
        plan.rank_spans(), plan.rank_block_count(), &sink,
        swbam::DefaultBamFilterOptions(), 1, &metrics) == 0;
    if (!swbam::mpi::AllRanksOk(ok)) return 1;
    unsigned long long records = static_cast<unsigned long long>(metrics.total_records);
    unsigned long long total_records = 0;
    MPI_Allreduce(&records, &total_records, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (total_records != expected || metrics.total_records != metrics.kept_records) return 1;
    std::printf("[codec-pack rank=%d] records=%llu body_bytes=%zu reserve=%zu\n",
                rank, records, body.writer.size, reserve);
    std::fflush(stdout);

    swbam::BgzfBlockBatch decoded, original, optimized;
    ok = !decoded.Allocate(64) && !original.Allocate(64) && !optimized.Allocate(64);
    CodecRatioPara *parameters = nullptr;
    ok = ok && posix_memalign(reinterpret_cast<void **>(&parameters), 64,
                             64 * sizeof(CodecRatioPara)) == 0;
    struct libdeflate_decompressor *decoder = libdeflate_alloc_decompressor();
    if (!swbam::mpi::AllRanksOk(ok && decoder)) {
        std::free(parameters); libdeflate_free_decompressor(decoder); return 1;
    }
    std::vector<unsigned char> scratch(BGZF_MAX_BLOCK_SIZE);
    unsigned long long counts[4] = {}; // original body, optimized body, raw bytes, blocks
    size_t offset = 0;
    while (offset < body.writer.size) {
        bam_block compressed[64] = {};
        std::memset(parameters, 0, 64 * sizeof(*parameters));
        size_t active = 0;
        while (active < 64 && offset < body.writer.size) {
            if (body.writer.size - offset < 26) { ok = 0; break; }
            unsigned char *data = reinterpret_cast<unsigned char *>(body.writer.data + offset);
            const size_t bytes = (size_t(data[16]) | (size_t(data[17]) << 8)) + 1;
            if (bytes < 26 || bytes > body.writer.size - offset) { ok = 0; break; }
            compressed[active].data = data;
            compressed[active].length = static_cast<unsigned int>(bytes);
            parameters[active] = {&compressed[active], &decoded.blocks()[active],
                                 &original.blocks()[active], &optimized.blocks()[active], 0};
            ++active; offset += bytes;
        }
        if (!ok) break;
        __real_athread_spawn(reinterpret_cast<void *>(slave_codec_ratio_probe), parameters, 1);
        athread_join();
        for (size_t i = 0; i < active; ++i) {
            if (parameters[i].status ||
                optimized.blocks()[i].length != compressed[i].length ||
                std::memcmp(optimized.blocks()[i].data, compressed[i].data, compressed[i].length) ||
                !Verify(decoder, original.blocks()[i], decoded.blocks()[i], &scratch) ||
                !Verify(decoder, optimized.blocks()[i], decoded.blocks()[i], &scratch)) {
                std::fprintf(stderr, "codec verification failed rank=%d block=%llu status=%d\n",
                             rank, counts[3], parameters[i].status);
                ok = 0; break;
            }
            counts[0] += original.blocks()[i].length;
            counts[1] += optimized.blocks()[i].length;
            counts[2] += decoded.blocks()[i].length;
            ++counts[3];
        }
        if (!ok) break;
    }
    libdeflate_free_decompressor(decoder);
    std::free(parameters);
    if (!swbam::mpi::AllRanksOk(ok && counts[1] == body.writer.size)) return 1;
    unsigned long long totals[4] = {};
    MPI_Reduce(counts, totals, 4, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    swbam::MemoryBamOutput prefix;
    if (rank == 0) ok = swbam_usage::BuildBamPrefix(
        const_cast<sam_hdr_t *>(input.backend()->header()), 1, &prefix) == 0;
    if (!swbam::mpi::AllRanksOk(ok)) return 1;
    if (rank == 0) {
        const unsigned long long original_bytes = totals[0] + prefix.size() + 28;
        const unsigned long long optimized_bytes = totals[1] + prefix.size() + 28;
        std::printf("{\"result\":\"PASS\",\"ranks\":%d,\"records\":%llu,"
            "\"input_bytes\":%zu,\"input_blocks\":%zu,\"output_blocks\":%llu,"
            "\"raw_body_bytes\":%llu,\"original_body_bytes\":%llu,\"optimized_body_bytes\":%llu,"
            "\"common_header_bytes\":%zu,\"common_eof_bytes\":28,"
            "\"original_bam_bytes\":%llu,\"optimized_bam_bytes\":%llu,"
            "\"increase_percent\":%.9f,\"optimized_byte_identical\":true,"
            "\"both_decoded_and_crc_verified\":true}\n",
            ranks, total_records, input.backend()->size(), plan.blocks.size(), totals[3],
            totals[2], totals[0], totals[1], prefix.size(), original_bytes, optimized_bytes,
            100.0 * (static_cast<double>(optimized_bytes) / original_bytes - 1.0));
    }
    return 0;
}
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    int status = 1;
    athread_init();
    if (argc == 3) {
        char *end = nullptr; errno = 0;
        const unsigned long long expected = std::strtoull(argv[2], &end, 10);
        if (!errno && end && !*end && expected) status = Run(argv[1], expected, rank, ranks);
    } else if (rank == 0) {
        std::fprintf(stderr, "Usage: %s input.bam expected_records\n", argv[0]);
    }
    athread_halt();
    MPI_Finalize();
    return status;
}
