#include "conversion_memory_chunks.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

namespace swbam_app {
namespace {

const size_t kCopyWindow = 64ull * 1024 * 1024;

class ChunkWorkspace {
public:
    ChunkWorkspace() : writer_(), sink_(&writer_), copy_(nullptr), guard_(0) {}
    ~ChunkWorkspace() { free(writer_.data); free(copy_); }

    int Initialize(ConversionChunkTiming *timing) {
        const double t0 = GetTime();
        copy_ = static_cast<char *>(malloc(kCopyWindow));
        timing->simulation_allocation += GetTime() - t0;
        return copy_ ? 0 : -1;
    }

    swbam::RankBodySink *Begin() { writer_.size = 0; return &sink_; }

    int Complete(ConversionChunkOutput *output, ConversionChunkTiming *timing) {
        const double copy_t0 = GetTime();
        for (size_t pos = 0; pos < writer_.size;) {
            const size_t count = std::min(kCopyWindow, writer_.size - pos);
            memcpy(copy_, writer_.data + pos, count);
            guard_ += static_cast<unsigned char>(copy_[0]);
            guard_ += static_cast<unsigned char>(copy_[count - 1]);
            pos += count;
        }
        timing->simulated_write += GetTime() - copy_t0;
        timing->copied_bytes += writer_.size;
        timing->max_output = std::max(timing->max_output, writer_.size);
        const double persist_t0 = GetTime();
        const int ret = output->Append(writer_.data, writer_.size);
        timing->persistence += GetTime() - persist_t0;
        ++timing->chunks;
        return ret;
    }

private:
    MemWriter writer_;
    swbam::MemoryRankBodySink sink_;
    char *copy_;
    volatile unsigned long long guard_;
};

template <typename Metrics, typename Value, size_t N>
void AddFields(Metrics *total, const Metrics &part,
               Value Metrics::*const (&fields)[N]) {
    for (size_t i = 0; i < N; ++i) total->*fields[i] += part.*fields[i];
}

void AddMetrics(swbam::operators::BamToSamMetrics *total,
                const swbam::operators::BamToSamMetrics &part) {
    typedef swbam::operators::BamToSamMetrics M;
    long long M::*const counts[] = {
        &M::input_blocks, &M::group_count, &M::total_records, &M::format_tiles};
    double M::*const times[] = {
        &M::t_decomp, &M::t_decomp_alloc, &M::t_decomp_inflate,
        &M::t_decomp_crc, &M::t_decomp_parse, &M::t_decomp_other,
        &M::t_format, &M::t_collect, &M::t_read, &M::t_write,
        &M::t_gather, &M::t_optimized_total, &M::t_alloc_init,
        &M::t_free_workspace};
    AddFields(total, part, counts);
    AddFields(total, part, times);
}

void AddMetrics(swbam::operators::SamToBamMetrics *total,
                const swbam::operators::SamToBamMetrics &part) {
    typedef swbam::operators::SamToBamMetrics M;
    long long M::*const counts[] = {
        &M::input_chunks, &M::chunk_groups, &M::total_records,
        &M::compress_groups, &M::bgzf_blocks,
        &M::parse_fast_records, &M::parse_fallback_records};
    double M::*const times[] = {
        &M::t_split, &M::t_copy_count, &M::t_parse, &M::t_parse_core,
        &M::t_parse_aux, &M::t_parse_cg, &M::t_parse_fallback,
        &M::t_parse_other, &M::t_pack, &M::t_compress,
        &M::t_compress_serialize, &M::t_compress_alloc,
        &M::t_compress_deflate, &M::t_compress_footer,
        &M::t_compress_other, &M::t_write, &M::t_gather,
        &M::t_optimized_total, &M::t_alloc_init, &M::t_setup_reset,
        &M::t_compress_setup, &M::t_status_check, &M::t_free_workspace};
    AddFields(total, part, counts);
    AddFields(total, part, times);
}

} // namespace

ConversionChunkTiming::ConversionChunkTiming()
    : core(0), simulated_write(0), preload(0), persistence(0),
      simulation_allocation(0), chunks(0), max_input(0), max_output(0),
      input_bytes(0), copied_bytes(0) {}

ConversionChunkOutput::ConversionChunkOutput() : keep_output_(false), size_(0) {}

int ConversionChunkOutput::Open(bool keep_output, const std::string &directory,
                                const std::string &prefix) {
    size_ = 0;
    keep_output_ = keep_output;
    return keep_output ? spool_.OpenTemporary(directory, prefix) : 0;
}

int ConversionChunkOutput::Close() { return spool_.Close(); }

int ConversionChunkOutput::Append(const void *data, size_t size) {
    if ((size && !data) || size > std::numeric_limits<uint64_t>::max() - size_) {
        errno = EOVERFLOW;
        return -1;
    }
    if (keep_output_ && spool_.Append(data, size) != 0) return -1;
    size_ += size;
    return 0;
}

int ConversionChunkOutput::Flush() { return keep_output_ ? spool_.Flush() : 0; }

int ConversionChunkOutput::ReadAt(uint64_t offset, void *data, size_t size) const {
    if (!keep_output_) { errno = ENOTSUP; return -1; }
    return spool_.ReadAt(offset, data, size);
}

int RunChunkedBamToSam(MemReader *reader, ConversionChunkOutput *output,
                       const sam_hdr_t *header, size_t chunk_bytes,
                       swbam::operators::BamToSamMetrics *metrics,
                       ConversionChunkTiming *timing) {
    if (!reader || !output || !header || !metrics || !timing ||
        !chunk_bytes || reader->pos > reader->size ||
        (reader->size && !reader->base)) return -1;
    *metrics = swbam::operators::BamToSamMetrics();
    ChunkWorkspace workspace;
    if (workspace.Initialize(timing) != 0) return -1;
    while (reader->pos < reader->size) {
        const double t0 = GetTime();
        const size_t begin = reader->pos;
        size_t end = begin;
        while (end < reader->size) {
            if (reader->size - end < BLOCK_HEADER_LENGTH) return -1;
            const unsigned char *src =
                reinterpret_cast<const unsigned char *>(reader->base + end);
            const size_t length = static_cast<size_t>(src[16] | src[17] << 8) + 1;
            if (length < BLOCK_HEADER_LENGTH || length > BGZF_MAX_BLOCK_SIZE ||
                length > reader->size - end) return -1;
            if (end > begin && length > chunk_bytes - (end - begin)) break;
            end += length;
            if (end - begin >= chunk_bytes) break;
        }
        MemReader part = {};
        part.base = reader->base + begin;
        part.size = end - begin;
        swbam::operators::BamToSamMetrics current = {};
        if (swbam::operators::RunBamToSamPipeline(
                &part, workspace.Begin(), header, &current) != 0 ||
            part.pos != part.size) return -1;
        AddMetrics(metrics, current);
        reader->pos = end;
        timing->core += GetTime() - t0;
        timing->input_bytes += part.size;
        timing->max_input = std::max(timing->max_input, part.size);
        if (workspace.Complete(output, timing) != 0) return -1;
    }
    metrics->t_optimized_total = timing->core;
    return 0;
}

int RunChunkedSamToBam(const swbam::SamInputBackend &input,
                       uint64_t offset, size_t size,
                       ConversionChunkOutput *output,
                       const sam_hdr_t *header, size_t chunk_bytes, int level,
                       swbam::operators::SamToBamMetrics *metrics,
                       ConversionChunkTiming *timing) {
    if (!output || !header || !metrics || !timing || !chunk_bytes ||
        offset > input.size() || size > input.size() - offset) return -1;
    *metrics = swbam::operators::SamToBamMetrics();
    ChunkWorkspace workspace;
    if (workspace.Initialize(timing) != 0) return -1;
    std::vector<char> storage;
    size_t processed = 0;
    while (processed < size) {
        const size_t requested = std::min(chunk_bytes, size - processed);
        const double preload_t0 = GetTime();
        if (input.ReadRange(offset + processed, requested, &storage) != 0 ||
            storage.size() != requested) return -1;
        timing->preload += GetTime() - preload_t0;
        const double t0 = GetTime();
        size_t count = requested;
        if (count < size - processed) {
            while (count && storage[count - 1] != '\n') --count;
            if (!count) {
                fprintf(stderr, "ERROR: SAM line exceeds conversion chunk size.\n");
                return -1;
            }
        }
        MemReader part = {};
        part.base = storage.data();
        part.size = count;
        swbam::operators::SamToBamMetrics current = {};
        if (swbam::operators::RunSamToBamPipeline(
                &part, workspace.Begin(), header, level, &current) != 0 ||
            part.pos != part.size) return -1;
        AddMetrics(metrics, current);
        processed += count;
        timing->core += GetTime() - t0;
        timing->input_bytes += count;
        timing->max_input = std::max(timing->max_input, requested);
        if (workspace.Complete(output, timing) != 0) return -1;
    }
    metrics->t_optimized_total = timing->core;
    return 0;
}

} // namespace swbam_app
