#include "conversion_memory_chunks.h"

#include <cassert>
#include <cstring>
#include <string>
#include <vector>

namespace {

const sam_hdr_t *TestHeader() { return reinterpret_cast<const sam_hdr_t *>(1); }

class TestSamInput : public swbam::SamInputBackend {
public:
    explicit TestSamInput(const std::string &text) : text_(text) {}
    int Open(const std::string &) { return 0; }
    void Close() {}
    const sam_hdr_t *header() const { return TestHeader(); }
    size_t size() const { return text_.size(); }
    uint64_t body_offset() const { return 0; }
    int SplitRanges(int, std::vector<long long> *, std::vector<long long> *) const {
        return -1;
    }
    int ReadRange(uint64_t offset, size_t size, std::vector<char> *data) const {
        if (offset > text_.size() || size > text_.size() - offset) return -1;
        data->assign(text_.begin() + offset, text_.begin() + offset + size);
        return 0;
    }
private:
    std::string text_;
};

std::vector<char> TestBam() {
    std::vector<char> bytes(24 * 7, 0);
    for (size_t i = 0; i < 7; ++i) {
        bytes[i * 24 + 16] = 23;
        bytes[i * 24 + 18] = 'a' + i;
    }
    return bytes;
}

} // namespace

// Fake kernels test the adapter's boundaries, ownership, and aggregation only.
namespace swbam { namespace operators {
int RunBamToSamPipeline(MemReader *reader, RankBodySink *sink,
                        const sam_hdr_t *, BamToSamMetrics *metrics) {
    *metrics = BamToSamMetrics();
    while (reader->pos < reader->size) {
        const size_t pos = reader->pos;
        const size_t length = static_cast<unsigned char>(reader->base[pos + 16]) + 1;
        assert(length == 24 && length <= reader->size - pos);
        if (sink->Append(reader->base + pos + 18, 1) != 0) return -1;
        reader->pos += length;
        ++metrics->total_records;
        ++metrics->input_blocks;
    }
    metrics->format_tiles = 2;
    metrics->t_format = 0.5;
    return 0;
}

int RunSamToBamPipeline(MemReader *reader, RankBodySink *sink,
                        const sam_hdr_t *, int, SamToBamMetrics *metrics) {
    *metrics = SamToBamMetrics();
    for (size_t i = 0; i < reader->size; ++i)
        if (reader->base[i] == '\n') ++metrics->total_records;
    if (reader->size && reader->base[reader->size - 1] != '\n')
        ++metrics->total_records;
    if (sink->Append(reader->base, reader->size) != 0) return -1;
    reader->pos = reader->size;
    metrics->t_compress = 0.25;
    return 0;
}
} }

int main() {
    std::vector<char> bam = TestBam();
    MemReader reader = {};
    reader.base = bam.data();
    reader.size = bam.size();
    swbam_app::ConversionChunkOutput output;
    assert(output.Open(true, "/tmp", "swbam-conversion-regression") == 0);
    swbam_app::ConversionChunkTiming timing;
    swbam::operators::BamToSamMetrics bam_metrics = {};
    assert(swbam_app::RunChunkedBamToSam(&reader, &output, TestHeader(), 50,
                                       &bam_metrics, &timing) == 0);
    assert(reader.pos == bam.size() && bam_metrics.total_records == 7);
    assert(bam_metrics.input_blocks == 7 && timing.chunks == 4);
    assert(bam_metrics.format_tiles == 8 && bam_metrics.t_format == 2.0);
    assert(timing.input_bytes == bam.size() && timing.copied_bytes == 7);
    char actual[7];
    assert(output.Flush() == 0 && output.ReadAt(0, actual, sizeof(actual)) == 0);
    assert(memcmp(actual, "abcdefg", 7) == 0);
    assert(output.Close() == 0);

    const std::string text = "HEADER\nfirst\nsecond-long\nlast";
    TestSamInput sam(text);
    assert(output.Open(true, "/tmp", "swbam-conversion-regression") == 0);
    timing = swbam_app::ConversionChunkTiming();
    swbam::operators::SamToBamMetrics sam_metrics = {};
    assert(swbam_app::RunChunkedSamToBam(sam, 7, text.size() - 7,
        &output, TestHeader(), 14, 1, &sam_metrics, &timing) == 0);
    assert(sam_metrics.total_records == 3 && timing.chunks == 3);
    assert(timing.input_bytes == text.size() - 7 &&
           timing.copied_bytes == text.size() - 7);
    assert(sam_metrics.t_compress == 0.75);
    std::vector<char> restored(output.size());
    assert(output.ReadAt(0, restored.data(), restored.size()) == 0);
    assert(std::string(restored.begin(), restored.end()) == text.substr(7));
    assert(output.Close() == 0);

    assert(output.Open(false, "/tmp", "unused") == 0);
    timing = swbam_app::ConversionChunkTiming();
    assert(swbam_app::RunChunkedSamToBam(sam, 7, text.size() - 7,
        &output, TestHeader(), 14, 1, &sam_metrics, &timing) == 0);
    assert(output.size() == text.size() - 7);
    assert(output.ReadAt(0, actual, 1) != 0);
    assert(output.Close() == 0);

    assert(output.Open(false, "/tmp", "unused") == 0);
    timing = swbam_app::ConversionChunkTiming();
    assert(swbam_app::RunChunkedSamToBam(sam, 7, text.size() - 7,
        &output, TestHeader(), 4, 1, &sam_metrics, &timing) != 0);
    reader.pos = 0;
    reader.size = 10;
    assert(swbam_app::RunChunkedBamToSam(&reader, &output, TestHeader(), 50,
                                       &bam_metrics, &timing) != 0);
    return 0;
}
