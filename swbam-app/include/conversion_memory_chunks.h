#ifndef SWBAM_APP_CONVERSION_MEMORY_CHUNKS_H
#define SWBAM_APP_CONVERSION_MEMORY_CHUNKS_H

#include "swbam/io.h"
#include "swbam/operators/bam_to_sam.h"
#include "swbam/operators/sam_to_bam.h"

namespace swbam_app {

struct ConversionChunkTiming {
    double core;
    double simulated_write;
    double preload;
    double persistence;
    double simulation_allocation;
    size_t chunks;
    size_t max_input;
    size_t max_output;
    uint64_t input_bytes;
    uint64_t copied_bytes;

    ConversionChunkTiming();
};

// The real body is persisted only after each pipeline has drained its CPEs.
// Without persistence, size/layout remain available but ReadAt is unsupported.
class ConversionChunkOutput : public swbam::RankBodySink {
public:
    ConversionChunkOutput();
    int Open(bool keep_output, const std::string &directory,
             const std::string &prefix);
    int Close();
    int Append(const void *data, size_t size);
    int Flush();
    uint64_t size() const { return size_; }
    int ReadAt(uint64_t offset, void *data, size_t size) const;
    bool is_memory() const { return false; }

private:
    swbam::SpoolRankBodySink spool_;
    bool keep_output_;
    uint64_t size_;
};

int RunChunkedBamToSam(MemReader *reader, ConversionChunkOutput *output,
                       const sam_hdr_t *header, size_t chunk_bytes,
                       swbam::operators::BamToSamMetrics *metrics,
                       ConversionChunkTiming *timing);

int RunChunkedSamToBam(const swbam::SamInputBackend &input,
                       uint64_t offset, size_t size,
                       ConversionChunkOutput *output,
                       const sam_hdr_t *header, size_t chunk_bytes, int level,
                       swbam::operators::SamToBamMetrics *metrics,
                       ConversionChunkTiming *timing);

} // namespace swbam_app

#endif
