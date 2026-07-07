#ifndef SWBAM_RAW_BAM_WRITER_H
#define SWBAM_RAW_BAM_WRITER_H

#include "swbam/io.h"
#include "swbam/raw_bam.h"

#include <cstddef>

namespace swbam {
namespace cpe {

struct RawBamWriterMetrics {
    long long records;
    long long raw_bytes;
    long long header_raw_bytes;
    long long header_blocks;
    long long packed_blocks;
    long long compressed_bytes;
    double pack;
    double source;
    double kernel;
    double consume;
    double total;

    RawBamWriterMetrics();
};

class RawBamWriter : public RawBamRecordConsumer {
public:
    RawBamWriter();
    ~RawBamWriter();

    int Initialize(BamOutputBackend *output,
                   int compression_level,
                   size_t chunk_blocks = 256);
    int InitializeBam(BamOutputBackend *output,
                      const sam_hdr_t *header,
                      int compression_level,
                      size_t chunk_blocks = 256);
    int ConsumeRaw(const RawBamRecordView *records, size_t count);
    int Finish(bool append_bgzf_eof = true);
    const RawBamWriterMetrics &metrics() const;

private:
    class Impl;
    Impl *impl_;

    RawBamWriter(const RawBamWriter &);
    RawBamWriter &operator=(const RawBamWriter &);
};

typedef RawBamWriter RawBamMemoryWriter;

} // namespace cpe
} // namespace swbam

#endif
