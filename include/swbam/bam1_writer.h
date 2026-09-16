#ifndef SWBAM_BAM1_WRITER_H
#define SWBAM_BAM1_WRITER_H

#include "swbam/bam1.h"
#include "swbam/raw_bam_writer.h"

#include <cstddef>

namespace swbam {
namespace cpe {

struct Bam1WriterMetrics {
    long long records;
    long long bam1_data_bytes;
    long long encoded_bytes;
    double encode;

    Bam1WriterMetrics();
};

class Bam1Writer : public Bam1RecordConsumer {
public:
    Bam1Writer();
    ~Bam1Writer();

    int Initialize(BamOutputBackend *output,
                   int compression_level,
                   size_t chunk_blocks = 256);
    int InitializeBam(BamOutputBackend *output,
                      const sam_hdr_t *header,
                      int compression_level,
                      size_t chunk_blocks = 256);
    int ConsumeBam1(const bam1_t *const *records, size_t count);
    int Finish(bool append_bgzf_eof = true);

    const Bam1WriterMetrics &metrics() const;
    const RawBamWriterMetrics &raw_metrics() const;

private:
    class Impl;
    Impl *impl_;

    Bam1Writer(const Bam1Writer &);
    Bam1Writer &operator=(const Bam1Writer &);
};

} // namespace cpe
} // namespace swbam

#endif
