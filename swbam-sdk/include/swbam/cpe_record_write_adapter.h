#ifndef SWBAM_CPE_RECORD_WRITE_ADAPTER_H
#define SWBAM_CPE_RECORD_WRITE_ADAPTER_H

#include "swbam/cpe_write_pipeline.h"

#include <cstddef>
#include <cstdint>

namespace swbam {
namespace cpe {

struct CpeRecordWriteMetrics {
    long long bgzf_blocks;
    double compress;
    double serialize;
    double alloc;
    double deflate;
    double footer;
    double other;
    double write;

    CpeRecordWriteMetrics();
};

// Adapts record plans and a RankBodySink to the shared BAM write pipeline.
// Submit joins the CPE before returning, so record pointers are not retained.
class CpeRecordWriteSession {
public:
    CpeRecordWriteSession();
    ~CpeRecordWriteSession();

    int Initialize(RankBodySink *sink, int compression_level);
    int Submit(const BamRecordPackBlock *blocks, size_t count);
    int Finish();
    const CpeRecordWriteMetrics &metrics() const;

private:
    class Impl;
    Impl *impl_;

    CpeRecordWriteSession(const CpeRecordWriteSession &);
    CpeRecordWriteSession &operator=(const CpeRecordWriteSession &);
};

} // namespace cpe
} // namespace swbam

#endif
