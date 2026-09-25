#ifndef SWBAM_CPE_SAM_WRITE_PIPELINE_H
#define SWBAM_CPE_SAM_WRITE_PIPELINE_H

#include "swbam/cpe_pipeline.h"

namespace swbam {
namespace cpe {

struct CpeSamWriteTiming {
    long long batches;
    double format;
    double write;

    CpeSamWriteTiming() : batches(0), format(0.0), write(0.0) {}
};

// The record slots are borrowed only until Submit returns. Formatted text is
// retained in two buffers so the caller can flush one batch during its next
// CPE read kernel; Submit also flushes before swapping if no overlap is used.
class CpeSamWriteSession {
public:
    CpeSamWriteSession();
    ~CpeSamWriteSession();

    int Initialize(RankBodySink *sink, const sam_hdr_t *header);
    bam1_t **RecordSlots();
    size_t capacity() const { return BATCH_SIZE; }
    int Submit(size_t count);
    int SubmitWithPrefetch(size_t count, CpeBatchPrefetch prefetch,
                           void *context);
    int FlushPrevious();
    int Finish();
    void Close();
    const CpeSamWriteTiming &timing() const { return timing_; }

private:
    CpeSamWriteSession(const CpeSamWriteSession &);
    CpeSamWriteSession &operator=(const CpeSamWriteSession &);

    int FlushPending();
    RankBodySink *sink_;
    SamFormatBatch *current_;
    SamFormatBatch *pending_;
    bool has_pending_;
    CpeSamWriteTiming timing_;
};

} // namespace cpe
} // namespace swbam

#endif
