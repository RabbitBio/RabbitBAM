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

struct SamWriteKernelSpec {
    // The entry accepts SamFormatBatch and writes per-core text/status.
    void *format_entry;
    size_t batch_capacity;
};

SamWriteKernelSpec DefaultSamWriteKernelSpec();

// The record slots are borrowed only until Submit returns. Formatted text is
// retained in two buffers so the caller can flush one batch during its next
// CPE read kernel; Submit also flushes before swapping if no overlap is used.
class CpeSamWriteSession {
public:
    CpeSamWriteSession();
    ~CpeSamWriteSession();

    int Initialize(RankBodySink *sink, const sam_hdr_t *header);
    int Initialize(RankBodySink *sink, const sam_hdr_t *header,
                   const SamWriteKernelSpec &kernel);
    bam1_t **RecordSlots();
    size_t capacity() const { return kernel_.batch_capacity; }
    // Convenience path for new applications. The existing RecordSlots path
    // avoids this pointer copy on performance-critical callers.
    int SubmitRecords(bam1_t *const *records, size_t count);
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
    SamWriteKernelSpec kernel_;
    CpeSamWriteTiming timing_;
};

} // namespace cpe
} // namespace swbam

#endif
