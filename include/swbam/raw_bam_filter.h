#ifndef SWBAM_RAW_BAM_FILTER_H
#define SWBAM_RAW_BAM_FILTER_H

#include "swbam/raw_bam.h"

#include <vector>

namespace swbam {
namespace cpe {

struct RawBamFilterMetrics {
    long long total_records;
    long long kept_records;
    long long dropped_records;
    long long total_bytes;
    long long kept_bytes;
    long long dropped_bytes;
    double filter;

    RawBamFilterMetrics();
};

class RawBamFilterConsumer : public RawBamRecordConsumer {
public:
    RawBamFilterConsumer(const BamFilterOptions &filter,
                         RawBamRecordConsumer *downstream);

    int ConsumeRaw(const RawBamRecordView *records, size_t count);
    const RawBamFilterMetrics &metrics() const { return metrics_; }

private:
    BamFilterOptions filter_;
    RawBamRecordConsumer *downstream_;
    bool noop_;
    RawBamFilterMetrics metrics_;
    std::vector<RawBamRecordView> kept_;
};

} // namespace cpe
} // namespace swbam

#endif
