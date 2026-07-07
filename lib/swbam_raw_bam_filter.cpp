#include "swbam/raw_bam_filter.h"

#include <cstring>

namespace swbam {
namespace cpe {
namespace {

uint16_t ReadLe16(const unsigned char *data) {
    return (uint16_t)data[0] | ((uint16_t)data[1] << 8);
}

int32_t ReadLe32Signed(const unsigned char *data) {
    uint32_t value = (uint32_t)data[0] |
                     ((uint32_t)data[1] << 8) |
                     ((uint32_t)data[2] << 16) |
                     ((uint32_t)data[3] << 24);
    int32_t result = 0;
    memcpy(&result, &value, sizeof(result));
    return result;
}

bool Matches(const RawBamRecordView &record,
             const BamFilterOptions &filter) {
    if (!record.encoded || record.encoded_size < 36 ||
        record.block_size + 4 != record.encoded_size) {
        return false;
    }

    const unsigned char *raw = record.encoded;
    const int32_t tid = ReadLe32Signed(raw + 4);
    const int mapq = raw[13];
    const uint16_t flag = ReadLe16(raw + 18);
    const int32_t read_length = ReadLe32Signed(raw + 20);

    if (filter.min_mapq >= 0 && mapq < filter.min_mapq) return false;
    if (filter.max_mapq >= 0 && mapq > filter.max_mapq) return false;
    if (filter.require_flag != 0 &&
        (flag & filter.require_flag) != filter.require_flag) return false;
    if (filter.exclude_flag != 0 &&
        (flag & filter.exclude_flag) != 0) return false;
    if (filter.ref_tid != -2 && tid != filter.ref_tid) return false;
    if (filter.min_read_len >= 0 &&
        read_length < filter.min_read_len) return false;
    if (filter.max_read_len >= 0 &&
        read_length > filter.max_read_len) return false;
    return true;
}

} // namespace

RawBamFilterMetrics::RawBamFilterMetrics()
    : total_records(0), kept_records(0), dropped_records(0),
      total_bytes(0), kept_bytes(0), dropped_bytes(0), filter(0.0) {}

RawBamFilterConsumer::RawBamFilterConsumer(
        const BamFilterOptions &filter,
        RawBamRecordConsumer *downstream)
    : filter_(filter), downstream_(downstream),
      noop_(bam_filter_is_noop(filter)) {
    kept_.reserve(64 * 1024);
}

int RawBamFilterConsumer::ConsumeRaw(
        const RawBamRecordView *records, size_t count) {
    if (!downstream_ || (!records && count != 0)) return -1;
    const double filter_t0 = GetTime();

    long long total_bytes = 0;
    for (size_t i = 0; i < count; ++i) {
        if (!records[i].encoded || records[i].encoded_size < 36 ||
            records[i].block_size + 4 != records[i].encoded_size) {
            return -1;
        }
        total_bytes += records[i].encoded_size;
    }

    if (noop_) {
        metrics_.total_records += (long long)count;
        metrics_.kept_records += (long long)count;
        metrics_.total_bytes += total_bytes;
        metrics_.kept_bytes += total_bytes;
        metrics_.filter += GetTime() - filter_t0;
        return downstream_->ConsumeRaw(records, count);
    }

    kept_.clear();
    long long kept_bytes = 0;
    for (size_t i = 0; i < count; ++i) {
        if (Matches(records[i], filter_)) {
            kept_.push_back(records[i]);
            kept_bytes += records[i].encoded_size;
        }
    }

    metrics_.total_records += (long long)count;
    metrics_.kept_records += (long long)kept_.size();
    metrics_.dropped_records +=
        (long long)count - (long long)kept_.size();
    metrics_.total_bytes += total_bytes;
    metrics_.kept_bytes += kept_bytes;
    metrics_.dropped_bytes += total_bytes - kept_bytes;
    metrics_.filter += GetTime() - filter_t0;

    return downstream_->ConsumeRaw(
        kept_.empty() ? nullptr : kept_.data(), kept_.size());
}

} // namespace cpe
} // namespace swbam
