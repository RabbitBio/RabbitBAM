#ifndef SWBAM_CPE_BAM_READ_STEPS_H
#define SWBAM_CPE_BAM_READ_STEPS_H

#include "swbam/cpe_bam_parser.h"
#include "swbam/cpe_codec.h"

#include <climits>
#include <cstring>

// These steps are composed inside one CPE launch. Policy calls are resolved
// at compile time; no per-record dispatch or intermediate batch is needed.
inline int swbam_cpe_decode_bam_block(
        bam_block *compressed, bam_block *decoded, int cpe_id,
        uint64_t *alloc_cycles, uint64_t *inflate_cycles,
        uint64_t *crc_cycles) {
    struct libdeflate_decompressor *decoder =
        swbam_cpe_get_decompressor(cpe_id, alloc_cycles);
    if (!decoder) return -1;
    return swbam_cpe_decode_bgzf(compressed, decoded, decoder,
                                 inflate_cycles, crc_cycles) == 0 ? 0 : -2;
}

struct SwbamCpeWalkResult {
    int count;
    int read_result;
    bam1_t *last_record;
};

template <typename Action>
inline SwbamCpeWalkResult swbam_cpe_walk_bam_records(
        bam_block *decoded, int capacity, Action &action) {
    SwbamCpeWalkResult result = {0, -1, nullptr};
    while (result.count < capacity) {
        bam1_t *record = action.Prepare(result.count);
        result.last_record = record;
        if (!record) {
            result.read_result = -5;
            break;
        }
        result.read_result = swbam_cpe_read_bam_record(decoded, record);
        if (result.read_result < 0) break;
        action.Process(record, result.count);
        ++result.count;
    }
    if (result.count == capacity && decoded->pos < decoded->length)
        result.read_result = -6;
    return result;
}

class SwbamCpeScratchRecord {
public:
    SwbamCpeScratchRecord(unsigned char *data, size_t capacity)
        : data_(data), capacity_(capacity) {
        std::memset(&record_, 0, sizeof(record_));
        record_.mempolicy = BAM_USER_OWNS_DATA;
    }

    bam1_t *Prepare(int) {
        record_.data = data_;
        record_.m_data = capacity_ > UINT32_MAX
            ? UINT32_MAX : static_cast<uint32_t>(capacity_);
        record_.l_data = 0;
        record_.mempolicy = BAM_USER_OWNS_DATA;
        return &record_;
    }

private:
    unsigned char *data_;
    size_t capacity_;
    bam1_t record_;
};

#endif
