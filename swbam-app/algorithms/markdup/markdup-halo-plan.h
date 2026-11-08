#ifndef SWBAM_MARKDUP_HALO_PLAN_H
#define SWBAM_MARKDUP_HALO_PLAN_H

#include "algorithm_types.h"
#include <algorithm>
#include <cstdint>
#include <vector>

namespace swbam_markdup {

struct HaloInterval {
    int32_t ref;
    int64_t low;
    int64_t high;
};

inline int64_t HaloOffset(int64_t coord, int64_t offset) {
    if (offset > 0 && coord > INT64_MAX - offset) return INT64_MAX;
    if (offset < 0 && coord < INT64_MIN - offset) return INT64_MIN;
    return coord + offset;
}

inline void MergeHalos(std::vector<HaloInterval> *halos) {
    std::sort(halos->begin(), halos->end(), [](const HaloInterval &a, const HaloInterval &b) {
        return a.ref < b.ref || (a.ref == b.ref && a.low < b.low);
    });
    size_t count = 0;
    for (const auto &halo : *halos) {
        if (count && (*halos)[count - 1].ref == halo.ref &&
            halo.low <= HaloOffset((*halos)[count - 1].high, 1)) {
            (*halos)[count - 1].high = std::max((*halos)[count - 1].high, halo.high);
        } else {
            (*halos)[count++] = halo;
        }
    }
    halos->resize(count);
}

inline bool HaloContainsPosition(const std::vector<HaloInterval> &halos,
                                 int32_t ref, int64_t pos) {
    for (const auto &halo : halos)
        if (halo.ref == ref && pos >= halo.low && pos <= halo.high) return true;
    return false;
}

// These roles match md_make_single_key/md_make_pair_key, including FF/RR pairs.
// Position support assumes the corresponding unclipped displacement <= L.
inline bool KeyUsesEnd(const MpiMarkdupKeyShared &key, bool *end) {
    if (key.single) {
        if (key.orientation != 2 && key.orientation != 3) return false;
        *end = key.orientation == 3;
        return true;
    }
    if (key.leftmost != 11 && key.leftmost != 13) return false;
    if (key.orientation == 5) *end = key.leftmost == 13;
    else if (key.orientation == 7) *end = key.leftmost == 11;
    else if (key.orientation == 2 || key.orientation == 3) *end = key.leftmost == 13;
    else return false;
    return true;
}

inline bool KeyPositionSupport(const MpiMarkdupKeyShared &key, int64_t length,
                               int64_t *low, int64_t *high) {
    bool end = false;
    if (length <= 0 || key.this_ref <= 0 || !KeyUsesEnd(key, &end)) return false;
    *low = end ? HaloOffset(key.this_coord, -length) : key.this_coord;
    *high = end ? key.this_coord : HaloOffset(key.this_coord, length);
    return true;
}

inline bool HaloCoversKey(const std::vector<HaloInterval> &halos,
                          const MpiMarkdupKeyShared &key, int64_t length) {
    int64_t low = 0, high = 0;
    if (!KeyPositionSupport(key, length, &low, &high)) return false;
    for (const auto &halo : halos)
        if (halo.ref == key.this_ref && halo.low <= low && halo.high >= high) return true;
    return false;
}

} // namespace swbam_markdup
#endif
