#ifndef SWBAM_MARKDUP_BOUNDARY_PROGRESS_H
#define SWBAM_MARKDUP_BOUNDARY_PROGRESS_H

#include <cstdint>
#include <vector>

namespace swbam_markdup {

inline int CompareBoundaryCoord(int32_t ar, int64_t ap, int32_t br, int64_t bp) {
    if (ar != br) return ar < br ? -1 : 1;
    if (ap != bp) return ap < bp ? -1 : 1;
    return 0;
}

// Future alignments beyond cut+2L cannot clip back into a cut+/-L boundary.
// Earlier ranks must advance or reach EOF: a long reference-consuming CIGAR
// can generate an end key far beyond that rank's aligned-start range.
// The clipping bound is the same longest-read contract as ordinary expiration.
template <typename Boundary>
bool BoundaryIsComplete(int32_t ref, int64_t coord, int max_read_length,
                        const std::vector<Boundary> &ranges,
                        const std::vector<long long> &progress) {
    if (ref <= 0) return true;
    if (max_read_length <= 0 || progress.size() != ranges.size() * 4) return false;
    const int64_t radius = (int64_t)max_read_length * 2;
    const int64_t upper = coord > INT64_MAX - radius ? INT64_MAX : coord + radius;
    for (size_t r = 0; r < ranges.size(); ++r) {
        if (!ranges[r].has_records) continue;
        if (CompareBoundaryCoord(ranges[r].first_ref, ranges[r].first_coord, ref, upper) > 0)
            continue;
        const long long *p = progress.data() + r * 4;
        if (!p[0]) return false;
        if (!p[1]) continue;
        if (p[2] <= 0 || CompareBoundaryCoord((int32_t)p[2], p[3], ref, upper) <= 0)
            return false;
    }
    return true;
}

} // namespace swbam_markdup

#endif
