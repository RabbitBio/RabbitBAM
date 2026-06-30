#include "swbam_mpi.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <new>
#include <vector>

namespace {

static uint64_t MdMix64(uint64_t value) {
    value += 0x9e3779b97f4a7c15ull;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ull;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebull;
    return value ^ (value >> 31);
}

static uint64_t MdRotl64(uint64_t value, int shift) {
    return (value << shift) | (value >> (64 - shift));
}

static uint64_t MdHashKeyFlat(const MpiMarkdupKeyShared &key) {
    uint64_t hash =
        (uint64_t)key.this_coord * 0x9e3779b97f4a7c15ull;
    hash ^= MdRotl64(
        (uint64_t)key.other_coord * 0xbf58476d1ce4e5b9ull, 23);
    hash ^= ((uint64_t)(uint32_t)key.this_ref << 32) |
            (uint64_t)(uint32_t)key.other_ref;
    hash ^= (uint64_t)(uint8_t)key.single << 8 |
            (uint64_t)(uint8_t)key.leftmost << 16 |
            (uint64_t)(uint8_t)key.orientation << 24;
    hash = MdMix64(hash);
    return hash ? hash : 1;
}

static bool MdKeyEqualFlat(const MpiMarkdupKeyShared &a,
                           const MpiMarkdupKeyShared &b) {
    return a.this_coord == b.this_coord &&
           a.other_coord == b.other_coord &&
           a.this_ref == b.this_ref &&
           a.other_ref == b.other_ref &&
           a.single == b.single &&
           a.leftmost == b.leftmost &&
           a.orientation == b.orientation;
}

static size_t MdNextTableCapacity(size_t expected) {
    if (expected == 0) return 0;
    if (expected > (SIZE_MAX / 2) - 16) return 0;
    size_t wanted = expected * 2 + 16;
    size_t cap = 16;
    while (cap < wanted) {
        if (cap > SIZE_MAX / 2) return 0;
        cap <<= 1;
    }
    return cap;
}

struct MdFlatBestTable {
    MpiMarkdupKeyShared *keys;
    uint64_t *hashes;
    size_t *best;
    unsigned char *used;
    size_t capacity;
    size_t mask;
    size_t filled;

    MdFlatBestTable()
        : keys(nullptr), hashes(nullptr), best(nullptr), used(nullptr),
          capacity(0), mask(0), filled(0) {}

    ~MdFlatBestTable() {
        Free();
    }

    int Init(size_t expected, double *alloc_time, double *init_time) {
        Free();
        double init_t0 = GetTime();
        capacity = MdNextTableCapacity(expected);
        if (expected == 0) {
            if (init_time) *init_time += GetTime() - init_t0;
            return 0;
        }
        if (capacity == 0) {
            if (init_time) *init_time += GetTime() - init_t0;
            return -1;
        }
        mask = capacity - 1;
        if (init_time) *init_time += GetTime() - init_t0;
        double alloc_t0 = GetTime();
        keys = (MpiMarkdupKeyShared *)aligned_alloc_custom(
            64, capacity * sizeof(MpiMarkdupKeyShared));
        hashes = (uint64_t *)aligned_alloc_custom(
            64, capacity * sizeof(uint64_t));
        best = (size_t *)aligned_alloc_custom(
            64, capacity * sizeof(size_t));
        used = aligned_alloc_custom(64, capacity);
        if (alloc_time) *alloc_time += GetTime() - alloc_t0;
        if (!keys || !hashes || !best || !used) return -1;
        init_t0 = GetTime();
        memset(used, 0, capacity);
        filled = 0;
        if (init_time) *init_time += GetTime() - init_t0;
        return 0;
    }

    int Prepare(size_t expected, double *alloc_time,
                double *init_time) {
        double init_t0 = GetTime();
        const size_t wanted = MdNextTableCapacity(expected);
        if (expected > 0 && wanted == 0) {
            if (init_time) *init_time += GetTime() - init_t0;
            return -1;
        }
        if (wanted > capacity) {
            if (init_time) *init_time += GetTime() - init_t0;
            Free();
            capacity = wanted;
            mask = capacity - 1;
            const double alloc_t0 = GetTime();
            keys = (MpiMarkdupKeyShared *)aligned_alloc_custom(
                64, capacity * sizeof(MpiMarkdupKeyShared));
            hashes = (uint64_t *)aligned_alloc_custom(
                64, capacity * sizeof(uint64_t));
            best = (size_t *)aligned_alloc_custom(
                64, capacity * sizeof(size_t));
            used = aligned_alloc_custom(64, capacity);
            if (alloc_time) *alloc_time += GetTime() - alloc_t0;
            if (!keys || !hashes || !best || !used) return -1;
            init_t0 = GetTime();
        }
        if (capacity > 0) memset(used, 0, capacity);
        filled = 0;
        if (init_time) *init_time += GetTime() - init_t0;
        return 0;
    }

    void Free() {
        if (keys) aligned_free_custom((unsigned char *)keys);
        if (hashes) aligned_free_custom((unsigned char *)hashes);
        if (best) aligned_free_custom((unsigned char *)best);
        if (used) aligned_free_custom(used);
        keys = nullptr;
        hashes = nullptr;
        best = nullptr;
        used = nullptr;
        capacity = 0;
        mask = 0;
        filled = 0;
    }

    int FindOrInsert(const MpiMarkdupKeyShared &key, uint64_t hash,
                     size_t initial_best, size_t *slot,
                     int *inserted) {
        if (capacity == 0) return -1;
        size_t pos = (size_t)hash & mask;
        for (;;) {
            if (!used[pos]) {
                if (filled + 1 >= capacity) return -1;
                used[pos] = 1;
                hashes[pos] = hash;
                keys[pos] = key;
                best[pos] = initial_best;
                filled++;
                *slot = pos;
                *inserted = 1;
                return 0;
            }
            if (hashes[pos] == hash && MdKeyEqualFlat(keys[pos], key)) {
                *slot = pos;
                *inserted = 0;
                return 0;
            }
            pos = (pos + 1) & mask;
        }
    }

    bool FindExisting(const MpiMarkdupKeyShared &key, uint64_t hash,
                      size_t *slot) const {
        if (capacity == 0) return false;
        size_t pos = (size_t)hash & mask;
        for (;;) {
            if (!used[pos]) return false;
            if (hashes[pos] == hash && MdKeyEqualFlat(keys[pos], key)) {
                *slot = pos;
                return true;
            }
            pos = (pos + 1) & mask;
        }
    }
};

static int MdCompareQnameStream(
        const MpiMarkdupCandidateShared &a,
        const MpiMarkdupCandidateShared &b,
        const std::vector<unsigned char> &qnames) {
    if (a.qname_offset > qnames.size() ||
        a.qname_len > qnames.size() - (size_t)a.qname_offset ||
        b.qname_offset > qnames.size() ||
        b.qname_len > qnames.size() - (size_t)b.qname_offset) {
        return 0;
    }
    size_t common = std::min((size_t)a.qname_len,
                             (size_t)b.qname_len);
    int cmp = common ? memcmp(
        qnames.data() + (size_t)a.qname_offset,
        qnames.data() + (size_t)b.qname_offset,
        common) : 0;
    if (cmp != 0) return cmp;
    if (a.qname_len < b.qname_len) return -1;
    if (a.qname_len > b.qname_len) return 1;
    return 0;
}

static bool MdPairBetterStream(
        const MpiMarkdupCandidateShared &a,
        const MpiMarkdupCandidateShared &b,
        const std::vector<unsigned char> &qnames) {
    if (a.qc_fail != b.qc_fail) return a.qc_fail < b.qc_fail;
    if (a.score != b.score) return a.score > b.score;
    int qcmp = MdCompareQnameStream(a, b, qnames);
    if (qcmp != 0) return qcmp < 0;
    return a.global_order < b.global_order;
}

static int MdPushDuplicateStream(
        const MpiMarkdupCandidateShared &candidate, int comm_size,
        std::vector<std::vector<uint64_t> > *duplicates_by_source) {
    const int source = candidate.source_rank;
    if (source < 0 || source >= comm_size) return -1;
    (*duplicates_by_source)[(size_t)source].push_back(candidate.ordinal);
    return 0;
}

} // namespace

struct MpiMarkdupStreamingWindow {
    int comm_size;
    int64_t max_read_length;
    int32_t first_ref;
    int64_t first_coord;
    int32_t last_ref;
    int64_t last_coord;
    size_t peak_working_bytes;
    std::vector<MpiMarkdupCandidateShared> carry;
    std::vector<unsigned char> carry_qnames;
    std::vector<MpiMarkdupCandidateShared> next_carry;
    std::vector<unsigned char> next_qnames;
    MdFlatBestTable pair_best;
    MdFlatBestTable single_best;
    std::vector<size_t> pair_slots;
    std::vector<size_t> single_slots;
    std::vector<size_t> marker_indices;
    std::vector<size_t> marker_carry_indices;

    MpiMarkdupStreamingWindow()
        : comm_size(1), max_read_length(300),
          first_ref(0), first_coord(0),
          last_ref(0), last_coord(0),
          peak_working_bytes(0) {}
};

namespace {

static bool MdWindowBoundaryKey(
        const MpiMarkdupStreamingWindow *window,
        const MpiMarkdupKeyShared &key) {
    if (window->first_ref > 0 && key.this_ref == window->first_ref &&
        key.this_coord >=
            window->first_coord - window->max_read_length &&
        key.this_coord <=
            window->first_coord + window->max_read_length) {
        return true;
    }
    return window->last_ref > 0 && key.this_ref == window->last_ref &&
           key.this_coord >=
               window->last_coord - window->max_read_length &&
           key.this_coord <=
               window->last_coord + window->max_read_length;
}

static bool MdWindowExpired(
        const MpiMarkdupStreamingWindow *window,
        const MpiMarkdupKeyShared &key,
        int progress_tid, int progress_pos) {
    if (progress_tid < 0) return true;
    const int32_t progress_ref = progress_tid + 1;
    const int64_t progress_coord = (int64_t)progress_pos + 1;
    if (key.this_ref < progress_ref) return true;
    if (key.this_ref > progress_ref) return false;
    return key.this_coord + window->max_read_length <= progress_coord;
}

static bool MdWindowKeepKey(
        const MpiMarkdupStreamingWindow *window,
        const MpiMarkdupKeyShared &key,
        int progress_tid, int progress_pos) {
    return MdWindowBoundaryKey(window, key) ||
           !MdWindowExpired(window, key, progress_tid, progress_pos);
}

static size_t MdFlatTableBytes(const MdFlatBestTable &table) {
    const size_t stride = sizeof(MpiMarkdupKeyShared) +
        sizeof(uint64_t) + sizeof(size_t) + 1;
    return table.capacity > SIZE_MAX / stride
        ? SIZE_MAX : table.capacity * stride;
}

} // namespace

MpiMarkdupStreamingWindow *MpiMarkdupStreamingWindowCreate(
        int comm_size, int max_read_length,
        int range_first_tid, int range_first_pos,
        int range_last_tid, int range_last_pos) {
    if (comm_size <= 0 || max_read_length <= 0) return nullptr;
    MpiMarkdupStreamingWindow *window =
        new (std::nothrow) MpiMarkdupStreamingWindow();
    if (!window) return nullptr;
    window->comm_size = comm_size;
    window->max_read_length = max_read_length;
    window->first_ref = range_first_tid >= 0
        ? range_first_tid + 1 : 0;
    window->first_coord = range_first_pos >= 0
        ? (int64_t)range_first_pos + 1 : 0;
    window->last_ref = range_last_tid >= 0
        ? range_last_tid + 1 : 0;
    window->last_coord = range_last_pos >= 0
        ? (int64_t)range_last_pos + 1 : 0;
    return window;
}

void MpiMarkdupStreamingWindowDestroy(
        MpiMarkdupStreamingWindow *window) {
    delete window;
}

int MpiMarkdupStreamingWindowProcess(
        MpiMarkdupStreamingWindow *window,
        const std::vector<MpiMarkdupCandidateShared> &owner_candidates,
        const std::vector<unsigned char> &owner_qnames,
        int progress_tid, int progress_pos,
        std::vector<std::vector<uint64_t> > *duplicates_by_source,
        MpiMarkdupStats *stats) {
    if (!window || !duplicates_by_source || !stats) return -1;
    if (duplicates_by_source->size() != (size_t)window->comm_size) {
        duplicates_by_source->assign(
            (size_t)window->comm_size, std::vector<uint64_t>());
    }
    const double group_t0 = GetTime();
    const size_t carry_count = window->carry.size();
    const size_t total = carry_count + owner_candidates.size();
    auto candidate_at = [&](size_t index)
        -> const MpiMarkdupCandidateShared & {
        return index < carry_count
            ? window->carry[index]
            : owner_candidates[index - carry_count];
    };
    auto qnames_at = [&](size_t index)
        -> const std::vector<unsigned char> & {
        return index < carry_count
            ? window->carry_qnames : owner_qnames;
    };
    auto pair_better = [&](size_t a_index, size_t b_index) -> bool {
        const MpiMarkdupCandidateShared &a = candidate_at(a_index);
        const MpiMarkdupCandidateShared &b = candidate_at(b_index);
        if (a.qc_fail != b.qc_fail) return a.qc_fail < b.qc_fail;
        if (a.score != b.score) return a.score > b.score;
        const std::vector<unsigned char> &a_qnames = qnames_at(a_index);
        const std::vector<unsigned char> &b_qnames = qnames_at(b_index);
        if (a.qname_offset > a_qnames.size() ||
            a.qname_len > a_qnames.size() - (size_t)a.qname_offset ||
            b.qname_offset > b_qnames.size() ||
            b.qname_len > b_qnames.size() - (size_t)b.qname_offset) {
            return a.global_order < b.global_order;
        }
        const size_t common = std::min(
            (size_t)a.qname_len, (size_t)b.qname_len);
        const int qcmp = common ? memcmp(
            a_qnames.data() + (size_t)a.qname_offset,
            b_qnames.data() + (size_t)b.qname_offset,
            common) : 0;
        if (qcmp != 0) return qcmp < 0;
        if (a.qname_len != b.qname_len) {
            return a.qname_len < b.qname_len;
        }
        return a.global_order < b.global_order;
    };

    size_t pair_expected = 0;
    size_t single_expected = 0;
    for (size_t i = 0; i < total; ++i) {
        const MpiMarkdupCandidateShared &candidate = candidate_at(i);
        if (!candidate.key.single) pair_expected++;
        else if (!candidate.paired_marker) single_expected++;
    }
    if (window->pair_best.Prepare(
            pair_expected, &stats->t_flat_alloc,
            &stats->t_flat_init) != 0 ||
        window->single_best.Prepare(
            single_expected, &stats->t_flat_alloc,
            &stats->t_flat_init) != 0) {
        return -1;
    }
    window->pair_slots.clear();
    window->single_slots.clear();
    window->marker_indices.clear();
    window->marker_carry_indices.clear();

    const double probe_t0 = GetTime();
    int status = 0;
    for (size_t i = 0; i < total; ++i) {
        const MpiMarkdupCandidateShared &candidate = candidate_at(i);
        const uint64_t hash = MdHashKeyFlat(candidate.key);
        size_t slot = 0;
        int inserted = 0;
        if (!candidate.key.single) {
            if (window->pair_best.FindOrInsert(
                    candidate.key, hash, i, &slot,
                    &inserted) != 0) {
                status = -1;
                break;
            }
            // Only winners close enough to cross a batch/rank boundary need
            // a second visit during carry construction.
            if (inserted && MdWindowKeepKey(
                    window, candidate.key,
                    progress_tid, progress_pos)) {
                window->pair_slots.push_back(slot);
            }
            if (inserted) continue;
            const size_t old_best = window->pair_best.best[slot];
            const MpiMarkdupCandidateShared &best =
                candidate_at(old_best);
            if (pair_better(i, old_best)) {
                if (MdPushDuplicateStream(best, window->comm_size,
                                          duplicates_by_source) != 0) {
                    status = -1;
                    break;
                }
                window->pair_best.best[slot] = i;
            } else if (MdPushDuplicateStream(
                           candidate, window->comm_size,
                           duplicates_by_source) != 0) {
                status = -1;
                break;
            }
            stats->pair_duplicates++;
            continue;
        }
        if (candidate.paired_marker) {
            window->marker_indices.push_back(i);
            if (MdWindowKeepKey(window, candidate.key,
                                progress_tid, progress_pos)) {
                window->marker_carry_indices.push_back(i);
            }
            continue;
        }
        if (window->single_best.FindOrInsert(
                candidate.key, hash, i, &slot,
                &inserted) != 0) {
            status = -1;
            break;
        }
        if (inserted && MdWindowKeepKey(
                window, candidate.key,
                progress_tid, progress_pos)) {
            window->single_slots.push_back(slot);
        }
        if (inserted) continue;
        const size_t old_best = window->single_best.best[slot];
        const MpiMarkdupCandidateShared &best = candidate_at(old_best);
        if (candidate.score > best.score ||
            (candidate.score == best.score &&
             candidate.global_order < best.global_order)) {
            if (MdPushDuplicateStream(best, window->comm_size,
                                      duplicates_by_source) != 0) {
                status = -1;
                break;
            }
            window->single_best.best[slot] = i;
        } else if (MdPushDuplicateStream(
                       candidate, window->comm_size,
                       duplicates_by_source) != 0) {
            status = -1;
            break;
        }
        stats->single_duplicates++;
    }
    const double probe_dt = GetTime() - probe_t0;
    stats->t_flat_probe += probe_dt;
    stats->t_group_sort += probe_dt;
    if (status != 0) return -1;

    const double finalize_t0 = GetTime();
    window->next_carry.clear();
    window->next_qnames.clear();
    auto keep = [&](size_t index) -> int {
        const MpiMarkdupCandidateShared &input = candidate_at(index);
        if (!MdWindowKeepKey(window, input.key,
                             progress_tid, progress_pos)) {
            return 0;
        }
        MpiMarkdupCandidateShared candidate = input;
        if (!candidate.key.single) {
            const std::vector<unsigned char> &qnames = qnames_at(index);
            if (candidate.qname_offset > qnames.size() ||
                candidate.qname_len > qnames.size() -
                    (size_t)candidate.qname_offset) {
                return -1;
            }
            const size_t qbase = window->next_qnames.size();
            if (candidate.qname_len > 0) {
                window->next_qnames.insert(
                    window->next_qnames.end(),
                    qnames.data() +
                        (size_t)candidate.qname_offset,
                    qnames.data() +
                        (size_t)candidate.qname_offset +
                        candidate.qname_len);
            }
            candidate.qname_offset = (uint64_t)qbase;
        }
        window->next_carry.push_back(candidate);
        return 0;
    };
    for (size_t m = 0; m < window->marker_indices.size(); ++m) {
        const size_t i = window->marker_indices[m];
        const MpiMarkdupCandidateShared &marker = candidate_at(i);
        size_t single_slot = 0;
        if (window->single_best.FindExisting(
                marker.key, MdHashKeyFlat(marker.key),
                &single_slot) &&
            window->single_best.best[single_slot] != SIZE_MAX) {
            const MpiMarkdupCandidateShared &candidate =
                candidate_at(window->single_best.best[single_slot]);
            if (MdPushDuplicateStream(
                    candidate, window->comm_size,
                    duplicates_by_source) != 0) {
                return -1;
            }
            window->single_best.best[single_slot] = SIZE_MAX;
            stats->single_duplicates++;
        }
    }
    for (size_t m = 0; m < window->marker_carry_indices.size(); ++m) {
        if (keep(window->marker_carry_indices[m]) != 0) return -1;
    }
    for (size_t i = 0; i < window->pair_slots.size(); ++i) {
        const size_t slot = window->pair_slots[i];
        if (keep(window->pair_best.best[slot]) != 0) return -1;
    }
    for (size_t i = 0; i < window->single_slots.size(); ++i) {
        const size_t slot = window->single_slots[i];
        if (window->single_best.best[slot] != SIZE_MAX &&
            keep(window->single_best.best[slot]) != 0) return -1;
    }
    size_t working_bytes =
        window->carry.capacity() * sizeof(MpiMarkdupCandidateShared) +
        window->carry_qnames.capacity() +
        window->next_carry.capacity() *
            sizeof(MpiMarkdupCandidateShared) +
        window->next_qnames.capacity() +
        window->pair_slots.capacity() * sizeof(size_t) +
        window->single_slots.capacity() * sizeof(size_t) +
        window->marker_indices.capacity() * sizeof(size_t) +
        window->marker_carry_indices.capacity() * sizeof(size_t);
    const size_t pair_bytes = MdFlatTableBytes(window->pair_best);
    const size_t single_bytes = MdFlatTableBytes(window->single_best);
    if (pair_bytes == SIZE_MAX || single_bytes == SIZE_MAX ||
        pair_bytes > SIZE_MAX - working_bytes ||
        single_bytes > SIZE_MAX - working_bytes - pair_bytes) {
        window->peak_working_bytes = SIZE_MAX;
    } else {
        working_bytes += pair_bytes + single_bytes;
        window->peak_working_bytes = std::max(
            window->peak_working_bytes, working_bytes);
    }
    window->carry.swap(window->next_carry);
    window->carry_qnames.swap(window->next_qnames);
    const double finalize_dt = GetTime() - finalize_t0;
    stats->t_flat_finalize += finalize_dt;
    stats->t_group_scan += finalize_dt;
    stats->t_group += GetTime() - group_t0;
    return 0;
}

size_t MpiMarkdupStreamingWindowMemory(
        const MpiMarkdupStreamingWindow *window) {
    if (!window) return 0;
    if (window->peak_working_bytes == SIZE_MAX) return SIZE_MAX;
    const size_t carry_bytes = window->carry.capacity() *
        sizeof(MpiMarkdupCandidateShared);
    if (window->carry_qnames.capacity() > SIZE_MAX - carry_bytes) {
        return SIZE_MAX;
    }
    return std::max(window->peak_working_bytes,
                    carry_bytes + window->carry_qnames.capacity());
}

int MpiMarkdupFindDuplicatesStreamingHash(
        std::vector<MpiMarkdupCandidateShared> *owner_candidates,
        const std::vector<unsigned char> &owner_qnames,
        int comm_size,
        std::vector<std::vector<uint64_t> > *duplicates_by_source,
        MpiMarkdupStats *stats) {
    const double group_t0 = GetTime();

    double alloc_t0 = GetTime();
    duplicates_by_source->assign(
        (size_t)comm_size, std::vector<uint64_t>());
    stats->t_flat_alloc += GetTime() - alloc_t0;

    double init_t0 = GetTime();
    const size_t total = owner_candidates->size();
    size_t pair_expected = 0;
    size_t single_expected = 0;
    for (size_t i = 0; i < total; ++i) {
        const MpiMarkdupCandidateShared &candidate =
            (*owner_candidates)[i];
        if (!candidate.key.single) {
            pair_expected++;
        } else if (!candidate.paired_marker) {
            single_expected++;
        }
    }
    stats->t_flat_init += GetTime() - init_t0;

    MdFlatBestTable pair_best;
    MdFlatBestTable single_best;
    int status = 0;
    double probe_t0 = 0.0;
    double finalize_t0 = 0.0;
    if (pair_best.Init(pair_expected, &stats->t_flat_alloc,
                       &stats->t_flat_init) != 0 ||
        single_best.Init(single_expected, &stats->t_flat_alloc,
                         &stats->t_flat_init) != 0) {
        status = -1;
        goto cleanup;
    }

    probe_t0 = GetTime();
    for (size_t i = 0; i < total; ++i) {
        const MpiMarkdupCandidateShared &candidate =
            (*owner_candidates)[i];
        if (!candidate.key.single) {
            const uint64_t hash = MdHashKeyFlat(candidate.key);
            size_t slot = 0;
            int inserted = 0;
            if (pair_best.FindOrInsert(candidate.key, hash, i, &slot,
                                      &inserted) != 0) {
                status = -1;
                goto cleanup_after_probe;
            }
            if (inserted) continue;

            size_t old_best = pair_best.best[slot];
            const MpiMarkdupCandidateShared &best =
                (*owner_candidates)[old_best];
            if (MdPairBetterStream(candidate, best, owner_qnames)) {
                if (MdPushDuplicateStream(best, comm_size,
                                          duplicates_by_source) != 0) {
                    status = -1;
                    goto cleanup_after_probe;
                }
                pair_best.best[slot] = i;
            } else {
                if (MdPushDuplicateStream(candidate, comm_size,
                                          duplicates_by_source) != 0) {
                    status = -1;
                    goto cleanup_after_probe;
                }
            }
            stats->pair_duplicates++;
            continue;
        }

        if (candidate.paired_marker) {
            continue;
        }

        const uint64_t hash = MdHashKeyFlat(candidate.key);
        size_t slot = 0;
        int inserted = 0;
        if (single_best.FindOrInsert(candidate.key, hash, i, &slot,
                                    &inserted) != 0) {
            status = -1;
            goto cleanup_after_probe;
        }
        if (inserted) continue;

        size_t old_best = single_best.best[slot];
        const MpiMarkdupCandidateShared &best =
            (*owner_candidates)[old_best];
        if (candidate.score > best.score ||
            (candidate.score == best.score &&
             candidate.global_order < best.global_order)) {
            if (MdPushDuplicateStream(best, comm_size,
                                      duplicates_by_source) != 0) {
                status = -1;
                goto cleanup_after_probe;
            }
            single_best.best[slot] = i;
        } else {
            if (MdPushDuplicateStream(candidate, comm_size,
                                      duplicates_by_source) != 0) {
                status = -1;
                goto cleanup_after_probe;
            }
        }
        stats->single_duplicates++;
    }
cleanup_after_probe:
    {
        const double probe_dt = GetTime() - probe_t0;
        stats->t_flat_probe += probe_dt;
        stats->t_group_sort += probe_dt;
    }
    if (status != 0) goto cleanup;

    finalize_t0 = GetTime();
    for (size_t i = 0; i < total; ++i) {
        const MpiMarkdupCandidateShared &marker =
            (*owner_candidates)[i];
        if (!marker.key.single || !marker.paired_marker) continue;
        const uint64_t hash = MdHashKeyFlat(marker.key);
        size_t slot = 0;
        if (!single_best.FindExisting(marker.key, hash, &slot)) continue;
        if (single_best.best[slot] == SIZE_MAX) continue;
        const MpiMarkdupCandidateShared &candidate =
            (*owner_candidates)[single_best.best[slot]];
        if (MdPushDuplicateStream(candidate, comm_size,
                                  duplicates_by_source) != 0) {
            status = -1;
            break;
        }
        single_best.best[slot] = SIZE_MAX;
        stats->single_duplicates++;
    }
    {
        const double finalize_dt = GetTime() - finalize_t0;
        stats->t_flat_finalize += finalize_dt;
        stats->t_group_scan += finalize_dt;
    }

cleanup:
    {
        double free_t0 = GetTime();
        pair_best.Free();
        single_best.Free();
        stats->t_flat_free += GetTime() - free_t0;
    }
    stats->t_group += GetTime() - group_t0;
    return status;
}
