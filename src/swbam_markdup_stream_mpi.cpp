#include "swbam_mpi.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
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
