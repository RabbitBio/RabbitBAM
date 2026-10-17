#include "swbam_mpi.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <string>
#include <utility>
#include <vector>

namespace {

const size_t kMdMarkerFilterBits = 1u << 16;
const size_t kMdMarkerFilterWords = kMdMarkerFilterBits / 64;

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

static size_t MdMarkerFilterIndex(const MpiMarkdupKeyShared &key) {
    uint64_t value = (uint64_t)key.this_coord;
    value ^= (uint64_t)(uint32_t)key.this_ref *
             0x9e3779b97f4a7c15ull;
    value ^= (uint64_t)(uint8_t)key.orientation *
             0xbf58476d1ce4e5b9ull;
    value ^= value >> 29;
    value *= 0x94d049bb133111ebull;
    value ^= value >> 32;
    return (size_t)value & (kMdMarkerFilterBits - 1);
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

    int ReservePersistent(size_t expected) {
        const size_t wanted = MdNextTableCapacity(expected);
        if (!wanted) return expected ? -1 : 0;
        if (wanted <= capacity) return 0;
        MdFlatBestTable next;
        if (next.Init(expected, nullptr, nullptr) != 0) return -1;
        for (size_t i = 0; i < capacity; ++i) {
            if (!used[i]) continue;
            size_t slot = 0;
            int inserted = 0;
            if (next.FindOrInsert(keys[i], hashes[i], best[i], &slot, &inserted) != 0)
                return -1;
        }
        std::swap(keys, next.keys);
        std::swap(hashes, next.hashes);
        std::swap(best, next.best);
        std::swap(used, next.used);
        std::swap(capacity, next.capacity);
        std::swap(mask, next.mask);
        std::swap(filled, next.filled);
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
        std::vector<std::vector<uint64_t> > *duplicates_by_source,
        std::vector<std::vector<MpiMarkdupDecisionShared> >
            *decisions_by_source) {
    const int source = candidate.source_rank;
    if (source < 0 || source >= comm_size) return -1;
    if (duplicates_by_source) {
        (*duplicates_by_source)[(size_t)source].push_back(
            candidate.ordinal);
    }
    if (decisions_by_source) {
        MpiMarkdupDecisionShared decision = {};
        decision.ordinal = candidate.ordinal;
        decision.duplicate = 1;
        (*decisions_by_source)[(size_t)source].push_back(decision);
    }
    return 0;
}

static int MdPushDuplicateStream(
        const MpiMarkdupCandidateShared &candidate, int comm_size,
        std::vector<std::vector<uint64_t> > *duplicates_by_source) {
    return MdPushDuplicateStream(
        candidate, comm_size, duplicates_by_source, nullptr);
}

static bool MdIsDecisionCandidate(
        const MpiMarkdupCandidateShared &candidate) {
    return !candidate.key.single || !candidate.paired_marker;
}

static int MdPushKeepDecision(
        const MpiMarkdupCandidateShared &candidate, int comm_size,
        std::vector<std::vector<MpiMarkdupDecisionShared> >
            *decisions_by_source) {
    if (!decisions_by_source || !MdIsDecisionCandidate(candidate)) {
        return 0;
    }
    const int source = candidate.source_rank;
    if (source < 0 || source >= comm_size) return -1;
    MpiMarkdupDecisionShared decision = {};
    decision.ordinal = candidate.ordinal;
    decision.duplicate = 0;
    (*decisions_by_source)[(size_t)source].push_back(decision);
    return 0;
}

static int MdEmitOptimisticDecision(
        const MpiMarkdupCandidateShared &candidate, int comm_size,
        int local_rank, int action,
        MpiMarkdupLocalDecisionFn apply_local, void *local_opaque,
        std::vector<std::vector<MpiMarkdupDecisionShared> >
            *remote_decisions_by_source) {
    if (!MdIsDecisionCandidate(candidate)) return 0;
    const int source = candidate.source_rank;
    if (source < 0 || source >= comm_size || !apply_local ||
        !remote_decisions_by_source) {
        return -1;
    }
    if (source == local_rank) {
        return apply_local(local_opaque, candidate.ordinal, action);
    }
    if (action == RB_MARKDUP_DECISION_PENDING) {
        // A remotely owned candidate is made pending by its source rank
        // before candidate exchange, so no reverse notification is needed.
        return 0;
    }
    MpiMarkdupDecisionShared decision = {};
    decision.ordinal = candidate.ordinal;
    decision.duplicate =
        action == RB_MARKDUP_DECISION_DUPLICATE ? 1 : 0;
    (*remote_decisions_by_source)[(size_t)source].push_back(decision);
    return 0;
}

struct MdBoundaryWinner {
    MpiMarkdupCandidateShared candidate = {};
    std::string qname;
    bool valid_qname = false;
    bool present = false;
    bool paired_marker = false;
};

// Only fallback boundaries use these persistent tables. Ordinary windows keep
// their existing batch-local tables and scheduling.
struct MdIncrementalBoundary {
    MdFlatBestTable pairs;
    MdFlatBestTable singles;
    std::vector<MdBoundaryWinner> pair_winners;
    std::vector<MdBoundaryWinner> single_winners;
    size_t qname_bytes = 0;
    uint64_t processed = 0;

    static bool ValidQname(const MpiMarkdupCandidateShared &candidate,
                           const std::vector<unsigned char> &qnames) {
        return candidate.qname_offset <= qnames.size() &&
            candidate.qname_len <= qnames.size() - (size_t)candidate.qname_offset;
    }

    bool PairBetter(const MpiMarkdupCandidateShared &candidate,
                    const std::vector<unsigned char> &qnames,
                    const MdBoundaryWinner &winner) const {
        const auto &old = winner.candidate;
        if (candidate.qc_fail != old.qc_fail) return candidate.qc_fail < old.qc_fail;
        if (candidate.score != old.score) return candidate.score > old.score;
        if (winner.valid_qname && ValidQname(candidate, qnames)) {
            const size_t common = std::min<size_t>(candidate.qname_len, winner.qname.size());
            const int cmp = common ? memcmp(qnames.data() + candidate.qname_offset,
                                            winner.qname.data(), common) : 0;
            if (cmp) return cmp < 0;
            if (candidate.qname_len != winner.qname.size())
                return candidate.qname_len < winner.qname.size();
        }
        return candidate.global_order < old.global_order;
    }

    void Store(MdBoundaryWinner *winner, const MpiMarkdupCandidateShared &candidate,
               const std::vector<unsigned char> &qnames) {
        qname_bytes -= winner->qname.capacity();
        winner->candidate = candidate;
        winner->present = true;
        winner->valid_qname = !candidate.key.single && ValidQname(candidate, qnames);
        if (winner->valid_qname && candidate.qname_len)
            winner->qname.assign((const char *)qnames.data() + candidate.qname_offset,
                                 candidate.qname_len);
        else
            winner->qname.clear();
        qname_bytes += winner->qname.capacity();
    }

    template <typename Duplicate, typename Pending>
    int Process(const MpiMarkdupCandidateShared &candidate,
                const std::vector<unsigned char> &qnames,
                Duplicate duplicate, Pending pending, MpiMarkdupStats *stats) {
        MdFlatBestTable &table = candidate.key.single ? singles : pairs;
        std::vector<MdBoundaryWinner> &winners =
            candidate.key.single ? single_winners : pair_winners;
        size_t slot = 0;
        const uint64_t hash = MdHashKeyFlat(candidate.key);
        if (!table.FindExisting(candidate.key, hash, &slot)) {
            if (table.ReservePersistent(table.filled + 1) != 0) return -1;
            int inserted = 0;
            if (table.FindOrInsert(candidate.key, hash, winners.size(), &slot, &inserted) != 0)
                return -1;
            winners.emplace_back();
            qname_bytes += winners.back().qname.capacity();
        }
        MdBoundaryWinner &winner = winners[table.best[slot]];
        ++processed;
        if (candidate.key.single && candidate.paired_marker) {
            if (winner.present) {
                if (duplicate(winner.candidate) != 0) return -1;
                ++stats->single_duplicates;
                winner.present = false;
            }
            winner.paired_marker = true;
            return 0;
        }
        if (winner.paired_marker) {
            if (duplicate(candidate) != 0) return -1;
            ++stats->single_duplicates;
            return 0;
        }
        if (!winner.present) {
            Store(&winner, candidate, qnames);
            return pending(candidate);
        }
        const bool better = candidate.key.single
            ? candidate.score > winner.candidate.score ||
              (candidate.score == winner.candidate.score &&
               candidate.global_order < winner.candidate.global_order)
            : PairBetter(candidate, qnames, winner);
        if (duplicate(better ? winner.candidate : candidate) != 0) return -1;
        if (candidate.key.single) ++stats->single_duplicates;
        else ++stats->pair_duplicates;
        if (better) {
            Store(&winner, candidate, qnames);
            if (pending(candidate) != 0) return -1;
        }
        return 0;
    }

    size_t Memory() const {
        const size_t table_stride = sizeof(MpiMarkdupKeyShared) + sizeof(uint64_t) + sizeof(size_t) + 1;
        return (pairs.capacity + singles.capacity) * table_stride +
            (pair_winners.capacity() + single_winners.capacity()) * sizeof(MdBoundaryWinner) +
            qname_bytes;
    }

    template <typename Keep>
    int Release(Keep keep) {
        for (const auto &winner : pair_winners)
            if (winner.present && keep(winner.candidate) != 0) return -1;
        for (const auto &winner : single_winners)
            if (winner.present && keep(winner.candidate) != 0) return -1;
        pairs.Free();
        singles.Free();
        std::vector<MdBoundaryWinner>().swap(pair_winners);
        std::vector<MdBoundaryWinner>().swap(single_winners);
        qname_bytes = 0;
        return 0;
    }
};

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
    std::vector<MpiMarkdupCandidateShared> boundary_carry;
    std::vector<unsigned char> boundary_qnames;
    std::vector<MpiMarkdupCandidateShared> next_boundary;
    std::vector<unsigned char> next_boundary_qnames;
    MdFlatBestTable pair_best;
    MdFlatBestTable single_best;
    MdFlatBestTable boundary_markers;
    MdIncrementalBoundary first_boundary;
    MdIncrementalBoundary last_boundary;
    std::vector<MpiMarkdupCandidateShared> interior_candidates;
    bool first_closed = false;
    bool last_closed = false;
    size_t peak_boundary_working_bytes = 0;
    uint64_t interior_marker_visits = 0;
    uint64_t interior_markers_retained = 0;
    size_t interior_marker_peak = 0;
    std::vector<size_t> pair_slots;
    std::vector<size_t> single_slots;
    std::vector<size_t> marker_indices;
    std::vector<size_t> marker_carry_indices;
    uint64_t marker_filter[kMdMarkerFilterWords];

    MpiMarkdupStreamingWindow()
        : comm_size(1), max_read_length(300),
          first_ref(0), first_coord(0),
          last_ref(0), last_coord(0),
          peak_working_bytes(0) {
        memset(marker_filter, 0, sizeof(marker_filter));
    }
};

namespace {

static bool MdWindowBoundaryKey(
        const MpiMarkdupStreamingWindow *window,
        const MpiMarkdupKeyShared &key) {
    if (!window->first_closed && window->first_ref > 0 && key.this_ref == window->first_ref &&
        key.this_coord >=
            window->first_coord - window->max_read_length &&
        key.this_coord <=
            window->first_coord + window->max_read_length) {
        return true;
    }
    return !window->last_closed && window->last_ref > 0 && key.this_ref == window->last_ref &&
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

static int MdStreamingWindowProcessImpl(
        MpiMarkdupStreamingWindow *window,
        const std::vector<MpiMarkdupCandidateShared> &owner_candidates,
        const std::vector<unsigned char> &owner_qnames,
        int progress_tid, int progress_pos,
        std::vector<std::vector<uint64_t> > *duplicates_by_source,
        std::vector<std::vector<MpiMarkdupDecisionShared> >
            *decisions_by_source,
        int optimistic_local_rank,
        MpiMarkdupLocalDecisionFn apply_local,
        void *local_opaque,
        bool optimistic,
        MpiMarkdupStats *stats) {
    if (!window || (!duplicates_by_source && !decisions_by_source) ||
        !stats) {
        return -1;
    }
    if (duplicates_by_source &&
        duplicates_by_source->size() != (size_t)window->comm_size) {
        duplicates_by_source->assign(
            (size_t)window->comm_size, std::vector<uint64_t>());
    }
    if (decisions_by_source &&
        decisions_by_source->size() != (size_t)window->comm_size) {
        decisions_by_source->assign(
            (size_t)window->comm_size,
            std::vector<MpiMarkdupDecisionShared>());
    }
    const bool emit_decisions = decisions_by_source != nullptr;
    const bool track_all_winners = emit_decisions && !optimistic;
    auto emit_duplicate = [&](const MpiMarkdupCandidateShared &candidate)
        -> int {
        return optimistic
            ? MdEmitOptimisticDecision(
                  candidate, window->comm_size, optimistic_local_rank,
                  RB_MARKDUP_DECISION_DUPLICATE, apply_local,
                  local_opaque, decisions_by_source)
            : MdPushDuplicateStream(
                  candidate, window->comm_size, duplicates_by_source,
                  decisions_by_source);
    };
    auto emit_keep = [&](const MpiMarkdupCandidateShared &candidate)
        -> int {
        return optimistic
            ? MdEmitOptimisticDecision(
                  candidate, window->comm_size, optimistic_local_rank,
                  RB_MARKDUP_DECISION_KEEP, apply_local,
                  local_opaque, decisions_by_source)
            : MdPushKeepDecision(
                  candidate, window->comm_size, decisions_by_source);
    };
    auto mark_pending = [&](const MpiMarkdupCandidateShared &candidate)
        -> int {
        return optimistic
            ? MdEmitOptimisticDecision(
                  candidate, window->comm_size, optimistic_local_rank,
                  RB_MARKDUP_DECISION_PENDING, apply_local,
                  local_opaque, decisions_by_source)
            : 0;
    };
    const double group_t0 = GetTime();
    const bool has_boundary = window->first_ref > 0 || window->last_ref > 0;
    const bool update_boundary = has_boundary && std::any_of(
        owner_candidates.begin(), owner_candidates.end(),
        [&](const MpiMarkdupCandidateShared &candidate) {
            return MdWindowBoundaryKey(window, candidate.key);
        });
    // Boundary keys cannot equal interior keys. Leave the deferred boundary
    // untouched until another boundary candidate arrives, or until EOF.
    const size_t ordinary_count = window->carry.size();
    const size_t carry_count = ordinary_count +
        (update_boundary ? window->boundary_carry.size() : 0);
    const size_t total = carry_count + owner_candidates.size();
    auto candidate_at = [&](size_t index)
        -> const MpiMarkdupCandidateShared & {
        if (index < ordinary_count) return window->carry[index];
        return index < carry_count
            ? window->boundary_carry[index - ordinary_count]
            : owner_candidates[index - carry_count];
    };
    auto qnames_at = [&](size_t index)
        -> const std::vector<unsigned char> & {
        if (index < ordinary_count) return window->carry_qnames;
        return index < carry_count ? window->boundary_qnames : owner_qnames;
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
    memset(window->marker_filter, 0, sizeof(window->marker_filter));
    for (size_t i = 0; i < total; ++i) {
        const MpiMarkdupCandidateShared &candidate = candidate_at(i);
        if (!candidate.key.single) pair_expected++;
        else if (!candidate.paired_marker) {
            single_expected++;
            const size_t bit = MdMarkerFilterIndex(candidate.key);
            window->marker_filter[bit >> 6] |=
                (uint64_t)1 << (bit & 63);
        }
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
            if (inserted && (i < carry_count || track_all_winners ||
                    (optimistic &&
                     candidate.source_rank != optimistic_local_rank) ||
                    MdWindowKeepKey(window, candidate.key,
                                    progress_tid, progress_pos))) {
                window->pair_slots.push_back(slot);
            }
            if (inserted) continue;
            const size_t old_best = window->pair_best.best[slot];
            const MpiMarkdupCandidateShared &best =
                candidate_at(old_best);
            if (pair_better(i, old_best)) {
                if (emit_duplicate(best) != 0) {
                    status = -1;
                    break;
                }
                window->pair_best.best[slot] = i;
                if (optimistic &&
                    candidate.source_rank != optimistic_local_rank &&
                    !MdWindowKeepKey(window, candidate.key,
                                     progress_tid, progress_pos)) {
                    window->pair_slots.push_back(slot);
                }
            } else if (emit_duplicate(candidate) != 0) {
                status = -1;
                break;
            }
            stats->pair_duplicates++;
            continue;
        }
        if (candidate.paired_marker) {
            const size_t bit = MdMarkerFilterIndex(candidate.key);
            if (window->marker_filter[bit >> 6] &
                    ((uint64_t)1 << (bit & 63))) {
                window->marker_indices.push_back(i);
            }
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
        if (inserted && (i < carry_count || track_all_winners ||
                (optimistic &&
                 candidate.source_rank != optimistic_local_rank) ||
                MdWindowKeepKey(window, candidate.key,
                                progress_tid, progress_pos))) {
            window->single_slots.push_back(slot);
        }
        if (inserted) continue;
        const size_t old_best = window->single_best.best[slot];
        const MpiMarkdupCandidateShared &best = candidate_at(old_best);
        if (candidate.score > best.score ||
            (candidate.score == best.score &&
             candidate.global_order < best.global_order)) {
            if (emit_duplicate(best) != 0) {
                status = -1;
                break;
            }
            window->single_best.best[slot] = i;
            if (optimistic &&
                candidate.source_rank != optimistic_local_rank &&
                !MdWindowKeepKey(window, candidate.key,
                                 progress_tid, progress_pos)) {
                window->single_slots.push_back(slot);
            }
        } else if (emit_duplicate(candidate) != 0) {
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
    if (update_boundary) {
        window->next_boundary.clear();
        window->next_boundary_qnames.clear();
    }
    auto keep = [&](size_t index) -> int {
        const MpiMarkdupCandidateShared &input = candidate_at(index);
        if (!MdWindowKeepKey(window, input.key,
                             progress_tid, progress_pos)) {
            return emit_keep(input);
        }
        // A retained fallback candidate was already marked pending when added.
        if (!(has_boundary && index < carry_count) && mark_pending(input) != 0)
            return -1;
        MpiMarkdupCandidateShared candidate = input;
        const bool boundary = has_boundary &&
            MdWindowBoundaryKey(window, input.key);
        std::vector<unsigned char> &next_qnames = boundary
            ? window->next_boundary_qnames : window->next_qnames;
        if (!candidate.key.single) {
            const std::vector<unsigned char> &qnames = qnames_at(index);
            if (candidate.qname_offset > qnames.size() ||
                candidate.qname_len > qnames.size() -
                    (size_t)candidate.qname_offset) {
                return -1;
            }
            const size_t qbase = next_qnames.size();
            if (candidate.qname_len > 0) {
                next_qnames.insert(
                    next_qnames.end(),
                    qnames.data() +
                        (size_t)candidate.qname_offset,
                    qnames.data() +
                        (size_t)candidate.qname_offset +
                        candidate.qname_len);
            }
            candidate.qname_offset = (uint64_t)qbase;
        }
        (boundary ? window->next_boundary : window->next_carry).push_back(candidate);
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
            if (emit_duplicate(candidate) != 0) {
                return -1;
            }
            window->single_best.best[single_slot] = SIZE_MAX;
            stats->single_duplicates++;
        }
    }
    // Fallback interior windows also need only one existence marker per key.
    // Dense coordinates may span many batches even away from a rank cut.
    const bool deduplicate_interior = optimistic && has_boundary;
    const bool deduplicate_markers = (update_boundary || deduplicate_interior) &&
        window->marker_carry_indices.size() > 1;
    if (deduplicate_markers && window->boundary_markers.Prepare(
            window->marker_carry_indices.size(), &stats->t_flat_alloc,
            &stats->t_flat_init) != 0) return -1;
    size_t retained_markers = 0;
    for (size_t m = 0; m < window->marker_carry_indices.size(); ++m) {
        const size_t index = window->marker_carry_indices[m];
        const auto &marker = candidate_at(index);
        if (deduplicate_markers &&
            (deduplicate_interior || MdWindowBoundaryKey(window, marker.key))) {
            size_t slot = 0;
            int inserted = 0;
            if (window->boundary_markers.FindOrInsert(marker.key,
                    MdHashKeyFlat(marker.key), index, &slot, &inserted) != 0) return -1;
            // A paired marker represents existence, not a record decision.
            if (!inserted) continue;
        }
        if (keep(index) != 0) return -1;
        ++retained_markers;
    }
    if (deduplicate_interior) {
        window->interior_marker_visits += window->marker_carry_indices.size();
        window->interior_markers_retained += retained_markers;
        window->interior_marker_peak = std::max(window->interior_marker_peak, retained_markers);
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
        (window->boundary_carry.capacity() + window->next_boundary.capacity()) *
            sizeof(MpiMarkdupCandidateShared) +
        window->boundary_qnames.capacity() + window->next_boundary_qnames.capacity() +
        window->pair_slots.capacity() * sizeof(size_t) +
        window->single_slots.capacity() * sizeof(size_t) +
        window->marker_indices.capacity() * sizeof(size_t) +
        window->marker_carry_indices.capacity() * sizeof(size_t) +
        sizeof(window->marker_filter);
    const size_t pair_bytes = MdFlatTableBytes(window->pair_best);
    const size_t single_bytes = MdFlatTableBytes(window->single_best);
    const size_t marker_bytes = MdFlatTableBytes(window->boundary_markers);
    if (pair_bytes == SIZE_MAX || single_bytes == SIZE_MAX || marker_bytes == SIZE_MAX ||
        pair_bytes > SIZE_MAX - working_bytes ||
        single_bytes > SIZE_MAX - working_bytes - pair_bytes ||
        marker_bytes > SIZE_MAX - working_bytes - pair_bytes - single_bytes) {
        window->peak_working_bytes = SIZE_MAX;
    } else {
        working_bytes += pair_bytes + single_bytes + marker_bytes;
        window->peak_working_bytes = std::max(
            window->peak_working_bytes, working_bytes);
    }
    window->carry.swap(window->next_carry);
    window->carry_qnames.swap(window->next_qnames);
    if (update_boundary) {
        window->boundary_carry.swap(window->next_boundary);
        window->boundary_qnames.swap(window->next_boundary_qnames);
    }
    static const bool pending_diagnostics = []() {
        const char *value = getenv("SWBAM_MARKDUP_PENDING_DIAG");
        return value && strcmp(value, "1") == 0;
    }();
    if (pending_diagnostics && optimistic &&
        (stats->group_count == 1 || stats->group_count % 256 == 0)) {
        const MpiMarkdupCandidateShared *first = nullptr;
        size_t local_pending = 0, remote_pending = 0;
        size_t paired_markers = 0, boundary_candidates = 0;
        const auto observe = [&](const std::vector<MpiMarkdupCandidateShared> &items) {
          for (const auto &candidate : items) {
            if (MdWindowBoundaryKey(window, candidate.key)) {
                ++boundary_candidates;
            }
            if (candidate.key.single && candidate.paired_marker) {
                ++paired_markers;
            }
            if (!MdIsDecisionCandidate(candidate)) continue;
            if (candidate.source_rank != optimistic_local_rank) {
                ++remote_pending;
                continue;
            }
            ++local_pending;
            if (!first || candidate.ordinal < first->ordinal) first = &candidate;
          }
        };
        observe(window->carry);
        observe(window->boundary_carry);
        fprintf(stderr, "[md-window rank=%d] batch=%lld progress=%d:%d "
                "incoming=%zu carry_in=%zu carry=%zu markers=%zu "
                "boundary_candidates=%zu boundary_update=%d local_pending=%zu remote_pending=%zu "
                "boundary_first=%d:%lld boundary_last=%d:%lld\n",
                optimistic_local_rank, stats->group_count, progress_tid, progress_pos,
                owner_candidates.size(), carry_count,
                window->carry.size() + window->boundary_carry.size(),
                paired_markers, boundary_candidates, (int)update_boundary,
                local_pending, remote_pending,
                window->first_ref, (long long)window->first_coord,
                window->last_ref, (long long)window->last_coord);
        if (first) {
            const char *reason = MdWindowBoundaryKey(window, first->key) ? "boundary" :
                (first->key.this_ref > progress_tid + 1 ? "future-ref" : "window");
            fprintf(stderr, "[md-window-first rank=%d] batch=%lld ordinal=%llu "
                    "global_block=%llu record=%u key=%d:%lld other=%d:%lld "
                    "single=%u leftmost=%u orientation=%u reason=%s\n",
                    optimistic_local_rank, stats->group_count,
                    (unsigned long long)first->ordinal,
                    (unsigned long long)(first->global_order >> 32),
                    (unsigned)(first->global_order & 0xffffffffu),
                    first->key.this_ref, (long long)first->key.this_coord,
                    first->key.other_ref, (long long)first->key.other_coord,
                    (unsigned)first->key.single, (unsigned)first->key.leftmost,
                    (unsigned)first->key.orientation, reason);
        }
        fflush(stderr);
    }
    const double finalize_dt = GetTime() - finalize_t0;
    stats->t_flat_finalize += finalize_dt;
    stats->t_group_scan += finalize_dt;
    stats->t_group += GetTime() - group_t0;
    return 0;
}

int MpiMarkdupStreamingWindowProcess(
        MpiMarkdupStreamingWindow *window,
        const std::vector<MpiMarkdupCandidateShared> &owner_candidates,
        const std::vector<unsigned char> &owner_qnames,
        int progress_tid, int progress_pos,
        std::vector<std::vector<uint64_t> > *duplicates_by_source,
        MpiMarkdupStats *stats) {
    return MdStreamingWindowProcessImpl(
        window, owner_candidates, owner_qnames,
        progress_tid, progress_pos, duplicates_by_source,
        nullptr, -1, nullptr, nullptr, false, stats);
}

int MpiMarkdupStreamingWindowProcessDecisions(
        MpiMarkdupStreamingWindow *window,
        const std::vector<MpiMarkdupCandidateShared> &owner_candidates,
        const std::vector<unsigned char> &owner_qnames,
        int progress_tid, int progress_pos,
        std::vector<std::vector<MpiMarkdupDecisionShared> >
            *decisions_by_source,
        MpiMarkdupStats *stats) {
    return MdStreamingWindowProcessImpl(
        window, owner_candidates, owner_qnames,
        progress_tid, progress_pos, nullptr,
        decisions_by_source, -1, nullptr, nullptr, false, stats);
}

int MpiMarkdupStreamingWindowProcessOptimistic(
        MpiMarkdupStreamingWindow *window,
        const std::vector<MpiMarkdupCandidateShared> &owner_candidates,
        const std::vector<unsigned char> &owner_qnames,
        int progress_tid, int progress_pos, int local_rank,
        MpiMarkdupLocalDecisionFn apply_local, void *local_opaque,
        std::vector<std::vector<MpiMarkdupDecisionShared> >
            *remote_decisions_by_source,
        MpiMarkdupStats *stats) {
    if (window && (window->first_ref > 0 || window->last_ref > 0)) {
        if (!stats || !apply_local || !remote_decisions_by_source ||
            remote_decisions_by_source->size() != (size_t)window->comm_size)
            return -1;
        const bool boundary_batch = std::any_of(owner_candidates.begin(), owner_candidates.end(),
            [&](const MpiMarkdupCandidateShared &candidate) {
                return MdWindowBoundaryKey(window, candidate.key);
            });
        if (boundary_batch) {
            const double t0 = GetTime();
            window->interior_candidates.clear();
            const auto emit = [&](const MpiMarkdupCandidateShared &candidate, int action) {
                return MdEmitOptimisticDecision(candidate, window->comm_size, local_rank,
                    action, apply_local, local_opaque, remote_decisions_by_source);
            };
            const auto duplicate = [&](const MpiMarkdupCandidateShared &candidate) {
                return emit(candidate, RB_MARKDUP_DECISION_DUPLICATE);
            };
            const auto pending = [&](const MpiMarkdupCandidateShared &candidate) {
                return emit(candidate, RB_MARKDUP_DECISION_PENDING);
            };
            for (const auto &candidate : owner_candidates) {
                if (!MdWindowBoundaryKey(window, candidate.key)) {
                    window->interior_candidates.push_back(candidate);
                    continue;
                }
                const bool first = !window->first_closed &&
                    candidate.key.this_ref == window->first_ref &&
                    candidate.key.this_coord >= window->first_coord - window->max_read_length &&
                    candidate.key.this_coord <= window->first_coord + window->max_read_length;
                MdIncrementalBoundary &boundary = first
                    ? window->first_boundary : window->last_boundary;
                if (boundary.Process(candidate, owner_qnames, duplicate, pending, stats) != 0)
                    return -1;
            }
            const double elapsed = GetTime() - t0;
            stats->t_group += elapsed;
            stats->t_flat_probe += elapsed;
            stats->t_group_sort += elapsed;
        }
        const int ret = MdStreamingWindowProcessImpl(
            window, boundary_batch ? window->interior_candidates : owner_candidates,
            owner_qnames, progress_tid, progress_pos, nullptr,
            remote_decisions_by_source, local_rank, apply_local, local_opaque, true, stats);
        const size_t boundary_bytes = window->first_boundary.Memory() + window->last_boundary.Memory() +
            window->interior_candidates.capacity() * sizeof(MpiMarkdupCandidateShared);
        window->peak_boundary_working_bytes = std::max(window->peak_boundary_working_bytes, boundary_bytes);
        const size_t ordinary_bytes = window->peak_working_bytes;
        // Boundary storage is accounted separately by Memory(), without scanning keys per batch.
        if (ordinary_bytes > SIZE_MAX - boundary_bytes) window->peak_working_bytes = SIZE_MAX;
        return ret;
    }
    return MdStreamingWindowProcessImpl(
        window, owner_candidates, owner_qnames,
        progress_tid, progress_pos, nullptr,
        remote_decisions_by_source, local_rank, apply_local,
        local_opaque, true, stats);
}

int MpiMarkdupStreamingWindowReleaseBoundaryOptimistic(
        MpiMarkdupStreamingWindow *window, bool first, bool last, int local_rank,
        MpiMarkdupLocalDecisionFn apply_local, void *local_opaque,
        std::vector<std::vector<MpiMarkdupDecisionShared> > *remote_decisions_by_source) {
    if (!window || !apply_local || !remote_decisions_by_source ||
        remote_decisions_by_source->size() != (size_t)window->comm_size) return -1;
    const auto keep = [&](const MpiMarkdupCandidateShared &candidate) {
        return MdEmitOptimisticDecision(candidate, window->comm_size, local_rank,
            RB_MARKDUP_DECISION_KEEP, apply_local, local_opaque, remote_decisions_by_source);
    };
    if (first && !window->first_closed) {
        if (window->first_boundary.Release(keep) != 0) return -1;
        window->first_closed = true;
    }
    if (last && !window->last_closed) {
        if (window->last_boundary.Release(keep) != 0) return -1;
        window->last_closed = true;
    }
    return 0;
}

int MpiMarkdupStreamingWindowFinalizeDecisions(
        MpiMarkdupStreamingWindow *window,
        std::vector<std::vector<MpiMarkdupDecisionShared> >
            *decisions_by_source) {
    if (!window || !decisions_by_source) return -1;
    if (decisions_by_source->size() != (size_t)window->comm_size) {
        decisions_by_source->assign(
            (size_t)window->comm_size,
            std::vector<MpiMarkdupDecisionShared>());
    }
    for (size_t i = 0; i < window->carry.size(); ++i) {
        if (MdPushKeepDecision(window->carry[i], window->comm_size,
                               decisions_by_source) != 0) {
            return -1;
        }
    }
    for (const auto &candidate : window->boundary_carry) {
        if (MdPushKeepDecision(candidate, window->comm_size, decisions_by_source) != 0)
            return -1;
    }
    window->boundary_carry.clear();
    window->boundary_qnames.clear();
    window->next_boundary.clear();
    window->next_boundary_qnames.clear();
    window->carry.clear();
    window->carry_qnames.clear();
    window->next_carry.clear();
    window->next_qnames.clear();
    return 0;
}

int MpiMarkdupStreamingWindowFinalizeOptimistic(
        MpiMarkdupStreamingWindow *window, int local_rank,
        MpiMarkdupLocalDecisionFn apply_local, void *local_opaque,
        std::vector<std::vector<MpiMarkdupDecisionShared> >
            *remote_decisions_by_source) {
    if (!window || !apply_local || !remote_decisions_by_source ||
        remote_decisions_by_source->size() !=
            (size_t)window->comm_size) {
        return -1;
    }
    if (window->first_ref > 0 || window->last_ref > 0) {
        fprintf(stderr, "[md-window-reuse rank=%d] marker_visits=%llu "
                "marker_retained=%llu marker_peak=%zu\n", local_rank,
                (unsigned long long)window->interior_marker_visits,
                (unsigned long long)window->interior_markers_retained,
                window->interior_marker_peak);
        fprintf(stderr, "[md-boundary-state rank=%d] mode=incremental candidates=%llu "
                "first_closed=%d last_closed=%d\n", local_rank,
                (unsigned long long)(window->first_boundary.processed + window->last_boundary.processed),
                (int)window->first_closed, (int)window->last_closed);
        if (MpiMarkdupStreamingWindowReleaseBoundaryOptimistic(window, true, true,
                local_rank, apply_local, local_opaque, remote_decisions_by_source) != 0) return -1;
    }
    for (size_t i = 0; i < window->carry.size(); ++i) {
        if (MdEmitOptimisticDecision(
                window->carry[i], window->comm_size, local_rank,
                RB_MARKDUP_DECISION_KEEP, apply_local, local_opaque,
                remote_decisions_by_source) != 0) {
            return -1;
        }
    }
    for (const auto &candidate : window->boundary_carry) {
        if (MdEmitOptimisticDecision(candidate, window->comm_size, local_rank,
                RB_MARKDUP_DECISION_KEEP, apply_local, local_opaque,
                remote_decisions_by_source) != 0) return -1;
    }
    window->boundary_carry.clear();
    window->boundary_qnames.clear();
    window->next_boundary.clear();
    window->next_boundary_qnames.clear();
    window->carry.clear();
    window->carry_qnames.clear();
    window->next_carry.clear();
    window->next_qnames.clear();
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
    const size_t ordinary = std::max(window->peak_working_bytes,
                                    carry_bytes + window->carry_qnames.capacity());
    const size_t boundary = window->peak_boundary_working_bytes;
    return boundary > SIZE_MAX - ordinary ? SIZE_MAX : ordinary + boundary;
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
