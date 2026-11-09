#include "../swbam-app/algorithms/markdup/swbam_markdup_stream_mpi.cpp"
#include "../swbam-app/algorithms/markdup/markdup-boundary-progress.h"

#include <cassert>
#include <map>

static MpiMarkdupCandidateShared Candidate(uint64_t order, bool marker,
                                           bool single = true) {
    MpiMarkdupCandidateShared candidate = {};
    candidate.key.this_ref = 1;
    candidate.key.this_coord = 1000;
    candidate.key.other_ref = single ? 0 : 1;
    candidate.key.other_coord = single ? 0 : 1200;
    candidate.key.single = single;
    candidate.paired_marker = marker;
    candidate.source_rank = order % 2;
    candidate.ordinal = order;
    candidate.global_order = order;
    candidate.score = (int)(order % 3);
    return candidate;
}

static void TestMarkerCarry(bool observe) {
    MpiMarkdupStreamingWindow *window =
        MpiMarkdupStreamingWindowCreate(2, 300, 0, 999, -1, -1);
    assert(window);
    MpiMarkdupStats stats = {};
    std::vector<std::vector<uint64_t> > duplicates;
    std::vector<unsigned char> qnames;
    for (unsigned batch = 0; batch < 32; ++batch) {
        std::vector<MpiMarkdupCandidateShared> candidates;
        for (unsigned i = 0; i < 256; ++i) {
            candidates.push_back(Candidate(batch * 256 + i, true));
        }
        assert(MpiMarkdupStreamingWindowProcess(window, candidates, qnames,
                   0, 10000, &duplicates, &stats) == 0);
        if (!observe) {
            assert(window->carry.empty());
            assert(window->boundary_carry.size() == 1);
        }
    }
    printf("marker_carry=%zu input_markers=8192 group_s=%.6f\n",
           window->carry.size() + window->boundary_carry.size(), stats.t_group);
    if (!observe) {
        const auto *parked = window->boundary_carry.data();
        std::vector<MpiMarkdupCandidateShared> empty;
        assert(MpiMarkdupStreamingWindowProcess(window, empty, qnames,
                   0, 10000, &duplicates, &stats) == 0);
        assert(window->boundary_carry.data() == parked);
        assert(window->boundary_carry.size() == 1);
    }
    assert(duplicates[0].empty() && duplicates[1].empty());
    MpiMarkdupStreamingWindowDestroy(window);
}

static void TestDuplicateDecisions() {
    MpiMarkdupStreamingWindow *window =
        MpiMarkdupStreamingWindowCreate(2, 300, 0, 999, -1, -1);
    assert(window);
    std::vector<MpiMarkdupCandidateShared> all;
    std::vector<unsigned char> qnames;
    std::vector<std::vector<uint64_t> > actual, expected;
    MpiMarkdupStats stats = {}, reference = {};
    for (unsigned batch = 0; batch < 4; ++batch) {
        std::vector<MpiMarkdupCandidateShared> candidates;
        for (unsigned i = 0; i < 10; ++i) {
            const uint64_t order = batch * 10 + i;
            auto candidate = Candidate(order, i < 4, i < 7);
            if (!candidate.key.single) {
                candidate.qname_offset = qnames.size();
                candidate.qname_len = 1;
                qnames.push_back((unsigned char)('a' + order % 7));
            }
            candidates.push_back(candidate);
            all.push_back(candidate);
        }
        assert(MpiMarkdupStreamingWindowProcess(window, candidates, qnames,
                   0, 10000, &actual, &stats) == 0);
    }
    assert(MpiMarkdupFindDuplicatesStreamingHash(&all, qnames, 2,
               &expected, &reference) == 0);
    for (size_t r = 0; r < 2; ++r) {
        std::sort(actual[r].begin(), actual[r].end());
        std::sort(expected[r].begin(), expected[r].end());
        assert(actual[r] == expected[r]);
        assert(std::adjacent_find(actual[r].begin(), actual[r].end()) == actual[r].end());
    }
    assert(stats.pair_duplicates == reference.pair_duplicates);
    assert(stats.single_duplicates == reference.single_duplicates);
    MpiMarkdupStreamingWindowDestroy(window);
}

static void TestBoundaryBatches() {
    for (uint32_t seed = 1; seed <= 40; ++seed) {
        uint32_t random = seed;
        MpiMarkdupStreamingWindow *window =
            MpiMarkdupStreamingWindowCreate(2, 300, 0, 999, 0, 1999);
        assert(window);
        std::vector<MpiMarkdupCandidateShared> all;
        std::vector<unsigned char> qnames;
        std::vector<std::vector<uint64_t> > actual, expected;
        MpiMarkdupStats stats = {}, reference = {};
        for (unsigned batch = 0; batch < 8; ++batch) {
            std::vector<MpiMarkdupCandidateShared> candidates;
            // Some batches have no boundary updates; others revisit old keys.
            if (batch % 3 != 1) {
                for (unsigned i = 0; i < 64; ++i) {
                    random = random * 1664525u + 1013904223u;
                    const uint64_t order = batch * 64 + i;
                    auto candidate = Candidate(order, random % 3 == 0,
                                               random % 4 != 0);
                    candidate.key.this_coord += (random >> 8) % 17;
                    if (batch >= 6) candidate.key.this_coord += 1000;
                    candidate.key.orientation = (random >> 16) % 2;
                    candidate.key.leftmost = (random >> 17) % 2;
                    candidate.score = (random >> 20) % 4;
                    candidate.qc_fail = (random >> 26) % 2;
                    if (!candidate.key.single) {
                        candidate.paired_marker = 0;
                        candidate.qname_offset = qnames.size();
                        candidate.qname_len = 1;
                        qnames.push_back((unsigned char)('a' + random % 5));
                    }
                    candidates.push_back(candidate);
                    all.push_back(candidate);
                }
            }
            assert(MpiMarkdupStreamingWindowProcess(window, candidates, qnames,
                       0, 10000, &actual, &stats) == 0);
        }
        assert(MpiMarkdupFindDuplicatesStreamingHash(&all, qnames, 2,
                   &expected, &reference) == 0);
        for (size_t r = 0; r < 2; ++r) {
            std::sort(actual[r].begin(), actual[r].end());
            std::sort(expected[r].begin(), expected[r].end());
            assert(actual[r] == expected[r]);
        }
        MpiMarkdupStreamingWindowDestroy(window);
    }
}

static int Apply(void *opaque, uint64_t ordinal, int action) {
    auto &states = *static_cast<std::map<uint64_t, int> *>(opaque);
    if (states.count(ordinal) && states[ordinal] == RB_MARKDUP_DECISION_DUPLICATE) {
        assert(action == RB_MARKDUP_DECISION_DUPLICATE);
    }
    states[ordinal] = action;
    return 0;
}

static void TestOptimistic(bool boundary) {
    MpiMarkdupStreamingWindow *window = MpiMarkdupStreamingWindowCreate(
        2, 300, boundary ? 0 : -1, boundary ? 999 : -1, -1, -1);
    assert(window);
    std::map<uint64_t, int> states;
    std::vector<MpiMarkdupCandidateShared> all;
    std::vector<unsigned char> qnames;
    std::vector<std::vector<MpiMarkdupDecisionShared> > decisions(2);
    MpiMarkdupStats stats = {}, reference = {};
    const auto route = [&]() {
        for (const auto &source : decisions) {
            for (const auto &decision : source) {
                Apply(&states, decision.ordinal, decision.duplicate);
            }
        }
        for (auto &source : decisions) source.clear();
    };
    for (unsigned batch = 0; batch < 6; ++batch) {
        std::vector<MpiMarkdupCandidateShared> candidates;
        if (batch % 3 != 1) {
            for (unsigned i = 0; i < 12; ++i) {
                auto candidate = Candidate(batch * 12 + i, i < 4, i < 8);
                if (!candidate.key.single) {
                    candidate.qname_offset = qnames.size();
                    candidate.qname_len = 1;
                    qnames.push_back((unsigned char)('a' + i % 3));
                }
                candidates.push_back(candidate);
                all.push_back(candidate);
                if (MdIsDecisionCandidate(candidate)) {
                    states[candidate.ordinal] = candidate.source_rank == 0
                        ? RB_MARKDUP_DECISION_KEEP : RB_MARKDUP_DECISION_PENDING;
                }
            }
        }
        assert(MpiMarkdupStreamingWindowProcessOptimistic(window, candidates,
                   qnames, 0, 1050, 0, Apply, &states, &decisions, &stats) == 0);
        route();
    }
    assert(MpiMarkdupStreamingWindowFinalizeOptimistic(window, 0, Apply,
               &states, &decisions) == 0);
    route();
    std::vector<std::vector<uint64_t> > duplicates;
    assert(MpiMarkdupFindDuplicatesStreamingHash(&all, qnames, 2,
               &duplicates, &reference) == 0);
    std::map<uint64_t, int> expected;
    for (const auto &candidate : all) {
        if (MdIsDecisionCandidate(candidate)) expected[candidate.ordinal] = 0;
    }
    for (const auto &source : duplicates) {
        for (uint64_t ordinal : source) expected[ordinal] = 1;
    }
    assert(states == expected);
    assert(window->carry.empty() && window->boundary_carry.empty());
    if (!boundary) {
        assert(window->boundary_markers.capacity > 0);
        assert(window->boundary_carry.capacity() == 0);
    }
    MpiMarkdupStreamingWindowDestroy(window);
}

struct PendingCalls {
    std::map<uint64_t, int> states;
    std::map<uint64_t, unsigned> pending_calls;
};

static int CountPending(void *opaque, uint64_t ordinal, int action) {
    auto &calls = *static_cast<PendingCalls *>(opaque);
    if (action == RB_MARKDUP_DECISION_PENDING) ++calls.pending_calls[ordinal];
    return Apply(&calls.states, ordinal, action);
}

static void TestPendingOnce(bool boundary) {
    auto *window = MpiMarkdupStreamingWindowCreate(2, 300,
        boundary ? 0 : -1, boundary ? 999 : -1, -1, -1);
    assert(window);
    PendingCalls calls;
    MpiMarkdupStats stats = {};
    std::vector<unsigned char> qnames;
    std::vector<std::vector<MpiMarkdupDecisionShared> > decisions(2);
    for (unsigned batch = 0; batch < 8; ++batch) {
        auto candidate = Candidate(batch * 2, false, false);
        candidate.score = batch == 0 ? 100 : 1;
        std::vector<MpiMarkdupCandidateShared> candidates(1, candidate);
        assert(MpiMarkdupStreamingWindowProcessOptimistic(window, candidates,
                   qnames, 0, 1050, 0, CountPending, &calls, &decisions, &stats) == 0);
    }
    assert(calls.pending_calls.size() == 1 && calls.pending_calls[0] == 1);
    assert(MpiMarkdupStreamingWindowFinalizeOptimistic(window, 0,
               CountPending, &calls, &decisions) == 0);
    assert(calls.states[0] == RB_MARKDUP_DECISION_KEEP);
    assert(stats.pair_duplicates == 7);
    MpiMarkdupStreamingWindowDestroy(window);
}

static void TestInteriorMarkers(bool observe, bool boundary = true) {
    auto *window = MpiMarkdupStreamingWindowCreate(2, 300,
        boundary ? 0 : -1, boundary ? 999 : -1,
        boundary ? 0 : -1, boundary ? 1999 : -1);
    PendingCalls calls;
    MpiMarkdupStats stats = {}, reference = {};
    std::vector<unsigned char> qnames;
    std::vector<MpiMarkdupCandidateShared> all;
    std::vector<std::vector<MpiMarkdupDecisionShared> > decisions(2);
    const auto route = [&]() {
        for (auto &source : decisions) {
            for (const auto &decision : source) Apply(&calls.states, decision.ordinal, decision.duplicate);
            source.clear();
        }
    };
    for (unsigned batch = 0; batch < 128; ++batch) {
        std::vector<MpiMarkdupCandidateShared> candidates;
        for (unsigned i = 0; i < 1024; ++i) {
            auto candidate = Candidate(batch * 2048 + i, true);
            candidate.key.this_coord = 10000 + i % 8;
            candidates.push_back(candidate);
            all.push_back(candidate);
        }
        if (batch % 3 == 0) {
            auto single = Candidate(batch * 2048 + 1024, false);
            single.key.this_coord = 10000 + batch % 16;
            calls.states[single.ordinal] = RB_MARKDUP_DECISION_KEEP;
            candidates.push_back(single);
            all.push_back(single);
        }
        assert(MpiMarkdupStreamingWindowProcessOptimistic(window, candidates, qnames,
            0, 10000 + batch / 16, 0, CountPending, &calls, &decisions, &stats) == 0);
        route();
    }
    const size_t retained = std::count_if(window->carry.begin(), window->carry.end(),
        [](const MpiMarkdupCandidateShared &candidate) { return candidate.paired_marker != 0; });
    printf("interior_marker_carry=%zu input_markers=131072 group_s=%.6f working_bytes=%zu\n",
        retained, stats.t_group, MpiMarkdupStreamingWindowMemory(window));
    if (!observe) assert(retained == 8);
    assert(MpiMarkdupStreamingWindowFinalizeOptimistic(window, 0, CountPending,
        &calls, &decisions) == 0);
    route();
    std::vector<std::vector<uint64_t> > duplicates;
    assert(MpiMarkdupFindDuplicatesStreamingHash(&all, qnames, 2, &duplicates, &reference) == 0);
    std::map<uint64_t, int> expected;
    for (const auto &candidate : all)
        if (MdIsDecisionCandidate(candidate)) expected[candidate.ordinal] = RB_MARKDUP_DECISION_KEEP;
    for (const auto &source : duplicates)
        for (uint64_t ordinal : source) expected[ordinal] = RB_MARKDUP_DECISION_DUPLICATE;
    assert(calls.states == expected);
    MpiMarkdupStreamingWindowDestroy(window);
}

static void TestInteriorMixedReference(bool boundary) {
    for (uint32_t seed = 1; seed <= 40; ++seed) {
        uint32_t random = seed;
        auto *window = MpiMarkdupStreamingWindowCreate(2, 300,
            boundary ? 0 : -1, boundary ? 999 : -1,
            boundary ? 0 : -1, boundary ? 1999 : -1);
        PendingCalls calls;
        MpiMarkdupStats stats = {}, reference = {};
        std::vector<MpiMarkdupCandidateShared> all;
        std::vector<unsigned char> full_qnames;
        std::vector<std::vector<MpiMarkdupDecisionShared> > decisions(2);
        const auto route = [&]() {
            for (auto &source : decisions) {
                for (const auto &decision : source) Apply(&calls.states, decision.ordinal, decision.duplicate);
                source.clear();
            }
        };
        for (unsigned batch = 0; batch < 20; ++batch) {
            std::vector<MpiMarkdupCandidateShared> candidates;
            std::vector<unsigned char> qnames;
            for (unsigned i = 0; i < 128; ++i) {
                random = random * 1664525u + 1013904223u;
                auto candidate = Candidate(batch * 128 + i, random % 3 == 0, random % 4 != 0);
                candidate.key.this_coord = 10000 + (random >> 8) % 23;
                candidate.key.orientation = (random >> 16) % 2;
                candidate.key.leftmost = (random >> 17) % 2;
                candidate.score = (random >> 20) % 7;
                candidate.qc_fail = (random >> 26) % 2;
                if (!candidate.key.single) {
                    candidate.paired_marker = 0;
                    candidate.qname_offset = qnames.size();
                    candidate.qname_len = 20 + random % 53;
                    qnames.insert(qnames.end(), candidate.qname_len, (unsigned char)('a' + random % 9));
                }
                candidates.push_back(candidate);
                auto original = candidate;
                original.qname_offset += full_qnames.size();
                all.push_back(original);
                if (MdIsDecisionCandidate(candidate)) calls.states[candidate.ordinal] =
                    candidate.source_rank == 0 ? RB_MARKDUP_DECISION_KEEP : RB_MARKDUP_DECISION_PENDING;
            }
            full_qnames.insert(full_qnames.end(), qnames.begin(), qnames.end());
            assert(MpiMarkdupStreamingWindowProcessOptimistic(window, candidates, qnames,
                0, 10000 + batch, 0, CountPending, &calls, &decisions, &stats) == 0);
            route();
        }
        std::vector<MpiMarkdupCandidateShared> empty;
        std::vector<unsigned char> empty_qnames;
        assert(MpiMarkdupStreamingWindowProcessOptimistic(window, empty, empty_qnames,
            0, 11000, 0, CountPending, &calls, &decisions, &stats) == 0);
        route();
        assert(window->carry.empty());
        assert(MpiMarkdupStreamingWindowFinalizeOptimistic(window, 0, CountPending,
            &calls, &decisions) == 0);
        route();
        std::vector<std::vector<uint64_t> > duplicates;
        assert(MpiMarkdupFindDuplicatesStreamingHash(&all, full_qnames, 2, &duplicates, &reference) == 0);
        std::map<uint64_t, int> expected;
        for (const auto &candidate : all)
            if (MdIsDecisionCandidate(candidate)) expected[candidate.ordinal] = RB_MARKDUP_DECISION_KEEP;
        for (const auto &source : duplicates)
            for (uint64_t ordinal : source) expected[ordinal] = RB_MARKDUP_DECISION_DUPLICATE;
        assert(calls.states == expected);
        assert(stats.pair_duplicates == reference.pair_duplicates);
        // Online singles can lose to another single before a later marker
        // arrives; final per-record decisions, not group counters, are compared.
        for (const auto &count : calls.pending_calls) assert(count.second == 1);
        MpiMarkdupStreamingWindowDestroy(window);
    }
}

static void TestIncrementalReference() {
    for (uint32_t seed = 1; seed <= 80; ++seed) {
        uint32_t random = seed;
        auto *window = MpiMarkdupStreamingWindowCreate(2, 300, 0, 999, 0, 1999);
        PendingCalls calls;
        std::vector<MpiMarkdupCandidateShared> all;
        std::vector<unsigned char> full_qnames;
        std::vector<std::vector<MpiMarkdupDecisionShared> > decisions(2);
        MpiMarkdupStats stats = {}, reference = {};
        for (unsigned batch = 0; batch < 12; ++batch) {
            std::vector<MpiMarkdupCandidateShared> candidates;
            std::vector<unsigned char> qnames;
            for (unsigned i = 0; i < 128; ++i) {
                random = random * 1664525u + 1013904223u;
                auto candidate = Candidate(batch * 128 + i, random % 5 == 0, random % 3 != 0);
                candidate.key.this_coord = (i % 3 == 0 ? 1000 : i % 3 == 1 ? 2000 :
                                               5000 + batch * 1000) + (random >> 8) % 11;
                candidate.key.orientation = (random >> 16) % 2;
                candidate.key.leftmost = (random >> 17) % 2;
                candidate.score = (random >> 20) % 7;
                candidate.qc_fail = (random >> 26) % 2;
                if (!candidate.key.single) {
                    candidate.paired_marker = 0;
                    candidate.qname_offset = qnames.size();
                    candidate.qname_len = 20 + random % 53;
                    qnames.insert(qnames.end(), candidate.qname_len, (unsigned char)('a' + random % 9));
                }
                auto original = candidate;
                original.qname_offset += full_qnames.size();
                all.push_back(original);
                candidates.push_back(candidate);
                if (MdIsDecisionCandidate(candidate)) calls.states[candidate.ordinal] =
                    candidate.source_rank == 0 ? RB_MARKDUP_DECISION_KEEP : RB_MARKDUP_DECISION_PENDING;
            }
            full_qnames.insert(full_qnames.end(), qnames.begin(), qnames.end());
            assert(MpiMarkdupStreamingWindowProcessOptimistic(window, candidates, qnames,
                0, 6000 + batch * 1000, 0, CountPending, &calls, &decisions, &stats) == 0);
            for (auto &source : decisions) {
                for (const auto &decision : source) Apply(&calls.states, decision.ordinal, decision.duplicate);
                source.clear();
            }
            assert(window->boundary_carry.empty() && window->boundary_qnames.empty());
        }
        assert(MpiMarkdupStreamingWindowFinalizeOptimistic(window, 0, CountPending,
            &calls, &decisions) == 0);
        for (const auto &source : decisions)
            for (const auto &decision : source) Apply(&calls.states, decision.ordinal, decision.duplicate);
        std::vector<std::vector<uint64_t> > duplicates;
        assert(MpiMarkdupFindDuplicatesStreamingHash(&all, full_qnames, 2, &duplicates, &reference) == 0);
        std::map<uint64_t, int> expected;
        for (const auto &candidate : all)
            if (MdIsDecisionCandidate(candidate)) expected[candidate.ordinal] = RB_MARKDUP_DECISION_KEEP;
        for (const auto &source : duplicates)
            for (uint64_t ordinal : source) expected[ordinal] = RB_MARKDUP_DECISION_DUPLICATE;
        assert(calls.states == expected);
        assert(stats.pair_duplicates == reference.pair_duplicates);
        assert(stats.single_duplicates == reference.single_duplicates);
        for (const auto &count : calls.pending_calls) assert(count.second == 1);
        assert(window->first_boundary.Memory() == 0 && window->last_boundary.Memory() == 0);
        MpiMarkdupStreamingWindowDestroy(window);
    }
}

static void TestIncrementalRelease() {
    auto *window = MpiMarkdupStreamingWindowCreate(2, 300, 0, 999, 0, 1999);
    PendingCalls calls;
    MpiMarkdupStats stats = {};
    std::vector<unsigned char> qnames;
    std::vector<std::vector<MpiMarkdupDecisionShared> > decisions(2);
    std::vector<MpiMarkdupCandidateShared> candidates;
    for (unsigned i = 0; i < 4096; ++i) {
        auto candidate = Candidate(i * 2, false, false);
        candidate.key.other_coord += i;
        candidates.push_back(candidate);
    }
    assert(MpiMarkdupStreamingWindowProcessOptimistic(window, candidates, qnames,
        0, 10000, 0, CountPending, &calls, &decisions, &stats) == 0);
    assert(window->first_boundary.pair_winners.size() == 4096);
    const size_t peak = MpiMarkdupStreamingWindowMemory(window);
    assert(MpiMarkdupStreamingWindowReleaseBoundaryOptimistic(window, true, false,
        0, CountPending, &calls, &decisions) == 0);
    assert(window->first_closed && !window->last_closed);
    assert(window->first_boundary.Memory() == 0);
    for (const auto &item : calls.states) assert(item.second == RB_MARKDUP_DECISION_KEEP);
    assert(MpiMarkdupStreamingWindowMemory(window) == peak);
    assert(MpiMarkdupStreamingWindowReleaseBoundaryOptimistic(window, true, false,
        0, CountPending, &calls, &decisions) == 0);
    MpiMarkdupStreamingWindowDestroy(window);
}

static void TestSafeProgress() {
    struct Boundary { bool has_records; int first_ref; int64_t first_coord; };
    const std::vector<Boundary> ranges{{true, 1, 100}, {true, 1, 1000},
                                      {false, 0, 0}, {true, 1, 5000}};
    std::vector<long long> progress{1, 1, 1, 900, 1, 1, 1, 1300,
                                   1, 0, 0, 0, 1, 1, 1, 6000};
    assert(!swbam_markdup::BoundaryIsComplete(1, 1000, 100, ranges, progress));
    progress[1] = 0; // The preceding input range has now reached EOF.
    assert(swbam_markdup::BoundaryIsComplete(1, 1000, 100, ranges, progress));
    progress[7] = 1200;
    assert(!swbam_markdup::BoundaryIsComplete(1, 1000, 100, ranges, progress));
    progress[7] = 1201;
    assert(swbam_markdup::BoundaryIsComplete(1, 1000, 100, ranges, progress));
    progress[4] = 0;
    assert(!swbam_markdup::BoundaryIsComplete(1, 1000, 100, ranges, progress));
    const std::vector<Boundary> distant{{true, 1, 100}, {true, 1, 200},
                                       {true, 1, 1000}, {true, 1, 5000}};
    progress = {1, 1, 1, 150, 1, 0, 1, 900, 1, 1, 1, 1300, 1, 1, 1, 6000};
    assert(!swbam_markdup::BoundaryIsComplete(1, 1000, 100, distant, progress));
    progress[1] = 0;
    assert(swbam_markdup::BoundaryIsComplete(1, 1000, 100, distant, progress));
}

static void TestBatchTableCollisions() {
    MdBatchBestTable table;
    std::vector<MpiMarkdupCandidateShared> candidates;
    for (unsigned i = 0; i < 48; ++i) {
        auto candidate = Candidate(i, false);
        candidate.key.this_coord += i;
        candidates.push_back(candidate);
    }
    const auto at = [&](size_t index) -> const MpiMarkdupCandidateShared & {
        return candidates[index];
    };
    assert(table.Prepare(candidates.size(), nullptr, nullptr) == 0);
    for (size_t i = 0; i < candidates.size(); ++i) {
        size_t slot = 0;
        int inserted = 0;
        // Force full hash collisions; equality must still inspect the key.
        assert(table.FindOrInsert(candidates[i].key, 7, i, &slot, &inserted, at) == 0);
        assert(inserted && table.slots[slot].best == i);
        if (i % 2 == 0) table.slots[slot].best = SIZE_MAX;
    }
    for (size_t i = 0; i < candidates.size(); ++i) {
        size_t slot = 0;
        assert(table.FindExisting(candidates[i].key, 7, &slot, at));
        assert(table.slots[slot].best == (i % 2 ? i : SIZE_MAX));
    }
    auto absent = candidates.back();
    ++absent.key.this_coord;
    size_t slot = 0;
    assert(!table.FindExisting(absent.key, 7, &slot, at));
    assert(table.Prepare(1, nullptr, nullptr) == 0);
    assert(!table.FindExisting(candidates[0].key, 7, &slot, at));
    assert(MdFlatTableBytes(table) == table.capacity * (sizeof(MdBatchBestTable::Slot) + 1));
    for (unsigned batch = 0; batch < 600; ++batch) {
        assert(table.Prepare(candidates.size(), nullptr, nullptr) == 0);
        assert(!table.FindExisting(candidates[0].key, 7, &slot, at));
        const size_t index = batch % candidates.size();
        int inserted = 0;
        assert(table.FindOrInsert(candidates[index].key, 7, index,
            &slot, &inserted, at) == 0 && inserted);
    }
}

static void TestPersistentTableGrowth() {
    MdFlatBestTable table;
    for (size_t i = 0; i < 128; ++i) {
        auto candidate = Candidate(i, false);
        candidate.key.this_coord += i;
        size_t slot = 0;
        int inserted = 0;
        assert(table.FindOrInsertPersistent(candidate.key, 7, i, &slot, &inserted) == 0);
        assert(inserted && table.best[slot] == i);
        const size_t capacity = table.capacity;
        assert(table.FindOrInsertPersistent(candidate.key, 7, 999, &slot, &inserted) == 0);
        assert(!inserted && table.best[slot] == i && table.capacity == capacity);
        assert(table.ReservePersistent(table.filled) == 0 && table.capacity == capacity);
    }
    for (size_t i = 0; i < 128; ++i) {
        auto candidate = Candidate(i, false);
        candidate.key.this_coord += i;
        size_t slot = 0;
        assert(table.FindExisting(candidate.key, 7, &slot) && table.best[slot] == i);
    }
}

static void TestQnameArenaCompaction() {
    auto *window = MpiMarkdupStreamingWindowCreate(1, 300, -1, -1, -1, -1);
    PendingCalls calls;
    MpiMarkdupStats stats = {}, reference = {};
    std::vector<std::vector<MpiMarkdupDecisionShared> > decisions(1);
    std::vector<MpiMarkdupCandidateShared> all;
    std::vector<unsigned char> all_names;
    size_t previous_bytes = 0;
    bool compacted = false;
    for (unsigned batch = 0; batch < 100; ++batch) {
        std::vector<MpiMarkdupCandidateShared> candidates;
        std::vector<unsigned char> names;
        const auto append = [&](MpiMarkdupCandidateShared candidate, unsigned char name) {
            candidate.source_rank = 0;
            candidate.qname_offset = names.size();
            candidate.qname_len = 200;
            names.insert(names.end(), candidate.qname_len, name);
            candidates.push_back(candidate);
            candidate.qname_offset += all_names.size();
            all.push_back(candidate);
            calls.states[candidate.ordinal] = RB_MARKDUP_DECISION_KEEP;
        };
        for (unsigned i = 0; i < 128; ++i) {
            auto candidate = Candidate(batch * 256 + i, false, false);
            candidate.key.this_coord = 1000 + batch * 400 + i;
            append(candidate, 'm');
        }
        if (batch == 0 || batch == 99) {
            auto candidate = Candidate(100000 + batch, false, false);
            candidate.key.this_coord = 200000;
            candidate.score = 100;
            append(candidate, batch == 0 ? 'z' : 'a');
        }
        all_names.insert(all_names.end(), names.begin(), names.end());
        assert(MpiMarkdupStreamingWindowProcessOptimistic(window, candidates, names,
            0, 1127 + batch * 400, 0, CountPending, &calls, &decisions, &stats) == 0);
        assert(decisions[0].empty());
        if (window->carry_qnames.size() < previous_bytes) compacted = true;
        previous_bytes = window->carry_qnames.size();
    }
    assert(compacted);
    assert(MpiMarkdupStreamingWindowFinalizeOptimistic(window, 0, CountPending,
        &calls, &decisions) == 0);
    std::vector<std::vector<uint64_t> > duplicates;
    assert(MpiMarkdupFindDuplicatesStreamingHash(&all, all_names, 1,
        &duplicates, &reference) == 0);
    assert((duplicates[0] == std::vector<uint64_t>{100000}));
    for (const auto &state : calls.states)
        assert(state.second == (state.first == 100000
            ? RB_MARKDUP_DECISION_DUPLICATE : RB_MARKDUP_DECISION_KEEP));
    for (const auto &count : calls.pending_calls) assert(count.second == 1);
    MpiMarkdupStreamingWindowDestroy(window);
}

static void TestEmptyQnameCarry() {
    auto *window = MpiMarkdupStreamingWindowCreate(1, 300, -1, -1, -1, -1);
    PendingCalls calls;
    MpiMarkdupStats stats = {};
    std::vector<std::vector<MpiMarkdupDecisionShared> > decisions(1);
    auto expired = Candidate(0, false, false);
    expired.source_rank = 0;
    expired.qname_len = 8;
    auto retained = Candidate(1, false, false);
    retained.source_rank = 0;
    retained.key.this_coord = 10000;
    retained.qname_offset = 8;
    std::vector<unsigned char> names(8, 'x');
    assert(MpiMarkdupStreamingWindowProcessOptimistic(window, {expired}, names,
        0, 1000, 0, CountPending, &calls, &decisions, &stats) == 0);
    assert(MpiMarkdupStreamingWindowProcessOptimistic(window, {retained}, names,
        0, 1400, 0, CountPending, &calls, &decisions, &stats) == 0);
    assert(window->carry_qnames.empty() && window->carry.size() == 1);
    assert(window->carry[0].qname_offset == 0);
    assert(MpiMarkdupStreamingWindowProcessOptimistic(window, {}, {},
        0, 1401, 0, CountPending, &calls, &decisions, &stats) == 0);
    assert(MpiMarkdupStreamingWindowFinalizeOptimistic(window, 0, CountPending,
        &calls, &decisions) == 0);
    MpiMarkdupStreamingWindowDestroy(window);
}

int main(int argc, char **argv) {
    if (argc == 2 && !strcmp(argv[1], "--interior-markers-observe")) {
        TestInteriorMarkers(true);
        return 0;
    }
    const bool observe = argc == 2 && !strcmp(argv[1], "--observe");
    TestMarkerCarry(observe);
    TestDuplicateDecisions();
    TestBoundaryBatches();
    TestOptimistic(true);
    TestOptimistic(false);
    TestPendingOnce(true);
    TestPendingOnce(false);
    TestIncrementalReference();
    TestIncrementalRelease();
    TestSafeProgress();
    TestBatchTableCollisions();
    TestPersistentTableGrowth();
    TestQnameArenaCompaction();
    TestEmptyQnameCarry();
    TestInteriorMarkers(false);
    TestInteriorMarkers(false, false);
    TestInteriorMixedReference(true);
    TestInteriorMixedReference(false);
    puts("markdup window decisions PASS");
}
