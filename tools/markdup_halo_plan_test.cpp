// Host metadata test: unused CPE paths are discarded, not emulated.
// Build with mpicxx, -ffunction-sections -fdata-sections and --gc-sections.
extern "C" int __real_athread_spawn(void *, void *, int);
extern "C" int athread_join();
#include "../swbam-app/algorithms/markdup/swbam_markdup_mpi.cpp"

#include <cassert>

static MpiMarkdupCandidateShared Candidate(uint64_t id, int64_t coord, int64_t score) {
    MpiMarkdupCandidateShared result = {};
    result.ordinal = result.global_order = id;
    result.source_rank = 0;
    result.key.this_ref = 1;
    result.key.this_coord = coord;
    result.key.single = 1;
    result.key.orientation = 2;
    result.score = score;
    return result;
}

static int Apply(void *opaque, uint64_t ordinal, int action) {
    auto &states = *static_cast<std::map<uint64_t, int> *>(opaque);
    if (states[ordinal] == RB_MARKDUP_DECISION_DUPLICATE)
        assert(action == RB_MARKDUP_DECISION_DUPLICATE);
    states[ordinal] = action;
    return 0;
}

static void TestSelectionAndPending() {
    MdOnePassBoundaryPlan plan;
    plan.halos.push_back({1, 9700, 10300});
    const auto a = Candidate(0, 9990, 4200);
    const auto u = Candidate(1, 10250, 3500); // POS 10250, 150M
    const auto v = Candidate(2, 10250, 4100); // POS 10320, 70S80M
    auto marker = Candidate(3, 10250, 100);
    marker.paired_marker = 1;
    std::vector<MpiMarkdupCandidateShared> selected;
    std::vector<unsigned char> names, selected_names;
    assert(MdAppendBoundaryCandidates({a, u, v, marker}, names,
        {9990, 10250, 10320, 10260}, 300, &plan, &selected, &selected_names) == 0);
    assert(selected.size() == 3); // v is outside the POS halo, despite its key.
    assert(plan.decisions.at(a.global_order) == kMdOnePassKeep);
    assert(plan.decisions.at(u.global_order) == kMdOnePassUnknown);
    assert(plan.candidates.empty()); // Selection alone never enables main-pass skipping.
    assert(!swbam_markdup::HaloCoversKey(plan.halos, marker.key, 300));

    std::map<uint64_t, int> states;
    states[u.ordinal] = states[v.ordinal] = RB_MARKDUP_DECISION_KEEP;
    auto *window = MpiMarkdupStreamingWindowCreate(1, 300, -1, -1, -1, -1);
    assert(window);
    MpiMarkdupStats stats = {};
    std::vector<std::vector<MpiMarkdupDecisionShared> > remote(1);
    assert(MpiMarkdupStreamingWindowProcessOptimistic(window, {u}, names,
        0, 10299, 0, Apply, &states, &remote, &stats) == 0);
    assert(states[u.ordinal] == RB_MARKDUP_DECISION_PENDING);
    assert(MpiMarkdupStreamingWindowProcessOptimistic(window, {v}, names,
        0, 10319, 0, Apply, &states, &remote, &stats) == 0);
    assert(states[u.ordinal] == RB_MARKDUP_DECISION_DUPLICATE);
    assert(states[v.ordinal] == RB_MARKDUP_DECISION_PENDING);
    assert(MpiMarkdupStreamingWindowProcessOptimistic(window, {}, names,
        0, 10599, 0, Apply, &states, &remote, &stats) == 0);
    assert(states[v.ordinal] == RB_MARKDUP_DECISION_KEEP);
    assert(remote[0].empty());
    MpiMarkdupStreamingWindowDestroy(window);
}

static void TestKeySupport() {
    auto key = Candidate(0, 10000, 0).key;
    int64_t low = 0, high = 0;
    assert(swbam_markdup::KeyPositionSupport(key, 300, &low, &high));
    assert(low == 10000 && high == 10300);
    key.orientation = 3;
    assert(swbam_markdup::KeyPositionSupport(key, 300, &low, &high));
    assert(low == 9700 && high == 10000);
    key.single = 0;
    for (int orientation : {2, 3, 5, 7}) {
        for (int leftmost : {11, 13}) {
            key.orientation = orientation;
            key.leftmost = leftmost;
            const bool end = orientation == 7 ? leftmost == 11 : leftmost == 13;
            assert(swbam_markdup::KeyPositionSupport(key, 300, &low, &high));
            assert(low == (end ? 9700 : 10000));
            assert(high == (end ? 10000 : 10300));
        }
    }
    key.leftmost = 0;
    assert(!swbam_markdup::KeyPositionSupport(key, 300, &low, &high));
    std::vector<swbam_markdup::HaloInterval> halos = {{1, 10, 20}, {1, 21, 30}, {2, 5, 8}};
    swbam_markdup::MergeHalos(&halos);
    assert(halos.size() == 2 && halos[0].low == 10 && halos[0].high == 30);
    assert(swbam_markdup::HaloOffset(INT64_MAX - 1, 300) == INT64_MAX);
    assert(swbam_markdup::HaloOffset(INT64_MIN + 1, -300) == INT64_MIN);
}

static void TestPartialSelection() {
    MdOnePassBoundaryPlan plan;
    plan.partial = true;
    plan.halos.push_back({1, 9700, 10300});
    const auto low = Candidate(0, 9990, 3500);
    const auto high = Candidate(1, 9990, 4100);
    auto marker = Candidate(2, 9990, 1);
    marker.paired_marker = 1;
    std::vector<MpiMarkdupCandidateShared> selected;
    std::vector<unsigned char> names;
    assert(swbam_markdup::HaloCoversKey(plan.halos, high.key, 300));
    assert(MdAppendBoundaryCandidates({low, high, marker}, {},
        {9990, 10020, 10030}, 300, &plan, &selected, &names) == 0);
    // Logical halo coverage is not proof that the fixed probe read all data.
    assert(plan.decisions.at(low.global_order) == kMdOnePassUnknown);
    assert(plan.decisions.at(high.global_order) == kMdOnePassUnknown);
    assert(plan.candidates.empty());
    MpiMarkdupStats stats = {};
    std::vector<std::vector<uint64_t> > duplicates;
    assert(MpiMarkdupFindDuplicatesStreamingHash(&selected, names, 1,
        &duplicates, &stats) == 0);
    // A paired marker beats even the higher-scoring singleton; it must survive
    // into the main pass to suppress additional, not-yet-probed singletons.
    std::sort(duplicates[0].begin(), duplicates[0].end());
    assert((duplicates[0] == std::vector<uint64_t>{0, 1}));
}

int main() {
    TestKeySupport();
    TestSelectionAndPending();
    TestPartialSelection();
    puts("PASS POS halo selection, partial DUP-only prepass, later winner, paired-key roles and retained marker");
}
