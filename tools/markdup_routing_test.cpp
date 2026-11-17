// Host MPI metadata regression; unused CPE and I/O paths are discarded.
extern "C" int __real_athread_spawn(void *, void *, int);
extern "C" int athread_join();
#include "../swbam-app/algorithms/markdup/swbam_markdup_mpi.cpp"

#include <cassert>

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    assert(ranks >= 2);
    for (int scenario = 0; scenario < 4; ++scenario) {
        const bool empty = scenario == 1 && rank == ranks - 1;
        const bool unplaced = scenario == 2 ||
            (scenario == 0 && rank == ranks - 1);
        const int tid = unplaced || empty ? -1 : 0;
        const int pos = unplaced || empty ? -1 :
            (scenario == 3 ? 100 : rank * 1000 + 100);
        std::vector<MdCoordBoundary> cuts;
        assert(MdBuildFirstCoordBoundaries(!empty, tid, pos, ranks, &cuts) == 0);
        assert(cuts.size() == (size_t)ranks);
        const int mapped_ranks = scenario == 2 ? 0 :
            scenario == 0 || scenario == 1 ? ranks - 1 : ranks;
        for (int r = 0; r < ranks; ++r)
            assert(cuts[(size_t)r].has_records == (r < mapped_ranks));
        if (mapped_ranks) {
            for (int64_t coord : {-100, 1, 101, 1101, 99999}) {
                MpiMarkdupKeyShared key = {};
                key.this_ref = 1;
                key.this_coord = coord;
                const int owner = MdCoordinateOwnerByFirst(key, cuts, 0);
                assert(owner >= 0 && owner < mapped_ranks);
                for (int source = 0; source < ranks; ++source)
                    assert(MdCoordinateOwnerByFirst(key, cuts, source) == owner);
            }
        }
        if (rank == 0)
            printf("PASS routing ranks=%d scenario=%d mapped_ranges=%d\n",
                   ranks, scenario, mapped_ranks);
    }
    std::vector<MdCoordBoundary> cuts;
    assert(MdBuildFirstCoordBoundaries(1, 0, (ranks - rank) * 100,
                                      ranks, &cuts) != 0);
    if (rank == 0) puts("PASS routing still rejects decreasing mapped coordinates");
    MPI_Finalize();
    return 0;
}
