#include "swbam/mpi_runtime.h"

#include <climits>
#include <vector>

namespace swbam {
namespace mpi {

int AllRanksOk(int local_ok, MPI_Comm comm) {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, comm);
    return global_ok;
}

double ReduceMaxCost(double local_cost, int root, MPI_Comm comm) {
    double max_cost = 0.0;
    MPI_Reduce(&local_cost, &max_cost, 1, MPI_DOUBLE, MPI_MAX, root, comm);
    return max_cost;
}

int PrepareMpiBamInputPlan(const BamInputBackend &input,
                           MpiBamInputPlan *plan,
                           int root, MPI_Comm comm) {
    if (!plan) return -1;

    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);
    if (root < 0 || root >= comm_size) return -1;

    int local_ok = 1;
    unsigned long long body_offset = 0;
    long long block_count = 0;
    if (rank == root) {
        plan->blocks.clear();
        if (input.ScanBlocks(&plan->blocks) != 0 ||
            plan->blocks.size() > static_cast<size_t>(INT_MAX)) {
            local_ok = 0;
        } else {
            body_offset = static_cast<unsigned long long>(input.body_offset());
            block_count = static_cast<long long>(plan->blocks.size());
        }
    }

    MPI_Bcast(&local_ok, 1, MPI_INT, root, comm);
    if (!local_ok) return -1;
    MPI_Bcast(&body_offset, 1, MPI_UNSIGNED_LONG_LONG, root, comm);
    MPI_Bcast(&block_count, 1, MPI_LONG_LONG, root, comm);
    if (block_count < 0 || block_count > INT_MAX) return -1;

    if (rank != root) plan->blocks.resize(static_cast<size_t>(block_count));
    std::vector<unsigned long long> offsets(static_cast<size_t>(block_count));
    std::vector<unsigned int> lengths(static_cast<size_t>(block_count));
    if (rank == root) {
        for (size_t i = 0; i < plan->blocks.size(); ++i) {
            offsets[i] = static_cast<unsigned long long>(plan->blocks[i].offset);
            lengths[i] = static_cast<unsigned int>(
                plan->blocks[i].compressed_size);
        }
    }
    if (block_count > 0) {
        MPI_Bcast(offsets.data(), static_cast<int>(block_count),
                  MPI_UNSIGNED_LONG_LONG, root, comm);
        MPI_Bcast(lengths.data(), static_cast<int>(block_count),
                  MPI_UNSIGNED, root, comm);
    }
    if (rank != root) {
        for (size_t i = 0; i < plan->blocks.size(); ++i) {
            plan->blocks[i].offset = static_cast<uint64_t>(offsets[i]);
            plan->blocks[i].compressed_size = static_cast<uint32_t>(lengths[i]);
        }
    }

    plan->body_offset = static_cast<uint64_t>(body_offset);
    size_t n = plan->blocks.size();
    plan->rank_begin = n * static_cast<size_t>(rank) /
                       static_cast<size_t>(comm_size);
    plan->rank_end = n * static_cast<size_t>(rank + 1) /
                     static_cast<size_t>(comm_size);
    return 0;
}

} // namespace mpi
} // namespace swbam
