#ifndef SWBAM_MPI_RUNTIME_H
#define SWBAM_MPI_RUNTIME_H

#include "swbam/io.h"

#include <mpi.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace swbam {
namespace mpi {

struct MpiBamInputPlan {
    std::vector<BgzfBlockSpan> blocks;
    size_t rank_begin;
    size_t rank_end;
    uint64_t body_offset;

    MpiBamInputPlan()
        : rank_begin(0), rank_end(0), body_offset(0) {}

    const BgzfBlockSpan *rank_spans() const {
        return rank_begin < rank_end ? blocks.data() + rank_begin : nullptr;
    }
    size_t rank_block_count() const { return rank_end - rank_begin; }
};

int AllRanksOk(int local_ok, MPI_Comm comm = MPI_COMM_WORLD);
double ReduceMaxCost(double local_cost, int root = 0,
                     MPI_Comm comm = MPI_COMM_WORLD);
int PrepareMpiBamInputPlan(const BamInputBackend &input,
                           MpiBamInputPlan *plan,
                           int root = 0,
                           MPI_Comm comm = MPI_COMM_WORLD);

} // namespace mpi
} // namespace swbam

#endif
