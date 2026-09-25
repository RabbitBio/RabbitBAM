#include "swbam/swbam.h"

#include <cstdio>
#include <mpi.h>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    if (argc != 2) {
        if (rank == 0) fprintf(stderr, "Usage: %s input.bam\n", argv[0]);
        MPI_Finalize();
        return 2;
    }
#ifdef PLATFORM_SUNWAY
    athread_init();
#endif
    int status = 1;
    {
        swbam::mpi::MpiBamInput input;
        swbam::mpi::MpiBamInputPlan plan;
        swbam::AdaptiveRankBodySink body;
        int ok = input.Open(argv[1], "memory", "", MPI_COMM_WORLD) == 0;
        if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
            ok = swbam::mpi::PrepareMpiBamInputPlan(
                *input.backend(), &plan, 0, MPI_COMM_WORLD) == 0;
            if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
                size_t reserve = 0;
                for (size_t i = plan.rank_begin; i < plan.rank_end; ++i)
                    reserve += plan.blocks[i].compressed_size;
                ok = body.Open("memory", 0, reserve,
                               "sdk-bam2bam-body", "") == 0;
                if (!ok) fprintf(stderr, "[rank %d] memory body allocation failed.\n", rank);
                if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
                    swbam::operators::BamTransformMetrics metrics = {};
                    const swbam::BamFilterOptions filter =
                        swbam::DefaultBamFilterOptions();
                    MPI_Barrier(MPI_COMM_WORLD);
                    const double start = MPI_Wtime();
                    ok = swbam::operators::RunBamTransformPipeline(
                        *input.backend(), plan.rank_spans(),
                        plan.rank_block_count(), &body, filter, 1,
                        &metrics) == 0 &&
                         metrics.kept_records == metrics.total_records;
                    const double seconds = MPI_Wtime() - start;
                    if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
                        long long records = 0;
                        unsigned long long body_bytes = 0;
                        double core = 0.0;
                        const unsigned long long local_bytes = body.size();
                        MPI_Reduce(&metrics.total_records, &records, 1,
                                   MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&local_bytes, &body_bytes, 1,
                                   MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0,
                                   MPI_COMM_WORLD);
                        MPI_Reduce(&seconds, &core, 1, MPI_DOUBLE,
                                   MPI_MAX, 0, MPI_COMM_WORLD);
                        if (rank == 0) {
                            printf("sdk_bam2bam ranks=%d records=%lld blocks=%zu "
                                   "input_bytes=%zu body_bytes=%llu core=%.6f\n",
                                   ranks, records, plan.blocks.size(),
                                   input.backend()->size(), body_bytes, core);
                        }
                        status = 0;
                    }
                }
            }
        }
        body.Close();
        input.Close();
    }
#ifdef PLATFORM_SUNWAY
    athread_halt();
#endif
    MPI_Finalize();
    return status;
}
