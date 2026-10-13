#include "swbam/swbam.h"
#ifdef SWBAM_ENABLE_RAM_STORAGE
#include "storage.h"
#else
#include "storage-native.h"
#endif

#include <cstdio>
#include <mpi.h>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

int main(int argc, char **argv) {
    swbam_usage::Session storage;
    if (storage.Prepare(&argc, &argv, true)) return 2;
    swbam_usage::BeginMpiInterval();
    const double e2e_start = swbam_usage::WallTime();
    MPI_Init(&argc, &argv);
    const double processing_start = swbam_usage::WallTime();
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    const bool layout_ok = swbam::mpi::AllRanksOk(storage.EnterMpi() == 0);
    if (argc != 2 || !layout_ok) {
        if (rank == 0) fprintf(stderr, "Usage: %s input.bam\n", argv[0]);
        const double processing_seconds = swbam_usage::WallTime() - processing_start;
        MPI_Finalize();
        const double elapsed = swbam_usage::WallTime() - e2e_start;
        swbam_usage::EndMpiInterval();
        return storage.Finish(2, elapsed, processing_seconds);
    }
#ifdef PLATFORM_SUNWAY
    athread_init();
#endif
    int status = 1;
    {
        swbam::mpi::MpiBamInput input;
        swbam::mpi::MpiBamInputPlan plan;
        int ok = swbam_usage::OpenInput(input, argv[1], "posix", "") == 0;
        if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
            ok = swbam::mpi::PrepareMpiBamInputPlan(
                *input.backend(), &plan, 0, MPI_COMM_WORLD) == 0;
            if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
                swbam::operators::RecordCountMetrics metrics = {};
                MPI_Barrier(MPI_COMM_WORLD);
                const double start = MPI_Wtime();
                ok = swbam::operators::RunRecordCountPipeline(
                    *input.backend(), plan.rank_spans(),
                    plan.rank_block_count(), &metrics) == 0;
                const double seconds = MPI_Wtime() - start;
                if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
                    long long records = 0;
                    double core = 0.0;
                    MPI_Reduce(&metrics.records, &records, 1, MPI_LONG_LONG,
                               MPI_SUM, 0, MPI_COMM_WORLD);
                    MPI_Reduce(&seconds, &core, 1, MPI_DOUBLE,
                               MPI_MAX, 0, MPI_COMM_WORLD);
                    if (rank == 0) {
                        printf("sdk_count ranks=%d records=%lld blocks=%zu "
                               "input_bytes=%zu core=%.6f\n",
                               ranks, records, plan.blocks.size(),
                               input.backend()->size(), core);
                    }
                    status = 0;
                }
            }
        }
        input.Close();
    }
#ifdef PLATFORM_SUNWAY
    storage.LeaveMpi();
    athread_halt();
#endif
    const double processing_seconds = swbam_usage::WallTime() - processing_start;
    MPI_Finalize();
    const double seconds = swbam_usage::WallTime() - e2e_start;
    swbam_usage::EndMpiInterval();
    return storage.Finish(status, seconds, processing_seconds);
}
