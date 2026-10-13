#include "swbam/swbam.h"
#ifdef SWBAM_ENABLE_RAM_STORAGE
#include "storage.h"
#else
#include "storage-native.h"
#endif

#include <cstdio>
#include <cstdlib>
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
    int rank = 0, ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    const bool layout_ok = swbam::mpi::AllRanksOk(storage.EnterMpi() == 0);
    if (argc != 2 || storage.output_path().empty() || !layout_ok) {
        if (rank == 0) fprintf(stderr, "Usage: %s input.bam -o output.bam [--storage ram ...]\n", argv[0]);
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
        swbam_usage::BodySink body;
        const auto run = [&]() -> int {
            int ok = swbam_usage::OpenInput(input, argv[1], "posix", "") == 0;
            if (!swbam::mpi::AllRanksOk(ok)) return 1;
            if (input.PreparePlan(&plan)) return 1;
            ok = body.Open("spool", 0, 8u * 1024u * 1024u, "sdk-body", "") == 0;
            if (!swbam::mpi::AllRanksOk(ok)) return 1;
            swbam::operators::BamTransformMetrics metrics = {};
            MPI_Barrier(MPI_COMM_WORLD);
            const double start = MPI_Wtime();
            ok = swbam::operators::RunBamTransformPipeline(
                *input.backend(), plan.rank_spans(), plan.rank_block_count(),
                &body, swbam::DefaultBamFilterOptions(), 1, &metrics) == 0 &&
                metrics.kept_records == metrics.total_records && body.Flush() == 0;
            const double seconds = MPI_Wtime() - start;
            if (!swbam::mpi::AllRanksOk(ok)) return 1;
            long long records = 0;
            unsigned long long bytes = 0, local_bytes = body.size();
            double core = 0.0;
            MPI_Reduce(&metrics.total_records, &records, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local_bytes, &bytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&seconds, &core, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            swbam::MemoryBamOutput prefix;
            if (rank == 0) {
                ok = swbam_usage::BuildBamPrefix(
                    const_cast<sam_hdr_t *>(input.backend()->header()), 1, &prefix) == 0;
                printf("sdk_bam2bam ranks=%d records=%lld blocks=%zu input_bytes=%zu body_bytes=%llu core=%.6f\n",
                       ranks, records, plan.blocks.size(), input.backend()->size(), bytes, core);
            }
            if (!swbam::mpi::AllRanksOk(ok)) return 1;
            const unsigned char eof[28] = {31,139,8,4,0,0,0,0,0,255,6,0,66,67,2,0,27,0,3,0,0,0,0,0,0,0,0,0};
            return swbam_usage::WriteOutput(storage.output_path(), body,
                prefix.size(), prefix.data(), sizeof(eof), eof) == 0 ? 0 : 1;
        };
        status = run();
        input.Close();
        if (body.Close()) status = 1;
    }
    storage.LeaveMpi();
#ifdef PLATFORM_SUNWAY
    athread_halt();
#endif
    const double processing_seconds = swbam_usage::WallTime() - processing_start;
    MPI_Finalize();
    const double seconds = swbam_usage::WallTime() - e2e_start;
    swbam_usage::EndMpiInterval();
    return storage.Finish(status, seconds, processing_seconds);
}
