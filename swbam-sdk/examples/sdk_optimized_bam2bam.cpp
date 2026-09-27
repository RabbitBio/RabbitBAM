#include "swbam/swbam.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include <unistd.h>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    const char *diagnostics = std::getenv("SWBAM_DIAGNOSTICS");
    const bool trace = diagnostics && std::strcmp(diagnostics, "0") != 0;
    const double diagnostic_start = trace ? MPI_Wtime() : 0.0;
    const auto mark = [&](const char *stage, int ok = 1, size_t value = 0) {
        if (!trace) return;
        fprintf(stderr, "[sdk-diag rank=%d pid=%ld] elapsed=%.6f stage=%s ok=%d value=%zu\n",
                rank, static_cast<long>(getpid()), MPI_Wtime() - diagnostic_start,
                stage, ok, value);
        fflush(stderr);
    };
    if (argc != 2) {
        if (rank == 0) fprintf(stderr, "Usage: %s input.bam\n", argv[0]);
        MPI_Finalize();
        return 2;
    }
#ifdef PLATFORM_SUNWAY
    mark("athread_init.begin");
    athread_init();
    mark("athread_init.done");
#endif
    int status = 1;
    {
        swbam::mpi::MpiBamInput input;
        swbam::mpi::MpiBamInputPlan plan;
        swbam::AdaptiveRankBodySink body;
        mark("open.begin");
        int ok = input.Open(argv[1], "memory", "", MPI_COMM_WORLD) == 0;
        mark("open.done", ok);
        if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
            mark("scan_plan.begin");
            ok = swbam::mpi::PrepareMpiBamInputPlan(
                *input.backend(), &plan, 0, MPI_COMM_WORLD) == 0;
            mark("scan_plan.done", ok, plan.rank_block_count());
            if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
                size_t reserve = 0;
                for (size_t i = plan.rank_begin; i < plan.rank_end; ++i)
                    reserve += plan.blocks[i].compressed_size;
                mark("body_reserve.begin", 1, reserve);
                ok = body.Open("memory", 0, reserve,
                               "sdk-bam2bam-body", "") == 0;
                mark("body_reserve.done", ok, reserve);
                if (!ok) fprintf(stderr, "[rank %d] memory body allocation failed.\n", rank);
                if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
                    swbam::operators::BamTransformMetrics metrics = {};
                    const swbam::BamFilterOptions filter =
                        swbam::DefaultBamFilterOptions();
                    mark("core_barrier.begin");
                    MPI_Barrier(MPI_COMM_WORLD);
                    mark("pipeline.begin");
                    const double start = MPI_Wtime();
                    ok = swbam::operators::RunBamTransformPipeline(
                        *input.backend(), plan.rank_spans(),
                        plan.rank_block_count(), &body, filter, 1,
                        &metrics) == 0 &&
                         metrics.kept_records == metrics.total_records;
                    const double seconds = MPI_Wtime() - start;
                    mark("pipeline.done", ok, body.size());
                    mark("result_collectives.begin", ok);
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
                        mark("result_collectives.done");
                    }
                }
            }
        }
        mark("body_close.begin");
        body.Close();
        mark("input_close.begin");
        input.Close();
        mark("input_close.done");
    }
#ifdef PLATFORM_SUNWAY
    mark("athread_halt.begin");
    athread_halt();
    mark("athread_halt.done");
#endif
    mark("mpi_finalize.begin");
    MPI_Finalize();
    return status;
}
