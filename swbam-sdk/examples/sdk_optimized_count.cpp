#include "swbam/swbam.h"
#include "sdk_chunked_bam.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    const char *progress_env = std::getenv("SWBAM_SDK_PROGRESS");
    const bool progress =
        (progress_env && std::strcmp(progress_env, "1") == 0) ||
        (argc == 3 && std::strcmp(argv[2], "--chunked-memory-progress") == 0);
    const double progress_start = MPI_Wtime();
    const auto trace = [&](const char *step, size_t chunk, size_t blocks) {
        if (!progress) return;
        fprintf(stderr, "[sdk-progress rank=%d] elapsed=%.3f step=%s "
                "chunk=%zu blocks=%zu\n", rank,
                MPI_Wtime() - progress_start, step, chunk, blocks);
        fflush(stderr);
    };
    const size_t chunk_limit = argc == 3
        ? sdk_bench::ChunkLimit(argv[2]) : 0;
    const bool chunked = chunk_limit != 0;
    if (argc != 2 && !chunked) {
        if (rank == 0) fprintf(stderr,
            "Usage: %s input.bam [--chunked-memory[=MiB]|--chunked-memory-progress]\n",
            argv[0]);
        MPI_Finalize();
        return 2;
    }
#ifdef PLATFORM_SUNWAY
    athread_init();
#endif
    int status = 1;
    {
        swbam::mpi::MpiBamInput input;
        sdk_bench::ChunkedBamInput chunk_input;
        swbam::mpi::MpiBamInputPlan plan;
        const double open_start = MPI_Wtime();
        int ok = chunked ? chunk_input.Open(argv[1]) == 0
                         : input.Open(argv[1], "memory", "", MPI_COMM_WORLD) == 0;
        const double open_wall = MPI_Wtime() - open_start;
        swbam::BamInputBackend *backend = chunked
            ? static_cast<swbam::BamInputBackend *>(&chunk_input)
            : input.backend();
        if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
            double index_prepare = 0.0;
            if (chunked && rank == 0) {
                trace("index_begin", 0, 0);
                const double start = MPI_Wtime();
                ok = chunk_input.PrepareIndex() == 0;
                index_prepare = MPI_Wtime() - start;
                trace("index_done", 0, 0);
            }
            if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
            trace("plan_begin", 0, 0);
            const double plan_start = MPI_Wtime();
            ok = swbam::mpi::PrepareMpiBamInputPlan(
                *backend, &plan, 0, MPI_COMM_WORLD) == 0;
            const double plan_wall = MPI_Wtime() - plan_start;
            trace("plan_done", 0, plan.rank_block_count());
            if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
                long long local_records = 0;
                double seconds = 0.0;
                double stage_read = 0.0, stage_prepare = 0.0;
                unsigned long long chunks = 0;
                const swbam::BgzfBlockSpan *spans = plan.rank_spans();
                const size_t count = plan.rank_block_count();
                if (!chunked) MPI_Barrier(MPI_COMM_WORLD);
                if (count == 0) {
                    swbam::operators::RecordCountMetrics metrics = {};
                    const double start = MPI_Wtime();
                    ok = swbam::operators::RunRecordCountPipeline(
                        *backend, nullptr, 0, &metrics) == 0;
                    seconds += MPI_Wtime() - start;
                }
                for (size_t begin = 0; begin < count && ok;) {
                    const size_t end = chunked
                        ? sdk_bench::NextChunkEnd(spans, count, begin,
                                                   chunk_limit)
                        : count;
                    trace("stage_begin", chunks + 1, end - begin);
                    if (chunked && chunk_input.Stage(
                            spans + begin, end - begin,
                            chunk_limit,
                            &stage_read, &stage_prepare) != 0) {
                        ok = 0;
                        break;
                    }
                    trace("stage_done", chunks + 1, end - begin);
                    swbam::operators::RecordCountMetrics metrics = {};
                    trace("core_begin", chunks + 1, end - begin);
                    const double start = MPI_Wtime();
                    ok = swbam::operators::RunRecordCountPipeline(
                        *backend, spans + begin, end - begin, &metrics) == 0;
                    seconds += MPI_Wtime() - start;
                    trace("core_done", chunks + 1, end - begin);
                    local_records += metrics.records;
                    begin = end;
                    ++chunks;
                }
                if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
                    long long records = 0;
                    double core = 0.0;
                    double max_open = 0.0, max_plan = 0.0;
                    double max_index_prepare = 0.0;
                    double max_stage_read = 0.0, max_stage_prepare = 0.0;
                    unsigned long long max_chunk_capacity = 0;
                    unsigned long long max_chunks = 0;
                    const unsigned long long local_capacity =
                        chunk_input.peak_capacity();
                    const double reduce_start = MPI_Wtime();
                    MPI_Reduce(&local_records, &records, 1, MPI_LONG_LONG,
                               MPI_SUM, 0, MPI_COMM_WORLD);
                    MPI_Reduce(&seconds, &core, 1, MPI_DOUBLE,
                               MPI_MAX, 0, MPI_COMM_WORLD);
                    MPI_Reduce(&open_wall, &max_open, 1, MPI_DOUBLE,
                               MPI_MAX, 0, MPI_COMM_WORLD);
                    MPI_Reduce(&plan_wall, &max_plan, 1, MPI_DOUBLE,
                               MPI_MAX, 0, MPI_COMM_WORLD);
                    MPI_Reduce(&index_prepare, &max_index_prepare, 1,
                               MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                    MPI_Reduce(&stage_read, &max_stage_read, 1, MPI_DOUBLE,
                               MPI_MAX, 0, MPI_COMM_WORLD);
                    MPI_Reduce(&stage_prepare, &max_stage_prepare, 1, MPI_DOUBLE,
                               MPI_MAX, 0, MPI_COMM_WORLD);
                    MPI_Reduce(&local_capacity, &max_chunk_capacity, 1,
                               MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
                    MPI_Reduce(&chunks, &max_chunks, 1,
                               MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
                    const double reduce_wall = MPI_Wtime() - reduce_start;
                    const double local_simulated =
                        open_wall + plan_wall + seconds + reduce_wall;
                    double simulated_total = 0.0;
                    MPI_Reduce(&local_simulated, &simulated_total, 1,
                               MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                    if (rank == 0) {
                        printf("sdk_count ranks=%d records=%lld blocks=%zu "
                               "input_bytes=%zu core=%.6f mode=%s "
                               "open_wall=%.6f scan_plan_wall=%.6f "
                               "index_prepare=%.6f "
                               "stage_read=%.6f stage_prepare=%.6f "
                               "peak_chunk_bytes=%llu max_chunks=%llu "
                               "simulated_total=%.6f\n",
                               ranks, records, plan.blocks.size(),
                               backend->size(), core,
                               chunked ? "chunked" : "full",
                               max_open, max_plan, max_index_prepare, max_stage_read,
                               max_stage_prepare, max_chunk_capacity,
                               max_chunks,
                               chunked ? simulated_total : 0.0);
                    }
                    status = 0;
                }
            }
            }
        }
        input.Close();
        chunk_input.Close();
    }
#ifdef PLATFORM_SUNWAY
    athread_halt();
#endif
    MPI_Finalize();
    return status;
}
