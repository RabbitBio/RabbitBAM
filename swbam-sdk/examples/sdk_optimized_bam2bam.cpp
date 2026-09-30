#include "swbam/swbam.h"
#include "sdk_chunked_bam.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
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
    const size_t chunk_limit = argc == 3
        ? sdk_bench::ChunkLimit(argv[2]) : 0;
    const bool chunked = chunk_limit != 0;
    const char *hash_env = std::getenv("SWBAM_SDK_BODY_HASH");
    const bool hash_body = hash_env && std::strcmp(hash_env, "1") == 0;
    if (argc != 2 && !chunked) {
        if (rank == 0) fprintf(stderr,
            "Usage: %s input.bam [--chunked-memory[=MiB]]\n", argv[0]);
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
        sdk_bench::ChunkedBamInput chunk_input;
        swbam::mpi::MpiBamInputPlan plan;
        swbam::AdaptiveRankBodySink body;
        sdk_bench::CircularBodySink circular_body;
        mark("open.begin");
        const double open_start = MPI_Wtime();
        int ok = chunked ? chunk_input.Open(argv[1]) == 0
                         : input.Open(argv[1], "memory", "", MPI_COMM_WORLD) == 0;
        const double open_wall = MPI_Wtime() - open_start;
        swbam::BamInputBackend *backend = chunked
            ? static_cast<swbam::BamInputBackend *>(&chunk_input)
            : input.backend();
        swbam::RankBodySink *sink = chunked && !hash_body
            ? static_cast<swbam::RankBodySink *>(&circular_body)
            : static_cast<swbam::RankBodySink *>(&body);
        mark("open.done", ok);
        if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
            double index_prepare = 0.0;
            if (chunked && rank == 0) {
                const double start = MPI_Wtime();
                ok = chunk_input.PrepareIndex() == 0;
                index_prepare = MPI_Wtime() - start;
            }
            if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
            mark("scan_plan.begin");
            const double plan_start = MPI_Wtime();
            ok = swbam::mpi::PrepareMpiBamInputPlan(
                *backend, &plan, 0, MPI_COMM_WORLD) == 0;
            const double plan_wall = MPI_Wtime() - plan_start;
            mark("scan_plan.done", ok, plan.rank_block_count());
            if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
                size_t reserve = 0;
                if (!chunked || hash_body) {
                    for (size_t i = plan.rank_begin; i < plan.rank_end; ++i)
                        reserve += plan.blocks[i].compressed_size;
                    // Leave room for level-1 output expansion without a large realloc.
                    if (reserve > std::numeric_limits<size_t>::max() - reserve / 4)
                        ok = 0;
                    else
                        reserve += reserve / 4;
                }
                mark("body_reserve.begin", 1, reserve);
                const double body_setup_start = MPI_Wtime();
                if (ok) ok = chunked && !hash_body
                    ? circular_body.Open(size_t(256) << 20) == 0
                    : body.Open("memory", 0, reserve,
                                "sdk-bam2bam-body", "") == 0;
                const double body_setup = MPI_Wtime() - body_setup_start;
                mark("body_reserve.done", ok, reserve);
                if (!ok) fprintf(stderr, "[rank %d] memory body allocation failed.\n", rank);
                if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
                    const swbam::BamFilterOptions filter =
                        swbam::DefaultBamFilterOptions();
                    mark("core_barrier.begin");
                    if (!chunked) MPI_Barrier(MPI_COMM_WORLD);
                    mark("pipeline.begin");
                    double seconds = 0.0, stage_read = 0.0;
                    double stage_prepare = 0.0;
                    long long local_records = 0;
                    unsigned long long chunks = 0;
                    const swbam::BgzfBlockSpan *spans = plan.rank_spans();
                    const size_t count = plan.rank_block_count();
                    if (count == 0) {
                        swbam::operators::BamTransformMetrics metrics = {};
                        const double start = MPI_Wtime();
                        ok = swbam::operators::RunBamTransformPipeline(
                            *backend, nullptr, 0, sink, filter, 1,
                            &metrics) == 0;
                        seconds += MPI_Wtime() - start;
                    }
                    for (size_t begin = 0; begin < count && ok;) {
                        const size_t end = chunked
                            ? sdk_bench::NextChunkEnd(spans, count, begin,
                                                       chunk_limit)
                            : count;
                        if (chunked && chunk_input.Stage(
                                spans + begin, end - begin,
                                chunk_limit,
                                &stage_read, &stage_prepare) != 0) {
                            ok = 0;
                            break;
                        }
                        swbam::operators::BamTransformMetrics metrics = {};
                        const double start = MPI_Wtime();
                        ok = swbam::operators::RunBamTransformPipeline(
                            *backend, spans + begin, end - begin,
                            sink, filter, 1, &metrics) == 0 &&
                             metrics.kept_records == metrics.total_records;
                        seconds += MPI_Wtime() - start;
                        local_records += metrics.total_records;
                        begin = end;
                        ++chunks;
                    }
                    mark("pipeline.done", ok, sink->size());
                    unsigned long long body_hash = 0;
                    if (ok && hash_body) {
                        const MemWriter *writer = body.memory_writer();
                        if (!writer || writer->size != sink->size()) {
                            ok = 0;
                        } else {
                            uint64_t hash = UINT64_C(14695981039346656037);
                            const unsigned char *data =
                                reinterpret_cast<const unsigned char *>(writer->data);
                            for (size_t i = 0; i < writer->size; ++i) {
                                hash ^= data[i];
                                hash *= UINT64_C(1099511628211);
                            }
                            body_hash = hash;
                        }
                    }
                    mark("result_collectives.begin", ok);
                    if (swbam::mpi::AllRanksOk(ok, MPI_COMM_WORLD)) {
                        long long records = 0;
                        unsigned long long body_bytes = 0;
                        unsigned long long global_hash = 0;
                        double core = 0.0;
                        double max_open = 0.0, max_plan = 0.0;
                        double max_index_prepare = 0.0;
                        double max_stage_read = 0.0, max_stage_prepare = 0.0;
                        double max_body_setup = 0.0;
                        unsigned long long max_chunk_capacity = 0;
                        unsigned long long max_chunks = 0;
                        const unsigned long long local_capacity =
                            chunk_input.peak_capacity();
                        const unsigned long long local_bytes = sink->size();
                        const double reduce_start = MPI_Wtime();
                        MPI_Reduce(&local_records, &records, 1,
                                   MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&local_bytes, &body_bytes, 1,
                                   MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0,
                                   MPI_COMM_WORLD);
                        MPI_Reduce(&body_hash, &global_hash, 1,
                                   MPI_UNSIGNED_LONG_LONG, MPI_BXOR, 0,
                                   MPI_COMM_WORLD);
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
                        MPI_Reduce(&body_setup, &max_body_setup, 1, MPI_DOUBLE,
                                   MPI_MAX, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&local_capacity, &max_chunk_capacity, 1,
                                   MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0,
                                   MPI_COMM_WORLD);
                        MPI_Reduce(&chunks, &max_chunks, 1,
                                   MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0,
                                   MPI_COMM_WORLD);
                        const double reduce_wall = MPI_Wtime() - reduce_start;
                        const double local_simulated =
                            open_wall + plan_wall + body_setup +
                            seconds + reduce_wall;
                        double simulated_total = 0.0;
                        MPI_Reduce(&local_simulated, &simulated_total, 1,
                                   MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                        if (rank == 0) {
                            printf("sdk_bam2bam ranks=%d records=%lld blocks=%zu "
                                   "input_bytes=%zu body_bytes=%llu core=%.6f "
                                   "mode=%s open_wall=%.6f scan_plan_wall=%.6f "
                                   "index_prepare=%.6f "
                                   "stage_read=%.6f stage_prepare=%.6f "
                                   "body_setup=%.6f peak_chunk_bytes=%llu "
                                   "max_chunks=%llu "
                                   "simulated_total=%.6f body_hash=%016llx\n",
                                   ranks, records, plan.blocks.size(),
                                   backend->size(), body_bytes, core,
                                   chunked ? "chunked" : "full",
                                   max_open, max_plan, max_index_prepare,
                                   max_stage_read,
                                   max_stage_prepare, max_body_setup,
                                   max_chunk_capacity, max_chunks,
                                   chunked ? simulated_total : 0.0,
                                   global_hash);
                        }
                        status = 0;
                        mark("result_collectives.done");
                    }
                }
            }
            }
        }
        mark("body_close.begin");
        body.Close();
        mark("input_close.begin");
        input.Close();
        chunk_input.Close();
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
