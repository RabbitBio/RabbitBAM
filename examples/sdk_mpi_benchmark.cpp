#include "swbam/swbam.h"

#include <cstdio>
#include <cstdlib>
#include <string>

#include <mpi.h>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

namespace {

class RawCountPostProcessor : public swbam::cpe::RawBamBatchPostProcessor {
public:
    RawCountPostProcessor() : records_(0), mapped_(0) {}

    int PostProcessRawBatch(
            const swbam::cpe::RawBamRecordView *records,
            size_t count) {
        if (!records && count != 0) return -1;
        for (size_t i = 0; i < count; ++i) {
            if (!records[i].encoded || records[i].encoded_size < 20) {
                return -1;
            }
            const unsigned char *data = records[i].encoded + 16;
            const unsigned int flag_nc =
                (unsigned int)data[0] |
                ((unsigned int)data[1] << 8) |
                ((unsigned int)data[2] << 16) |
                ((unsigned int)data[3] << 24);
            ++records_;
            if (((flag_nc >> 16) & BAM_FUNMAP) == 0) ++mapped_;
        }
        return 0;
    }

    long long records() const { return records_; }
    long long mapped() const { return mapped_; }

private:
    long long records_;
    long long mapped_;
};

class Bam1CountPostProcessor : public swbam::cpe::Bam1BatchPostProcessor {
public:
    Bam1CountPostProcessor() : records_(0), mapped_(0) {}

    int PostProcessBam1Batch(
            const bam1_t *const *records, size_t count) {
        if (!records && count != 0) return -1;
        for (size_t i = 0; i < count; ++i) {
            if (!records[i]) return -1;
            ++records_;
            if ((records[i]->core.flag & BAM_FUNMAP) == 0) ++mapped_;
        }
        return 0;
    }

    long long records() const { return records_; }
    long long mapped() const { return mapped_; }

private:
    long long records_;
    long long mapped_;
};

double ReduceMax(double value, int rank) {
    double result = 0.0;
    MPI_Reduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, 0,
               MPI_COMM_WORLD);
    return rank == 0 ? result : 0.0;
}

long long ReduceSum(long long value, int rank) {
    long long result = 0;
    MPI_Reduce(&value, &result, 1, MPI_LONG_LONG, MPI_SUM, 0,
               MPI_COMM_WORLD);
    return rank == 0 ? result : 0;
}

unsigned long long ReduceSumBytes(unsigned long long value, int rank) {
    unsigned long long result = 0;
    MPI_Reduce(&value, &result, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0,
               MPI_COMM_WORLD);
    return rank == 0 ? result : 0;
}

int AllRanksOk(int local_ok) {
    return swbam::mpi::AllRanksOk(local_ok, MPI_COMM_WORLD);
}

size_t OutputReserveBytes(size_t input_size, int comm_size) {
    const size_t reserve_slack = 128u * 1024u * 1024u;
    const size_t rank_input = input_size / (size_t)comm_size;
    if (rank_input > (static_cast<size_t>(-1) - reserve_slack) / 3u * 2u) {
        return 0;
    }
    return rank_input + rank_input / 2u + reserve_slack;
}

int RunRawRead(const swbam::BamInputBackend &input,
               const swbam::mpi::MpiBamInputPlan &plan,
               int rank, int comm_size) {
    RawCountPostProcessor post_processor;
    swbam::cpe::CpeReadPipelineTiming timing;
    swbam::cpe::ComposableDecodeMetrics decode;
    swbam::cpe::ComposableRawBamMetrics raw;

    MPI_Barrier(MPI_COMM_WORLD);
    const double wall_t0 = MPI_Wtime();
    int local_ok = swbam::cpe::RunComposableRawBamPipeline(
        input, plan.rank_spans(), plan.rank_block_count(),
        &post_processor, &timing, &decode, &raw) == 0;
    const double local_core = MPI_Wtime() - wall_t0;
    if (!AllRanksOk(local_ok)) return -1;

    const long long records = ReduceSum(post_processor.records(), rank);
    const long long mapped = ReduceSum(post_processor.mapped(), rank);
    const double core = ReduceMax(local_core, rank);
    const double read = ReduceMax(timing.read, rank);
    const double cpe = ReduceMax(timing.kernel, rank);
    const double post = ReduceMax(timing.post_process, rank);
    if (rank == 0) {
        printf("sdk_mpi mode=raw_read ranks=%d records=%lld mapped=%lld "
               "core=%.6f read=%.6f cpe=%.6f post=%.6f\n",
               comm_size, records, mapped, core, read, cpe, post);
    }
    return 0;
}

int RunBam1Read(const swbam::BamInputBackend &input,
                const swbam::mpi::MpiBamInputPlan &plan,
                int rank, int comm_size) {
    Bam1CountPostProcessor post_processor;
    swbam::cpe::CpeReadPipelineTiming timing;
    swbam::cpe::ComposableDecodeMetrics decode;
    swbam::cpe::ComposableBam1Metrics bam1;

    MPI_Barrier(MPI_COMM_WORLD);
    const double wall_t0 = MPI_Wtime();
    int local_ok = swbam::cpe::RunComposableBam1Pipeline(
        input, plan.rank_spans(), plan.rank_block_count(),
        &post_processor, &timing, &decode, &bam1) == 0;
    const double local_core = MPI_Wtime() - wall_t0;
    if (!AllRanksOk(local_ok)) return -1;

    const long long records = ReduceSum(post_processor.records(), rank);
    const long long mapped = ReduceSum(post_processor.mapped(), rank);
    const double core = ReduceMax(local_core, rank);
    const double read = ReduceMax(timing.read, rank);
    const double cpe = ReduceMax(timing.kernel, rank);
    const double pipeline_post = ReduceMax(timing.post_process, rank);
    const double materialize = ReduceMax(bam1.materialize, rank);
    const double post = ReduceMax(bam1.post_process, rank);
    if (rank == 0) {
        printf("sdk_mpi mode=bam1_read ranks=%d records=%lld mapped=%lld "
               "core=%.6f read=%.6f cpe=%.6f pipeline_post=%.6f "
               "materialize=%.6f post=%.6f\n",
               comm_size, records, mapped, core, read, cpe,
               pipeline_post, materialize, post);
    }
    return 0;
}

int RunRawReadWrite(const swbam::BamInputBackend &input,
                    const swbam::mpi::MpiBamInputPlan &plan,
                    size_t output_reserve,
                    int rank, int comm_size) {
    swbam::MemoryBamOutput output;
    swbam::cpe::RawBamWriter writer;
    swbam::cpe::CpeReadPipelineTiming timing;
    swbam::cpe::ComposableDecodeMetrics decode;
    swbam::cpe::ComposableRawBamMetrics raw;
    int local_ok = output_reserve > 0 &&
                   output.Reserve(output_reserve) == 0;
    if (!AllRanksOk(local_ok)) return -1;

    MPI_Barrier(MPI_COMM_WORLD);
    const double wall_t0 = MPI_Wtime();
    if (writer.InitializeBam(&output, input.header(), 1) != 0 ||
        swbam::cpe::RunComposableRawBamPipeline(
            input, plan.rank_spans(), plan.rank_block_count(),
            &writer, &timing, &decode, &raw) != 0 ||
        writer.Finish() != 0) {
        local_ok = 0;
    }
    const double local_core = MPI_Wtime() - wall_t0;
    if (!AllRanksOk(local_ok)) return -1;

    const swbam::cpe::RawBamWriterMetrics &writer_metrics =
        writer.metrics();
    const long long records = ReduceSum(writer_metrics.records, rank);
    const unsigned long long output_bytes = ReduceSumBytes(
        (unsigned long long)output.bytes_written(), rank);
    const double core = ReduceMax(local_core, rank);
    const double pipeline = ReduceMax(timing.total, rank);
    const double read = ReduceMax(timing.read, rank);
    const double cpe = ReduceMax(timing.kernel, rank);
    const double pipeline_post = ReduceMax(timing.post_process, rank);
    const double pack = ReduceMax(writer_metrics.pack, rank);
    const double compress = ReduceMax(writer_metrics.kernel, rank);
    const double writer_total = ReduceMax(writer_metrics.total, rank);
    if (rank == 0) {
        printf("sdk_mpi mode=raw_rw ranks=%d records=%lld "
               "output_bytes=%llu core=%.6f pipeline=%.6f read=%.6f "
               "cpe=%.6f pipeline_post=%.6f pack=%.6f compress=%.6f "
               "writer=%.6f\n",
               comm_size, records, output_bytes, core, pipeline, read,
               cpe, pipeline_post, pack, compress, writer_total);
    }
    return 0;
}

int RunBam1ReadWrite(const swbam::BamInputBackend &input,
                     const swbam::mpi::MpiBamInputPlan &plan,
                     size_t output_reserve,
                     int rank, int comm_size) {
    swbam::MemoryBamOutput output;
    swbam::cpe::Bam1Writer writer;
    swbam::cpe::CpeReadPipelineTiming timing;
    swbam::cpe::ComposableDecodeMetrics decode;
    swbam::cpe::ComposableBam1Metrics bam1;
    int local_ok = output_reserve > 0 &&
                   output.Reserve(output_reserve) == 0;
    if (!AllRanksOk(local_ok)) return -1;

    MPI_Barrier(MPI_COMM_WORLD);
    const double wall_t0 = MPI_Wtime();
    if (writer.InitializeBam(&output, input.header(), 1) != 0 ||
        swbam::cpe::RunComposableBam1Pipeline(
            input, plan.rank_spans(), plan.rank_block_count(),
            &writer, &timing, &decode, &bam1) != 0 ||
        writer.Finish() != 0) {
        local_ok = 0;
    }
    const double local_core = MPI_Wtime() - wall_t0;
    if (!AllRanksOk(local_ok)) return -1;

    const swbam::cpe::Bam1WriterMetrics &writer_metrics =
        writer.metrics();
    const swbam::cpe::RawBamWriterMetrics &raw_writer =
        writer.raw_metrics();
    const long long records = ReduceSum(writer_metrics.records, rank);
    const unsigned long long output_bytes = ReduceSumBytes(
        (unsigned long long)output.bytes_written(), rank);
    const double core = ReduceMax(local_core, rank);
    const double pipeline = ReduceMax(timing.total, rank);
    const double read = ReduceMax(timing.read, rank);
    const double cpe = ReduceMax(timing.kernel, rank);
    const double pipeline_post = ReduceMax(timing.post_process, rank);
    const double materialize = ReduceMax(bam1.materialize, rank);
    const double encode = ReduceMax(writer_metrics.encode, rank);
    const double pack = ReduceMax(raw_writer.pack, rank);
    const double compress = ReduceMax(raw_writer.kernel, rank);
    if (rank == 0) {
        printf("sdk_mpi mode=bam1_rw ranks=%d records=%lld "
               "output_bytes=%llu core=%.6f pipeline=%.6f read=%.6f "
               "cpe=%.6f pipeline_post=%.6f materialize=%.6f "
               "encode=%.6f pack=%.6f compress=%.6f\n",
               comm_size, records, output_bytes, core, pipeline, read,
               cpe, pipeline_post, materialize, encode, pack, compress);
    }
    return 0;
}

} // namespace

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    if (argc != 2) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s input.bam\n", argv[0]);
        }
        MPI_Finalize();
        return 2;
    }

#ifdef PLATFORM_SUNWAY
    athread_init();
#endif

    int exit_code = 1;
    swbam::mpi::MpiBamInput mpi_input;
    swbam::mpi::MpiBamInputPlan plan;
    const double open_t0 = MPI_Wtime();
    int local_ok = mpi_input.Open(
        argv[1], "memory", "", MPI_COMM_WORLD) == 0;
    const double local_open = MPI_Wtime() - open_t0;
    if (!AllRanksOk(local_ok)) goto cleanup;

    {
        const double plan_t0 = MPI_Wtime();
        local_ok = swbam::mpi::PrepareMpiBamInputPlan(
            *mpi_input.backend(), &plan, 0, MPI_COMM_WORLD) == 0;
        const double local_plan = MPI_Wtime() - plan_t0;
        if (!AllRanksOk(local_ok)) goto cleanup;

        const double open_max = ReduceMax(local_open, rank);
        const double plan_max = ReduceMax(local_plan, rank);
        if (rank == 0) {
            printf("sdk_mpi_setup ranks=%d blocks=%zu input_bytes=%zu "
                   "open=%.6f plan=%.6f\n",
                   comm_size, plan.blocks.size(),
                   mpi_input.backend()->size(), open_max, plan_max);
        }
    }

    {
        const size_t output_reserve = OutputReserveBytes(
            mpi_input.backend()->size(), comm_size);
        if (RunRawRead(*mpi_input.backend(), plan, rank, comm_size) != 0 ||
            RunRawReadWrite(*mpi_input.backend(), plan, output_reserve,
                            rank, comm_size) != 0 ||
            RunBam1Read(*mpi_input.backend(), plan, rank, comm_size) != 0 ||
            RunBam1ReadWrite(*mpi_input.backend(), plan, output_reserve,
                             rank, comm_size) != 0) {
            goto cleanup;
        }
    }

    exit_code = 0;

cleanup:
    mpi_input.Close();
#ifdef PLATFORM_SUNWAY
    athread_halt();
#endif
    MPI_Finalize();
    return exit_code;
}
