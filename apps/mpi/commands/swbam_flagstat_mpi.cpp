#include "swbam_mpi.h"
#include "swbam/io.h"
#include "swbam/mpi_runtime.h"
#include "swbam/operators/flagstat.h"

#include <cstdio>

#include <mpi.h>

namespace {

using swbam::operators::FlagstatCounterId;
using swbam::operators::FlagstatCounts;
using swbam::operators::FlagstatMetrics;

void MpiFlagstatFormatPct(char *buf, size_t buf_size, long long value, long long total) {
    if (total > 0) {
        snprintf(buf, buf_size, "%.2f%%", 100.0 * (double)value / (double)total);
    } else {
        snprintf(buf, buf_size, "N/A");
    }
}

void MpiFlagstatPrintSimple(const FlagstatCounts &counts,
                            FlagstatCounterId id,
                            const char *label) {
    printf("%lld + %lld %s\n", counts.values[id][0], counts.values[id][1], label);
}

void MpiFlagstatPrintPct(const FlagstatCounts &counts,
                         FlagstatCounterId id,
                         FlagstatCounterId denom_id,
                         const char *label) {
    char pass_pct[32];
    char fail_pct[32];
    MpiFlagstatFormatPct(pass_pct, sizeof(pass_pct),
                         counts.values[id][0], counts.values[denom_id][0]);
    MpiFlagstatFormatPct(fail_pct, sizeof(fail_pct),
                         counts.values[id][1], counts.values[denom_id][1]);
    printf("%lld + %lld %s (%s : %s)\n",
           counts.values[id][0],
           counts.values[id][1],
           label,
           pass_pct,
           fail_pct);
}

void MpiFlagstatPrintSamtoolsStyle(const FlagstatCounts &counts) {
    MpiFlagstatPrintSimple(counts, swbam::operators::kFlagstatTotal,
                           "in total (QC-passed reads + QC-failed reads)");
    MpiFlagstatPrintSimple(counts, swbam::operators::kFlagstatPrimary, "primary");
    MpiFlagstatPrintSimple(counts, swbam::operators::kFlagstatSecondary, "secondary");
    MpiFlagstatPrintSimple(counts, swbam::operators::kFlagstatSupplementary, "supplementary");
    MpiFlagstatPrintSimple(counts, swbam::operators::kFlagstatDuplicates, "duplicates");
    MpiFlagstatPrintSimple(counts, swbam::operators::kFlagstatPrimaryDuplicates,
                           "primary duplicates");
    MpiFlagstatPrintPct(counts, swbam::operators::kFlagstatMapped,
                       swbam::operators::kFlagstatTotal, "mapped");
    MpiFlagstatPrintPct(counts, swbam::operators::kFlagstatPrimaryMapped,
                       swbam::operators::kFlagstatPrimary, "primary mapped");
    MpiFlagstatPrintSimple(counts, swbam::operators::kFlagstatPaired,
                           "paired in sequencing");
    MpiFlagstatPrintSimple(counts, swbam::operators::kFlagstatRead1, "read1");
    MpiFlagstatPrintSimple(counts, swbam::operators::kFlagstatRead2, "read2");
    MpiFlagstatPrintPct(counts, swbam::operators::kFlagstatProperlyPaired,
                       swbam::operators::kFlagstatPaired, "properly paired");
    MpiFlagstatPrintSimple(counts, swbam::operators::kFlagstatPairMapped,
                           "with itself and mate mapped");
    MpiFlagstatPrintPct(counts, swbam::operators::kFlagstatSingletons,
                       swbam::operators::kFlagstatPaired, "singletons");
    MpiFlagstatPrintSimple(counts, swbam::operators::kFlagstatDiffChr,
                           "with mate mapped to a different chr");
    MpiFlagstatPrintSimple(counts, swbam::operators::kFlagstatDiffChrMapq5,
                           "with mate mapped to a different chr (mapQ>=5)");
}

void MpiFlagstatReduceStats(const FlagstatMetrics &local_stats,
                            FlagstatMetrics *global_stats) {
    long long local_long[3] = {
        local_stats.input_blocks,
        local_stats.batch_count,
        local_stats.total_records
    };
    long long global_long[3] = {};
    double local_double[9] = {
        local_stats.read,
        local_stats.kernel,
        local_stats.decomp_alloc,
        local_stats.decomp_inflate,
        local_stats.decomp_crc,
        local_stats.decomp_parse,
        local_stats.decomp_other,
        local_stats.post_process,
        local_stats.total
    };
    double global_double[9] = {};

    MPI_Reduce(local_long, global_long, 3, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_double, global_double, 9, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (global_stats) {
        global_stats->input_blocks = global_long[0];
        global_stats->batch_count = global_long[1];
        global_stats->total_records = global_long[2];
        global_stats->read = global_double[0];
        global_stats->kernel = global_double[1];
        global_stats->decomp_alloc = global_double[2];
        global_stats->decomp_inflate = global_double[3];
        global_stats->decomp_crc = global_double[4];
        global_stats->decomp_parse = global_double[5];
        global_stats->decomp_other = global_double[6];
        global_stats->post_process = global_double[7];
        global_stats->total = global_double[8];
    }
}

} // namespace

int ProcessFlagstatMPI(CmdInfo *cmd_info) {
    double t_init = GetTime();

    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int exit_code = 1;
    int local_ok = 1;
    swbam::mpi::MpiBamInput input_handle;
    swbam::BamInputBackend *input = nullptr;
    swbam::mpi::MpiBamInputPlan input_plan;
    FlagstatCounts local_counts;
    FlagstatCounts global_counts;
    FlagstatMetrics local_stats;
    FlagstatMetrics global_stats;

    double init_cost = GetTime() - t_init;
    double init_cost_max = swbam::mpi::ReduceMaxCost(init_cost);
    if (rank == 0) {
        printf("111Complete the initialization cost %lf-----\n", init_cost_max);
    }

    if (!cmd_info || input_handle.Open(
            cmd_info->in_file_name_, cmd_info->io_backend_,
            cmd_info->io_memory_limit_) != 0) {
        fprintf(stderr, "[rank %d] ERROR: cannot open BAM input backend for %s\n",
                rank, cmd_info ? cmd_info->in_file_name_.c_str() : "(null)");
        local_ok = 0;
    }
    input = input_handle.backend();
    if (local_ok && (!input || input->format() != bam)) local_ok = 0;

    {
        const double preload_cost_max = swbam::mpi::ReduceMaxCost(
            input_handle.data_open_cost());
        const double header_cost_max = swbam::mpi::ReduceMaxCost(
            input_handle.header_open_cost());
        if (rank == 0 && local_ok) {
            printf("MPI BAM input backend=%s auto_memory_budget=%llu\n",
                   input_handle.selected_backend().c_str(),
                   (unsigned long long)input_handle.auto_memory_budget());
            printf("222Complete the input open cost %lf--\n", preload_cost_max);
            printf("333Complete the head cost %lf---\n", header_cost_max);
        }
    }
    if (!swbam::mpi::AllRanksOk(local_ok)) goto cleanup;

    {
        double body_total_t0 = GetTime();

        double stage41_t0 = GetTime();
        if (rank == 0) {
            printf("Enable MPI FLAGSTAT mode (%d MPE + %d CPEs)!!!\n",
                   comm_size, comm_size * 64);
        }
        if (swbam::mpi::PrepareMpiBamInputPlan(
                *input, &input_plan) != 0) {
            if (rank == 0) {
                fprintf(stderr, "ERROR: failed to prepare BGZF input plan for flagstat.\n");
            }
            local_ok = 0;
        }
        if (rank == 0 && local_ok) {
            printf("MPI BAM scan complete. data_blocks=%lld body_start=%lld header_end=%lld\n",
                   (long long)input_plan.blocks.size(),
                   (long long)input_plan.body_offset,
                   (long long)input->body_offset());
        }
        double stage41_cost = GetTime() - stage41_t0;
        double stage41_cost_max = swbam::mpi::ReduceMaxCost(stage41_cost);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.1 scan/split/broadcast/select cost %lf\n", stage41_cost_max);
        }
        if (!swbam::mpi::AllRanksOk(local_ok)) goto cleanup;

        double stage42_t0 = GetTime();
        double stage42_cost = GetTime() - stage42_t0;
        double stage42_cost_max = swbam::mpi::ReduceMaxCost(stage42_cost);
        if (rank == 0) {
            printf("Complete the 4.2 init reader cost %lf\n", stage42_cost_max);
        }

        double stage43_t0 = GetTime();
        if (swbam::operators::RunFlagstatPipeline(
                *input, input_plan.rank_spans(),
                input_plan.rank_block_count(),
                &local_counts, &local_stats) != 0) {
            local_ok = 0;
        }
        double stage43_cost = GetTime() - stage43_t0;
        int global_ok = swbam::mpi::AllRanksOk(local_ok);
        double stage43_cost_max = swbam::mpi::ReduceMaxCost(stage43_cost);
        if (rank == 0 && global_ok) {
            printf("Complete the 4.3 OptimizedFlagstatMPI cost %lf\n", stage43_cost_max);
        }
        if (!global_ok) goto cleanup;

        double stage44_t0 = GetTime();
        MPI_Reduce(&local_counts.values[0][0], &global_counts.values[0][0],
                   swbam::operators::kFlagstatCounterCount * 2,
                   MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MpiFlagstatReduceStats(local_stats, rank == 0 ? &global_stats : nullptr);
        double stage44_cost = GetTime() - stage44_t0;
        double stage44_cost_max = swbam::mpi::ReduceMaxCost(stage44_cost);
        if (rank == 0) {
            MpiFlagstatPrintSamtoolsStyle(global_counts);
            printf("OptimizedFlagstatMPI finished. ranks=%d in_blocks=%lld groups=%lld total_records=%lld\n",
                   comm_size,
                   global_stats.input_blocks,
                   global_stats.batch_count,
                   global_stats.total_records);
            printf("  read_sum=%.3f  decomp_sum=%.3f  merge_sum=%.3f  optimized_total_sum=%.3f\n",
                   global_stats.read,
                   global_stats.kernel,
                   global_stats.post_process,
                   global_stats.total);
            printf("  decomp_detail_sum alloc=%.3f  inflate=%.3f  crc=%.3f  parse=%.3f  other=%.3f\n",
                   global_stats.decomp_alloc,
                   global_stats.decomp_inflate,
                   global_stats.decomp_crc,
                   global_stats.decomp_parse,
                   global_stats.decomp_other);
            printf("Complete the 4.4 flagstat reduce/print cost %lf\n", stage44_cost_max);
        }

        double body_total_cost = GetTime() - body_total_t0;
        double body_total_cost_max = swbam::mpi::ReduceMaxCost(body_total_cost);
        if (rank == 0) {
            printf("444Complete the total body cost %lf\n", body_total_cost_max);
        }
    }

    exit_code = 0;

cleanup:
    {
        double close_t0 = GetTime();
        input_handle.Close();
        double close_cost = GetTime() - close_t0;
        double close_cost_max = swbam::mpi::ReduceMaxCost(close_cost);
        if (rank == 0) printf("666close the files cost %lf-----\n", close_cost_max);
    }

    return exit_code;
}
