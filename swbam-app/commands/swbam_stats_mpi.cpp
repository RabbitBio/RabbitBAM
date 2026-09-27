#include "swbam_mpi.h"
#include "swbam/io.h"
#include "swbam/mpi_runtime.h"
#include "swbam/operators/stats_basic.h"

#include <cmath>
#include <cstdio>
#include <vector>

#include <mpi.h>

namespace {

using MpiBasicStatsCounts = swbam::operators::StatsBasicCounts;
using MpiStatsSortState = swbam::operators::StatsSortState;
using MpiStatsPerf = swbam::operators::StatsBasicMetrics;

const int kStatsInsertBins = swbam::operators::kStatsInsertBins;
const double kStatsInsertMainBulk = 0.99;

enum StatsOutputValueId {
    STATS_NREADS_1ST = swbam::operators::kStatsFirstFragments,
    STATS_NREADS_2ND = swbam::operators::kStatsLastFragments,
    STATS_NREADS_OTHER = swbam::operators::kStatsOtherFragments,
    STATS_NREADS_FILTERED = swbam::operators::kStatsFilteredSequences,
    STATS_NREADS_DUP = swbam::operators::kStatsDuplicateReads,
    STATS_NREADS_UNMAPPED = swbam::operators::kStatsUnmappedReads,
    STATS_NREADS_SINGLE_MAPPED = swbam::operators::kStatsSingleMappedReads,
    STATS_NREADS_PAIRED_AND_MAPPED =
        swbam::operators::kStatsPairedAndMappedReads,
    STATS_NREADS_PROPERLY_PAIRED =
        swbam::operators::kStatsProperlyPairedReads,
    STATS_NREADS_PAIRED_TECH = swbam::operators::kStatsPairedReads,
    STATS_NREADS_ANOMALOUS = swbam::operators::kStatsAnomalousReads,
    STATS_NREADS_MQ0 = swbam::operators::kStatsMq0Reads,
    STATS_NREADS_QCFAILED = swbam::operators::kStatsQcFailedReads,
    STATS_NREADS_SECONDARY = swbam::operators::kStatsSecondaryReads,
    STATS_NREADS_SUPPLEMENTARY =
        swbam::operators::kStatsSupplementaryReads,
    STATS_TOTAL_LEN = swbam::operators::kStatsTotalLength,
    STATS_TOTAL_LEN_1ST = swbam::operators::kStatsFirstTotalLength,
    STATS_TOTAL_LEN_2ND = swbam::operators::kStatsLastTotalLength,
    STATS_TOTAL_LEN_DUP = swbam::operators::kStatsDuplicateTotalLength,
    STATS_NBASES_MAPPED = swbam::operators::kStatsMappedBases,
    STATS_NBASES_MAPPED_CIGAR = swbam::operators::kStatsMappedCigarBases,
    STATS_NBASES_TRIMMED = swbam::operators::kStatsTrimmedBases,
    STATS_NMISMATCHES = swbam::operators::kStatsMismatches,
    STATS_MAX_LEN = swbam::operators::kStatsMaxLength,
    STATS_MAX_LEN_1ST = swbam::operators::kStatsFirstMaxLength,
    STATS_MAX_LEN_2ND = swbam::operators::kStatsLastMaxLength,
    STATS_SUM_QUAL = swbam::operators::kStatsQualitySum,
    STATS_LONG_COUNT = swbam::operators::kStatsValueCount
};

enum StatsOrientationOutputId {
    ORIENT_DIAG_REF_IN = swbam::operators::kStatsOrientRefIn,
    ORIENT_DIAG_REF_OUT = swbam::operators::kStatsOrientRefOut,
    ORIENT_DIAG_REF_OTHER = swbam::operators::kStatsOrientRefOther,
    ORIENT_DIAG_DOC_IN = swbam::operators::kStatsOrientDocIn,
    ORIENT_DIAG_DOC_OUT = swbam::operators::kStatsOrientDocOut,
    ORIENT_DIAG_DOC_OTHER = swbam::operators::kStatsOrientDocOther,
    ORIENT_DIAG_REF_IN_LT_READ =
        swbam::operators::kStatsOrientRefInLtRead,
    ORIENT_DIAG_REF_IN_LT_2READ =
        swbam::operators::kStatsOrientRefInLt2Read,
    ORIENT_DIAG_REF_IN_GE_2READ =
        swbam::operators::kStatsOrientRefInGe2Read,
    ORIENT_DIAG_REF_OUT_LT_READ =
        swbam::operators::kStatsOrientRefOutLtRead,
    ORIENT_DIAG_REF_OUT_LT_2READ =
        swbam::operators::kStatsOrientRefOutLt2Read,
    ORIENT_DIAG_REF_OUT_GE_2READ =
        swbam::operators::kStatsOrientRefOutGe2Read,
    ORIENT_DIAG_POS_NEG_STRAND_NEG =
        swbam::operators::kStatsOrientPosNegStrandNeg,
    ORIENT_DIAG_POS_NEG_STRAND_POS =
        swbam::operators::kStatsOrientPosNegStrandPos,
    ORIENT_DIAG_POS_ZERO_STRAND_NEG =
        swbam::operators::kStatsOrientPosZeroStrandNeg,
    ORIENT_DIAG_POS_ZERO_STRAND_POS =
        swbam::operators::kStatsOrientPosZeroStrandPos,
    ORIENT_DIAG_POS_POS_STRAND_NEG =
        swbam::operators::kStatsOrientPosPosStrandNeg,
    ORIENT_DIAG_POS_POS_STRAND_POS =
        swbam::operators::kStatsOrientPosPosStrandPos,
    ORIENT_DIAG_READ1_LEFT = swbam::operators::kStatsOrientRead1Left,
    ORIENT_DIAG_READ1_RIGHT = swbam::operators::kStatsOrientRead1Right,
    ORIENT_DIAG_READ1_SAME_POS =
        swbam::operators::kStatsOrientRead1SamePos,
    ORIENT_DIAG_READ2_LEFT = swbam::operators::kStatsOrientRead2Left,
    ORIENT_DIAG_READ2_RIGHT = swbam::operators::kStatsOrientRead2Right,
    ORIENT_DIAG_READ2_SAME_POS =
        swbam::operators::kStatsOrientRead2SamePos,
    ORIENT_DIAG_ISIZE_NEG = swbam::operators::kStatsOrientIsizeNeg,
    ORIENT_DIAG_ISIZE_ZERO = swbam::operators::kStatsOrientIsizeZero,
    ORIENT_DIAG_ISIZE_POS = swbam::operators::kStatsOrientIsizePos,
    ORIENT_DIAG_COUNT =
        swbam::operators::kStatsOrientationDiagnosticCount
};

void MpiStatsReducePerf(const MpiStatsPerf &local_stats, MpiStatsPerf *global_stats) {
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

int MpiStatsIsSortedGlobal(int comm_size, const std::vector<long long> &sort_values) {
    int sorted = 1;
    int have_prev = 0;
    long long prev_tid = 0;
    long long prev_pos = 0;
    for (int r = 0; r < comm_size; ++r) {
        const long long *state = &sort_values[(size_t)r * 6];
        const long long has_coord = state[0];
        const long long local_sorted = state[1];
        const long long first_tid = state[2];
        const long long first_pos = state[3];
        const long long last_tid = state[4];
        const long long last_pos = state[5];
        if (!local_sorted) sorted = 0;
        if (!has_coord) continue;
        if (have_prev &&
            (first_tid < prev_tid ||
             (first_tid == prev_tid && first_pos < prev_pos))) {
            sorted = 0;
        }
        have_prev = 1;
        prev_tid = last_tid;
        prev_pos = last_pos;
    }
    return sorted;
}

void MpiStatsComputeInsertSummary(const MpiBasicStatsCounts &counts,
                                  double *avg_isize,
                                  double *sd_isize,
                                  long long *n_inward,
                                  long long *n_outward,
                                  long long *n_other) {
    long long total_pairs = 0;
    double bulk = 0.0;
    double avg = 0.0;
    int ibulk = 0;
    long long full_inward = 0;
    long long full_outward = 0;
    long long full_other = 0;

    for (int isize = 0; isize < kStatsInsertBins; ++isize) {
        const long long inward = counts.isize_inward[isize] / 2;
        const long long outward = counts.isize_outward[isize] / 2;
        const long long other = counts.isize_other[isize] / 2;
        full_inward += inward;
        full_outward += outward;
        full_other += other;
        total_pairs += inward + outward + other;
    }

    long long main_pairs = total_pairs;
    for (int isize = 0; isize < kStatsInsertBins; ++isize) {
        const long long num = counts.isize_inward[isize] / 2 +
                              counts.isize_outward[isize] / 2 +
                              counts.isize_other[isize] / 2;
        if (num > 0) ibulk = isize + 1;
        bulk += (double)num;
        avg += (double)isize * (double)num;
        if (total_pairs > 0 && bulk / (double)total_pairs > kStatsInsertMainBulk) {
            ibulk = isize + 1;
            main_pairs = (long long)bulk;
            break;
        }
    }

    avg /= main_pairs ? (double)main_pairs : 1.0;
    double sd = 0.0;
    for (int isize = 1; isize < ibulk; ++isize) {
        const long long num = counts.isize_inward[isize] / 2 +
                              counts.isize_outward[isize] / 2 +
                              counts.isize_other[isize] / 2;
        const double diff = (double)isize - avg;
        sd += (double)num * diff * diff / (main_pairs ? (double)main_pairs : 1.0);
    }

    *avg_isize = avg;
    *sd_isize = sqrt(sd);
    *n_inward = full_inward;
    *n_outward = full_outward;
    *n_other = full_other;
}

void MpiStatsPrintBasic(const MpiBasicStatsCounts &counts, int is_sorted) {
    const long long nreads_1st = counts.values[STATS_NREADS_1ST];
    const long long nreads_2nd = counts.values[STATS_NREADS_2ND];
    const long long nreads_other = counts.values[STATS_NREADS_OTHER];
    const long long sequences = nreads_1st + nreads_2nd + nreads_other;
    const long long raw_total = counts.values[STATS_NREADS_FILTERED] + sequences;
    const long long reads_mapped = counts.values[STATS_NREADS_PAIRED_AND_MAPPED] +
                                   counts.values[STATS_NREADS_SINGLE_MAPPED];
    const double error_rate = counts.values[STATS_NBASES_MAPPED_CIGAR] > 0
        ? (double)counts.values[STATS_NMISMATCHES] / (double)counts.values[STATS_NBASES_MAPPED_CIGAR]
        : 0.0;
    const double average_length = sequences > 0
        ? (double)counts.values[STATS_TOTAL_LEN] / (double)sequences
        : 0.0;
    const double average_first_length = nreads_1st > 0
        ? (double)counts.values[STATS_TOTAL_LEN_1ST] / (double)nreads_1st
        : 0.0;
    const double average_last_length = nreads_2nd > 0
        ? (double)counts.values[STATS_TOTAL_LEN_2ND] / (double)nreads_2nd
        : 0.0;
    const double average_quality = counts.values[STATS_TOTAL_LEN] > 0
        ? (double)counts.values[STATS_SUM_QUAL] / (double)counts.values[STATS_TOTAL_LEN]
        : 0.0;
    const double properly_paired_pct = sequences > 0
        ? 100.0 * (double)counts.values[STATS_NREADS_PROPERLY_PAIRED] / (double)sequences
        : 0.0;

    double avg_isize = 0.0;
    double sd_isize = 0.0;
    long long n_inward = 0;
    long long n_outward = 0;
    long long n_other = 0;
    MpiStatsComputeInsertSummary(counts, &avg_isize, &sd_isize,
                                 &n_inward, &n_outward, &n_other);

    printf("# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.\n");
    printf("SN\traw total sequences:\t%lld\t# excluding supplementary and secondary reads\n", raw_total);
    printf("SN\tfiltered sequences:\t%lld\n", counts.values[STATS_NREADS_FILTERED]);
    printf("SN\tsequences:\t%lld\n", sequences);
    printf("SN\tis sorted:\t%d\t# %s by coordinate\n",
           is_sorted ? 1 : 0, is_sorted ? "sorted" : "not sorted");
    printf("SN\t1st fragments:\t%lld\n", nreads_1st);
    printf("SN\tlast fragments:\t%lld\n", nreads_2nd);
    printf("SN\treads mapped:\t%lld\n", reads_mapped);
    printf("SN\treads mapped and paired:\t%lld\t# paired-end technology bit set + both mates mapped\n",
           counts.values[STATS_NREADS_PAIRED_AND_MAPPED]);
    printf("SN\treads unmapped:\t%lld\n", counts.values[STATS_NREADS_UNMAPPED]);
    printf("SN\treads properly paired:\t%lld\t# proper-pair bit set\n",
           counts.values[STATS_NREADS_PROPERLY_PAIRED]);
    printf("SN\treads paired:\t%lld\t# paired-end technology bit set\n",
           counts.values[STATS_NREADS_PAIRED_TECH]);
    printf("SN\treads duplicated:\t%lld\t# PCR or optical duplicate bit set\n",
           counts.values[STATS_NREADS_DUP]);
    printf("SN\treads MQ0:\t%lld\t# mapped and MQ=0\n", counts.values[STATS_NREADS_MQ0]);
    printf("SN\treads QC failed:\t%lld\n", counts.values[STATS_NREADS_QCFAILED]);
    printf("SN\tnon-primary alignments:\t%lld\n", counts.values[STATS_NREADS_SECONDARY]);
    printf("SN\tsupplementary alignments:\t%lld\n", counts.values[STATS_NREADS_SUPPLEMENTARY]);
    printf("SN\ttotal length:\t%lld\t# ignores clipping\n", counts.values[STATS_TOTAL_LEN]);
    printf("SN\ttotal first fragment length:\t%lld\t# ignores clipping\n",
           counts.values[STATS_TOTAL_LEN_1ST]);
    printf("SN\ttotal last fragment length:\t%lld\t# ignores clipping\n",
           counts.values[STATS_TOTAL_LEN_2ND]);
    printf("SN\tbases mapped:\t%lld\t# ignores clipping\n", counts.values[STATS_NBASES_MAPPED]);
    printf("SN\tbases mapped (cigar):\t%lld\t# more accurate\n",
           counts.values[STATS_NBASES_MAPPED_CIGAR]);
    printf("SN\tbases trimmed:\t%lld\n", counts.values[STATS_NBASES_TRIMMED]);
    printf("SN\tbases duplicated:\t%lld\n", counts.values[STATS_TOTAL_LEN_DUP]);
    printf("SN\tmismatches:\t%lld\t# from NM fields\n", counts.values[STATS_NMISMATCHES]);
    printf("SN\terror rate:\t%e\t# mismatches / bases mapped (cigar)\n", error_rate);
    printf("SN\taverage length:\t%.0f\n", average_length);
    printf("SN\taverage first fragment length:\t%.0f\n", average_first_length);
    printf("SN\taverage last fragment length:\t%.0f\n", average_last_length);
    printf("SN\tmaximum length:\t%lld\n", counts.values[STATS_MAX_LEN]);
    printf("SN\tmaximum first fragment length:\t%lld\n", counts.values[STATS_MAX_LEN_1ST]);
    printf("SN\tmaximum last fragment length:\t%lld\n", counts.values[STATS_MAX_LEN_2ND]);
    printf("SN\taverage quality:\t%.1f\n", average_quality);
    printf("SN\tinsert size average:\t%.1f\n", avg_isize);
    printf("SN\tinsert size standard deviation:\t%.1f\n", sd_isize);
    printf("SN\tinward oriented pairs:\t%lld\n", n_inward);
    printf("SN\toutward oriented pairs:\t%lld\n", n_outward);
    printf("SN\tpairs with other orientation:\t%lld\n", n_other);
    printf("SN\tpairs on different chromosomes:\t%lld\n",
           counts.values[STATS_NREADS_ANOMALOUS] / 2);
    printf("SN\tpercentage of properly paired reads (%%):\t%.1f\n", properly_paired_pct);
}

void MpiStatsPrintOrientationDiag(const MpiBasicStatsCounts &counts) {
    long long active_in = 0;
    long long active_out = 0;
    long long active_other = 0;
    for (int isize = 0; isize < kStatsInsertBins; ++isize) {
        active_in += counts.isize_inward[isize] / 2;
        active_out += counts.isize_outward[isize] / 2;
        active_other += counts.isize_other[isize] / 2;
    }

    printf("ORIENT_DIAG\tactive_pairs\tinward=%lld\toutward=%lld\tother=%lld\n",
           active_in, active_out, active_other);
    printf("ORIENT_DIAG\tref_formula_pairs\tinward=%lld\toutward=%lld\tother=%lld\n",
           counts.orientation_diagnostics[ORIENT_DIAG_REF_IN] / 2,
           counts.orientation_diagnostics[ORIENT_DIAG_REF_OUT] / 2,
           counts.orientation_diagnostics[ORIENT_DIAG_REF_OTHER] / 2);
    printf("ORIENT_DIAG\tread_strand_pairs\tinward=%lld\toutward=%lld\tother=%lld\n",
           counts.orientation_diagnostics[ORIENT_DIAG_DOC_IN] / 2,
           counts.orientation_diagnostics[ORIENT_DIAG_DOC_OUT] / 2,
           counts.orientation_diagnostics[ORIENT_DIAG_DOC_OTHER] / 2);
    printf("ORIENT_DIAG\tref_in_by_isize_pairs\tlt_read=%lld\tlt_2read=%lld\tge_2read=%lld\n",
           counts.orientation_diagnostics[ORIENT_DIAG_REF_IN_LT_READ] / 2,
           counts.orientation_diagnostics[ORIENT_DIAG_REF_IN_LT_2READ] / 2,
           counts.orientation_diagnostics[ORIENT_DIAG_REF_IN_GE_2READ] / 2);
    printf("ORIENT_DIAG\tref_out_by_isize_pairs\tlt_read=%lld\tlt_2read=%lld\tge_2read=%lld\n",
           counts.orientation_diagnostics[ORIENT_DIAG_REF_OUT_LT_READ] / 2,
           counts.orientation_diagnostics[ORIENT_DIAG_REF_OUT_LT_2READ] / 2,
           counts.orientation_diagnostics[ORIENT_DIAG_REF_OUT_GE_2READ] / 2);
    printf("ORIENT_DIAG\tpos_strand_records\tpos_neg_strand_neg=%lld\tpos_neg_strand_pos=%lld\tpos_zero_strand_neg=%lld\tpos_zero_strand_pos=%lld\tpos_pos_strand_neg=%lld\tpos_pos_strand_pos=%lld\n",
           counts.orientation_diagnostics[ORIENT_DIAG_POS_NEG_STRAND_NEG],
           counts.orientation_diagnostics[ORIENT_DIAG_POS_NEG_STRAND_POS],
           counts.orientation_diagnostics[ORIENT_DIAG_POS_ZERO_STRAND_NEG],
           counts.orientation_diagnostics[ORIENT_DIAG_POS_ZERO_STRAND_POS],
           counts.orientation_diagnostics[ORIENT_DIAG_POS_POS_STRAND_NEG],
           counts.orientation_diagnostics[ORIENT_DIAG_POS_POS_STRAND_POS]);
    printf("ORIENT_DIAG\tread_order_position_records\tread1_left=%lld\tread1_right=%lld\tread1_same=%lld\tread2_left=%lld\tread2_right=%lld\tread2_same=%lld\n",
           counts.orientation_diagnostics[ORIENT_DIAG_READ1_LEFT],
           counts.orientation_diagnostics[ORIENT_DIAG_READ1_RIGHT],
           counts.orientation_diagnostics[ORIENT_DIAG_READ1_SAME_POS],
           counts.orientation_diagnostics[ORIENT_DIAG_READ2_LEFT],
           counts.orientation_diagnostics[ORIENT_DIAG_READ2_RIGHT],
           counts.orientation_diagnostics[ORIENT_DIAG_READ2_SAME_POS]);
    printf("ORIENT_DIAG\tisize_sign_records\tneg=%lld\tzero=%lld\tpos=%lld\n",
           counts.orientation_diagnostics[ORIENT_DIAG_ISIZE_NEG],
           counts.orientation_diagnostics[ORIENT_DIAG_ISIZE_ZERO],
           counts.orientation_diagnostics[ORIENT_DIAG_ISIZE_POS]);
}


} // namespace

int ProcessStatsMPI(CmdInfo *cmd_info) {
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
    MpiBasicStatsCounts local_counts = {};
    MpiBasicStatsCounts global_counts = {};
    MpiStatsSortState local_sort = {};
    MpiStatsPerf local_perf = {};
    MpiStatsPerf global_perf = {};

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
            printf("Enable MPI STATS --basic mode (%d MPE + %d CPEs)!!!\n",
                   comm_size, comm_size * 64);
        }
        if (swbam::mpi::PrepareMpiBamInputPlan(
                *input, &input_plan) != 0) {
            if (rank == 0) {
                fprintf(stderr, "ERROR: failed to prepare BGZF input plan for stats.\n");
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
        const int collect_diag = cmd_info->verbose_ ? 1 : 0;
        if (swbam::operators::RunStatsBasicPipeline(
                *input, input_plan.rank_spans(),
                input_plan.rank_block_count(),
                &local_counts, &local_sort, &local_perf,
                collect_diag != 0) != 0) {
            local_ok = 0;
        }
        double stage43_cost = GetTime() - stage43_t0;
        int global_ok = swbam::mpi::AllRanksOk(local_ok);
        double stage43_cost_max = swbam::mpi::ReduceMaxCost(stage43_cost);
        if (rank == 0 && global_ok) {
            printf("Complete the 4.3 OptimizedStatsMPI cost %lf\n", stage43_cost_max);
        }
        if (!global_ok) goto cleanup;

        double stage44_t0 = GetTime();
        MPI_Reduce(local_counts.values, global_counts.values,
                   STATS_LONG_COUNT, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        long long local_max_values[3] = {
            local_counts.values[STATS_MAX_LEN],
            local_counts.values[STATS_MAX_LEN_1ST],
            local_counts.values[STATS_MAX_LEN_2ND]
        };
        long long global_max_values[3] = {};
        MPI_Reduce(local_max_values, global_max_values,
                   3, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            global_counts.values[STATS_MAX_LEN] = global_max_values[0];
            global_counts.values[STATS_MAX_LEN_1ST] = global_max_values[1];
            global_counts.values[STATS_MAX_LEN_2ND] = global_max_values[2];
        }
        MPI_Reduce(local_counts.isize_inward, global_counts.isize_inward,
                   kStatsInsertBins, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(local_counts.isize_outward, global_counts.isize_outward,
                   kStatsInsertBins, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(local_counts.isize_other, global_counts.isize_other,
                   kStatsInsertBins, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        if (cmd_info->verbose_) {
            MPI_Reduce(local_counts.orientation_diagnostics,
                       global_counts.orientation_diagnostics,
                       ORIENT_DIAG_COUNT, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        MpiStatsReducePerf(local_perf, rank == 0 ? &global_perf : nullptr);

        long long local_sort_values[6] = {
            local_sort.has_coord,
            local_sort.sorted,
            local_sort.first_tid,
            local_sort.first_pos,
            local_sort.last_tid,
            local_sort.last_pos
        };
        std::vector<long long> global_sort_values;
        if (rank == 0) global_sort_values.resize((size_t)comm_size * 6);
        MPI_Gather(local_sort_values, 6, MPI_LONG_LONG,
                   rank == 0 ? global_sort_values.data() : nullptr,
                   6, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

        double stage44_cost = GetTime() - stage44_t0;
        double stage44_cost_max = swbam::mpi::ReduceMaxCost(stage44_cost);
        if (rank == 0) {
            int is_sorted = MpiStatsIsSortedGlobal(comm_size, global_sort_values);
            MpiStatsPrintBasic(global_counts, is_sorted);
            if (cmd_info->verbose_) MpiStatsPrintOrientationDiag(global_counts);
            printf("OptimizedStatsMPI finished. ranks=%d in_blocks=%lld groups=%lld total_records=%lld\n",
                   comm_size,
                   global_perf.input_blocks,
                   global_perf.batch_count,
                   global_perf.total_records);
            printf("  read_sum=%.3f  decomp_sum=%.3f  merge_sum=%.3f  optimized_total_sum=%.3f\n",
                   global_perf.read,
                   global_perf.kernel,
                   global_perf.post_process,
                   global_perf.total);
            printf("  decomp_detail_sum alloc=%.3f  inflate=%.3f  crc=%.3f  parse=%.3f  other=%.3f\n",
                   global_perf.decomp_alloc,
                   global_perf.decomp_inflate,
                   global_perf.decomp_crc,
                   global_perf.decomp_parse,
                   global_perf.decomp_other);
            printf("Complete the 4.4 stats reduce/print cost %lf\n", stage44_cost_max);
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
