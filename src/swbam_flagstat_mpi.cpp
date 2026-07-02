#include "swbam_mpi.h"
#include "swbam/cpe_pipeline.h"
#include "swbam/io.h"
#include "swbam/mpi_runtime.h"

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <vector>

#include <mpi.h>

extern "C" {
    void slave_mpi_flagstat_count();
}

namespace {

const int kFlagstatNB = 64;

enum FlagstatCounterId {
    FLAGSTAT_TOTAL = 0,
    FLAGSTAT_PRIMARY,
    FLAGSTAT_SECONDARY,
    FLAGSTAT_SUPPLEMENTARY,
    FLAGSTAT_DUPLICATES,
    FLAGSTAT_PRIMARY_DUPLICATES,
    FLAGSTAT_MAPPED,
    FLAGSTAT_PRIMARY_MAPPED,
    FLAGSTAT_PAIRED,
    FLAGSTAT_READ1,
    FLAGSTAT_READ2,
    FLAGSTAT_PROPERLY_PAIRED,
    FLAGSTAT_PAIR_MAPPED,
    FLAGSTAT_SINGLETONS,
    FLAGSTAT_DIFF_CHR,
    FLAGSTAT_DIFF_CHR_MAPQ5,
    FLAGSTAT_COUNTER_COUNT
};

struct MpiFlagstatCounts {
    long long values[FLAGSTAT_COUNTER_COUNT][2];
};

struct MpiFlagstatStats {
    long long input_blocks;
    long long group_count;
    long long total_records;
    double t_read;
    double t_decomp;
    double t_decomp_alloc;
    double t_decomp_inflate;
    double t_decomp_crc;
    double t_decomp_parse;
    double t_decomp_other;
    double t_count;
    double t_fused_total;
};

BamFilterOptions MpiFlagstatNoFilter() {
    BamFilterOptions filter;
    filter.min_mapq = -1;
    filter.max_mapq = -1;
    filter.require_flag = 0;
    filter.exclude_flag = 0;
    filter.ref_tid = -2;
    filter.min_read_len = -1;
    filter.max_read_len = -1;
    return filter;
}

void MpiFlagstatAccumulateDecompDetail(const MpiFlagstatCountPara *paras,
                                       int active_blocks,
                                       double decomp_wall,
                                       MpiFlagstatStats *stats) {
    if (!stats || active_blocks <= 0 || decomp_wall <= 0.0) return;

    const MpiFlagstatCountPara *critical = nullptr;
    uint64_t critical_total = 0;
    for (int b = 0; b < active_blocks; ++b) {
        if (paras[b].decomp_total_cycles >= critical_total) {
            critical_total = paras[b].decomp_total_cycles;
            critical = &paras[b];
        }
    }
    if (!critical || critical_total == 0) {
        stats->t_decomp_other += decomp_wall;
        return;
    }

    const double scale = decomp_wall / (double)critical_total;
    const double alloc_time = scale * (double)critical->decomp_alloc_cycles;
    const double inflate_time = scale * (double)critical->decomp_inflate_cycles;
    const double crc_time = scale * (double)critical->decomp_crc_cycles;
    const double parse_time = scale * (double)critical->decomp_parse_cycles;
    double other_time = decomp_wall - alloc_time - inflate_time - crc_time - parse_time;
    if (other_time < 0.0) other_time = 0.0;

    stats->t_decomp_alloc += alloc_time;
    stats->t_decomp_inflate += inflate_time;
    stats->t_decomp_crc += crc_time;
    stats->t_decomp_parse += parse_time;
    stats->t_decomp_other += other_time;
}

void MpiFlagstatAddRecord(const bam1_t *record, MpiFlagstatCounts *counts) {
    const uint16_t flag = record->core.flag;
    const int bucket = (flag & BAM_FQCFAIL) ? 1 : 0;
    const bool is_secondary = (flag & BAM_FSECONDARY) != 0;
    const bool is_supplementary = (flag & BAM_FSUPPLEMENTARY) != 0;
    const bool is_primary = !is_secondary && !is_supplementary;
    const bool is_paired = (flag & BAM_FPAIRED) != 0;
    const bool is_mapped = (flag & BAM_FUNMAP) == 0;
    const bool mate_mapped = (flag & BAM_FMUNMAP) == 0;

    counts->values[FLAGSTAT_TOTAL][bucket]++;
    if (is_primary) counts->values[FLAGSTAT_PRIMARY][bucket]++;
    if (is_secondary) counts->values[FLAGSTAT_SECONDARY][bucket]++;
    if (is_supplementary) counts->values[FLAGSTAT_SUPPLEMENTARY][bucket]++;
    if (flag & BAM_FDUP) counts->values[FLAGSTAT_DUPLICATES][bucket]++;
    if (is_primary && (flag & BAM_FDUP)) counts->values[FLAGSTAT_PRIMARY_DUPLICATES][bucket]++;
    if (is_mapped) counts->values[FLAGSTAT_MAPPED][bucket]++;
    if (is_primary && is_mapped) counts->values[FLAGSTAT_PRIMARY_MAPPED][bucket]++;

    if (is_primary && is_paired) {
        counts->values[FLAGSTAT_PAIRED][bucket]++;
        if (flag & BAM_FREAD1) counts->values[FLAGSTAT_READ1][bucket]++;
        if (flag & BAM_FREAD2) counts->values[FLAGSTAT_READ2][bucket]++;
        if ((flag & BAM_FPROPER_PAIR) && is_mapped) {
            counts->values[FLAGSTAT_PROPERLY_PAIRED][bucket]++;
        }
        if (is_mapped && mate_mapped) {
            counts->values[FLAGSTAT_PAIR_MAPPED][bucket]++;
            if (record->core.tid != record->core.mtid) {
                counts->values[FLAGSTAT_DIFF_CHR][bucket]++;
                if (record->core.qual >= 5) counts->values[FLAGSTAT_DIFF_CHR_MAPQ5][bucket]++;
            }
        }
        if (is_mapped && !mate_mapped) {
            counts->values[FLAGSTAT_SINGLETONS][bucket]++;
        }
    }
}

void MpiFlagstatMergeSlice(const MpiFlagstatCountSlice &slice,
                           MpiFlagstatCounts *counts) {
    for (int i = 0; i < FLAGSTAT_COUNTER_COUNT; ++i) {
        counts->values[i][0] += slice.values[i][0];
        counts->values[i][1] += slice.values[i][1];
    }
}

void MpiFlagstatFormatPct(char *buf, size_t buf_size, long long value, long long total) {
    if (total > 0) {
        snprintf(buf, buf_size, "%.2f%%", 100.0 * (double)value / (double)total);
    } else {
        snprintf(buf, buf_size, "N/A");
    }
}

void MpiFlagstatPrintSimple(const MpiFlagstatCounts &counts,
                            FlagstatCounterId id,
                            const char *label) {
    printf("%lld + %lld %s\n", counts.values[id][0], counts.values[id][1], label);
}

void MpiFlagstatPrintPct(const MpiFlagstatCounts &counts,
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

void MpiFlagstatPrintSamtoolsStyle(const MpiFlagstatCounts &counts) {
    MpiFlagstatPrintSimple(counts, FLAGSTAT_TOTAL,
                           "in total (QC-passed reads + QC-failed reads)");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_PRIMARY, "primary");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_SECONDARY, "secondary");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_SUPPLEMENTARY, "supplementary");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_DUPLICATES, "duplicates");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_PRIMARY_DUPLICATES, "primary duplicates");
    MpiFlagstatPrintPct(counts, FLAGSTAT_MAPPED, FLAGSTAT_TOTAL, "mapped");
    MpiFlagstatPrintPct(counts, FLAGSTAT_PRIMARY_MAPPED, FLAGSTAT_PRIMARY, "primary mapped");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_PAIRED, "paired in sequencing");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_READ1, "read1");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_READ2, "read2");
    MpiFlagstatPrintPct(counts, FLAGSTAT_PROPERLY_PAIRED, FLAGSTAT_PAIRED, "properly paired");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_PAIR_MAPPED, "with itself and mate mapped");
    MpiFlagstatPrintPct(counts, FLAGSTAT_SINGLETONS, FLAGSTAT_PAIRED, "singletons");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_DIFF_CHR, "with mate mapped to a different chr");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_DIFF_CHR_MAPQ5,
                           "with mate mapped to a different chr (mapQ>=5)");
}

void MpiFlagstatReduceStats(const MpiFlagstatStats &local_stats,
                            MpiFlagstatStats *global_stats) {
    long long local_long[3] = {
        local_stats.input_blocks,
        local_stats.group_count,
        local_stats.total_records
    };
    long long global_long[3] = {};
    double local_double[9] = {
        local_stats.t_read,
        local_stats.t_decomp,
        local_stats.t_decomp_alloc,
        local_stats.t_decomp_inflate,
        local_stats.t_decomp_crc,
        local_stats.t_decomp_parse,
        local_stats.t_decomp_other,
        local_stats.t_count,
        local_stats.t_fused_total
    };
    double global_double[9] = {};

    MPI_Reduce(local_long, global_long, 3, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_double, global_double, 9, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (global_stats) {
        global_stats->input_blocks = global_long[0];
        global_stats->group_count = global_long[1];
        global_stats->total_records = global_long[2];
        global_stats->t_read = global_double[0];
        global_stats->t_decomp = global_double[1];
        global_stats->t_decomp_alloc = global_double[2];
        global_stats->t_decomp_inflate = global_double[3];
        global_stats->t_decomp_crc = global_double[4];
        global_stats->t_decomp_parse = global_double[5];
        global_stats->t_decomp_other = global_double[6];
        global_stats->t_count = global_double[7];
        global_stats->t_fused_total = global_double[8];
    }
}

class FlagstatCpeOperator : public swbam::cpe::CpeBatchOperator {
public:
    FlagstatCpeOperator(MpiFlagstatCounts *counts, MpiFlagstatStats *stats)
        : counts_(counts), stats_(stats), count_slices_(nullptr),
          scratch_data_(nullptr) {
        memset(paras_, 0, sizeof(paras_));
    }

    const char *name() const { return "flagstat"; }
    size_t batch_capacity() const { return kFlagstatNB; }

    int Initialize() {
        if (!counts_) return -1;
        memset(counts_, 0, sizeof(*counts_));
        if (stats_) memset(stats_, 0, sizeof(*stats_));
        count_slices_ = (MpiFlagstatCountSlice *)aligned_alloc_custom(
            64, (size_t)kFlagstatNB * sizeof(MpiFlagstatCountSlice));
        scratch_data_ = aligned_alloc_custom(
            64, (size_t)kFlagstatNB * MPI_BAM_BLOCK_ARENA_SIZE);
        return count_slices_ && scratch_data_ ? 0 : -1;
    }

    void Shutdown() {
        if (count_slices_) {
            aligned_free_custom((unsigned char *)count_slices_);
            count_slices_ = nullptr;
        }
        if (scratch_data_) {
            aligned_free_custom(scratch_data_);
            scratch_data_ = nullptr;
        }
    }

    int Prepare(const swbam::BgzfBlockBatch &compressed,
                swbam::BgzfBlockBatch *decoded,
                size_t active_blocks) {
        memset(count_slices_, 0,
               (size_t)kFlagstatNB * sizeof(MpiFlagstatCountSlice));
        for (int b = 0; b < kFlagstatNB; ++b) {
            MpiFlagstatCountPara &para = paras_[b];
            para.block_id = b;
            para.scratch_data = scratch_data_ +
                (size_t)b * MPI_BAM_BLOCK_ARENA_SIZE;
            para.scratch_capacity = MPI_BAM_BLOCK_ARENA_SIZE;
            para.counts = &count_slices_[b];
            para.n_total_records = 0;
            para.record_index = 0;
            para.actual_value = 0;
            para.limit_value = 0;
            para.limit_id = BOUNDS_LIMIT_NONE;
            para.decomp_alloc_cycles = 0;
            para.decomp_inflate_cycles = 0;
            para.decomp_crc_cycles = 0;
            para.decomp_parse_cycles = 0;
            para.decomp_total_cycles = 0;
            if ((size_t)b < active_blocks) {
                para.input_block = const_cast<bam_block *>(
                    &compressed.blocks()[b]);
                para.un_comp_block = &decoded->blocks()[b];
                para.status = 0;
            } else {
                para.input_block = nullptr;
                para.un_comp_block = nullptr;
                para.status = -1;
            }
        }
        return 0;
    }

    void *kernel_entry() const {
        return (void *)slave_mpi_flagstat_count;
    }
    void *kernel_arguments() { return paras_; }

    void ObserveKernel(double wall_seconds, size_t active_blocks) {
        if (!stats_) return;
        MpiFlagstatAccumulateDecompDetail(
            paras_, (int)active_blocks, wall_seconds, stats_);
    }

    int Validate(size_t active_blocks) const {
        for (size_t b = 0; b < active_blocks; ++b) {
            const MpiFlagstatCountPara &para = paras_[b];
            if (para.status == 0) continue;
            if (para.status == -3) {
                fprintf(stderr,
                        "ERROR: MPI flagstat capacity exceeded on input block %zu. limit_id=%d limit=%lld actual=%lld record=%d.\n",
                        b, para.limit_id, para.limit_value,
                        para.actual_value, para.record_index);
            } else {
                fprintf(stderr,
                        "ERROR: MPI flagstat count failed on input block %zu with status %d.\n",
                        b, para.status);
            }
            return -1;
        }
        return 0;
    }

    int Consume(size_t active_blocks, long long *records_processed) {
        long long records = 0;
        for (size_t b = 0; b < active_blocks; ++b) {
            records += paras_[b].n_total_records;
            MpiFlagstatMergeSlice(count_slices_[b], counts_);
        }
        if (records_processed) *records_processed = records;
        return 0;
    }

    int Finish() { return 0; }

private:
    MpiFlagstatCounts *counts_;
    MpiFlagstatStats *stats_;
    MpiFlagstatCountPara paras_[kFlagstatNB];
    MpiFlagstatCountSlice *count_slices_;
    unsigned char *scratch_data_;
};

int FusedFlagstatMPI(const swbam::BamInputBackend &input,
                     const swbam::BgzfBlockSpan *spans,
                     size_t span_count,
                     MpiFlagstatCounts *counts,
                     MpiFlagstatStats *stats) {
    FlagstatCpeOperator op(counts, stats);
    swbam::cpe::CpeReadPipelineTiming timing;
    const int ret = swbam::cpe::RunCpeReadPipeline(
        input, spans, span_count, &op, &timing);
    if (stats) {
        stats->input_blocks = timing.input_blocks;
        stats->group_count = timing.batch_count;
        stats->total_records = timing.total_records;
        stats->t_read = timing.read;
        stats->t_decomp = timing.kernel;
        stats->t_count = timing.consume;
        stats->t_fused_total = timing.total;
    }
    return ret;
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
    swbam::MemoryBamInput input;
    swbam::mpi::MpiBamInputPlan input_plan;
    MpiFlagstatCounts local_counts = {};
    MpiFlagstatCounts global_counts = {};
    MpiFlagstatStats local_stats = {};
    MpiFlagstatStats global_stats = {};

    double init_cost = GetTime() - t_init;
    double init_cost_max = swbam::mpi::ReduceMaxCost(init_cost);
    if (rank == 0) {
        printf("111Complete the initialization cost %lf-----\n", init_cost_max);
    }

    {
        double preload_t0 = GetTime();
        if (input.Load(cmd_info->in_file_name_) != 0) {
            fprintf(stderr, "[rank %d] ERROR: cannot preload input %s into memory\n",
                    rank, cmd_info->in_file_name_.c_str());
            local_ok = 0;
        }
        double preload_cost = GetTime() - preload_t0;
        double preload_cost_max = swbam::mpi::ReduceMaxCost(preload_cost);
        if (rank == 0 && local_ok) printf("222Complete the memory cost %lf--\n", preload_cost_max);
    }
    if (!swbam::mpi::AllRanksOk(local_ok)) goto cleanup;

    {
        double header_t0 = GetTime();

        if (input.ParseHeader() != 0) {
            fprintf(stderr, "[rank %d] ERROR: cannot read BAM header from %s\n",
                    rank, cmd_info->in_file_name_.c_str());
            local_ok = 0;
        }
        if (local_ok && input.format() != bam) {
                if (rank == 0) {
                    fprintf(stderr, "ERROR: RabbitBAM-MPI flagstat only supports BAM input in v1.\n");
                }
                local_ok = 0;
        }

        double header_cost = GetTime() - header_t0;
        double header_cost_max = swbam::mpi::ReduceMaxCost(header_cost);
        if (rank == 0 && local_ok) {
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
                input, &input_plan) != 0) {
            if (rank == 0) {
                fprintf(stderr, "ERROR: failed to prepare BGZF input plan for flagstat.\n");
            }
            local_ok = 0;
        }
        if (rank == 0 && local_ok) {
            printf("MPI BAM scan complete. data_blocks=%lld body_start=%lld header_end=%lld\n",
                   (long long)input_plan.blocks.size(),
                   (long long)input_plan.body_offset,
                   (long long)input.body_offset());
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
        if (FusedFlagstatMPI(input, input_plan.rank_spans(),
                             input_plan.rank_block_count(),
                             &local_counts, &local_stats) != 0) {
            local_ok = 0;
        }
        double stage43_cost = GetTime() - stage43_t0;
        int global_ok = swbam::mpi::AllRanksOk(local_ok);
        double stage43_cost_max = swbam::mpi::ReduceMaxCost(stage43_cost);
        if (rank == 0 && global_ok) {
            printf("Complete the 4.3 FusedFlagstatMPI cost %lf\n", stage43_cost_max);
        }
        if (!global_ok) goto cleanup;

        double stage44_t0 = GetTime();
        MPI_Reduce(&local_counts.values[0][0], &global_counts.values[0][0],
                   FLAGSTAT_COUNTER_COUNT * 2, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MpiFlagstatReduceStats(local_stats, rank == 0 ? &global_stats : nullptr);
        double stage44_cost = GetTime() - stage44_t0;
        double stage44_cost_max = swbam::mpi::ReduceMaxCost(stage44_cost);
        if (rank == 0) {
            MpiFlagstatPrintSamtoolsStyle(global_counts);
            printf("FusedFlagstatMPI finished. ranks=%d in_blocks=%lld groups=%lld total_records=%lld\n",
                   comm_size,
                   global_stats.input_blocks,
                   global_stats.group_count,
                   global_stats.total_records);
            printf("  read_sum=%.3f  decomp_sum=%.3f  merge_sum=%.3f  fused_total_sum=%.3f\n",
                   global_stats.t_read,
                   global_stats.t_decomp,
                   global_stats.t_count,
                   global_stats.t_fused_total);
            printf("  decomp_detail_sum alloc=%.3f  inflate=%.3f  crc=%.3f  parse=%.3f  other=%.3f\n",
                   global_stats.t_decomp_alloc,
                   global_stats.t_decomp_inflate,
                   global_stats.t_decomp_crc,
                   global_stats.t_decomp_parse,
                   global_stats.t_decomp_other);
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
        input.Close();
        double close_cost = GetTime() - close_t0;
        double close_cost_max = swbam::mpi::ReduceMaxCost(close_cost);
        if (rank == 0) printf("666close the files cost %lf-----\n", close_cost_max);
    }

    return exit_code;
}
