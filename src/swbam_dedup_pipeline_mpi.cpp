#include "CmdInfo.h"
#include "swbam_mpi.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <mpi.h>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
extern "C" {
    void slave_mpi_collate_release_caches(void *);
    void slave_mpi_sort_release_caches(void *);
    void slave_mpi_common_release_caches(void *);
}
#endif

namespace {

double PipelineReduceMax(double local_cost) {
    double max_cost = 0.0;
    MPI_Reduce(&local_cost, &max_cost, 1, MPI_DOUBLE,
               MPI_MAX, 0, MPI_COMM_WORLD);
    return max_cost;
}

int PipelineAllRanksOk(int local_ok) {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    return global_ok;
}

int PipelineBroadcastStage(const char *name, MpiMemoryBam *bam,
                           int rank, double *broadcast_cost) {
    double t0 = GetTime();
    int local_ok = MpiBroadcastMemoryBam(bam, 0) == 0 ? 1 : 0;
    if (local_ok && bam && bam->size > 0) {
        char *refreshed = (char *)malloc(bam->size);
        if (!refreshed) {
            local_ok = 0;
        } else {
            memcpy(refreshed, bam->data, bam->size);
            free(bam->data);
            bam->data = refreshed;
        }
    }
    *broadcast_cost = PipelineReduceMax(GetTime() - t0);
    if (rank == 0 && local_ok) {
        printf("%s broadcast complete. size=%zu cost=%lf\n",
               name, bam->size, *broadcast_cost);
    }
    return PipelineAllRanksOk(local_ok) ? 0 : -1;
}

std::string PipelineStageDumpPath(const char *prefix,
                                  const char *stage) {
    std::string path = prefix ? prefix : "";
    path += ".";
    path += stage ? stage : "stage";
    path += ".bam";
    return path;
}

int PipelineDumpStage(const char *stage, const MpiMemoryBam &bam,
                      const std::string &path, int rank) {
    double t0 = GetTime();
    int local_ok = 1;
    if (rank == 0) {
        if (!bam.data || bam.size == 0 ||
            MpiCommonDumpMemoryToFile(path, bam.data,
                                      bam.size) != 0) {
            fprintf(stderr,
                    "ERROR: failed to dump pipeline %s BAM to %s.\n",
                    stage ? stage : "stage", path.c_str());
            local_ok = 0;
        } else {
            printf("pipeline debug dump %s size=%zu path=%s cost=%lf\n",
                   stage ? stage : "stage", bam.size,
                   path.c_str(), GetTime() - t0);
        }
    }
    return PipelineAllRanksOk(local_ok) ? 0 : -1;
}

int PipelineMaybeDebugDump(const char *stage, const MpiMemoryBam &bam,
                           int rank) {
    const char *prefix =
        getenv("RABBITBAM_PIPELINE_DEBUG_PREFIX");
    if (!prefix || !*prefix) return 0;
    return PipelineDumpStage(stage, bam,
                             PipelineStageDumpPath(prefix, stage),
                             rank);
}

void PipelineDumpFailureInput(const char *stage,
                              const MpiMemoryBam &bam,
                              const std::string &out_path,
                              int rank) {
    std::string path = out_path;
    path += ".";
    path += stage ? stage : "failed";
    path += "_input.bam";
    (void)PipelineDumpStage(stage, bam, path, rank);
}

void PipelineReleaseSlaveCaches(const char *stage, int rank) {
#ifdef PLATFORM_SUNWAY
    double t0 = GetTime();
    __real_athread_spawn(
        (void *)slave_mpi_collate_release_caches,
        nullptr, 1);
    athread_join();
    __real_athread_spawn(
        (void *)slave_mpi_sort_release_caches,
        nullptr, 1);
    athread_join();
    __real_athread_spawn(
        (void *)slave_mpi_common_release_caches,
        nullptr, 1);
    athread_join();
    double cost = PipelineReduceMax(GetTime() - t0);
    if (rank == 0) {
        printf("%s CPE cache release cost %lf\n",
               stage, cost);
    }
#else
    (void)stage;
    (void)rank;
#endif
}

void PipelineResetSlaveRuntime(const char *stage, int rank) {
#ifdef PLATFORM_SUNWAY
    double t0 = GetTime();
    athread_halt();
    athread_init();
    double cost = PipelineReduceMax(GetTime() - t0);
    if (rank == 0) {
        printf("%s CPE runtime reset cost %lf\n",
               stage, cost);
    }
#else
    (void)stage;
    (void)rank;
#endif
}

} // namespace

int ProcessDedupPipelineMPI(CmdInfo *cmd_info) {
    double pipeline_t0 = GetTime();
    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int local_ok = 1;
    int exit_code = 1;
    MpiMemoryBam current;
    MpiMemoryBam next;
    MpiMemoryBam final_bam;
    double pipeline_core = 0.0;
    double collate_core = 0.0;
    double fixmate_core = 0.0;
    double sort_core = 0.0;
    double markdup_core = 0.0;
    double p0_cost = 0.0;
    double p1_cost = 0.0;
    double p2_cost = 0.0;
    double p3_cost = 0.0;
    double p4_cost = 0.0;
    double p0_load_cost = 0.0;
    double p0_broadcast_cost = 0.0;
    double broadcast_cost = 0.0;

    if (cmd_info->markdup_remove_dups_ ||
        cmd_info->markdup_clear_) {
        if (rank == 0) {
            fprintf(stderr,
                    "ERROR: dedup-pipeline v1 supports default markdup only; -r/-c are not supported.\n");
        }
        local_ok = 0;
    }
    if (!PipelineAllRanksOk(local_ok)) goto cleanup;

    if (rank == 0) {
        printf("Enable MPI BAM DEDUP-PIPELINE mode "
               "(%d MPE + %d CPEs)!!!\n",
               comm_size, comm_size * 64);
        printf("dedup-pipeline stages: collate -> fixmate -m -> sort -> markdup\n");
        printf("dedup-pipeline compression=%d memory='%s' bins=%d\n",
               cmd_info->compress_level_,
               cmd_info->pipeline_memory_.c_str(),
               cmd_info->collate_bins_);
    }

    {
        double t0 = GetTime();
        double load_t0 = GetTime();
        if (rank == 0 &&
            MpiCommonLoadFileToMemory(cmd_info->in_file_name_,
                                      &current.data,
                                      &current.size) != 0) {
            fprintf(stderr,
                    "ERROR: failed to load dedup-pipeline input %s.\n",
                    cmd_info->in_file_name_.c_str());
            local_ok = 0;
        }
        int load_ok = PipelineAllRanksOk(local_ok);
        p0_load_cost = PipelineReduceMax(GetTime() - load_t0);
        if (!load_ok) goto cleanup;
        double broadcast_t0 = GetTime();
        if (MpiBroadcastMemoryBam(&current, 0) != 0) {
            local_ok = 0;
        }
        p0_broadcast_cost =
            PipelineReduceMax(GetTime() - broadcast_t0);
        p0_cost = PipelineReduceMax(GetTime() - t0);
        int p0_ok = PipelineAllRanksOk(local_ok);
        if (rank == 0 && p0_ok) {
            printf("P0 load input memory cost %lf broadcast %lf total %lf size=%zu\n",
                   p0_load_cost, p0_broadcast_cost,
                   p0_cost, current.size);
        }
        if (!p0_ok) goto cleanup;
    }

    {
        CmdInfo stage = *cmd_info;
        stage.collate_memory_ = cmd_info->pipeline_memory_;
        stage.collate_temp_prefix_.clear();
        double t0 = GetTime();
        if (MpiCollateMemoryToMemory(&stage, current.data,
                                     current.size, &next,
                                     &collate_core) != 0) {
            PipelineDumpFailureInput(
                "collate", current,
                cmd_info->out_file_name_, rank);
            local_ok = 0;
        }
        if (!PipelineAllRanksOk(local_ok)) goto cleanup;
        PipelineReleaseSlaveCaches("P1 collate", rank);
        PipelineResetSlaveRuntime("P1 collate", rank);
        if (PipelineBroadcastStage("P1 collate output",
                                   &next, rank,
                                   &broadcast_cost) != 0) {
            local_ok = 0;
        }
        if (local_ok &&
            PipelineMaybeDebugDump("collate", next,
                                   rank) != 0) {
            local_ok = 0;
        }
        p1_cost = PipelineReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("P1 collate memory output/broadcast cost %lf "
                   "core=%lf broadcast=%lf\n",
                   p1_cost, collate_core, broadcast_cost);
        }
        MpiMemoryBamFree(&current);
        current = next;
        next.data = nullptr;
        next.size = 0;
        pipeline_core += collate_core;
    }
    if (!PipelineAllRanksOk(local_ok)) goto cleanup;

    {
        CmdInfo stage = *cmd_info;
        stage.fixmate_mate_score_ = true;
        double t0 = GetTime();
        if (MpiFixmateMemoryToMemory(&stage, current.data,
                                     current.size, &next,
                                     &fixmate_core) != 0) {
            PipelineDumpFailureInput(
                "fixmate", current,
                cmd_info->out_file_name_, rank);
            local_ok = 0;
        }
        if (!PipelineAllRanksOk(local_ok)) goto cleanup;
        PipelineReleaseSlaveCaches("P2 fixmate", rank);
        PipelineResetSlaveRuntime("P2 fixmate", rank);
        if (PipelineBroadcastStage("P2 fixmate output",
                                   &next, rank,
                                   &broadcast_cost) != 0) {
            local_ok = 0;
        }
        if (local_ok &&
            PipelineMaybeDebugDump("fixmate", next,
                                   rank) != 0) {
            local_ok = 0;
        }
        p2_cost = PipelineReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("P2 fixmate memory output/broadcast cost %lf "
                   "core=%lf broadcast=%lf\n",
                   p2_cost, fixmate_core, broadcast_cost);
        }
        MpiMemoryBamFree(&current);
        current = next;
        next.data = nullptr;
        next.size = 0;
        pipeline_core += fixmate_core;
    }
    if (!PipelineAllRanksOk(local_ok)) goto cleanup;

    {
        CmdInfo stage = *cmd_info;
        stage.sort_memory_ = cmd_info->pipeline_memory_;
        stage.sort_temp_prefix_.clear();
        double t0 = GetTime();
        if (MpiSortMemoryToMemory(&stage, current.data,
                                  current.size, &next,
                                  &sort_core) != 0) {
            PipelineDumpFailureInput(
                "sort", current,
                cmd_info->out_file_name_, rank);
            local_ok = 0;
        }
        if (!PipelineAllRanksOk(local_ok)) goto cleanup;
        PipelineReleaseSlaveCaches("P3 sort", rank);
        PipelineResetSlaveRuntime("P3 sort", rank);
        if (PipelineBroadcastStage("P3 sort output",
                                   &next, rank,
                                   &broadcast_cost) != 0) {
            local_ok = 0;
        }
        if (local_ok &&
            PipelineMaybeDebugDump("sort", next,
                                   rank) != 0) {
            local_ok = 0;
        }
        p3_cost = PipelineReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("P3 sort memory output/broadcast cost %lf "
                   "core=%lf broadcast=%lf\n",
                   p3_cost, sort_core, broadcast_cost);
        }
        MpiMemoryBamFree(&current);
        current = next;
        next.data = nullptr;
        next.size = 0;
        pipeline_core += sort_core;
    }
    if (!PipelineAllRanksOk(local_ok)) goto cleanup;

    {
        CmdInfo stage = *cmd_info;
        stage.markdup_memory_ = cmd_info->pipeline_memory_;
        stage.markdup_remove_dups_ = false;
        stage.markdup_clear_ = false;
        double t0 = GetTime();
        if (MpiMarkdupMemoryToMemory(&stage, current.data,
                                     current.size, &final_bam,
                                     &markdup_core) != 0) {
            PipelineDumpFailureInput(
                "markdup", current,
                cmd_info->out_file_name_, rank);
            local_ok = 0;
        }
        if (!PipelineAllRanksOk(local_ok)) goto cleanup;
        if (PipelineMaybeDebugDump("markdup",
                                   final_bam, rank) != 0) {
            local_ok = 0;
        }
        if (!PipelineAllRanksOk(local_ok)) goto cleanup;
        double dump_t0 = GetTime();
        if (rank == 0 &&
            MpiCommonDumpMemoryToFile(cmd_info->out_file_name_,
                                      final_bam.data,
                                      final_bam.size) != 0) {
            fprintf(stderr,
                    "ERROR: failed to dump dedup-pipeline output %s.\n",
                    cmd_info->out_file_name_.c_str());
            local_ok = 0;
        }
        double dump_cost = PipelineReduceMax(GetTime() - dump_t0);
        p4_cost = PipelineReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("P4 markdup final dump cost %lf core=%lf dump=%lf size=%zu\n",
                   p4_cost, markdup_core, dump_cost,
                   final_bam.size);
        }
        pipeline_core += markdup_core;
    }
    if (!PipelineAllRanksOk(local_ok)) goto cleanup;

    if (rank == 0) {
        printf("dedup-pipeline pipeline_core=%lf "
               "collate=%lf fixmate=%lf sort=%lf markdup=%lf\n",
               pipeline_core, collate_core, fixmate_core,
               sort_core, markdup_core);
        printf("dedup-pipeline pipeline_total=%lf "
               "P0=%lf P1=%lf P2=%lf P3=%lf P4=%lf\n",
               GetTime() - pipeline_t0,
               p0_cost, p1_cost, p2_cost, p3_cost, p4_cost);
    }
    exit_code = 0;

cleanup:
    MpiMemoryBamFree(&current);
    MpiMemoryBamFree(&next);
    MpiMemoryBamFree(&final_bam);
    return exit_code;
}
