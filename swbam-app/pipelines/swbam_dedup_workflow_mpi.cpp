#include "CmdInfo.h"
#include "swbam_mpi.h"

#include <cerrno>
#include <cstdio>
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

typedef int (*WorkflowStageFn)(CmdInfo *);

int WorkflowAllRanksOk(int local_ok) {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    return global_ok;
}

double WorkflowReduceMax(double local_seconds) {
    double max_seconds = 0.0;
    MPI_Reduce(&local_seconds, &max_seconds, 1, MPI_DOUBLE,
               MPI_MAX, 0, MPI_COMM_WORLD);
    return max_seconds;
}

int WorkflowRunStage(const char *name, WorkflowStageFn function,
                     CmdInfo *stage, int rank, double *wall_seconds) {
    MPI_Barrier(MPI_COMM_WORLD);
    const double t0 = GetTime();
    const int local_ok = function && stage && function(stage) == 0;
    const int global_ok = WorkflowAllRanksOk(local_ok);
    const double stage_wall = WorkflowReduceMax(GetTime() - t0);
    if (wall_seconds) *wall_seconds = stage_wall;
    if (rank == 0) {
        printf("dedup-workflow stage=%s status=%s wall=%lf input=%s output=%s\n",
               name ? name : "unknown", global_ok ? "ok" : "failed",
               stage_wall, stage ? stage->in_file_name_.c_str() : "",
               stage ? stage->out_file_name_.c_str() : "");
    }
    return global_ok ? 0 : -1;
}

void WorkflowResetCpeRuntime(const char *after_stage, int rank) {
#ifdef PLATFORM_SUNWAY
    const double t0 = GetTime();
    __real_athread_spawn(
        (void *)slave_mpi_collate_release_caches, nullptr, 1);
    athread_join();
    __real_athread_spawn(
        (void *)slave_mpi_sort_release_caches, nullptr, 1);
    athread_join();
    __real_athread_spawn(
        (void *)slave_mpi_common_release_caches, nullptr, 1);
    athread_join();
    athread_halt();
    athread_init();
    const double reset_wall = WorkflowReduceMax(GetTime() - t0);
    if (rank == 0) {
        printf("dedup-workflow reset_after=%s wall=%lf\n",
               after_stage ? after_stage : "stage", reset_wall);
    }
#else
    (void)after_stage;
    (void)rank;
#endif
}

void WorkflowRemoveIntermediate(const std::string &path, int rank) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0 && std::remove(path.c_str()) != 0 && errno != ENOENT) {
        fprintf(stderr,
                "WARNING: failed to remove dedup-workflow intermediate "
                "%s: %s\n",
                path.c_str(), strerror(errno));
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

bool WorkflowPathsValid(const CmdInfo &cmd_info,
                        const std::string &collate_path,
                        const std::string &fixmate_path,
                        const std::string &sort_path) {
    if (cmd_info.in_file_name_.empty() ||
        cmd_info.out_file_name_.empty() ||
        cmd_info.workflow_temp_prefix_.empty()) {
        return false;
    }
    const std::string &input = cmd_info.in_file_name_;
    const std::string &output = cmd_info.out_file_name_;
    return input != output &&
           input != collate_path && input != fixmate_path &&
           input != sort_path && output != collate_path &&
           output != fixmate_path && output != sort_path &&
           collate_path != fixmate_path && collate_path != sort_path &&
           fixmate_path != sort_path;
}

} // namespace

int ProcessDedupWorkflowMPI(CmdInfo *cmd_info) {
    if (!cmd_info) return 1;

    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    const double workflow_t0 = GetTime();
    const std::string collate_path =
        cmd_info->workflow_temp_prefix_ + ".collate.bam";
    const std::string fixmate_path =
        cmd_info->workflow_temp_prefix_ + ".fixmate.bam";
    const std::string sort_path =
        cmd_info->workflow_temp_prefix_ + ".sort.bam";
    int local_ok = WorkflowPathsValid(
        *cmd_info, collate_path, fixmate_path, sort_path) ? 1 : 0;
    if (!WorkflowAllRanksOk(local_ok)) {
        if (rank == 0) {
            fprintf(stderr,
                    "ERROR: dedup-workflow requires distinct input, output "
                    "and non-empty -T intermediate paths.\n");
        }
        return 1;
    }

    if (rank == 0) {
        printf("Enable MPI BAM DEDUP-WORKFLOW mode "
               "(%d MPE + %d CPEs)!!!\n",
               comm_size, comm_size * 64);
        printf("dedup-workflow stages: collate -> fixmate -m -> "
               "sort -> markdup coordinate-stream-two-pass\n");
        printf("dedup-workflow temp_prefix=%s keep_intermediates=%d "
               "memory='%s' bins=%d compression=%d\n",
               cmd_info->workflow_temp_prefix_.c_str(),
               cmd_info->workflow_keep_intermediates_ ? 1 : 0,
               cmd_info->pipeline_memory_.c_str(),
               cmd_info->collate_bins_, cmd_info->compress_level_);
    }

    double collate_wall = 0.0;
    double fixmate_wall = 0.0;
    double sort_wall = 0.0;
    double markdup_wall = 0.0;

    CmdInfo collate = *cmd_info;
    collate.in_file_name_ = cmd_info->in_file_name_;
    collate.out_file_name_ = collate_path;
    collate.collate_memory_ = cmd_info->pipeline_memory_;
    collate.collate_temp_prefix_ =
        cmd_info->workflow_temp_prefix_ + ".collate-runs";
    if (WorkflowRunStage("collate", ProcessCollateMPI, &collate,
                         rank, &collate_wall) != 0) {
        return 1;
    }
    WorkflowResetCpeRuntime("collate", rank);

    CmdInfo fixmate = *cmd_info;
    fixmate.in_file_name_ = collate_path;
    fixmate.out_file_name_ = fixmate_path;
    fixmate.fixmate_mate_score_ = true;
    if (WorkflowRunStage("fixmate", ProcessFixmateMPI, &fixmate,
                         rank, &fixmate_wall) != 0) {
        return 1;
    }
    WorkflowResetCpeRuntime("fixmate", rank);
    if (!cmd_info->workflow_keep_intermediates_) {
        WorkflowRemoveIntermediate(collate_path, rank);
    }

    CmdInfo sort = *cmd_info;
    sort.in_file_name_ = fixmate_path;
    sort.out_file_name_ = sort_path;
    sort.sort_memory_ = cmd_info->pipeline_memory_;
    sort.sort_temp_prefix_ =
        cmd_info->workflow_temp_prefix_ + ".sort-runs";
    if (WorkflowRunStage("sort", ProcessSortMPI, &sort,
                         rank, &sort_wall) != 0) {
        return 1;
    }
    WorkflowResetCpeRuntime("sort", rank);
    if (!cmd_info->workflow_keep_intermediates_) {
        WorkflowRemoveIntermediate(fixmate_path, rank);
    }

    CmdInfo markdup = *cmd_info;
    markdup.in_file_name_ = sort_path;
    markdup.out_file_name_ = cmd_info->out_file_name_;
    markdup.markdup_memory_ = cmd_info->pipeline_memory_;
    markdup.markdup_remove_dups_ = false;
    markdup.markdup_clear_ = false;
    // Reuse the validated coordinate-stream implementation used by the
    // in-memory pipeline without enabling the newer one-pass experiment.
    markdup.markdup_streaming_ = true;
    markdup.markdup_force_two_pass_ = true;
    if (WorkflowRunStage("markdup", ProcessMarkdupMPI, &markdup,
                         rank, &markdup_wall) != 0) {
        return 1;
    }
    if (!cmd_info->workflow_keep_intermediates_) {
        WorkflowRemoveIntermediate(sort_path, rank);
    }

    const double workflow_wall =
        WorkflowReduceMax(GetTime() - workflow_t0);
    if (rank == 0) {
        printf("dedup-workflow complete wall=%lf collate=%lf "
               "fixmate=%lf sort=%lf markdup=%lf\n",
               workflow_wall, collate_wall, fixmate_wall,
               sort_wall, markdup_wall);
        if (cmd_info->workflow_keep_intermediates_) {
            printf("dedup-workflow intermediates: %s %s %s\n",
                   collate_path.c_str(), fixmate_path.c_str(),
                   sort_path.c_str());
        }
    }
    return 0;
}
