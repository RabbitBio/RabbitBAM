#include "CLI11.hpp"
#include "CmdInfo.h"
#include "Globals.h"
#include "swbam_mpi.h"

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

#include <cstdio>
#include <mpi.h>
#include <string>
#include <vector>

int main(int argc, char **argv) {
    double t0 = GetTime();
    double t1 = 0.0;

    double mpi_init_t0 = GetTime();
    MPI_Init(&argc, &argv);
    double mpi_init_cost = GetTime() - mpi_init_t0;

    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

#ifndef PLATFORM_SUNWAY
    if (my_rank == 0) {
        fprintf(stderr, "RabbitBAM-MPI is only supported on the Sunway platform.\n");
    }
    MPI_Finalize();
    return 1;
#else
    athread_init();
#endif

    CmdInfo cmd_info;
    bool stats_basic = false;
    CLI::App app("RabbitBAM-MPI");

    CLI::App *run_all = app.add_subcommand("run_all", "Run MPI conversion (BAM -> BAM, BAM -> SAM, and SAM -> BAM implemented)");
    run_all->add_option("-i,--inFile", cmd_info.in_file_name_, "input sam/bam name")->required()->check(CLI::ExistingFile);
    run_all->add_option("-o,--outFile", cmd_info.out_file_name_, "output sam/bam name")->required();
    run_all->add_option("--io-backend", cmd_info.io_backend_,
                        "Input backend: memory, posix, mpiio, or auto")
        ->default_val("memory")
        ->check(CLI::IsMember(std::vector<std::string>{"memory", "posix", "mpiio", "auto"}));
    run_all->add_option("--io-memory-limit", cmd_info.io_memory_limit_,
                        "Per-rank memory input limit used by --io-backend auto")
        ->default_val("8G");
    run_all->add_option("--io-output-backend", cmd_info.io_output_backend_,
                        "Output backend: memory or mpiio")
        ->default_val("memory")
        ->check(CLI::IsMember(std::vector<std::string>{"memory", "mpiio"}));
    run_all->add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);
    run_all->add_flag("--validate-bounds", cmd_info.validate_bounds_, "Boundary validation is not supported by RabbitBAM-MPI")->default_val(false);
    run_all->add_option("--min-mapq", cmd_info.min_mapq_, "Keep reads with MAPQ >= value");
    run_all->add_option("--max-mapq", cmd_info.max_mapq_, "Keep reads with MAPQ <= value");
    run_all->add_option("--require-flag", cmd_info.require_flag_, "Keep reads whose FLAG contains all bits in value");
    run_all->add_option("--exclude-flag", cmd_info.exclude_flag_, "Drop reads whose FLAG contains any bit in value");
    run_all->add_option("--ref-name", cmd_info.ref_name_, "Keep reads mapped to the given reference name");
    run_all->add_option("--min-read-len", cmd_info.min_read_len_, "Keep reads with read length >= value");
    run_all->add_option("--max-read-len", cmd_info.max_read_len_, "Keep reads with read length <= value");
    run_all->add_option("--compress-level", cmd_info.compress_level_, "MPI BAM output compression level: 0, 1, or 6")->default_val(1);

    CLI::App *flagstat = app.add_subcommand("flagstat", "Run MPI BAM flagstat");
    flagstat->add_option("-i,--inFile", cmd_info.in_file_name_, "input bam name")->required()->check(CLI::ExistingFile);
    flagstat->add_option("--io-backend", cmd_info.io_backend_,
                         "Input backend: memory, posix, mpiio, or auto")
        ->default_val("memory")
        ->check(CLI::IsMember(std::vector<std::string>{"memory", "posix", "mpiio", "auto"}));
    flagstat->add_option("--io-memory-limit", cmd_info.io_memory_limit_,
                         "Per-rank memory input limit used by --io-backend auto")
        ->default_val("8G");
    flagstat->add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);

    CLI::App *stats = app.add_subcommand("stats", "Run MPI BAM stats");
    stats->add_option("-i,--inFile", cmd_info.in_file_name_, "input bam name")->required()->check(CLI::ExistingFile);
    stats->add_option("--io-backend", cmd_info.io_backend_,
                      "Input backend: memory, posix, mpiio, or auto")
        ->default_val("memory")
        ->check(CLI::IsMember(std::vector<std::string>{"memory", "posix", "mpiio", "auto"}));
    stats->add_option("--io-memory-limit", cmd_info.io_memory_limit_,
                      "Per-rank memory input limit used by --io-backend auto")
        ->default_val("8G");
    stats->add_flag("--basic", stats_basic, "Only output samtools stats SN summary numbers")->default_val(false);
    stats->add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);

    CLI::App *io_check = app.add_subcommand(
        "io-check", "Validate the generic MPI BGZF decode pipeline");
    io_check->add_option("-i,--inFile", cmd_info.in_file_name_,
                         "input bam name")->required()->check(CLI::ExistingFile);

    CLI::App *sort = app.add_subcommand("sort", "Run MPI BAM coordinate sort");
    sort->add_option("-i,--inFile", cmd_info.in_file_name_, "input bam name")->required()->check(CLI::ExistingFile);
    sort->add_option("-o,--outFile", cmd_info.out_file_name_, "output bam name")->required();
    sort->add_option("-m,--memory", cmd_info.sort_memory_, "Sort memory limit per MPI rank, e.g. 4G or 4096M");
    sort->add_option("-T,--temp-prefix", cmd_info.sort_temp_prefix_,
                     "External sort temporary prefix or directory");
    sort->add_option("--compress-level", cmd_info.compress_level_, "MPI BAM output compression level: 0, 1, or 6")->default_val(1);
    sort->add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);

    CLI::App *collate = app.add_subcommand("collate", "Run MPI BAM name collation");
    collate->add_option("-i,--inFile", cmd_info.in_file_name_, "input BAM name")->required()->check(CLI::ExistingFile);
    collate->add_option("-o,--outFile", cmd_info.out_file_name_, "output BAM name")->required();
    CLI::Option *collate_bins_option =
        collate->add_option("-n,--bins", cmd_info.collate_bins_, "Number of logical QNAME hash bins")->default_val(64);
    collate->add_option("-m,--memory", cmd_info.collate_memory_, "Collate memory limit per MPI rank, e.g. 4G or 4096M");
    collate->add_option("-T,--temp-prefix", cmd_info.collate_temp_prefix_, "External collate temporary prefix or directory");
    collate->add_option("--compress-level", cmd_info.compress_level_, "MPI BAM output compression level: 0, 1, or 6")->default_val(1);
    collate->add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);

    CLI::App *markdup = app.add_subcommand("markdup", "Run MPI BAM duplicate marking");
    markdup->add_option("-i,--inFile", cmd_info.in_file_name_, "coordinate-sorted fixmate BAM name")->required()->check(CLI::ExistingFile);
    markdup->add_option("-o,--outFile", cmd_info.out_file_name_, "output bam name")->required();
    markdup->add_flag("-r,--remove-dups", cmd_info.markdup_remove_dups_, "Remove duplicate records")->default_val(false);
    markdup->add_flag("-c,--clear", cmd_info.markdup_clear_, "Clear existing duplicate flags and dt/do tags first")->default_val(false);
    markdup->add_flag("--include-fails", cmd_info.markdup_include_fails_, "Include QC-fail records in duplicate detection")->default_val(false);
    markdup->add_flag("--stream", cmd_info.markdup_streaming_, "Use bounded-memory coordinate streaming markdup")->default_val(false);
    markdup->add_option("-l,--max-read-length", cmd_info.markdup_max_read_length_, "Expected maximum read length for streaming eviction")->default_val(300)->check(CLI::PositiveNumber);
    markdup->add_option("-m,--memory", cmd_info.markdup_memory_, "Markdup candidate memory limit per MPI rank, e.g. 4G or 4096M");
    markdup->add_option("--compress-level", cmd_info.compress_level_, "MPI BAM output compression level: 0, 1, or 6")->default_val(1);
    markdup->add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);

    CLI::App *fixmate = app.add_subcommand("fixmate", "Run MPI BAM fixmate");
    fixmate->add_option("-i,--inFile", cmd_info.in_file_name_, "name-collated or queryname-sorted BAM name")->required()->check(CLI::ExistingFile);
    fixmate->add_option("-o,--outFile", cmd_info.out_file_name_, "output bam name")->required();
    fixmate->add_flag("-m", cmd_info.fixmate_mate_score_, "Add mate score ms tags")->required();
    fixmate->add_option("--compress-level", cmd_info.compress_level_, "MPI BAM output compression level: 0, 1, or 6")->default_val(1);
    fixmate->add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);

    CLI::App *dedup_pipeline = app.add_subcommand("dedup-pipeline", "Run MPI BAM in-memory dedup pipeline: collate -> fixmate -m -> sort -> markdup");
    dedup_pipeline->add_option("-i,--inFile", cmd_info.in_file_name_, "input BAM name")->required()->check(CLI::ExistingFile);
    dedup_pipeline->add_option("-o,--outFile", cmd_info.out_file_name_, "final output BAM name")->required();
    CLI::Option *pipeline_bins_option =
        dedup_pipeline->add_option("-n,--bins", cmd_info.collate_bins_, "Number of logical QNAME hash bins for collate")->default_val(64);
    dedup_pipeline->add_option("-m,--memory", cmd_info.pipeline_memory_, "Pipeline memory limit per MPI rank, e.g. 8G or 8192M");
    dedup_pipeline->add_option("-l,--max-read-length", cmd_info.markdup_max_read_length_, "Expected maximum read length for streaming markdup eviction")->default_val(300)->check(CLI::PositiveNumber);
    dedup_pipeline->add_option("--compress-level", cmd_info.compress_level_, "MPI BAM output compression level: 0, 1, or 6")->default_val(1);
    dedup_pipeline->add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);

    CLI11_PARSE(app, argc, argv);

    int exit_code = 1;
    CLI::App *selected_command = nullptr;
    bool is_run_all = false;
    bool is_flagstat = false;
    bool is_stats = false;
    bool is_io_check = false;
    bool is_sort = false;
    bool is_collate = false;
    bool is_markdup = false;
    bool is_fixmate = false;
    bool is_dedup_pipeline = false;
    if (app.get_subcommands().empty()) {
        if (my_rank == 0) fprintf(stderr, "ERROR: You should input one command.\n");
        goto cleanup;
    }
    if (app.get_subcommands().size() > 1) {
        if (my_rank == 0) fprintf(stderr, "ERROR: You should input one command.\n");
        goto cleanup;
    }
    selected_command = app.get_subcommands()[0];
    is_run_all = selected_command->get_name() == "run_all";
    is_flagstat = selected_command->get_name() == "flagstat";
    is_stats = selected_command->get_name() == "stats";
    is_io_check = selected_command->get_name() == "io-check";
    is_sort = selected_command->get_name() == "sort";
    is_collate = selected_command->get_name() == "collate";
    is_markdup = selected_command->get_name() == "markdup";
    is_fixmate = selected_command->get_name() == "fixmate";
    is_dedup_pipeline = selected_command->get_name() == "dedup-pipeline";
    cmd_info.collate_bins_explicit_ =
        (is_collate && collate_bins_option &&
         collate_bins_option->count() > 0) ||
        (is_dedup_pipeline && pipeline_bins_option &&
         pipeline_bins_option->count() > 0);
    if (!is_run_all && !is_flagstat && !is_stats && !is_io_check && !is_sort && !is_collate && !is_markdup && !is_fixmate && !is_dedup_pipeline) {
        if (my_rank == 0) fprintf(stderr, "ERROR: RabbitBAM-MPI only supports run_all, flagstat, stats, io-check, sort, collate, markdup, fixmate, and dedup-pipeline.\n");
        goto cleanup;
    }
    if (is_stats && !stats_basic) {
        if (my_rank == 0) fprintf(stderr, "ERROR: RabbitBAM-MPI stats v1 requires --basic.\n");
        goto cleanup;
    }
    if ((is_run_all || is_sort || is_collate || is_markdup || is_fixmate || is_dedup_pipeline) &&
        cmd_info.compress_level_ != 0 &&
        cmd_info.compress_level_ != 1 &&
        cmd_info.compress_level_ != 6) {
        if (my_rank == 0) {
            fprintf(stderr, "ERROR: --compress-level only supports 0, 1, or 6 in RabbitBAM-MPI.\n");
        }
        goto cleanup;
    }

    t1 = GetTime();
    exit_code = is_io_check ? ProcessIoCheckMPI(&cmd_info)
                         : (is_dedup_pipeline ? ProcessDedupPipelineMPI(&cmd_info)
                         : (is_sort ? ProcessSortMPI(&cmd_info)
                         : (is_collate ? ProcessCollateMPI(&cmd_info)
                         : (is_markdup ? ProcessMarkdupMPI(&cmd_info)
                                      : (is_fixmate ? ProcessFixmateMPI(&cmd_info)
                                                    : (is_stats ? ProcessStatsMPI(&cmd_info)
                                                                : (is_flagstat ? ProcessFlagstatMPI(&cmd_info) : ProcessSwBamMPI(&cmd_info))))))));
    if (my_rank == 0) {
        printf("%s rank0 time is %lf--\n",
               is_io_check ? "ProcessIoCheckMPI" :
               (is_dedup_pipeline ? "ProcessDedupPipelineMPI" :
               (is_sort ? "ProcessSortMPI" :
               (is_collate ? "ProcessCollateMPI" :
                (is_markdup ? "ProcessMarkdupMPI" :
                (is_fixmate ? "ProcessFixmateMPI" :
                 (is_stats ? "ProcessStatsMPI" : (is_flagstat ? "ProcessFlagstatMPI" : "ProcessSwBamMPI"))))))),
               GetTime() - t1);
    }


cleanup:
#ifdef PLATFORM_SUNWAY
    athread_halt();
#endif

    double finalize_t0 = GetTime();
    MPI_Finalize();
    double finalize_cost = GetTime() - finalize_t0;

    if (my_rank == 0) {
        printf("MPI_Init rank0 time %lf\n", mpi_init_cost);
        printf("MPI_Finalize rank0 time %lf\n", finalize_cost);
        printf("TOT TIME %lf\n", GetTime() - t0);
    }
    return exit_code;
}
