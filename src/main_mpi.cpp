#include "CLI11.hpp"
#include "CmdInfo.h"
#include "Globals.h"
#include "swbam_mpi.h"

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

#include <cstdio>
#include <mpi.h>

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
    CLI::App app("RabbitBAM-MPI");

    CLI::App *run_all = app.add_subcommand("run_all", "Run MPI conversion (BAM -> BAM, BAM -> SAM, and SAM -> BAM implemented)");
    run_all->add_option("-i,--inFile", cmd_info.in_file_name_, "input sam/bam name")->required()->check(CLI::ExistingFile);
    run_all->add_option("-o,--outFile", cmd_info.out_file_name_, "output sam/bam name")->required();
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

    CLI11_PARSE(app, argc, argv);

    int exit_code = 1;
    if (app.get_subcommands().empty()) {
        if (my_rank == 0) fprintf(stderr, "ERROR: You should input one command.\n");
        goto cleanup;
    }
    if (app.get_subcommands().size() > 1) {
        if (my_rank == 0) fprintf(stderr, "ERROR: You should input one command.\n");
        goto cleanup;
    }
    if (app.get_subcommands()[0]->get_name() != "run_all") {
        if (my_rank == 0) fprintf(stderr, "ERROR: RabbitBAM-MPI only supports run_all.\n");
        goto cleanup;
    }
    if (cmd_info.compress_level_ != 0 &&
        cmd_info.compress_level_ != 1 &&
        cmd_info.compress_level_ != 6) {
        if (my_rank == 0) {
            fprintf(stderr, "ERROR: --compress-level only supports 0, 1, or 6 in RabbitBAM-MPI.\n");
        }
        goto cleanup;
    }

    t1 = GetTime();
    exit_code = ProcessSwBamMPI(&cmd_info);
    if (my_rank == 0) {
        printf("ProcessSwBamMPI rank0 time is %lf--\n", GetTime() - t1);
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
