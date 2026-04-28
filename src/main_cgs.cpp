#include "CLI11.hpp"
#include "CmdInfo.h"
#include "Globals.h"
#include "swbam.h"

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

#include <cstring>
#include <cstdio>

int main(int argc, char **argv) {
    double t0 = GetTime();

#ifndef PLATFORM_SUNWAY
    fprintf(stderr, "RabbitBAM-CGS is only supported on the Sunway platform.\n");
    return 1;
#else
    athread_init_cgs();
#endif

    CmdInfo cmd_info;
    CLI::App app("RabbitBAM-CGS");

    CLI::App *run_all = app.add_subcommand("run_all", "Run BAM/SAM conversion in CGS mode");
    run_all->add_option("-i,--inFile", cmd_info.in_file_name_, "input SAM/BAM name")->required()->check(CLI::ExistingFile);
    run_all->add_option("-o,--outFile", cmd_info.out_file_name_, "output SAM/BAM name")->required();
    run_all->add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);
    run_all->add_flag("--validate-bounds", cmd_info.validate_bounds_, "Boundary validation is not supported by RabbitBAM-CGS")->default_val(false);
    run_all->add_option("--min-mapq", cmd_info.min_mapq_, "Keep reads with MAPQ >= value");
    run_all->add_option("--max-mapq", cmd_info.max_mapq_, "Keep reads with MAPQ <= value");
    run_all->add_option("--require-flag", cmd_info.require_flag_, "Keep reads whose FLAG contains all bits in value");
    run_all->add_option("--exclude-flag", cmd_info.exclude_flag_, "Drop reads whose FLAG contains any bit in value");
    run_all->add_option("--ref-name", cmd_info.ref_name_, "Keep reads mapped to the given reference name");
    run_all->add_option("--min-read-len", cmd_info.min_read_len_, "Keep reads with read length >= value");
    run_all->add_option("--max-read-len", cmd_info.max_read_len_, "Keep reads with read length <= value");

    CLI11_PARSE(app, argc, argv);

    int exit_code = 1;
    if (app.get_subcommands().empty()) {
        fprintf(stderr, "ERROR: You should input one command.\n");
        goto cleanup;
    }
    if (app.get_subcommands().size() > 1) {
        fprintf(stderr, "ERROR: You should input one command.\n");
        goto cleanup;
    }
    if (app.get_subcommands()[0]->get_name() != "run_all") {
        fprintf(stderr, "ERROR: RabbitBAM-CGS only supports run_all.\n");
        goto cleanup;
    }

    exit_code = ProcessSwBamCGS(&cmd_info);

cleanup:
#ifdef PLATFORM_SUNWAY
    athread_halt();
#endif
    printf("TOT TIME %lf\n", GetTime() - t0);
    return exit_code;
}
