#include "CLI11.hpp"
#include "CmdInfo.h"
#include "Globals.h"
#include "swbam.h"

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#ifdef USE_SWLU
#include "swlu.h"
#endif
#endif

int main(int argc, char **argv) {


#ifdef PLATFORM_SUNWAY
    // Initialize athread for Sunway platform
    athread_init();
#endif

#if 0
    volatile int aa = 0;
#ifdef DEBUG
    fprintf(stderr, "DEBUG aa: %d\n", aa);
#endif
    while(aa == 0) {
        
    }
#ifdef DEBUG
    fprintf(stderr, "DEBUG aa: %d\n", aa);
#endif
#endif

#ifdef USE_SWLU
    swlu_debug_init();
    swlu_prof_init();
#endif

    double ttt = GetTime();

    CmdInfo cmd_info;
    CLI::App app("SWBAM");

    CLI::App *run_all = app.add_subcommand("run_all", "Run full function (sam/bam processing)");
    run_all->add_option("-i,--inFile", cmd_info.in_file_name_, "input sam/bam name")->required()->check(CLI::ExistingFile);
    run_all->add_option("-o,--outFile", cmd_info.out_file_name_, "output sam/bam name")->required();
    run_all->add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);

    CLI::App *test_api = app.add_subcommand("test_api", "Test internal API functions");
    test_api->add_option("-i,--inFile", cmd_info.in_file_name_, "input sam/bam name")->required()->check(CLI::ExistingFile);
    test_api->add_option("-o,--outFile", cmd_info.out_file_name_, "output sam/bam name")->required();
    test_api->add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);

    CLI11_PARSE(app, argc, argv);

    if (app.get_subcommands().size() > 1) {
        fprintf(stderr, "ERROR: You should input one command!\n");
        return 0;
    }

    //完整程序-------------------------------------------------------------------------------------
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "run_all") == 0) {
        printf("Input file: %s\n", cmd_info.in_file_name_.c_str());
        printf("Output file: %s\n", cmd_info.out_file_name_.c_str());

        if(cmd_info.verbose_){
            printf("Enable verbose logging\n");
        }

        SwBam *swbam = new SwBam(&cmd_info);
        swbam->ProcessSwBam();
        delete swbam;


    }


    //测试部分接口函数-------------------------------------------------------------------------------
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "test_api") == 0) {
        printf("Input file: %s\n", cmd_info.in_file_name_.c_str());
        printf("Output file: %s\n", cmd_info.out_file_name_.c_str());

        if(cmd_info.verbose_){
            printf("Enable verbose logging\n");
        }


        
    }

    printf("TOT TIME %lf\n", GetTime() - ttt);


#ifdef USE_SWLU
    swlu_prof_print();
#endif

#ifdef PLATFORM_SUNWAY
    // Halt athread for Sunway platform
    athread_halt();
#endif

    return 0;
}