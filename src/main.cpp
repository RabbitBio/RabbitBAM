#include "CLI11.hpp"
#include "CmdInfo.h"
#include "Globals.h"
#include "swbam.h"
#include "BamTools.h"

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#ifdef USE_SWLU
#include "swlu.h"
#endif
#endif

int main(int argc, char **argv) {


    double ttt = GetTime();


#ifdef PLATFORM_SUNWAY
    // Initialize athread for Sunway platform
    athread_init();
#endif


#define DEBUG
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

    CmdInfo cmd_info;
    CLI::App app("SWBAM");

    CLI::App *run_all = app.add_subcommand("run_all", "Run full function (sam/bam processing)");
    run_all->add_option("-i,--inFile", cmd_info.in_file_name_, "input sam/bam name")->required()->check(CLI::ExistingFile);
    run_all->add_option("-o,--outFile", cmd_info.out_file_name_, "output sam/bam name")->required();
    run_all->add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);

    CLI::App *check_cross = app.add_subcommand("check_cross_block",
        "Check whether a BAM file contains records that span across BGZF block boundaries. "
        "The bam2sam fast-path assumes no cross-block records; use this command to verify.");
    check_cross->add_option("-i,--inFile", cmd_info.in_file_name_, "input BAM file")->required()->check(CLI::ExistingFile);

    CLI::App *test_api = app.add_subcommand("test_api", "Test internal API functions");
    // test_api->add_option("-i,--inFile", cmd_info.in_file_name_, "input sam/bam name")->required()->check(CLI::ExistingFile);
    // test_api->add_option("-o,--outFile", cmd_info.out_file_name_, "output sam/bam name")->required();
    test_api->add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);

    CLI11_PARSE(app, argc, argv);

    if (app.get_subcommands().size() > 1) {
        fprintf(stderr, "ERROR: You should input one command!\n");
        return 0;
    }

    //完整程序-------------------------------------------------------------------------------------
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "run_all") == 0) {
        //printf("Input file: %s\n", cmd_info.in_file_name_.c_str());
        //printf("Output file: %s\n", cmd_info.out_file_name_.c_str());

        if(cmd_info.verbose_){
            printf("Enable verbose logging\n");
        }

        SwBam *swbam = new SwBam(&cmd_info);
        swbam->ProcessSwBam();
        delete swbam;


    }

    //检测 BAM 记录是否跨 BGZF 块----------------------------------------------------------------------
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "check_cross_block") == 0) {
        // ret == 0: 无跨块，安全；ret == 1: 有跨块，不安全；ret == -1: 文件打开/读取错误
        int ret = check_bam_cross_block(cmd_info.in_file_name_.c_str());
        if (ret < 0) {
            fprintf(stderr, "ERROR: check_bam_cross_block failed\n");
            return 1;  
        }
    
    }


    //测试部分接口函数-------------------------------------------------------------------------------
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "test_api") == 0) {
        // printf("Input file: %s\n", cmd_info.in_file_name_.c_str());
        // printf("Output file: %s\n", cmd_info.out_file_name_.c_str());

        // 探测当前环境下 malloc 的极限
        size_t GB = 1024LL * 1024LL * 1024LL;
        size_t test_size = 1 * GB; 

        while (true) {
            void* ptr = malloc(test_size);
            if (ptr) {
                printf("Successfully allocated: %zu GB\n", test_size / GB);
                free(ptr);
                test_size += 1 * GB;
            } else {
                printf("FAILED to allocate: %zu GB. Your limit is around here.\n", test_size / GB);
                break;
            }
        }

        if(cmd_info.verbose_){
            printf("Enable verbose logging\n");
        }

    }



#ifdef USE_SWLU
    swlu_prof_print();
#endif

#ifdef PLATFORM_SUNWAY
    // Halt athread for Sunway platform
    athread_halt();
#endif

    printf("TOT TIME %lf\n", GetTime() - ttt);

    return 0;
}