#include "sunway/CLI11.hpp"
#include "sunway/CmdInfo.h"
#include "sunway/Globals.h"
#include "swbam.h"

int main(int argc, char **argv) {

    double ttt = GetTime();

    CmdInfo cmd_info;
    CLI::App app("SWBAM");

    app.add_option("-i,--inFile", cmd_info.in_file_name_, "input sam/bam name")->required();
    app.add_option("-o,--outFile", cmd_info.out_file_name_, "output sam/bam name")->required();
    app.add_flag("--verbose", cmd_info.verbose_, "Enable verbose logging")->default_val(false);

    CLI11_PARSE(app, argc, argv);

    if(cmd_info.verbose_){
        printf("Input file: %s\n", cmd_info.in_file_name_.c_str());
        printf("Output file: %s\n", cmd_info.out_file_name_.c_str());
    }

    SwBam *swbam = new SwBam(&cmd_info);
    swbam->ProcessSwBam();
    delete swbam;

    printf("TOT TIME %lf\n", GetTime() - ttt);

    return 0;
}