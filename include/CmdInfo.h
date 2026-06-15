
#ifndef CMDINFO_H
#define CMDINFO_H

#include "Globals.h"
#include <cstdint>
#include <string>


class CmdInfo {
public:
    CmdInfo();

public:
    std::string in_file_name_;   
    std::string out_file_name_;  
    bool verbose_;       
    bool validate_bounds_;
    int min_mapq_;
    int max_mapq_;
    uint32_t require_flag_;
    uint32_t exclude_flag_;
    std::string ref_name_;
    int min_read_len_;
    int max_read_len_;
    int compress_level_;
    std::string sort_memory_;
    std::string sort_temp_prefix_;
    int collate_bins_;
    std::string collate_memory_;
    std::string collate_temp_prefix_;
    bool markdup_remove_dups_;
    bool markdup_clear_;
    bool markdup_include_fails_;
    std::string markdup_memory_;
    bool fixmate_mate_score_;

};

#endif //CMDINFO_H
