
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
         
};

#endif //CMDINFO_H
