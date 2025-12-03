
#ifndef CMDINFO_H
#define CMDINFO_H

#include "Globals.h"
#include <string>


class CmdInfo {
public:
    CmdInfo();

public:
    std::string in_file_name_;   
    std::string out_file_name_;  
    bool verbose_;       
         
};

#endif //CMDINFO_H