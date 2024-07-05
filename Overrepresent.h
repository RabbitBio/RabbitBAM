
//
// Created by 赵展 on 2021/4/12.
//

#ifndef BAMSTATUS_OVERREPRESENT_H
#define BAMSTATUS_OVERREPRESENT_H
#include "config.h"
#include <htslib/sam.h>
#include <algorithm>
using namespace std;
class Overrepresent {
public:
    Overrepresent(int Total,double Center=0.001);
    void insert(bam1_t* b);
    void status();
public:
    char **sequence;
    int Total;
    int Pos;
    double Center;
    double OverrepresentDate=-1;
};


#endif //BAMSTATUS_OVERREPRESENT_H
