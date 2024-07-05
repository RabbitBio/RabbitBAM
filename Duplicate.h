//
// Created by 赵展 on 2021/4/7.
//

#ifndef BAMSTATUS_DUPLICATE_H
#define BAMSTATUS_DUPLICATE_H

#include <htslib/sam.h>
#include <cstring>
#include <memory>
#include <cmath>
#include <mutex>
#include "config.h"
using namespace std;

class Duplicate {
public:
    Duplicate();
    ~Duplicate();
    void statusSeq(bam1_t *b);
    uint64 seq2int(uint8 *b,int start,int keybit,bool &valid);
    void addRecord(uint32 key,unsigned long long kmer32,int gc);
    double statAll(int* hist, double* meanGC, int histSize);

    void add(Duplicate* b);
public:
    int *Counts;
    unsigned long long *Dup;
    int *GC;
    int CountNum;
    int KeyBase;
    int KeyBit;
    mutex mtx;
};


#endif //BAMSTATUS_DUPLICATE_H
