#ifndef BAMSATUS_BAMSTATUS_H
#define BAMSATUS_BAMSTATUS_H

#include <htslib/sam.h>
#include <cstdint>
#include <fstream>
#include <set>
#include "config.h"
#include "Duplicate.h"
#include "Overrepresent.h"

using namespace std;

//const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
class BamStatus {
public:
    BamStatus();

    BamStatus(string filename);

    ~BamStatus();

    void statusbam(bam1_t *b);

    void statusAll();

    void contentstatus();

    void add(BamStatus *b);

    void print();

    void reportHTML(ofstream *fout);

    void reportHTML(ofstream *fout, Duplicate *duplicate, Overrepresent *overrepresent);

    void reportHTML(ofstream *fout, Duplicate *duplicate, Overrepresent *overrepresent, sam_hdr_t *hdr);

public:
    int **NumberList;
    int **Qualitylist;
    int **QualityPositonList;
    int **Content;
    double **DoubleContent;
    int *LengthSequence;
    int *QualitySequence;
    int **Kmer;
    int ChooseKmerNum;
    int *ChooseKmerPos;
    int *ChooseKmerKey;
    int KmerBase = 5;
    int KmerBit;
    int ChromosomeNumber = 2000;
    int *Chromosome;
    int total_number = 0;
    int total_aligen_number = 0;
    int max_len = 0;
    int min_len = -1;

    int ContentLen;
    int QulityLen;
    set<int> Chr;
    string filename;
};


#endif //BAMSATUS_BAMSTATUS_H
