
#ifndef SWBAM_H
#define SWBAM_H

#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>
#include <functional>


// #include <htslib/sam.h>
// #include <htslib/hts.h>
#include "BamRead.h"
#include "BamComplete.h"
#include "BamWrite.h"
#include "BamWriteComplete.h"
#include "sunway/BamTools.h"
#include "sunway/CmdInfo.h"

class SwBam {
public:
    SwBam(CmdInfo *cmd_info);

    SwBam();

    ~SwBam();

    void ProcessSwBam();

private:
    void ConsumerSwBamTask(BamRead *read, BamComplete *complete);
    void ProducerSwBamTask(BGZF *fp, BamRead *read);
    //void WriteSwBamTask();
    int writeBam1_t(samFile *fp, const sam_hdr_t *h, const bam1_t *b);

    void ConsumerSwBamTask2(BamWrite *write, BamWriteComplete *complete);
    void ProducerSwBamTask2(samFile *fp, BamWrite *write);
    int SwBam:: writeBlockTobam(BGZF *fp, bam_block *block);

private:
    CmdInfo *cmd_info_;
    sam_hdr_t *hdr;
    samFile *sin;
    samFile *sout;

    //bam2sam------------------
    BamRead *read;
    BamComplete *complete;

    //sam2bam------------------
    BamWrite *write;
    BamWriteComplete *writeComplete;

};




#endif // SWBAM_H