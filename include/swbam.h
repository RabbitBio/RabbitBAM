
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
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <algorithm>


// #include <htslib/sam.h>
// #include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <htslib/hfile.h>

#include "BamRead.h"
#include "BamComplete.h"
#include "BamWrite.h"
#include "BamWriteComplete.h"
#include "BamTools.h"
#include "CmdInfo.h"

class SwBam {
public:
    SwBam(CmdInfo *cmd_info);

    SwBam();

    ~SwBam();

    void ProcessSwBam();

private:
    void ConsumerSwBamTask(BamRead *read, BamComplete *complete);
    void ProducerSwBamTask(BGZF *fp, BamRead *read);
    void ProducerSwBamTask_memory(BGZF *fp, BamRead *read , char *bam_mem, size_t bam_size);
    int writeBam1_tToSam(samFile *fp, const sam_hdr_t *h, const bam1_t *b);

    void ConsumerSwBamTask2(BamWrite *write, BamWriteComplete *complete);
    void ProducerSwBamTask2(samFile *fp, BamWrite *write, sam_hdr_t *h);
    void ProducerSwBamTask2_parallel(samFile *fp, BamWrite *write, sam_hdr_t *h);
    void ProducerSwBamTask2_parallel_memory(BamWrite *write, sam_hdr_t *h , MemReader reader);
    void ProducerSwBamTask2_parallel_memory_OP( BamWrite *write, sam_hdr_t *h , MemReader reader);
    int writeBlockTobam(BGZF *fp, bam_block *block);

    static inline const char *get_sam_open_mode(const std::string &out_name) {
        size_t len = out_name.size();
        if (len >= 4 && out_name.substr(len - 4) == ".bam") {
            return "wb";   // BAM (BGZF compressed)
        } else {
            return "w";    // SAM text
        }
    }


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