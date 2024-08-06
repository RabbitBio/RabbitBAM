#ifndef BAMSTATUS_BAMREADER_H
#define BAMSTATUS_BAMREADER_H

#include "BamTools.h"
#include "BamRead.h"
#include "BamCompleteBlock.h"
#include "BamCompress.h"
#include "Buffer.h"
#include <thread>


class BamReader {

public:

    BamReader(std::string file_name, int n_thread = 8, bool is_tgs = false);

    BamReader(std::string file_name, int read_block, int compress_block, int compress_complete_block, int n_thread = 8, bool is_tgs = false);

    ~BamReader() {
        read_thread->join();
        for (int i = 0; i < n_thread; i++) compress_thread[i]->join();
        assign_thread->join();
    }

    sam_hdr_t *getHeader();

    bool getBam1_t(bam1_t *b);

    std::vector<bam1_t *> getBam1_t_parallel(std::vector<bam1_t *> b_vec[THREAD_NUM_P]);

    bam_complete_block *getBamCompleteClock();

    void backBamCompleteBlock(bam_complete_block *un_comp);


private:
    BamRead *read;
    BamCompress *compress;
    BamCompleteBlock *completeBlock;
    samFile *sin;
    sam_hdr_t *hdr;
    samFile *output;
    std::thread *read_thread;
    std::thread **compress_thread;
    std::thread *assign_thread;
    bam_complete_block *un_comp;
    int n_thread;
};


#endif //BAMSTATUS_BAMREADER_H
