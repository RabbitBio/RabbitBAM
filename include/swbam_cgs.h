#ifndef SWBAM_CGS_H
#define SWBAM_CGS_H

#include "swbam.h"

#include <cstddef>
#include <vector>

bool CgsIsBamLikeFormat(int format);
bool CgsIsSamLikeFormat(int format);
bool CgsHasBamFilterRequest(const CmdInfo *cmd_info);
BamFilterOptions CgsBuildBamFilterOptions(const CmdInfo *cmd_info);
const char *CgsOutputOpenMode(const std::string &out_name);

void InitMemWriterCGS(MemWriter &w, size_t cap);
int EnsureMemWriterCapacityCGS(MemWriter &w, size_t add);
int WriteBlockToMemCGS(MemWriter &w, bam_block *block);
int WriteBytesToMemCGS(MemWriter &w, const char *data, size_t len);
int MemReadBlockCGS(char *base, size_t size, size_t &pos, bam_block *block);
bool MemGetlineCGS(MemReader &r, kstring_t *ks);
void InitEmptyCompParaCGS(Comp_Para *para, int block_id);
int CheckCgsResources();

class CgsCrossAllocator {
public:
    CgsCrossAllocator();
    ~CgsCrossAllocator();
    void *allocate(size_t alignment, size_t size);

private:
    std::vector<void *> allocations_;
};

struct CgsBlockSet {
    bam_block *blocks;
    unsigned char *data;
    int n;
};

struct CgsRecordSet {
    bam1_t *records;
    unsigned char *data;
    bam1_t **ptrs;
    uint32_t *bam_lens;
    int capacity;
};

struct CgsPackedPlan {
    bam1_t **records;
    int n_records;
    uint32_t total_len;
};

struct CgsPackWorkspace {
    bam1_t **records;
    CgsPackedPlan *plans;
    int capacity;
    int max_plans;
    int active_blocks;
    int total_records;
    bam1_t **current_begin;
    int current_records;
    uint32_t current_len;
};

int AllocBlockSet(CgsCrossAllocator &alloc, CgsBlockSet *set, int n);
int AllocRecordSet(CgsCrossAllocator &alloc, CgsRecordSet *set, int total_records);
void ResetRecordPointersForBlock(CgsRecordSet *set, int block_id, int records_per_block);
void ResetRecordRange(CgsRecordSet *set, int begin, int count);
int AllocPackWorkspace(CgsCrossAllocator &alloc, CgsPackWorkspace *workspace, int capacity, int n_plans);
void ResetPackWorkspace(CgsPackWorkspace *workspace);
int SealCurrentPackBlock(CgsPackWorkspace *workspace);
int AppendRecordToPackWorkspace(CgsPackWorkspace *workspace, bam1_t *record, uint32_t bam_len);

int FusedBamToBamCGS(MemReader &reader, MemWriter &mem_writer, const BamFilterOptions &filter);
int FusedBamToSamCGS(MemReader &reader, MemWriter &mem_writer, sam_hdr_t *hdr);
int FusedSamToBamCGS(MemReader reader, MemWriter &mem_writer, sam_hdr_t *hdr);

#endif
