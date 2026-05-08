#ifndef SWBAM_MPI_H
#define SWBAM_MPI_H

#include "swbam.h"

#include <cstddef>

struct MpiBamToBamStats {
    long long input_blocks;
    long long group_count;
    long long total_records;
    long long kept_records;
    long long dropped_records;
    long long bgzf_blocks;
    long long pack_records;
    double t_decomp_filter;
    double t_pack;
    double t_compress;
    double t_read;
    double t_write;
    double t_mpi_write;
    double t_fused_total;
    double t_alloc_init;
    double t_initial_read;
    double t_prepare_decomp;
    double t_decomp_check;
    double t_compress_prepare;
    double t_compress_check;
    double t_empty_group_path;
    double t_final_flush;
    double t_free_workspace;
};

struct MpiBamToSamStats {
    long long input_blocks;
    long long group_count;
    long long total_records;
    long long format_tiles;
    double t_decomp;
    double t_format;
    double t_collect;
    double t_read;
    double t_write;
    double t_gather;
    double t_fused_total;
    double t_alloc_init;
    double t_free_workspace;
};

struct MpiSamToBamStats {
    long long input_chunks;
    long long chunk_groups;
    long long total_records;
    long long compress_groups;
    long long bgzf_blocks;
    double t_split;
    double t_copy_count;
    double t_parse;
    double t_pack;
    double t_compress;
    double t_write;
    double t_gather;
    double t_fused_total;
    double t_alloc_init;
    double t_setup_reset;
    double t_compress_setup;
    double t_status_check;
    double t_free_workspace;
};

int ProcessSwBamMPI(CmdInfo *cmd_info);
int ProcessBamToBamMPI(CmdInfo *cmd_info);

int FusedBamToBamMPI(MemReader &reader, MemWriter &mem_writer,
                     const BamFilterOptions &filter,
                     MpiBamToBamStats *stats);
int FusedBamToSamMPI(MemReader &reader, MemWriter &mem_writer,
                     sam_hdr_t *hdr,
                     MpiBamToSamStats *stats);
int FusedSamToBamMPI(MemReader &reader, MemWriter &mem_writer,
                     sam_hdr_t *hdr,
                     MpiSamToBamStats *stats);

int MpiWriteBlockToMem(MemWriter &w, bam_block *block);
int MpiWriteBytesToMem(MemWriter &w, const char *data, size_t len);

#endif
