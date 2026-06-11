#ifndef SWBAM_MPI_H
#define SWBAM_MPI_H

#include "swbam.h"

#include <cstddef>

#define ENABLE_MASKING

struct MpiBamToBamStats {
    long long input_blocks;
    long long group_count;
    long long total_records;
    long long kept_records;
    long long dropped_records;
    long long bgzf_blocks;
    long long pack_records;
    double t_decomp_filter;
    double t_decomp_alloc;
    double t_decomp_inflate;
    double t_decomp_crc;
    double t_decomp_parse;
    double t_decomp_other;
    double t_pack;
    double t_compress;
    double t_compress_serialize;
    double t_compress_alloc;
    double t_compress_deflate;
    double t_compress_footer;
    double t_compress_other;
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
    double t_decomp_alloc;
    double t_decomp_inflate;
    double t_decomp_crc;
    double t_decomp_parse;
    double t_decomp_other;
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
    long long parse_fast_records;
    long long parse_fallback_records;
    double t_split;
    double t_copy_count;
    double t_parse;
    double t_parse_core;
    double t_parse_aux;
    double t_parse_cg;
    double t_parse_fallback;
    double t_parse_other;
    double t_pack;
    double t_compress;
    double t_compress_serialize;
    double t_compress_alloc;
    double t_compress_deflate;
    double t_compress_footer;
    double t_compress_other;
    double t_write;
    double t_gather;
    double t_fused_total;
    double t_alloc_init;
    double t_setup_reset;
    double t_compress_setup;
    double t_status_check;
    double t_free_workspace;
};

struct MpiSortStats {
    long long sort_mode;
    long long input_blocks;
    long long local_records;
    long long received_records;
    long long sample_records;
    long long bgzf_blocks;
    long long external_runs;
    long long external_segments;
    long long bucket_self_records;
    long long bucket_remote_records;
    long long bucket_self_raw_bytes;
    long long bucket_remote_raw_bytes;
    long long temp_read_bytes;
    long long temp_write_bytes;
    long long tracked_peak_bytes;
    long long run_arena_bytes;
    long long merge_fan_in;
    long long consolidation_passes;
    long long cpe_calibration_cycles;
    double t_setup;
    double t_extract;
    double t_extract_read;
    double t_extract_read_unhidden;
    double t_extract_prepare;
    double t_extract_merge;
    double t_extract_alloc;
    double t_extract_inflate;
    double t_extract_crc;
    double t_extract_parse;
    double t_extract_other;
    double t_local_sort;
    double t_sample;
    double t_partition;
    double t_bucket_count;
    double t_bucket_pack;
    double t_exchange;
    double t_mpi_exchange;
    double t_mpi_simulated;
    double t_offset_fix;
    double t_final_sort;
    double t_temp_write_sim;
    double t_temp_read_sim;
    double t_temp_write_actual;
    double t_temp_read_actual;
    double t_temp_open_actual;
    double t_run_sort;
    double t_run_bucket;
    double t_run_exchange;
    double t_run_merge;
    double t_merge_unhidden;
    double t_merge_simulated;
    double t_merge_temp_read_sim;
    double t_consolidation_simulated;
    double t_consolidation_temp_read_sim;
    double t_consolidation_temp_write_sim;
    double t_compress;
    double t_compress_simulated;
    double t_compress_pack;
    double t_compress_alloc;
    double t_compress_deflate;
    double t_compress_footer;
    double t_compress_other;
    double t_write;
    double t_status_check;
    double t_cleanup;
    double t_fused_total;
    double t_fused_actual;
    double t_cpe_calibration_wall;
};

struct MpiMarkdupStats {
    long long input_blocks;
    long long group_count;
    long long total_records;
    long long examined_records;
    long long excluded_records;
    long long pair_candidates;
    long long single_candidates;
    long long owner_candidates;
    long long pair_duplicates;
    long long single_duplicates;
    long long marked_records;
    long long cleared_records;
    long long removed_records;
    long long qname_bytes;
    long long mpi_candidate_bytes;
    long long mpi_result_bytes;
    long long bgzf_blocks;
    long long tracked_peak_bytes;
    double t_candidate_decomp;
    double t_candidate_extract;
    double t_candidate_exchange;
    double t_group;
    double t_result_exchange;
    double t_rewrite_decomp;
    double t_rewrite;
    double t_pack;
    double t_compress;
    double t_read;
    double t_write;
    double t_fused_total;
};

struct MpiFixmateStats {
    long long input_blocks;
    long long group_count;
    long long total_records;
    long long paired_groups;
    long long singleton_groups;
    long long secondary_records;
    long long supplementary_records;
    long long boundary_groups;
    long long boundary_bytes;
    long long mq_updates;
    long long mc_updates;
    long long ms_updates;
    long long bgzf_blocks;
    double t_read;
    double t_decompress;
    double t_group_scan;
    double t_plan;
    double t_boundary_exchange;
    double t_rewrite;
    double t_pack;
    double t_compress;
    double t_write;
    double t_fused_total;
};

int ProcessSwBamMPI(CmdInfo *cmd_info);
int ProcessBamToBamMPI(CmdInfo *cmd_info);
int ProcessFlagstatMPI(CmdInfo *cmd_info);
int ProcessStatsMPI(CmdInfo *cmd_info);
int ProcessSortMPI(CmdInfo *cmd_info);
int ProcessMarkdupMPI(CmdInfo *cmd_info);
int ProcessFixmateMPI(CmdInfo *cmd_info);

int FusedBamToBamMPI(MemReader &reader, MemWriter &mem_writer,
                     const BamFilterOptions &filter,
                     int compress_level,
                     MpiBamToBamStats *stats);
int FusedBamToSamMPI(MemReader &reader, MemWriter &mem_writer,
                     sam_hdr_t *hdr,
                     MpiBamToSamStats *stats);
int FusedSamToBamMPI(MemReader &reader, MemWriter &mem_writer,
                     sam_hdr_t *hdr,
                     int compress_level,
                     MpiSamToBamStats *stats);
int FusedBamSortMPI(MemReader &reader, MemWriter &mem_writer,
                    long long global_block_begin,
                    int rank,
                    int comm_size,
                    int compress_level,
                    size_t memory_limit,
                    MpiSortStats *stats);
int FusedBamExternalSortMPI(MemReader &reader, MemWriter &mem_writer,
                            long long global_block_begin,
                            int rank,
                            int comm_size,
                            int compress_level,
                            size_t memory_limit,
                            const char *temp_prefix,
                            MpiSortStats *stats);

int MpiWriteBlockToMem(MemWriter &w, bam_block *block);
int MpiWriteBytesToMem(MemWriter &w, const char *data, size_t len);

#endif
