#ifndef SWBAM_MPI_H
#define SWBAM_MPI_H

#include "swbam.h"

#include <cstddef>
#include <string>
#include <vector>

#define ENABLE_MASKING

struct MpiMemoryBam {
    char *data;
    size_t size;

    MpiMemoryBam() : data(nullptr), size(0) {}
};

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
    long long resident_run_records;
    long long resident_run_raw_bytes;
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

struct MpiCollateStats {
    long long mode;
    long long input_blocks;
    long long total_records;
    long long received_records;
    long long qname_groups;
    long long hash_collision_groups;
    long long bgzf_blocks;
    long long runs;
    long long segments;
    long long resident_run_records;
    long long resident_run_raw_bytes;
    long long temp_read_bytes;
    long long temp_write_bytes;
    long long tracked_peak_bytes;
    long long run_arena_bytes;
    long long merge_fan_in;
    long long consolidation_passes;
    long long exchange_self_records;
    long long exchange_remote_records;
    long long exchange_self_raw_bytes;
    long long exchange_remote_raw_bytes;
    long long exchange_chunks;
    double t_extract;
    double t_extract_alloc;
    double t_extract_read;
    double t_extract_raw_resize;
    double t_extract_setup;
    double t_extract_status;
    double t_extract_merge;
    double t_extract_free;
    double t_local_sort;
    double t_exchange;
    double t_mpi;
    double t_exchange_bound;
    double t_exchange_self_pack;
    double t_exchange_self_append;
    double t_exchange_remote_pack;
    double t_exchange_recv_append;
    double t_memory_track;
    double t_merge;
    double t_compress;
    double t_compress_fill;
    double t_write;
    double t_temp_read_actual;
    double t_temp_write_actual;
    double t_temp_read_sim;
    double t_temp_write_sim;
    double t_fused;
    double t_actual;
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
    double t_extract_status;
    double t_extract_merge;
    double t_extract_memcheck;
    double t_candidate_exchange;
    double t_candidate_pack;
    double t_candidate_exchange_setup;
    double t_candidate_exchange_mpi;
    double t_candidate_cleanup;
    double t_boundary_check;
    double t_group_prepare;
    double t_group_sort;
    double t_group_scan;
    double t_group;
    double t_flat_alloc;
    double t_flat_init;
    double t_flat_probe;
    double t_flat_finalize;
    double t_flat_free;
    double t_owner_cleanup;
    double t_result_exchange;
    double t_bitmap_merge;
    double t_result_cleanup;
    double t_rewrite_decomp;
    double t_rewrite;
    double t_pack;
    double t_compress;
    double t_read;
    double t_write;
    double t_rank_sync;
    double t_cpe_launch;
    double t_cpe_sync;
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
    double t_record_collect;
    double t_group_scan;
    double t_process_setup;
    double t_process_alloc;
    double t_plan;
    double t_output_index;
    double t_boundary_exchange;
    double t_rewrite;
    double t_pack;
    double t_compress_setup;
    double t_compress;
    double t_write;
    double t_rank_sync;
    double t_workspace_free;
    double t_fused_total;
};

int ProcessSwBamMPI(CmdInfo *cmd_info);
int ProcessBamToBamMPI(CmdInfo *cmd_info);
int ProcessFlagstatMPI(CmdInfo *cmd_info);
int ProcessStatsMPI(CmdInfo *cmd_info);
int ProcessSortMPI(CmdInfo *cmd_info);
int ProcessMarkdupMPI(CmdInfo *cmd_info);
int ProcessFixmateMPI(CmdInfo *cmd_info);
int ProcessCollateMPI(CmdInfo *cmd_info);
int ProcessDedupPipelineMPI(CmdInfo *cmd_info);

void MpiMemoryBamFree(MpiMemoryBam *bam);
int MpiBroadcastMemoryBam(MpiMemoryBam *bam, int root);
int MpiCollateMemoryToMemory(CmdInfo *cmd_info,
                             const char *input_memory,
                             size_t input_size,
                             MpiMemoryBam *output_bam,
                             double *core_cost);
int MpiFixmateMemoryToMemory(CmdInfo *cmd_info,
                             const char *input_memory,
                             size_t input_size,
                             MpiMemoryBam *output_bam,
                             double *core_cost);
int MpiSortMemoryToMemory(CmdInfo *cmd_info,
                          const char *input_memory,
                          size_t input_size,
                          MpiMemoryBam *output_bam,
                          double *core_cost);
int MpiMarkdupMemoryToMemory(CmdInfo *cmd_info,
                             const char *input_memory,
                             size_t input_size,
                             MpiMemoryBam *output_bam,
                             double *core_cost);

struct MpiMarkdupStreamingWindow;
MpiMarkdupStreamingWindow *MpiMarkdupStreamingWindowCreate(
    int comm_size, int max_read_length,
    int range_first_tid, int range_first_pos,
    int range_last_tid, int range_last_pos);
void MpiMarkdupStreamingWindowDestroy(
    MpiMarkdupStreamingWindow *window);
int MpiMarkdupStreamingWindowProcess(
    MpiMarkdupStreamingWindow *window,
    const std::vector<MpiMarkdupCandidateShared> &owner_candidates,
    const std::vector<unsigned char> &owner_qnames,
    int progress_tid, int progress_pos,
    std::vector<std::vector<uint64_t> > *duplicates_by_source,
    MpiMarkdupStats *stats);
size_t MpiMarkdupStreamingWindowMemory(
    const MpiMarkdupStreamingWindow *window);

int MpiCommonLoadFileToMemory(const std::string &path,
                              char **data, size_t *size);
int MpiCommonScanBgzfBlocksInMemory(
    const char *base, size_t size, long long body_start,
    std::vector<long long> *offsets,
    std::vector<long long> *lengths);
int MpiCommonSelectBlockRangeFromMemory(
    char *base, size_t input_size,
    const std::vector<long long> &offsets,
    const std::vector<long long> &lengths,
    long long begin, long long end,
    char **data, size_t *size);
int MpiCommonDumpMemoryToFile(const std::string &path,
                              const char *data, size_t size);
int MpiCommonInitMemWriter(MemWriter &writer, size_t capacity);
int MpiCommonBuildBamHeaderMemory(sam_hdr_t *header,
                                  int compress_level,
                                  char **data, size_t *size);
int MpiCommonReadBamHeaderFromMemory(const char *data, size_t size,
                                     sam_hdr_t **header,
                                     long long *body_start);

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
