#ifndef SWBAM_APP_ALGORITHM_TYPES_H
#define SWBAM_APP_ALGORITHM_TYPES_H

#include "swbam/bam_types.h"

// Shared verbatim by the application's MPE code and algorithm CPE kernels.
struct MpiSortRecordMetaShared {
    int32_t tid;
    int32_t pos;
    uint16_t flag;
    uint16_t pad;
    uint32_t raw_len;
    uint32_t pad2;
    uint64_t raw_offset;
    uint64_t global_order;
};

struct MpiSortExtractPara {
    int block_id;
    bam_block *input_block;
    bam_block *un_comp_block;
    unsigned char *raw_arena;
    size_t raw_capacity;
    size_t raw_used;
    uint64_t raw_base_offset;
    MpiSortRecordMetaShared *records;
    int record_capacity;
    int n_records;
    long long global_block_index;
    int status;
    int record_index;
    long long actual_value;
    long long limit_value;
    int limit_id;
    uint64_t decomp_alloc_cycles;
    uint64_t decomp_inflate_cycles;
    uint64_t decomp_crc_cycles;
    uint64_t decomp_parse_cycles;
    uint64_t decomp_total_cycles;
};

struct MpiCollateRecordMetaShared {
    uint32_t hash;
    uint32_t bin;
    uint16_t qname_len;
    uint8_t flag_order;
    uint8_t pad;
    uint32_t raw_len;
    uint64_t raw_offset;
    uint64_t global_order;
};

struct MpiCollateExtractPara {
    int block_id;
    bam_block *input_block;
    bam_block *un_comp_block;
    unsigned char *raw_arena;
    size_t raw_capacity;
    size_t raw_used;
    uint64_t raw_base_offset;
    MpiCollateRecordMetaShared *records;
    int record_capacity;
    int n_records;
    int n_bins;
    long long global_block_index;
    int status;
    int record_index;
    long long actual_value;
    long long limit_value;
    int limit_id;
    uint64_t decomp_alloc_cycles;
    uint64_t decomp_inflate_cycles;
    uint64_t decomp_crc_cycles;
    uint64_t decomp_parse_cycles;
    uint64_t decomp_total_cycles;
};

struct MpiSortRawCompressPara {
    int block_id;
    bam_block *un_comp_block;
    int un_comp_size;
    bam_block *output_block;
    int output_size;
    int status;
    int compress_level;
    uint64_t compress_pack_cycles;
    uint64_t compress_alloc_cycles;
    uint64_t compress_deflate_cycles;
    uint64_t compress_footer_cycles;
    uint64_t compress_total_cycles;
};

struct MpiSortBucketPackPara {
    int core_id;
    const MpiSortRecordMetaShared *local_records;
    const unsigned char *local_raw;
    size_t local_raw_size;
    const int *bucket_ids;
    const uint64_t *record_meta_in_bucket;
    const uint64_t *record_raw_in_bucket;
    const long long *send_meta_displs;
    const long long *send_raw_displs;
    MpiSortRecordMetaShared *send_meta;
    unsigned char *send_raw;
    uint64_t send_meta_capacity;
    uint64_t send_raw_capacity;
    size_t record_begin;
    size_t record_end;
    int status;
    int record_index;
    long long actual_value;
    long long limit_value;
    int limit_id;
    uint64_t pack_cycles;
    uint64_t total_cycles;
};

struct MpiSortRangePackPara {
    int core_id;
    const MpiSortRecordMetaShared *local_records;
    const unsigned char *local_raw;
    size_t local_raw_size;
    const uint64_t *record_raw_offsets;
    MpiSortRecordMetaShared *send_meta;
    unsigned char *send_raw;
    uint64_t send_meta_capacity;
    uint64_t send_raw_capacity;
    size_t record_begin;
    size_t record_end;
    size_t output_record_begin;
    int status;
    int record_index;
    long long actual_value;
    long long limit_value;
    int limit_id;
    uint64_t pack_cycles;
    uint64_t total_cycles;
};

struct MpiMarkdupKeyShared {
    int64_t this_coord;
    int64_t other_coord;
    int32_t this_ref;
    int32_t other_ref;
    int8_t single;
    int8_t leftmost;
    int8_t orientation;
    int8_t pad;
};

struct MpiMarkdupCandidateShared {
    MpiMarkdupKeyShared key;
    uint64_t ordinal;
    uint64_t global_order;
    int64_t score;
    uint64_t qname_offset;
    uint16_t qname_len;
    uint8_t paired_marker;
    uint8_t qc_fail;
    int32_t source_rank;
};

struct MpiMarkdupDecisionShared {
    uint64_t ordinal;
    uint8_t duplicate;
    uint8_t pad[7];
};

struct MpiMarkdupExtractPara {
    int block_id;
    bam1_t **records;
    int n_records;
    uint64_t ordinal_base;
    uint64_t global_block_index;
    int source_rank;
    int include_fails;
    MpiMarkdupCandidateShared *candidates;
    int candidate_capacity;
    int n_candidates;
    unsigned char *qname_arena;
    size_t qname_capacity;
    size_t qname_used;
    long long examined;
    long long excluded;
    long long pair_candidates;
    long long single_candidates;
    int first_tid;
    int first_pos;
    int last_tid;
    int last_pos;
    int has_records;
    int status;
    int record_index;
    long long actual_value;
    long long limit_value;
    int limit_id;
};

struct MpiMarkdupRewritePara {
    int block_id;
    bam1_t **records;
    uint32_t *bam_lens;
    int n_records;
    uint64_t ordinal_base;
    const uint8_t *duplicate_bitmap;
    size_t duplicate_bitmap_bytes;
    int clear_old;
    int remove_dups;
    int n_kept_records;
    uint32_t kept_total_len;
    long long marked_records;
    long long cleared_records;
    long long removed_records;
    int status;
    int record_index;
};

struct MpiFixmateGroupShared {
    int32_t begin;
    int32_t count;
};

struct MpiFixmateRecordPlanShared {
    bam1_core_t core;
    uint32_t output_data_len;
    int32_t mate_index;
    uint8_t rewrite_tags;
    uint8_t update_mq;
    uint8_t update_mc;
    uint8_t update_ms;
};

struct MpiFixmatePlanPara {
    bam1_t **records;
    MpiFixmateRecordPlanShared *plans;
    const MpiFixmateGroupShared *groups;
    int group_begin;
    int group_end;
    long long paired_groups;
    long long singleton_groups;
    long long secondary_records;
    long long supplementary_records;
    long long mq_updates;
    long long mc_updates;
    long long ms_updates;
    int status;
    int record_index;
};

struct MpiFixmateRewritePara {
    bam1_t **records;
    bam1_t *output_records;
    const MpiFixmateRecordPlanShared *plans;
    const uint64_t *output_offsets;
    unsigned char *output_data;
    size_t output_capacity;
    int record_begin;
    int record_end;
    int status;
    int record_index;
};

#endif // SWBAM_APP_ALGORITHM_TYPES_H
