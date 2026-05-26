#include "swbam_mpi.h"

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

#include <mpi.h>

extern "C" {
    void slave_mpi_decompress_bam2bam_passthrough();
}

namespace {

const int kFlagstatNB = 64;

const unsigned char kFlagstatBgzfEofBlock[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
};

enum FlagstatCounterId {
    FLAGSTAT_TOTAL = 0,
    FLAGSTAT_PRIMARY,
    FLAGSTAT_SECONDARY,
    FLAGSTAT_SUPPLEMENTARY,
    FLAGSTAT_DUPLICATES,
    FLAGSTAT_PRIMARY_DUPLICATES,
    FLAGSTAT_MAPPED,
    FLAGSTAT_PRIMARY_MAPPED,
    FLAGSTAT_PAIRED,
    FLAGSTAT_READ1,
    FLAGSTAT_READ2,
    FLAGSTAT_PROPERLY_PAIRED,
    FLAGSTAT_PAIR_MAPPED,
    FLAGSTAT_SINGLETONS,
    FLAGSTAT_DIFF_CHR,
    FLAGSTAT_DIFF_CHR_MAPQ5,
    FLAGSTAT_COUNTER_COUNT
};

struct MpiFlagstatCounts {
    long long values[FLAGSTAT_COUNTER_COUNT][2];
};

struct MpiFlagstatStats {
    long long input_blocks;
    long long group_count;
    long long total_records;
    double t_read;
    double t_decomp;
    double t_decomp_alloc;
    double t_decomp_inflate;
    double t_decomp_crc;
    double t_decomp_parse;
    double t_decomp_other;
    double t_count;
    double t_fused_total;
};

struct MpiFlagstatBlockSet {
    bam_block *blocks;
    unsigned char *data;
    int n;
};

struct MpiFlagstatRecordSet {
    bam1_t *records;
    unsigned char *data;
    bam1_t **ptrs;
    uint32_t *bam_lens;
    int capacity;
};

bool MpiFlagstatIsBamLikeFormat(int format) {
    return format == bam || format == binary_format;
}

bool MpiFlagstatIsSamLikeFormat(int format) {
    return format == sam || format == text_format;
}

int MpiFlagstatNormalizeFormat(int format) {
    if (MpiFlagstatIsBamLikeFormat(format)) return bam;
    if (MpiFlagstatIsSamLikeFormat(format)) return sam;
    return format;
}

int MpiFlagstatAllRanksOk(int local_ok) {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return global_ok;
}

double MpiFlagstatReduceMaxCost(double local_cost) {
    double max_cost = 0.0;
    MPI_Reduce(&local_cost, &max_cost, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return max_cost;
}

int MpiFlagstatLoadFileToMemory(const std::string &path, char **data, size_t *size) {
    *data = nullptr;
    *size = 0;

    FILE *fp = fopen(path.c_str(), "rb");
    if (!fp) return -1;
    if (fseeko(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        return -1;
    }

    off_t end = ftello(fp);
    if (end < 0) {
        fclose(fp);
        return -1;
    }
    if (fseeko(fp, 0, SEEK_SET) != 0) {
        fclose(fp);
        return -1;
    }
    if ((unsigned long long)end > (unsigned long long)SIZE_MAX) {
        fclose(fp);
        return -1;
    }

    *size = (size_t)end;
    *data = *size ? (char *)malloc(*size) : nullptr;
    if (*size > 0 && !*data) {
        fclose(fp);
        return -1;
    }
    if (*size > 0 && fread(*data, 1, *size, fp) != *size) {
        fclose(fp);
        free(*data);
        *data = nullptr;
        *size = 0;
        return -1;
    }

    fclose(fp);
    return 0;
}

int MpiFlagstatScanBgzfBlocksInMemory(const char *base, size_t size, long long body_start,
                                      std::vector<long long> *offsets,
                                      std::vector<long long> *lengths) {
    if (!base || body_start < 0 ||
        (unsigned long long)body_start > (unsigned long long)size) {
        return -1;
    }

    long long pos = body_start;
    while ((unsigned long long)pos < (unsigned long long)size) {
        if ((unsigned long long)pos + BLOCK_HEADER_LENGTH > (unsigned long long)size) return -1;
        const unsigned char *header = (const unsigned char *)(base + pos);
        int block_len = (int)header[16] | ((int)header[17] << 8);
        block_len += 1;
        if (block_len <= 0 ||
            (unsigned long long)pos + (unsigned long long)block_len > (unsigned long long)size) {
            return -1;
        }

        bool is_eof = block_len == (int)sizeof(kFlagstatBgzfEofBlock) &&
                      memcmp(base + pos, kFlagstatBgzfEofBlock, sizeof(kFlagstatBgzfEofBlock)) == 0;
        if (is_eof) break;

        offsets->push_back(pos);
        lengths->push_back(block_len);
        pos += block_len;
    }

    return 0;
}

int MpiFlagstatSelectBlockRangeFromMemory(char *base, size_t input_size,
                                          const std::vector<long long> &offsets,
                                          const std::vector<long long> &lengths,
                                          long long begin,
                                          long long end,
                                          char **data,
                                          size_t *size) {
    *data = nullptr;
    *size = 0;
    if (begin >= end) return 0;

    long long start = offsets[(size_t)begin];
    long long stop = offsets[(size_t)(end - 1)] + lengths[(size_t)(end - 1)];
    if (start < 0 || stop < start ||
        (unsigned long long)stop > (unsigned long long)input_size) {
        return -1;
    }
    long long total = stop - start;
    if ((unsigned long long)total > (unsigned long long)SIZE_MAX) return -1;

    *data = base + start;
    *size = (size_t)total;
    return 0;
}

int MpiFlagstatAllocateBlockSet(MpiFlagstatBlockSet *set, int n) {
    set->n = n;
    set->blocks = (bam_block *)aligned_alloc_custom(64, (size_t)n * sizeof(bam_block));
    set->data = aligned_alloc_custom(64, (size_t)n * BGZF_MAX_BLOCK_SIZE);
    if (!set->blocks || !set->data) return -1;
    memset(set->blocks, 0, (size_t)n * sizeof(bam_block));
    for (int i = 0; i < n; ++i) {
        set->blocks[i].data = set->data + (size_t)i * BGZF_MAX_BLOCK_SIZE;
        set->blocks[i].block_id = i;
    }
    return 0;
}

void MpiFlagstatFreeBlockSet(MpiFlagstatBlockSet *set) {
    if (set->blocks) aligned_free_custom((unsigned char *)set->blocks);
    if (set->data) aligned_free_custom(set->data);
    set->blocks = nullptr;
    set->data = nullptr;
    set->n = 0;
}

int MpiFlagstatAllocateRecordSet(MpiFlagstatRecordSet *set, int total_records) {
    set->capacity = total_records;
    set->records = (bam1_t *)aligned_alloc_custom(64, (size_t)total_records * sizeof(bam1_t));
    set->data = aligned_alloc_custom(64, (size_t)total_records * INIT_DATA_SIZE);
    set->ptrs = (bam1_t **)aligned_alloc_custom(64, (size_t)total_records * sizeof(bam1_t *));
    set->bam_lens = (uint32_t *)aligned_alloc_custom(64, (size_t)total_records * sizeof(uint32_t));
    if (!set->records || !set->data || !set->ptrs || !set->bam_lens) return -1;
    memset(set->records, 0, (size_t)total_records * sizeof(bam1_t));
    memset(set->bam_lens, 0, (size_t)total_records * sizeof(uint32_t));
    for (int i = 0; i < total_records; ++i) {
        set->records[i].data = set->data + (size_t)i * INIT_DATA_SIZE;
        set->records[i].m_data = INIT_DATA_SIZE;
        set->records[i].l_data = 0;
        set->records[i].mempolicy = BAM_USER_OWNS_DATA;
        set->ptrs[i] = &set->records[i];
    }
    return 0;
}

void MpiFlagstatFreeRecordSet(MpiFlagstatRecordSet *set) {
    if (set->records) aligned_free_custom((unsigned char *)set->records);
    if (set->data) aligned_free_custom(set->data);
    if (set->ptrs) aligned_free_custom((unsigned char *)set->ptrs);
    if (set->bam_lens) aligned_free_custom((unsigned char *)set->bam_lens);
    set->records = nullptr;
    set->data = nullptr;
    set->ptrs = nullptr;
    set->bam_lens = nullptr;
    set->capacity = 0;
}

int MpiFlagstatMemReadBlock(char *base, size_t size, size_t &pos, bam_block *block) {
    if (pos >= size) return -1;
    if (pos + BLOCK_HEADER_LENGTH > size) return -1;
    int bsize = (int)(unsigned char)base[pos + 16] |
                ((int)(unsigned char)base[pos + 17] << 8);
    bsize += 1;
    if (bsize <= 0 || pos + (size_t)bsize > size) return -1;
    memcpy(block->data, base + pos, (size_t)bsize);
    block->length = bsize;
    block->pos = 0;
    block->errcode = 0;
    block->block_address = (int64_t)pos;
    pos += (size_t)bsize;
    return bsize;
}

BamFilterOptions MpiFlagstatNoFilter() {
    BamFilterOptions filter;
    filter.min_mapq = -1;
    filter.max_mapq = -1;
    filter.require_flag = 0;
    filter.exclude_flag = 0;
    filter.ref_tid = -2;
    filter.min_read_len = -1;
    filter.max_read_len = -1;
    return filter;
}

void MpiFlagstatAccumulateDecompDetail(const Bam2BamPara *paras,
                                       int active_blocks,
                                       double decomp_wall,
                                       MpiFlagstatStats *stats) {
    if (!stats || active_blocks <= 0 || decomp_wall <= 0.0) return;

    const Bam2BamPara *critical = nullptr;
    uint64_t critical_total = 0;
    for (int b = 0; b < active_blocks; ++b) {
        if (paras[b].decomp_total_cycles >= critical_total) {
            critical_total = paras[b].decomp_total_cycles;
            critical = &paras[b];
        }
    }
    if (!critical || critical_total == 0) {
        stats->t_decomp_other += decomp_wall;
        return;
    }

    const double scale = decomp_wall / (double)critical_total;
    const double alloc_time = scale * (double)critical->decomp_alloc_cycles;
    const double inflate_time = scale * (double)critical->decomp_inflate_cycles;
    const double crc_time = scale * (double)critical->decomp_crc_cycles;
    const double parse_time = scale * (double)critical->decomp_parse_cycles;
    double other_time = decomp_wall - alloc_time - inflate_time - crc_time - parse_time;
    if (other_time < 0.0) other_time = 0.0;

    stats->t_decomp_alloc += alloc_time;
    stats->t_decomp_inflate += inflate_time;
    stats->t_decomp_crc += crc_time;
    stats->t_decomp_parse += parse_time;
    stats->t_decomp_other += other_time;
}

void MpiFlagstatAddRecord(const bam1_t *record, MpiFlagstatCounts *counts) {
    const uint16_t flag = record->core.flag;
    const int bucket = (flag & BAM_FQCFAIL) ? 1 : 0;
    const bool is_secondary = (flag & BAM_FSECONDARY) != 0;
    const bool is_supplementary = (flag & BAM_FSUPPLEMENTARY) != 0;
    const bool is_primary = !is_secondary && !is_supplementary;
    const bool is_paired = (flag & BAM_FPAIRED) != 0;
    const bool is_mapped = (flag & BAM_FUNMAP) == 0;
    const bool mate_mapped = (flag & BAM_FMUNMAP) == 0;

    counts->values[FLAGSTAT_TOTAL][bucket]++;
    if (is_primary) counts->values[FLAGSTAT_PRIMARY][bucket]++;
    if (is_secondary) counts->values[FLAGSTAT_SECONDARY][bucket]++;
    if (is_supplementary) counts->values[FLAGSTAT_SUPPLEMENTARY][bucket]++;
    if (flag & BAM_FDUP) counts->values[FLAGSTAT_DUPLICATES][bucket]++;
    if (is_primary && (flag & BAM_FDUP)) counts->values[FLAGSTAT_PRIMARY_DUPLICATES][bucket]++;
    if (is_mapped) counts->values[FLAGSTAT_MAPPED][bucket]++;
    if (is_primary && is_mapped) counts->values[FLAGSTAT_PRIMARY_MAPPED][bucket]++;

    if (is_primary && is_paired) {
        counts->values[FLAGSTAT_PAIRED][bucket]++;
        if (flag & BAM_FREAD1) counts->values[FLAGSTAT_READ1][bucket]++;
        if (flag & BAM_FREAD2) counts->values[FLAGSTAT_READ2][bucket]++;
        if ((flag & BAM_FPROPER_PAIR) && is_mapped) {
            counts->values[FLAGSTAT_PROPERLY_PAIRED][bucket]++;
        }
        if (is_mapped && mate_mapped) {
            counts->values[FLAGSTAT_PAIR_MAPPED][bucket]++;
            if (record->core.tid != record->core.mtid) {
                counts->values[FLAGSTAT_DIFF_CHR][bucket]++;
                if (record->core.qual >= 5) counts->values[FLAGSTAT_DIFF_CHR_MAPQ5][bucket]++;
            }
        }
        if (is_mapped && !mate_mapped) {
            counts->values[FLAGSTAT_SINGLETONS][bucket]++;
        }
    }
}

void MpiFlagstatFormatPct(char *buf, size_t buf_size, long long value, long long total) {
    if (total > 0) {
        snprintf(buf, buf_size, "%.2f%%", 100.0 * (double)value / (double)total);
    } else {
        snprintf(buf, buf_size, "N/A");
    }
}

void MpiFlagstatPrintSimple(const MpiFlagstatCounts &counts,
                            FlagstatCounterId id,
                            const char *label) {
    printf("%lld + %lld %s\n", counts.values[id][0], counts.values[id][1], label);
}

void MpiFlagstatPrintPct(const MpiFlagstatCounts &counts,
                         FlagstatCounterId id,
                         FlagstatCounterId denom_id,
                         const char *label) {
    char pass_pct[32];
    char fail_pct[32];
    MpiFlagstatFormatPct(pass_pct, sizeof(pass_pct),
                         counts.values[id][0], counts.values[denom_id][0]);
    MpiFlagstatFormatPct(fail_pct, sizeof(fail_pct),
                         counts.values[id][1], counts.values[denom_id][1]);
    printf("%lld + %lld %s (%s : %s)\n",
           counts.values[id][0],
           counts.values[id][1],
           label,
           pass_pct,
           fail_pct);
}

void MpiFlagstatPrintSamtoolsStyle(const MpiFlagstatCounts &counts) {
    MpiFlagstatPrintSimple(counts, FLAGSTAT_TOTAL,
                           "in total (QC-passed reads + QC-failed reads)");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_PRIMARY, "primary");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_SECONDARY, "secondary");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_SUPPLEMENTARY, "supplementary");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_DUPLICATES, "duplicates");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_PRIMARY_DUPLICATES, "primary duplicates");
    MpiFlagstatPrintPct(counts, FLAGSTAT_MAPPED, FLAGSTAT_TOTAL, "mapped");
    MpiFlagstatPrintPct(counts, FLAGSTAT_PRIMARY_MAPPED, FLAGSTAT_PRIMARY, "primary mapped");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_PAIRED, "paired in sequencing");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_READ1, "read1");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_READ2, "read2");
    MpiFlagstatPrintPct(counts, FLAGSTAT_PROPERLY_PAIRED, FLAGSTAT_PAIRED, "properly paired");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_PAIR_MAPPED, "with itself and mate mapped");
    MpiFlagstatPrintPct(counts, FLAGSTAT_SINGLETONS, FLAGSTAT_PAIRED, "singletons");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_DIFF_CHR, "with mate mapped to a different chr");
    MpiFlagstatPrintSimple(counts, FLAGSTAT_DIFF_CHR_MAPQ5,
                           "with mate mapped to a different chr (mapQ>=5)");
}

void MpiFlagstatReduceStats(const MpiFlagstatStats &local_stats,
                            MpiFlagstatStats *global_stats) {
    long long local_long[3] = {
        local_stats.input_blocks,
        local_stats.group_count,
        local_stats.total_records
    };
    long long global_long[3] = {};
    double local_double[9] = {
        local_stats.t_read,
        local_stats.t_decomp,
        local_stats.t_decomp_alloc,
        local_stats.t_decomp_inflate,
        local_stats.t_decomp_crc,
        local_stats.t_decomp_parse,
        local_stats.t_decomp_other,
        local_stats.t_count,
        local_stats.t_fused_total
    };
    double global_double[9] = {};

    MPI_Reduce(local_long, global_long, 3, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_double, global_double, 9, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (global_stats) {
        global_stats->input_blocks = global_long[0];
        global_stats->group_count = global_long[1];
        global_stats->total_records = global_long[2];
        global_stats->t_read = global_double[0];
        global_stats->t_decomp = global_double[1];
        global_stats->t_decomp_alloc = global_double[2];
        global_stats->t_decomp_inflate = global_double[3];
        global_stats->t_decomp_crc = global_double[4];
        global_stats->t_decomp_parse = global_double[5];
        global_stats->t_decomp_other = global_double[6];
        global_stats->t_count = global_double[7];
        global_stats->t_fused_total = global_double[8];
    }
}

int FusedFlagstatMPI(MemReader &reader,
                     MpiFlagstatCounts *counts,
                     MpiFlagstatStats *stats) {
    double fused_t0 = GetTime();
    int ret = -1;
    Bam2BamPara paras[kFlagstatNB];
    MpiFlagstatBlockSet input_blocks = {};
    MpiFlagstatBlockSet un_blocks = {};
    MpiFlagstatRecordSet record_set = {};
    BamFilterOptions no_filter = MpiFlagstatNoFilter();
    int next_n_blocks = 0;

    auto do_read_group = [&](int *n_blocks) -> int {
        double read_t0 = GetTime();
        int count = 0;
        for (int b = 0; b < kFlagstatNB; ++b) {
            bam_block *blk = &input_blocks.blocks[b];
            int read_ret = MpiFlagstatMemReadBlock(reader.base, reader.size, reader.pos, blk);
            if (read_ret < 0 || blk->length == 28) break;
            blk->block_id = b;
            blk->pos = 0;
            count++;
        }
        *n_blocks = count;
        if (stats) stats->t_read += GetTime() - read_t0;
        return 0;
    };

    if (counts) memset(counts, 0, sizeof(*counts));
    if (stats) memset(stats, 0, sizeof(*stats));

    if (MpiFlagstatAllocateBlockSet(&input_blocks, kFlagstatNB) != 0 ||
        MpiFlagstatAllocateBlockSet(&un_blocks, kFlagstatNB) != 0 ||
        MpiFlagstatAllocateRecordSet(&record_set,
                                     kFlagstatNB * (int)MAX_RECORDS_PER_BLOCK) != 0) {
        fprintf(stderr, "ERROR: failed to allocate MPI flagstat workspace.\n");
        goto cleanup;
    }

    if (do_read_group(&next_n_blocks) != 0) goto cleanup;
    while (next_n_blocks > 0) {
        int n_blocks = next_n_blocks;
        next_n_blocks = 0;
        if (stats) {
            stats->input_blocks += n_blocks;
            stats->group_count++;
        }

        for (int b = 0; b < kFlagstatNB; ++b) {
            paras[b].filter = no_filter;
            paras[b].block_id = b;
            paras[b].output_records = record_set.ptrs + (size_t)b * MAX_RECORDS_PER_BLOCK;
            paras[b].record_base = record_set.records + (size_t)b * MAX_RECORDS_PER_BLOCK;
            paras[b].bam_lens = record_set.bam_lens + (size_t)b * MAX_RECORDS_PER_BLOCK;
            paras[b].n_total_records = 0;
            paras[b].n_kept_records = 0;
            paras[b].kept_total_len = 0;
            paras[b].decomp_alloc_cycles = 0;
            paras[b].decomp_inflate_cycles = 0;
            paras[b].decomp_crc_cycles = 0;
            paras[b].decomp_parse_cycles = 0;
            paras[b].decomp_total_cycles = 0;
            if (b < n_blocks) {
                paras[b].input_block = &input_blocks.blocks[b];
                paras[b].un_comp_block = &un_blocks.blocks[b];
                paras[b].status = 0;
            } else {
                paras[b].input_block = nullptr;
                paras[b].un_comp_block = nullptr;
                paras[b].status = -1;
            }
        }

        double decomp_t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_decompress_bam2bam_passthrough, paras, 1);
        athread_join();
        double decomp_wall = GetTime() - decomp_t0;
        if (stats) {
            stats->t_decomp += decomp_wall;
            MpiFlagstatAccumulateDecompDetail(paras, n_blocks, decomp_wall, stats);
        }

        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].status != 0) {
                fprintf(stderr, "ERROR: MPI flagstat decompress failed on input block %d with status %d.\n",
                        b, paras[b].status);
                goto cleanup;
            }
        }

        double count_t0 = GetTime();
        for (int b = 0; b < n_blocks; ++b) {
            if (stats) stats->total_records += paras[b].n_total_records;
            for (int r = 0; r < paras[b].n_total_records; ++r) {
                MpiFlagstatAddRecord(paras[b].output_records[r], counts);
            }
        }
        if (stats) stats->t_count += GetTime() - count_t0;

        if (do_read_group(&next_n_blocks) != 0) goto cleanup;
    }

    ret = 0;

cleanup:
    if (stats) stats->t_fused_total = GetTime() - fused_t0;
    MpiFlagstatFreeRecordSet(&record_set);
    MpiFlagstatFreeBlockSet(&input_blocks);
    MpiFlagstatFreeBlockSet(&un_blocks);
    return ret;
}

} // namespace

int ProcessFlagstatMPI(CmdInfo *cmd_info) {
    double t_init = GetTime();

    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int exit_code = 1;
    int local_ok = 1;
    samFile *sin = nullptr;
    sam_hdr_t *hdr = nullptr;
    hFILE *input_mem_hfile = nullptr;
    char *input_file_mem = nullptr;
    size_t input_file_size = 0;
    char *rank_input_mem = nullptr;
    size_t rank_input_size = 0;
    long long body_start = 0;
    long long n_blocks = 0;
    int input_format = -1;
    std::vector<long long> block_offsets;
    std::vector<long long> block_lengths;
    MemReader reader = {};
    MpiFlagstatCounts local_counts = {};
    MpiFlagstatCounts global_counts = {};
    MpiFlagstatStats local_stats = {};
    MpiFlagstatStats global_stats = {};

    double init_cost = GetTime() - t_init;
    double init_cost_max = MpiFlagstatReduceMaxCost(init_cost);
    if (rank == 0) {
        printf("111Complete the initialization cost %lf-----\n", init_cost_max);
    }

    {
        double preload_t0 = GetTime();
        if (MpiFlagstatLoadFileToMemory(cmd_info->in_file_name_, &input_file_mem, &input_file_size) != 0) {
            fprintf(stderr, "[rank %d] ERROR: cannot preload input %s into memory\n",
                    rank, cmd_info->in_file_name_.c_str());
            local_ok = 0;
        }
        double preload_cost = GetTime() - preload_t0;
        double preload_cost_max = MpiFlagstatReduceMaxCost(preload_cost);
        if (rank == 0 && local_ok) printf("222Complete the memory cost %lf--\n", preload_cost_max);
    }
    if (!MpiFlagstatAllRanksOk(local_ok)) goto cleanup;

    {
        double header_t0 = GetTime();

        input_mem_hfile = hopen("mem:", "rb:", input_file_mem, input_file_size);
        if (!input_mem_hfile) {
            fprintf(stderr, "[rank %d] ERROR: cannot open preloaded BAM memory for %s\n",
                    rank, cmd_info->in_file_name_.c_str());
            local_ok = 0;
        }
        if (local_ok) {
            sin = (samFile *)hts_hopen(input_mem_hfile, "data", "rb");
            if (!sin) {
                fprintf(stderr, "[rank %d] ERROR: cannot create HTS input handle from memory\n", rank);
                if (hclose(input_mem_hfile) != 0) {
                    fprintf(stderr, "[rank %d] ERROR: closing failed HTS memory handle failed.\n", rank);
                }
                input_mem_hfile = nullptr;
                input_file_mem = nullptr;
                input_file_size = 0;
                local_ok = 0;
            } else {
                input_mem_hfile = nullptr;
            }
        }

        if (local_ok) {
            hdr = sam_hdr_read(sin);
            if (!hdr) {
                fprintf(stderr, "[rank %d] ERROR: cannot read header from %s\n",
                        rank, cmd_info->in_file_name_.c_str());
                local_ok = 0;
            }
        }
        if (local_ok) {
            input_format = MpiFlagstatNormalizeFormat(sin->format.format);
            if (input_format != bam) {
                if (rank == 0) {
                    fprintf(stderr, "ERROR: RabbitBAM-MPI flagstat only supports BAM input in v1.\n");
                }
                local_ok = 0;
            }
        }
        if (local_ok) {
            body_start = (long long)sin->fp.bgzf->block_address;
            if (body_start < 0 ||
                (unsigned long long)body_start > (unsigned long long)input_file_size) {
                fprintf(stderr, "[rank %d] ERROR: invalid BAM body start offset %lld.\n", rank, body_start);
                local_ok = 0;
            }
        }

        double header_cost = GetTime() - header_t0;
        double header_cost_max = MpiFlagstatReduceMaxCost(header_cost);
        if (rank == 0 && local_ok) {
            printf("333Complete the head cost %lf---\n", header_cost_max);
        }
    }
    if (!MpiFlagstatAllRanksOk(local_ok)) goto cleanup;

    {
        double body_total_t0 = GetTime();

        double stage41_t0 = GetTime();
        if (rank == 0) {
            printf("Enable MPI FLAGSTAT mode (%d MPE + %d CPEs)!!!\n",
                   comm_size, comm_size * 64);
            if (MpiFlagstatScanBgzfBlocksInMemory(input_file_mem, input_file_size, body_start,
                                                  &block_offsets, &block_lengths) != 0) {
                fprintf(stderr, "ERROR: failed to scan input BGZF blocks for flagstat.\n");
                local_ok = 0;
            }
            n_blocks = (long long)block_offsets.size();
            if (local_ok && n_blocks > (long long)INT_MAX) {
                fprintf(stderr, "ERROR: too many BGZF blocks for MPI_Bcast in RabbitBAM-MPI flagstat v1.\n");
                local_ok = 0;
            }
            if (local_ok) {
                printf("MPI BAM scan complete. data_blocks=%lld body_start=%lld header_end=%lld\n",
                       n_blocks, body_start, (long long)sin->fp.bgzf->block_address);
            }
        }

        MPI_Bcast(&local_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (!local_ok) goto cleanup;
        MPI_Bcast(&body_start, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&n_blocks, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        if (rank != 0) {
            block_offsets.resize((size_t)n_blocks);
            block_lengths.resize((size_t)n_blocks);
        }
        if (n_blocks > 0) {
            MPI_Bcast(block_offsets.data(), (int)n_blocks, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(block_lengths.data(), (int)n_blocks, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        }

        long long begin = n_blocks * rank / comm_size;
        long long end = n_blocks * (rank + 1) / comm_size;
        if (MpiFlagstatSelectBlockRangeFromMemory(input_file_mem, input_file_size,
                                                  block_offsets, block_lengths,
                                                  begin, end,
                                                  &rank_input_mem, &rank_input_size) != 0) {
            fprintf(stderr, "[rank %d] ERROR: failed to select assigned BGZF block range [%lld, %lld) from memory.\n",
                    rank, begin, end);
            local_ok = 0;
        }
        double stage41_cost = GetTime() - stage41_t0;
        double stage41_cost_max = MpiFlagstatReduceMaxCost(stage41_cost);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.1 scan/split/broadcast/select cost %lf\n", stage41_cost_max);
        }
        if (!MpiFlagstatAllRanksOk(local_ok)) goto cleanup;

        double stage42_t0 = GetTime();
        reader.base = rank_input_mem;
        reader.size = rank_input_size;
        reader.pos = 0;
        double stage42_cost = GetTime() - stage42_t0;
        double stage42_cost_max = MpiFlagstatReduceMaxCost(stage42_cost);
        if (rank == 0) {
            printf("Complete the 4.2 init reader cost %lf\n", stage42_cost_max);
        }

        double stage43_t0 = GetTime();
        if (FusedFlagstatMPI(reader, &local_counts, &local_stats) != 0) {
            local_ok = 0;
        }
        double stage43_cost = GetTime() - stage43_t0;
        int global_ok = MpiFlagstatAllRanksOk(local_ok);
        double stage43_cost_max = MpiFlagstatReduceMaxCost(stage43_cost);
        if (rank == 0 && global_ok) {
            printf("Complete the 4.3 FusedFlagstatMPI cost %lf\n", stage43_cost_max);
        }
        if (!global_ok) goto cleanup;

        double stage44_t0 = GetTime();
        MPI_Reduce(&local_counts.values[0][0], &global_counts.values[0][0],
                   FLAGSTAT_COUNTER_COUNT * 2, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MpiFlagstatReduceStats(local_stats, rank == 0 ? &global_stats : nullptr);
        double stage44_cost = GetTime() - stage44_t0;
        double stage44_cost_max = MpiFlagstatReduceMaxCost(stage44_cost);
        if (rank == 0) {
            MpiFlagstatPrintSamtoolsStyle(global_counts);
            printf("FusedFlagstatMPI finished. ranks=%d in_blocks=%lld groups=%lld total_records=%lld\n",
                   comm_size,
                   global_stats.input_blocks,
                   global_stats.group_count,
                   global_stats.total_records);
            printf("  read_sum=%.3f  decomp_sum=%.3f  count_sum=%.3f  fused_total_sum=%.3f\n",
                   global_stats.t_read,
                   global_stats.t_decomp,
                   global_stats.t_count,
                   global_stats.t_fused_total);
            printf("  decomp_detail_sum alloc=%.3f  inflate=%.3f  crc=%.3f  parse=%.3f  other=%.3f\n",
                   global_stats.t_decomp_alloc,
                   global_stats.t_decomp_inflate,
                   global_stats.t_decomp_crc,
                   global_stats.t_decomp_parse,
                   global_stats.t_decomp_other);
            printf("Complete the 4.4 flagstat reduce/print cost %lf\n", stage44_cost_max);
        }

        double body_total_cost = GetTime() - body_total_t0;
        double body_total_cost_max = MpiFlagstatReduceMaxCost(body_total_cost);
        if (rank == 0) {
            printf("444Complete the total body cost %lf\n", body_total_cost_max);
        }
    }

    exit_code = 0;

cleanup:
    {
        double close_t0 = GetTime();
        if (hdr) sam_hdr_destroy(hdr);
        if (sin) {
            int ret = hts_close(sin);
            if (ret < 0) fprintf(stderr, "[rank %d] ERROR: closing input failed.\n", rank);
            input_file_mem = nullptr;
        } else if (input_mem_hfile) {
            if (hclose(input_mem_hfile) != 0) {
                fprintf(stderr, "[rank %d] ERROR: closing memory hFILE failed.\n", rank);
            }
            input_file_mem = nullptr;
        } else if (input_file_mem) {
            free(input_file_mem);
            input_file_mem = nullptr;
        }
        double close_cost = GetTime() - close_t0;
        double close_cost_max = MpiFlagstatReduceMaxCost(close_cost);
        if (rank == 0) printf("666close the files cost %lf-----\n", close_cost_max);
    }

    return exit_code;
}
