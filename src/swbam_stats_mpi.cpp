#include "swbam_mpi.h"

#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

#include <mpi.h>

extern "C" {
    void slave_mpi_stats_basic_count();
}

namespace {

const int kStatsNB = 64;
const int kStatsMaxInsertSize = 8000;
const int kStatsInsertBins = kStatsMaxInsertSize + 1;
const double kStatsInsertMainBulk = 0.99;

const unsigned char kStatsBgzfEofBlock[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
};

enum StatsLongId {
    STATS_NREADS_1ST = 0,
    STATS_NREADS_2ND,
    STATS_NREADS_OTHER,
    STATS_NREADS_FILTERED,
    STATS_NREADS_DUP,
    STATS_NREADS_UNMAPPED,
    STATS_NREADS_SINGLE_MAPPED,
    STATS_NREADS_PAIRED_AND_MAPPED,
    STATS_NREADS_PROPERLY_PAIRED,
    STATS_NREADS_PAIRED_TECH,
    STATS_NREADS_ANOMALOUS,
    STATS_NREADS_MQ0,
    STATS_NREADS_QCFAILED,
    STATS_NREADS_SECONDARY,
    STATS_NREADS_SUPPLEMENTARY,
    STATS_TOTAL_LEN,
    STATS_TOTAL_LEN_1ST,
    STATS_TOTAL_LEN_2ND,
    STATS_TOTAL_LEN_DUP,
    STATS_NBASES_MAPPED,
    STATS_NBASES_MAPPED_CIGAR,
    STATS_NBASES_TRIMMED,
    STATS_NMISMATCHES,
    STATS_MAX_LEN,
    STATS_MAX_LEN_1ST,
    STATS_MAX_LEN_2ND,
    STATS_SUM_QUAL,
    STATS_LONG_COUNT
};

enum StatsOrientDiagId {
    ORIENT_DIAG_REF_IN = 0,
    ORIENT_DIAG_REF_OUT,
    ORIENT_DIAG_REF_OTHER,
    ORIENT_DIAG_DOC_IN,
    ORIENT_DIAG_DOC_OUT,
    ORIENT_DIAG_DOC_OTHER,
    ORIENT_DIAG_REF_IN_LT_READ,
    ORIENT_DIAG_REF_IN_LT_2READ,
    ORIENT_DIAG_REF_IN_GE_2READ,
    ORIENT_DIAG_REF_OUT_LT_READ,
    ORIENT_DIAG_REF_OUT_LT_2READ,
    ORIENT_DIAG_REF_OUT_GE_2READ,
    ORIENT_DIAG_POS_NEG_STRAND_NEG,
    ORIENT_DIAG_POS_NEG_STRAND_POS,
    ORIENT_DIAG_POS_ZERO_STRAND_NEG,
    ORIENT_DIAG_POS_ZERO_STRAND_POS,
    ORIENT_DIAG_POS_POS_STRAND_NEG,
    ORIENT_DIAG_POS_POS_STRAND_POS,
    ORIENT_DIAG_READ1_LEFT,
    ORIENT_DIAG_READ1_RIGHT,
    ORIENT_DIAG_READ1_SAME_POS,
    ORIENT_DIAG_READ2_LEFT,
    ORIENT_DIAG_READ2_RIGHT,
    ORIENT_DIAG_READ2_SAME_POS,
    ORIENT_DIAG_ISIZE_NEG,
    ORIENT_DIAG_ISIZE_ZERO,
    ORIENT_DIAG_ISIZE_POS,
    ORIENT_DIAG_COUNT
};

struct MpiBasicStatsCounts {
    long long values[STATS_LONG_COUNT];
    long long isize_inward[kStatsInsertBins];
    long long isize_outward[kStatsInsertBins];
    long long isize_other[kStatsInsertBins];
    long long orient_diag[ORIENT_DIAG_COUNT];
};

struct MpiStatsSortState {
    long long has_coord;
    long long sorted;
    long long first_tid;
    long long first_pos;
    long long last_tid;
    long long last_pos;
};

struct MpiStatsPerf {
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

struct MpiStatsBlockSet {
    bam_block *blocks;
    unsigned char *data;
    int n;
};

bool MpiStatsIsBamLikeFormat(int format) {
    return format == bam || format == binary_format;
}

bool MpiStatsIsSamLikeFormat(int format) {
    return format == sam || format == text_format;
}

int MpiStatsNormalizeFormat(int format) {
    if (MpiStatsIsBamLikeFormat(format)) return bam;
    if (MpiStatsIsSamLikeFormat(format)) return sam;
    return format;
}

int MpiStatsAllRanksOk(int local_ok) {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return global_ok;
}

double MpiStatsReduceMaxCost(double local_cost) {
    double max_cost = 0.0;
    MPI_Reduce(&local_cost, &max_cost, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return max_cost;
}

int MpiStatsLoadFileToMemory(const std::string &path, char **data, size_t *size) {
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

int MpiStatsScanBgzfBlocksInMemory(const char *base, size_t size, long long body_start,
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

        bool is_eof = block_len == (int)sizeof(kStatsBgzfEofBlock) &&
                      memcmp(base + pos, kStatsBgzfEofBlock, sizeof(kStatsBgzfEofBlock)) == 0;
        if (is_eof) break;

        offsets->push_back(pos);
        lengths->push_back(block_len);
        pos += block_len;
    }

    return 0;
}

int MpiStatsSelectBlockRangeFromMemory(char *base, size_t input_size,
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

int MpiStatsAllocateBlockSet(MpiStatsBlockSet *set, int n) {
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

void MpiStatsFreeBlockSet(MpiStatsBlockSet *set) {
    if (set->blocks) aligned_free_custom((unsigned char *)set->blocks);
    if (set->data) aligned_free_custom(set->data);
    set->blocks = nullptr;
    set->data = nullptr;
    set->n = 0;
}

int MpiStatsMemReadBlock(char *base, size_t size, size_t &pos, bam_block *block) {
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

BamFilterOptions MpiStatsNoFilter() {
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

void MpiStatsAccumulateDecompDetail(const MpiStatsBasicCountPara *paras,
                                    int active_blocks,
                                    double decomp_wall,
                                    MpiStatsPerf *stats) {
    if (!stats || active_blocks <= 0 || decomp_wall <= 0.0) return;

    const MpiStatsBasicCountPara *critical = nullptr;
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

int MpiStatsUnclippedLength(const bam1_t *record) {
    int read_len = record->core.l_qseq;
    const uint32_t *cigar = bam_get_cigar(record);
    for (int i = 0; i < record->core.n_cigar; ++i) {
        if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
            read_len += bam_cigar_oplen(cigar[i]);
        }
    }
    return read_len;
}

long long MpiStatsMappedCigarBases(const bam1_t *record) {
    long long total = 0;
    const uint32_t *cigar = bam_get_cigar(record);
    for (int i = 0; i < record->core.n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CEQUAL || op == BAM_CDIFF) {
            total += bam_cigar_oplen(cigar[i]);
        }
    }
    return total;
}

uint32_t MpiStatsReadOrder(const bam1_t *record) {
    const uint16_t flag = record->core.flag;
    if ((flag & BAM_FPAIRED) == 0) return 1;
    return ((flag & BAM_FREAD1) ? 1u : 0u) + ((flag & BAM_FREAD2) ? 2u : 0u);
}

void MpiStatsUpdateSort(const bam1_t *record, MpiStatsSortState *sort_state) {
    const long long tid = record->core.tid;
    const long long pos = record->core.pos;
    if (!sort_state->has_coord) {
        sort_state->has_coord = 1;
        sort_state->first_tid = tid;
        sort_state->first_pos = pos;
        sort_state->last_tid = tid;
        sort_state->last_pos = pos;
        return;
    }
    if (tid < sort_state->last_tid ||
        (tid == sort_state->last_tid && pos < sort_state->last_pos)) {
        sort_state->sorted = 0;
    }
    sort_state->last_tid = tid;
    sort_state->last_pos = pos;
}

enum MpiStatsOrientClass {
    MPI_STATS_ORIENT_IN = 0,
    MPI_STATS_ORIENT_OUT = 1,
    MPI_STATS_ORIENT_OTHER = 2
};

MpiStatsOrientClass MpiStatsReferenceOrientClass(const bam1_t *record) {
    const uint16_t flag = record->core.flag;
    const long long pos_fst = record->core.mpos - record->core.pos;
    const int is_fst = (flag & BAM_FREAD1) ? 1 : -1;
    const int is_fwd = (flag & BAM_FREVERSE) ? -1 : 1;
    const int is_mfwd = (flag & BAM_FMREVERSE) ? -1 : 1;

    if (is_fwd * is_mfwd > 0) return MPI_STATS_ORIENT_OTHER;
    if (is_fst * pos_fst >= 0) {
        return (is_fst * is_fwd > 0) ? MPI_STATS_ORIENT_IN : MPI_STATS_ORIENT_OUT;
    }
    return (is_fst * is_fwd > 0) ? MPI_STATS_ORIENT_OUT : MPI_STATS_ORIENT_IN;
}

MpiStatsOrientClass MpiStatsReadStrandOrientClass(const bam1_t *record) {
    const uint16_t flag = record->core.flag;
    const int is_fst = (flag & BAM_FREAD1) ? 1 : -1;
    const int is_fwd = (flag & BAM_FREVERSE) ? -1 : 1;
    const int is_mfwd = (flag & BAM_FMREVERSE) ? -1 : 1;

    if (is_fwd * is_mfwd > 0) return MPI_STATS_ORIENT_OTHER;
    return (is_fst * is_fwd > 0) ? MPI_STATS_ORIENT_IN : MPI_STATS_ORIENT_OUT;
}

void MpiStatsAddInsertDiag(const bam1_t *record,
                           int isize,
                           long long signed_isize,
                           MpiBasicStatsCounts *counts) {
    const uint16_t flag = record->core.flag;
    const MpiStatsOrientClass ref_class = MpiStatsReferenceOrientClass(record);
    const MpiStatsOrientClass doc_class = MpiStatsReadStrandOrientClass(record);

    if (ref_class == MPI_STATS_ORIENT_IN) counts->orient_diag[ORIENT_DIAG_REF_IN]++;
    else if (ref_class == MPI_STATS_ORIENT_OUT) counts->orient_diag[ORIENT_DIAG_REF_OUT]++;
    else counts->orient_diag[ORIENT_DIAG_REF_OTHER]++;

    if (doc_class == MPI_STATS_ORIENT_IN) counts->orient_diag[ORIENT_DIAG_DOC_IN]++;
    else if (doc_class == MPI_STATS_ORIENT_OUT) counts->orient_diag[ORIENT_DIAG_DOC_OUT]++;
    else counts->orient_diag[ORIENT_DIAG_DOC_OTHER]++;

    const int read_len = record->core.l_qseq > 0 ? record->core.l_qseq : 1;
    int size_bucket = 2;
    if (isize < read_len) size_bucket = 0;
    else if (isize < 2 * read_len) size_bucket = 1;

    if (ref_class == MPI_STATS_ORIENT_IN) {
        counts->orient_diag[ORIENT_DIAG_REF_IN_LT_READ + size_bucket]++;
    } else if (ref_class == MPI_STATS_ORIENT_OUT) {
        counts->orient_diag[ORIENT_DIAG_REF_OUT_LT_READ + size_bucket]++;
    }

    const long long pos_fst = record->core.mpos - record->core.pos;
    const int is_fst = (flag & BAM_FREAD1) ? 1 : -1;
    const int is_fwd = (flag & BAM_FREVERSE) ? -1 : 1;
    const long long pos_term = is_fst * pos_fst;
    const int strand_term = is_fst * is_fwd;
    int sign_base = ORIENT_DIAG_POS_ZERO_STRAND_NEG;
    if (pos_term < 0) sign_base = ORIENT_DIAG_POS_NEG_STRAND_NEG;
    else if (pos_term > 0) sign_base = ORIENT_DIAG_POS_POS_STRAND_NEG;
    counts->orient_diag[sign_base + (strand_term > 0 ? 1 : 0)]++;

    if (flag & BAM_FREAD1) {
        if (pos_fst > 0) counts->orient_diag[ORIENT_DIAG_READ1_LEFT]++;
        else if (pos_fst < 0) counts->orient_diag[ORIENT_DIAG_READ1_RIGHT]++;
        else counts->orient_diag[ORIENT_DIAG_READ1_SAME_POS]++;
    } else if (flag & BAM_FREAD2) {
        if (pos_fst > 0) counts->orient_diag[ORIENT_DIAG_READ2_LEFT]++;
        else if (pos_fst < 0) counts->orient_diag[ORIENT_DIAG_READ2_RIGHT]++;
        else counts->orient_diag[ORIENT_DIAG_READ2_SAME_POS]++;
    }

    if (signed_isize < 0) counts->orient_diag[ORIENT_DIAG_ISIZE_NEG]++;
    else if (signed_isize > 0) counts->orient_diag[ORIENT_DIAG_ISIZE_POS]++;
    else counts->orient_diag[ORIENT_DIAG_ISIZE_ZERO]++;
}

void MpiStatsAddInsert(const bam1_t *record, MpiBasicStatsCounts *counts, int collect_diag) {
    const uint16_t flag = record->core.flag;
    if ((flag & BAM_FPAIRED) == 0 || (flag & BAM_FUNMAP) || (flag & BAM_FMUNMAP)) return;
    if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) return;

    const long long signed_isize = record->core.isize;
    long long isize_ll = signed_isize;
    if (isize_ll < 0) isize_ll = -isize_ll;
    if (isize_ll > kStatsMaxInsertSize) isize_ll = kStatsMaxInsertSize;
    const int isize = (int)isize_ll;
    if (isize == 0 && record->core.tid != record->core.mtid) return;

    if (collect_diag) MpiStatsAddInsertDiag(record, isize, signed_isize, counts);

    const MpiStatsOrientClass orient_class = MpiStatsReferenceOrientClass(record);
    if (orient_class == MPI_STATS_ORIENT_OTHER) {
        counts->isize_other[isize]++;
    } else if (orient_class == MPI_STATS_ORIENT_OUT) {
        counts->isize_outward[isize]++;
    } else {
        counts->isize_inward[isize]++;
    }
}

void MpiStatsAddRecord(bam1_t *record,
                       MpiBasicStatsCounts *counts,
                       MpiStatsSortState *sort_state,
                       int collect_diag) {
    const uint16_t flag = record->core.flag;
    const bool is_secondary = (flag & BAM_FSECONDARY) != 0;
    const bool is_supplementary = (flag & BAM_FSUPPLEMENTARY) != 0;
    const bool is_original = !is_secondary && !is_supplementary;

    if (is_secondary) {
        counts->values[STATS_NREADS_SECONDARY]++;
        return;
    }
    if (is_supplementary) {
        counts->values[STATS_NREADS_SUPPLEMENTARY]++;
    }

    const int seq_len = record->core.l_qseq;
    if (seq_len == 0) return;

    if (flag & BAM_FDUP) {
        counts->values[STATS_TOTAL_LEN_DUP] += seq_len;
        counts->values[STATS_NREADS_DUP]++;
    }

    const uint32_t order = MpiStatsReadOrder(record);
    const int read_len = MpiStatsUnclippedLength(record);
    if (read_len > counts->values[STATS_MAX_LEN]) counts->values[STATS_MAX_LEN] = read_len;
    if (order == 1 && read_len > counts->values[STATS_MAX_LEN_1ST]) {
        counts->values[STATS_MAX_LEN_1ST] = read_len;
    }
    if (order == 2 && read_len > counts->values[STATS_MAX_LEN_2ND]) {
        counts->values[STATS_MAX_LEN_2ND] = read_len;
    }

    if (is_original) {
        counts->values[STATS_TOTAL_LEN] += seq_len;
        if (flag & BAM_FQCFAIL) counts->values[STATS_NREADS_QCFAILED]++;
        if (flag & BAM_FPAIRED) counts->values[STATS_NREADS_PAIRED_TECH]++;

        if (order == 1) {
            counts->values[STATS_NREADS_1ST]++;
            counts->values[STATS_TOTAL_LEN_1ST] += seq_len;
        } else if (order == 2) {
            counts->values[STATS_NREADS_2ND]++;
            counts->values[STATS_TOTAL_LEN_2ND] += seq_len;
        } else {
            counts->values[STATS_NREADS_OTHER]++;
        }

        if (order == 1 || order == 2) {
            uint8_t *qual = bam_get_qual(record);
            for (int i = 0; i < seq_len; ++i) {
                counts->values[STATS_SUM_QUAL] += qual[i];
            }
        }

        if (flag & BAM_FUNMAP) {
            counts->values[STATS_NREADS_UNMAPPED]++;
        } else {
            counts->values[STATS_NBASES_MAPPED] += seq_len;
            if (record->core.qual == 0) counts->values[STATS_NREADS_MQ0]++;
            if ((flag & BAM_FPAIRED) && !(flag & BAM_FMUNMAP)) {
                counts->values[STATS_NREADS_PAIRED_AND_MAPPED]++;
                if (((flag & (BAM_FPAIRED | BAM_FPROPER_PAIR)) == (BAM_FPAIRED | BAM_FPROPER_PAIR))) {
                    counts->values[STATS_NREADS_PROPERLY_PAIRED]++;
                }
                if (record->core.tid != record->core.mtid) counts->values[STATS_NREADS_ANOMALOUS]++;
            } else {
                counts->values[STATS_NREADS_SINGLE_MAPPED]++;
            }
        }
    }

    if (flag & BAM_FUNMAP) return;

    MpiStatsAddInsert(record, counts, collect_diag);

    uint8_t *nm = bam_aux_get(record, "NM");
    if (nm) counts->values[STATS_NMISMATCHES] += bam_aux2i(nm);
    counts->values[STATS_NBASES_MAPPED_CIGAR] += MpiStatsMappedCigarBases(record);
    MpiStatsUpdateSort(record, sort_state);
}

void MpiStatsMergeSlice(const MpiStatsBasicCountSlice &slice,
                        MpiBasicStatsCounts *counts) {
    for (int i = 0; i < STATS_LONG_COUNT; ++i) {
        if (i == STATS_MAX_LEN || i == STATS_MAX_LEN_1ST || i == STATS_MAX_LEN_2ND) {
            if (slice.values[i] > counts->values[i]) counts->values[i] = slice.values[i];
            continue;
        }
        counts->values[i] += slice.values[i];
    }
    for (int i = 0; i < kStatsInsertBins; ++i) {
        counts->isize_inward[i] += slice.isize_inward[i];
        counts->isize_outward[i] += slice.isize_outward[i];
        counts->isize_other[i] += slice.isize_other[i];
    }
    for (int i = 0; i < ORIENT_DIAG_COUNT; ++i) {
        counts->orient_diag[i] += slice.orient_diag[i];
    }
}

void MpiStatsMergeBlockSort(const MpiStatsBlockSortState &block_sort,
                            MpiStatsSortState *rank_sort) {
    if (!block_sort.has_coord) return;
    if (!rank_sort->has_coord) {
        rank_sort->has_coord = 1;
        rank_sort->first_tid = block_sort.first_tid;
        rank_sort->first_pos = block_sort.first_pos;
        rank_sort->last_tid = block_sort.last_tid;
        rank_sort->last_pos = block_sort.last_pos;
        if (!block_sort.sorted) rank_sort->sorted = 0;
        return;
    }
    if (block_sort.first_tid < rank_sort->last_tid ||
        (block_sort.first_tid == rank_sort->last_tid &&
         block_sort.first_pos < rank_sort->last_pos)) {
        rank_sort->sorted = 0;
    }
    if (!block_sort.sorted) rank_sort->sorted = 0;
    rank_sort->last_tid = block_sort.last_tid;
    rank_sort->last_pos = block_sort.last_pos;
}

void MpiStatsReducePerf(const MpiStatsPerf &local_stats, MpiStatsPerf *global_stats) {
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

int MpiStatsIsSortedGlobal(int comm_size, const std::vector<long long> &sort_values) {
    int sorted = 1;
    int have_prev = 0;
    long long prev_tid = 0;
    long long prev_pos = 0;
    for (int r = 0; r < comm_size; ++r) {
        const long long *state = &sort_values[(size_t)r * 6];
        const long long has_coord = state[0];
        const long long local_sorted = state[1];
        const long long first_tid = state[2];
        const long long first_pos = state[3];
        const long long last_tid = state[4];
        const long long last_pos = state[5];
        if (!local_sorted) sorted = 0;
        if (!has_coord) continue;
        if (have_prev &&
            (first_tid < prev_tid ||
             (first_tid == prev_tid && first_pos < prev_pos))) {
            sorted = 0;
        }
        have_prev = 1;
        prev_tid = last_tid;
        prev_pos = last_pos;
    }
    return sorted;
}

void MpiStatsComputeInsertSummary(const MpiBasicStatsCounts &counts,
                                  double *avg_isize,
                                  double *sd_isize,
                                  long long *n_inward,
                                  long long *n_outward,
                                  long long *n_other) {
    long long total_pairs = 0;
    double bulk = 0.0;
    double avg = 0.0;
    int ibulk = 0;
    long long full_inward = 0;
    long long full_outward = 0;
    long long full_other = 0;

    for (int isize = 0; isize < kStatsInsertBins; ++isize) {
        const long long inward = counts.isize_inward[isize] / 2;
        const long long outward = counts.isize_outward[isize] / 2;
        const long long other = counts.isize_other[isize] / 2;
        full_inward += inward;
        full_outward += outward;
        full_other += other;
        total_pairs += inward + outward + other;
    }

    long long main_pairs = total_pairs;
    for (int isize = 0; isize < kStatsInsertBins; ++isize) {
        const long long num = counts.isize_inward[isize] / 2 +
                              counts.isize_outward[isize] / 2 +
                              counts.isize_other[isize] / 2;
        if (num > 0) ibulk = isize + 1;
        bulk += (double)num;
        avg += (double)isize * (double)num;
        if (total_pairs > 0 && bulk / (double)total_pairs > kStatsInsertMainBulk) {
            ibulk = isize + 1;
            main_pairs = (long long)bulk;
            break;
        }
    }

    avg /= main_pairs ? (double)main_pairs : 1.0;
    double sd = 0.0;
    for (int isize = 1; isize < ibulk; ++isize) {
        const long long num = counts.isize_inward[isize] / 2 +
                              counts.isize_outward[isize] / 2 +
                              counts.isize_other[isize] / 2;
        const double diff = (double)isize - avg;
        sd += (double)num * diff * diff / (main_pairs ? (double)main_pairs : 1.0);
    }

    *avg_isize = avg;
    *sd_isize = sqrt(sd);
    *n_inward = full_inward;
    *n_outward = full_outward;
    *n_other = full_other;
}

void MpiStatsPrintBasic(const MpiBasicStatsCounts &counts, int is_sorted) {
    const long long nreads_1st = counts.values[STATS_NREADS_1ST];
    const long long nreads_2nd = counts.values[STATS_NREADS_2ND];
    const long long nreads_other = counts.values[STATS_NREADS_OTHER];
    const long long sequences = nreads_1st + nreads_2nd + nreads_other;
    const long long raw_total = counts.values[STATS_NREADS_FILTERED] + sequences;
    const long long reads_mapped = counts.values[STATS_NREADS_PAIRED_AND_MAPPED] +
                                   counts.values[STATS_NREADS_SINGLE_MAPPED];
    const double error_rate = counts.values[STATS_NBASES_MAPPED_CIGAR] > 0
        ? (double)counts.values[STATS_NMISMATCHES] / (double)counts.values[STATS_NBASES_MAPPED_CIGAR]
        : 0.0;
    const double average_length = sequences > 0
        ? (double)counts.values[STATS_TOTAL_LEN] / (double)sequences
        : 0.0;
    const double average_first_length = nreads_1st > 0
        ? (double)counts.values[STATS_TOTAL_LEN_1ST] / (double)nreads_1st
        : 0.0;
    const double average_last_length = nreads_2nd > 0
        ? (double)counts.values[STATS_TOTAL_LEN_2ND] / (double)nreads_2nd
        : 0.0;
    const double average_quality = counts.values[STATS_TOTAL_LEN] > 0
        ? (double)counts.values[STATS_SUM_QUAL] / (double)counts.values[STATS_TOTAL_LEN]
        : 0.0;
    const double properly_paired_pct = sequences > 0
        ? 100.0 * (double)counts.values[STATS_NREADS_PROPERLY_PAIRED] / (double)sequences
        : 0.0;

    double avg_isize = 0.0;
    double sd_isize = 0.0;
    long long n_inward = 0;
    long long n_outward = 0;
    long long n_other = 0;
    MpiStatsComputeInsertSummary(counts, &avg_isize, &sd_isize,
                                 &n_inward, &n_outward, &n_other);

    printf("# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.\n");
    printf("SN\traw total sequences:\t%lld\t# excluding supplementary and secondary reads\n", raw_total);
    printf("SN\tfiltered sequences:\t%lld\n", counts.values[STATS_NREADS_FILTERED]);
    printf("SN\tsequences:\t%lld\n", sequences);
    printf("SN\tis sorted:\t%d\t# %s by coordinate\n",
           is_sorted ? 1 : 0, is_sorted ? "sorted" : "not sorted");
    printf("SN\t1st fragments:\t%lld\n", nreads_1st);
    printf("SN\tlast fragments:\t%lld\n", nreads_2nd);
    printf("SN\treads mapped:\t%lld\n", reads_mapped);
    printf("SN\treads mapped and paired:\t%lld\t# paired-end technology bit set + both mates mapped\n",
           counts.values[STATS_NREADS_PAIRED_AND_MAPPED]);
    printf("SN\treads unmapped:\t%lld\n", counts.values[STATS_NREADS_UNMAPPED]);
    printf("SN\treads properly paired:\t%lld\t# proper-pair bit set\n",
           counts.values[STATS_NREADS_PROPERLY_PAIRED]);
    printf("SN\treads paired:\t%lld\t# paired-end technology bit set\n",
           counts.values[STATS_NREADS_PAIRED_TECH]);
    printf("SN\treads duplicated:\t%lld\t# PCR or optical duplicate bit set\n",
           counts.values[STATS_NREADS_DUP]);
    printf("SN\treads MQ0:\t%lld\t# mapped and MQ=0\n", counts.values[STATS_NREADS_MQ0]);
    printf("SN\treads QC failed:\t%lld\n", counts.values[STATS_NREADS_QCFAILED]);
    printf("SN\tnon-primary alignments:\t%lld\n", counts.values[STATS_NREADS_SECONDARY]);
    printf("SN\tsupplementary alignments:\t%lld\n", counts.values[STATS_NREADS_SUPPLEMENTARY]);
    printf("SN\ttotal length:\t%lld\t# ignores clipping\n", counts.values[STATS_TOTAL_LEN]);
    printf("SN\ttotal first fragment length:\t%lld\t# ignores clipping\n",
           counts.values[STATS_TOTAL_LEN_1ST]);
    printf("SN\ttotal last fragment length:\t%lld\t# ignores clipping\n",
           counts.values[STATS_TOTAL_LEN_2ND]);
    printf("SN\tbases mapped:\t%lld\t# ignores clipping\n", counts.values[STATS_NBASES_MAPPED]);
    printf("SN\tbases mapped (cigar):\t%lld\t# more accurate\n",
           counts.values[STATS_NBASES_MAPPED_CIGAR]);
    printf("SN\tbases trimmed:\t%lld\n", counts.values[STATS_NBASES_TRIMMED]);
    printf("SN\tbases duplicated:\t%lld\n", counts.values[STATS_TOTAL_LEN_DUP]);
    printf("SN\tmismatches:\t%lld\t# from NM fields\n", counts.values[STATS_NMISMATCHES]);
    printf("SN\terror rate:\t%e\t# mismatches / bases mapped (cigar)\n", error_rate);
    printf("SN\taverage length:\t%.0f\n", average_length);
    printf("SN\taverage first fragment length:\t%.0f\n", average_first_length);
    printf("SN\taverage last fragment length:\t%.0f\n", average_last_length);
    printf("SN\tmaximum length:\t%lld\n", counts.values[STATS_MAX_LEN]);
    printf("SN\tmaximum first fragment length:\t%lld\n", counts.values[STATS_MAX_LEN_1ST]);
    printf("SN\tmaximum last fragment length:\t%lld\n", counts.values[STATS_MAX_LEN_2ND]);
    printf("SN\taverage quality:\t%.1f\n", average_quality);
    printf("SN\tinsert size average:\t%.1f\n", avg_isize);
    printf("SN\tinsert size standard deviation:\t%.1f\n", sd_isize);
    printf("SN\tinward oriented pairs:\t%lld\n", n_inward);
    printf("SN\toutward oriented pairs:\t%lld\n", n_outward);
    printf("SN\tpairs with other orientation:\t%lld\n", n_other);
    printf("SN\tpairs on different chromosomes:\t%lld\n",
           counts.values[STATS_NREADS_ANOMALOUS] / 2);
    printf("SN\tpercentage of properly paired reads (%%):\t%.1f\n", properly_paired_pct);
}

void MpiStatsPrintOrientationDiag(const MpiBasicStatsCounts &counts) {
    long long active_in = 0;
    long long active_out = 0;
    long long active_other = 0;
    for (int isize = 0; isize < kStatsInsertBins; ++isize) {
        active_in += counts.isize_inward[isize] / 2;
        active_out += counts.isize_outward[isize] / 2;
        active_other += counts.isize_other[isize] / 2;
    }

    printf("ORIENT_DIAG\tactive_pairs\tinward=%lld\toutward=%lld\tother=%lld\n",
           active_in, active_out, active_other);
    printf("ORIENT_DIAG\tref_formula_pairs\tinward=%lld\toutward=%lld\tother=%lld\n",
           counts.orient_diag[ORIENT_DIAG_REF_IN] / 2,
           counts.orient_diag[ORIENT_DIAG_REF_OUT] / 2,
           counts.orient_diag[ORIENT_DIAG_REF_OTHER] / 2);
    printf("ORIENT_DIAG\tread_strand_pairs\tinward=%lld\toutward=%lld\tother=%lld\n",
           counts.orient_diag[ORIENT_DIAG_DOC_IN] / 2,
           counts.orient_diag[ORIENT_DIAG_DOC_OUT] / 2,
           counts.orient_diag[ORIENT_DIAG_DOC_OTHER] / 2);
    printf("ORIENT_DIAG\tref_in_by_isize_pairs\tlt_read=%lld\tlt_2read=%lld\tge_2read=%lld\n",
           counts.orient_diag[ORIENT_DIAG_REF_IN_LT_READ] / 2,
           counts.orient_diag[ORIENT_DIAG_REF_IN_LT_2READ] / 2,
           counts.orient_diag[ORIENT_DIAG_REF_IN_GE_2READ] / 2);
    printf("ORIENT_DIAG\tref_out_by_isize_pairs\tlt_read=%lld\tlt_2read=%lld\tge_2read=%lld\n",
           counts.orient_diag[ORIENT_DIAG_REF_OUT_LT_READ] / 2,
           counts.orient_diag[ORIENT_DIAG_REF_OUT_LT_2READ] / 2,
           counts.orient_diag[ORIENT_DIAG_REF_OUT_GE_2READ] / 2);
    printf("ORIENT_DIAG\tpos_strand_records\tpos_neg_strand_neg=%lld\tpos_neg_strand_pos=%lld\tpos_zero_strand_neg=%lld\tpos_zero_strand_pos=%lld\tpos_pos_strand_neg=%lld\tpos_pos_strand_pos=%lld\n",
           counts.orient_diag[ORIENT_DIAG_POS_NEG_STRAND_NEG],
           counts.orient_diag[ORIENT_DIAG_POS_NEG_STRAND_POS],
           counts.orient_diag[ORIENT_DIAG_POS_ZERO_STRAND_NEG],
           counts.orient_diag[ORIENT_DIAG_POS_ZERO_STRAND_POS],
           counts.orient_diag[ORIENT_DIAG_POS_POS_STRAND_NEG],
           counts.orient_diag[ORIENT_DIAG_POS_POS_STRAND_POS]);
    printf("ORIENT_DIAG\tread_order_position_records\tread1_left=%lld\tread1_right=%lld\tread1_same=%lld\tread2_left=%lld\tread2_right=%lld\tread2_same=%lld\n",
           counts.orient_diag[ORIENT_DIAG_READ1_LEFT],
           counts.orient_diag[ORIENT_DIAG_READ1_RIGHT],
           counts.orient_diag[ORIENT_DIAG_READ1_SAME_POS],
           counts.orient_diag[ORIENT_DIAG_READ2_LEFT],
           counts.orient_diag[ORIENT_DIAG_READ2_RIGHT],
           counts.orient_diag[ORIENT_DIAG_READ2_SAME_POS]);
    printf("ORIENT_DIAG\tisize_sign_records\tneg=%lld\tzero=%lld\tpos=%lld\n",
           counts.orient_diag[ORIENT_DIAG_ISIZE_NEG],
           counts.orient_diag[ORIENT_DIAG_ISIZE_ZERO],
           counts.orient_diag[ORIENT_DIAG_ISIZE_POS]);
}

int FusedStatsMPI(MemReader &reader,
                  MpiBasicStatsCounts *counts,
                  MpiStatsSortState *sort_state,
                  MpiStatsPerf *perf,
                  int collect_diag) {
    double fused_t0 = GetTime();
    int ret = -1;
    MpiStatsBasicCountPara paras[kStatsNB];
    MpiStatsBlockSet input_a = {};
    MpiStatsBlockSet input_b = {};
    MpiStatsBlockSet un_a = {};
    MpiStatsBlockSet un_b = {};
    MpiStatsBasicCountSlice *count_slices = nullptr;
    unsigned char *scratch_data = nullptr;
    const size_t scratch_stride = MPI_BAM_BLOCK_ARENA_SIZE;
    MpiStatsBlockSet *cur_input = nullptr;
    MpiStatsBlockSet *next_input = nullptr;
    MpiStatsBlockSet *cur_un = nullptr;
    MpiStatsBlockSet *next_un = nullptr;
    int n_blocks = 0;

    auto do_read_group = [&](MpiStatsBlockSet *input_blocks, int *n_blocks) -> int {
        double read_t0 = GetTime();
        int count = 0;
        for (int b = 0; b < kStatsNB; ++b) {
            bam_block *blk = &input_blocks->blocks[b];
            int read_ret = MpiStatsMemReadBlock(reader.base, reader.size, reader.pos, blk);
            if (read_ret < 0 || blk->length == 28) break;
            blk->block_id = b;
            blk->pos = 0;
            count++;
        }
        *n_blocks = count;
        if (perf) perf->t_read += GetTime() - read_t0;
        return 0;
    };

    if (counts) memset(counts, 0, sizeof(*counts));
    if (perf) memset(perf, 0, sizeof(*perf));
    if (sort_state) {
        memset(sort_state, 0, sizeof(*sort_state));
        sort_state->sorted = 1;
        sort_state->first_tid = -1;
        sort_state->last_tid = -1;
    }

    count_slices = (MpiStatsBasicCountSlice *)aligned_alloc_custom(
        64, (size_t)kStatsNB * sizeof(MpiStatsBasicCountSlice));
    scratch_data = aligned_alloc_custom(64, (size_t)kStatsNB * scratch_stride);
    if (MpiStatsAllocateBlockSet(&input_a, kStatsNB) != 0 ||
        MpiStatsAllocateBlockSet(&input_b, kStatsNB) != 0 ||
        MpiStatsAllocateBlockSet(&un_a, kStatsNB) != 0 ||
        MpiStatsAllocateBlockSet(&un_b, kStatsNB) != 0 ||
        !count_slices || !scratch_data) {
        fprintf(stderr, "ERROR: failed to allocate MPI stats workspace.\n");
        goto cleanup;
    }
    memset(count_slices, 0, (size_t)kStatsNB * sizeof(MpiStatsBasicCountSlice));

    cur_input = &input_a;
    next_input = &input_b;
    cur_un = &un_a;
    next_un = &un_b;
    if (do_read_group(cur_input, &n_blocks) != 0) goto cleanup;
    while (n_blocks > 0) {
        if (perf) {
            perf->input_blocks += n_blocks;
            perf->group_count++;
        }

        for (int b = 0; b < kStatsNB; ++b) {
            paras[b].block_id = b;
            paras[b].scratch_data = scratch_data + (size_t)b * scratch_stride;
            paras[b].scratch_capacity = scratch_stride;
            paras[b].counts = &count_slices[b];
            paras[b].collect_diag = collect_diag;
            paras[b].n_total_records = 0;
            paras[b].record_index = 0;
            paras[b].actual_value = 0;
            paras[b].limit_value = 0;
            paras[b].limit_id = BOUNDS_LIMIT_NONE;
            memset(&paras[b].sort_state, 0, sizeof(paras[b].sort_state));
            paras[b].sort_state.sorted = 1;
            paras[b].sort_state.first_tid = -1;
            paras[b].sort_state.last_tid = -1;
            paras[b].decomp_alloc_cycles = 0;
            paras[b].decomp_inflate_cycles = 0;
            paras[b].decomp_crc_cycles = 0;
            paras[b].decomp_parse_cycles = 0;
            paras[b].decomp_total_cycles = 0;
            if (b < n_blocks) {
                paras[b].input_block = &cur_input->blocks[b];
                paras[b].un_comp_block = &cur_un->blocks[b];
                paras[b].status = 0;
            } else {
                paras[b].input_block = nullptr;
                paras[b].un_comp_block = nullptr;
                paras[b].status = -1;
            }
        }

        double decomp_t0 = GetTime();
        __real_athread_spawn((void *)slave_mpi_stats_basic_count, paras, 1);
        int next_n_blocks = 0;
        #ifdef ENABLE_MASKING
        if (do_read_group(next_input, &next_n_blocks) != 0) goto cleanup;
        #endif
        athread_join();
        double decomp_wall = GetTime() - decomp_t0;
        #ifndef ENABLE_MASKING
        if (do_read_group(next_input, &next_n_blocks) != 0) goto cleanup;
        #endif
        if (perf) {
            perf->t_decomp += decomp_wall;
            MpiStatsAccumulateDecompDetail(paras, n_blocks, decomp_wall, perf);
        }

        for (int b = 0; b < n_blocks; ++b) {
            if (paras[b].status != 0) {
                if (paras[b].status == -3) {
                    fprintf(stderr,
                            "ERROR: MPI stats capacity exceeded on input block %d. limit_id=%d limit=%lld actual=%lld record=%d.\n",
                            b, paras[b].limit_id, paras[b].limit_value,
                            paras[b].actual_value, paras[b].record_index);
                } else {
                    fprintf(stderr,
                            "ERROR: MPI stats count failed on input block %d "
                            "with status %d. record=%d actual=%lld limit_id=%d.\n",
                            b, paras[b].status, paras[b].record_index,
                            paras[b].actual_value, paras[b].limit_id);
                }
                goto cleanup;
            }
        }

        double count_t0 = GetTime();
        for (int b = 0; b < n_blocks; ++b) {
            if (perf) perf->total_records += paras[b].n_total_records;
            if (sort_state) MpiStatsMergeBlockSort(paras[b].sort_state, sort_state);
        }
        if (perf) perf->t_count += GetTime() - count_t0;

        std::swap(cur_input, next_input);
        std::swap(cur_un, next_un);
        n_blocks = next_n_blocks;
    }

    if (counts) {
        double count_t0 = GetTime();
        for (int b = 0; b < kStatsNB; ++b) {
            MpiStatsMergeSlice(count_slices[b], counts);
        }
        if (perf) perf->t_count += GetTime() - count_t0;
    }

    ret = 0;

cleanup:
    if (perf) perf->t_fused_total = GetTime() - fused_t0;
    if (count_slices) aligned_free_custom((unsigned char *)count_slices);
    if (scratch_data) aligned_free_custom(scratch_data);
    MpiStatsFreeBlockSet(&input_a);
    MpiStatsFreeBlockSet(&input_b);
    MpiStatsFreeBlockSet(&un_a);
    MpiStatsFreeBlockSet(&un_b);
    return ret;
}

} // namespace

int ProcessStatsMPI(CmdInfo *cmd_info) {
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
    MpiBasicStatsCounts local_counts = {};
    MpiBasicStatsCounts global_counts = {};
    MpiStatsSortState local_sort = {};
    MpiStatsPerf local_perf = {};
    MpiStatsPerf global_perf = {};

    double init_cost = GetTime() - t_init;
    double init_cost_max = MpiStatsReduceMaxCost(init_cost);
    if (rank == 0) {
        printf("111Complete the initialization cost %lf-----\n", init_cost_max);
    }

    {
        double preload_t0 = GetTime();
        if (MpiStatsLoadFileToMemory(cmd_info->in_file_name_, &input_file_mem, &input_file_size) != 0) {
            fprintf(stderr, "[rank %d] ERROR: cannot preload input %s into memory\n",
                    rank, cmd_info->in_file_name_.c_str());
            local_ok = 0;
        }
        double preload_cost = GetTime() - preload_t0;
        double preload_cost_max = MpiStatsReduceMaxCost(preload_cost);
        if (rank == 0 && local_ok) printf("222Complete the memory cost %lf--\n", preload_cost_max);
    }
    if (!MpiStatsAllRanksOk(local_ok)) goto cleanup;

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
            input_format = MpiStatsNormalizeFormat(sin->format.format);
            if (input_format != bam) {
                if (rank == 0) {
                    fprintf(stderr, "ERROR: RabbitBAM-MPI stats --basic only supports BAM input in v1.\n");
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
        double header_cost_max = MpiStatsReduceMaxCost(header_cost);
        if (rank == 0 && local_ok) {
            printf("333Complete the head cost %lf---\n", header_cost_max);
        }
    }
    if (!MpiStatsAllRanksOk(local_ok)) goto cleanup;

    {
        double body_total_t0 = GetTime();

        double stage41_t0 = GetTime();
        if (rank == 0) {
            printf("Enable MPI STATS --basic mode (%d MPE + %d CPEs)!!!\n",
                   comm_size, comm_size * 64);
            if (MpiStatsScanBgzfBlocksInMemory(input_file_mem, input_file_size, body_start,
                                               &block_offsets, &block_lengths) != 0) {
                fprintf(stderr, "ERROR: failed to scan input BGZF blocks for stats.\n");
                local_ok = 0;
            }
            n_blocks = (long long)block_offsets.size();
            if (local_ok && n_blocks > (long long)INT_MAX) {
                fprintf(stderr, "ERROR: too many BGZF blocks for MPI_Bcast in RabbitBAM-MPI stats v1.\n");
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
        if (MpiStatsSelectBlockRangeFromMemory(input_file_mem, input_file_size,
                                               block_offsets, block_lengths,
                                               begin, end,
                                               &rank_input_mem, &rank_input_size) != 0) {
            fprintf(stderr, "[rank %d] ERROR: failed to select assigned BGZF block range [%lld, %lld) from memory.\n",
                    rank, begin, end);
            local_ok = 0;
        }
        double stage41_cost = GetTime() - stage41_t0;
        double stage41_cost_max = MpiStatsReduceMaxCost(stage41_cost);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.1 scan/split/broadcast/select cost %lf\n", stage41_cost_max);
        }
        if (!MpiStatsAllRanksOk(local_ok)) goto cleanup;

        double stage42_t0 = GetTime();
        reader.base = rank_input_mem;
        reader.size = rank_input_size;
        reader.pos = 0;
        double stage42_cost = GetTime() - stage42_t0;
        double stage42_cost_max = MpiStatsReduceMaxCost(stage42_cost);
        if (rank == 0) {
            printf("Complete the 4.2 init reader cost %lf\n", stage42_cost_max);
        }

        double stage43_t0 = GetTime();
        const int collect_diag = cmd_info->verbose_ ? 1 : 0;
        if (FusedStatsMPI(reader, &local_counts, &local_sort, &local_perf, collect_diag) != 0) {
            local_ok = 0;
        }
        double stage43_cost = GetTime() - stage43_t0;
        int global_ok = MpiStatsAllRanksOk(local_ok);
        double stage43_cost_max = MpiStatsReduceMaxCost(stage43_cost);
        if (rank == 0 && global_ok) {
            printf("Complete the 4.3 FusedStatsMPI cost %lf\n", stage43_cost_max);
        }
        if (!global_ok) goto cleanup;

        double stage44_t0 = GetTime();
        MPI_Reduce(local_counts.values, global_counts.values,
                   STATS_LONG_COUNT, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        long long local_max_values[3] = {
            local_counts.values[STATS_MAX_LEN],
            local_counts.values[STATS_MAX_LEN_1ST],
            local_counts.values[STATS_MAX_LEN_2ND]
        };
        long long global_max_values[3] = {};
        MPI_Reduce(local_max_values, global_max_values,
                   3, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            global_counts.values[STATS_MAX_LEN] = global_max_values[0];
            global_counts.values[STATS_MAX_LEN_1ST] = global_max_values[1];
            global_counts.values[STATS_MAX_LEN_2ND] = global_max_values[2];
        }
        MPI_Reduce(local_counts.isize_inward, global_counts.isize_inward,
                   kStatsInsertBins, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(local_counts.isize_outward, global_counts.isize_outward,
                   kStatsInsertBins, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(local_counts.isize_other, global_counts.isize_other,
                   kStatsInsertBins, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        if (cmd_info->verbose_) {
            MPI_Reduce(local_counts.orient_diag, global_counts.orient_diag,
                       ORIENT_DIAG_COUNT, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        MpiStatsReducePerf(local_perf, rank == 0 ? &global_perf : nullptr);

        long long local_sort_values[6] = {
            local_sort.has_coord,
            local_sort.sorted,
            local_sort.first_tid,
            local_sort.first_pos,
            local_sort.last_tid,
            local_sort.last_pos
        };
        std::vector<long long> global_sort_values;
        if (rank == 0) global_sort_values.resize((size_t)comm_size * 6);
        MPI_Gather(local_sort_values, 6, MPI_LONG_LONG,
                   rank == 0 ? global_sort_values.data() : nullptr,
                   6, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

        double stage44_cost = GetTime() - stage44_t0;
        double stage44_cost_max = MpiStatsReduceMaxCost(stage44_cost);
        if (rank == 0) {
            int is_sorted = MpiStatsIsSortedGlobal(comm_size, global_sort_values);
            MpiStatsPrintBasic(global_counts, is_sorted);
            if (cmd_info->verbose_) MpiStatsPrintOrientationDiag(global_counts);
            printf("FusedStatsMPI finished. ranks=%d in_blocks=%lld groups=%lld total_records=%lld\n",
                   comm_size,
                   global_perf.input_blocks,
                   global_perf.group_count,
                   global_perf.total_records);
            printf("  read_sum=%.3f  decomp_sum=%.3f  merge_sum=%.3f  fused_total_sum=%.3f\n",
                   global_perf.t_read,
                   global_perf.t_decomp,
                   global_perf.t_count,
                   global_perf.t_fused_total);
            printf("  decomp_detail_sum alloc=%.3f  inflate=%.3f  crc=%.3f  parse=%.3f  other=%.3f\n",
                   global_perf.t_decomp_alloc,
                   global_perf.t_decomp_inflate,
                   global_perf.t_decomp_crc,
                   global_perf.t_decomp_parse,
                   global_perf.t_decomp_other);
            printf("Complete the 4.4 stats reduce/print cost %lf\n", stage44_cost_max);
        }

        double body_total_cost = GetTime() - body_total_t0;
        double body_total_cost_max = MpiStatsReduceMaxCost(body_total_cost);
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
        double close_cost_max = MpiStatsReduceMaxCost(close_cost);
        if (rank == 0) printf("666close the files cost %lf-----\n", close_cost_max);
    }

    return exit_code;
}
