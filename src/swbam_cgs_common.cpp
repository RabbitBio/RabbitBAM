#include "swbam_cgs.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

bool CgsIsBamLikeFormat(int format) {
    return format == bam || format == binary_format;
}

bool CgsIsSamLikeFormat(int format) {
    return format == sam || format == text_format;
}

bool CgsHasBamFilterRequest(const CmdInfo *cmd_info) {
    return cmd_info->min_mapq_ >= 0 ||
           cmd_info->max_mapq_ >= 0 ||
           cmd_info->require_flag_ != 0 ||
           cmd_info->exclude_flag_ != 0 ||
           !cmd_info->ref_name_.empty() ||
           cmd_info->min_read_len_ >= 0 ||
           cmd_info->max_read_len_ >= 0;
}

BamFilterOptions CgsBuildBamFilterOptions(const CmdInfo *cmd_info) {
    BamFilterOptions filter;
    filter.min_mapq = cmd_info->min_mapq_;
    filter.max_mapq = cmd_info->max_mapq_;
    filter.require_flag = cmd_info->require_flag_;
    filter.exclude_flag = cmd_info->exclude_flag_;
    filter.ref_tid = -2;
    filter.min_read_len = cmd_info->min_read_len_;
    filter.max_read_len = cmd_info->max_read_len_;
    return filter;
}

const char *CgsOutputOpenMode(const std::string &out_name) {
    size_t len = out_name.size();
    if (len >= 4 && out_name.substr(len - 4) == ".bam") return "wb";
    return "w";
}

void InitMemWriterCGS(MemWriter &w, size_t cap) {
    if (cap == 0) cap = 64 * 1024 * 1024;
    w.data = (char *)malloc(cap);
    w.capacity = w.data ? cap : 0;
    w.size = 0;
}

int EnsureMemWriterCapacityCGS(MemWriter &w, size_t add) {
    if (w.size + add <= w.capacity) return 0;
    size_t next = w.capacity ? w.capacity : 64 * 1024 * 1024;
    while (w.size + add > next) {
        size_t prev = next;
        next *= 2;
        if (next < prev) return -1;
    }
    char *new_data = (char *)realloc(w.data, next);
    if (!new_data) return -1;
    w.data = new_data;
    w.capacity = next;
    return 0;
}

int WriteBlockToMemCGS(MemWriter &w, bam_block *block) {
    if (EnsureMemWriterCapacityCGS(w, block->length) != 0) return -1;
    memcpy(w.data + w.size, block->data, block->length);
    w.size += block->length;
    return 0;
}

int WriteBytesToMemCGS(MemWriter &w, const char *data, size_t len) {
    if (len == 0) return 0;
    if (EnsureMemWriterCapacityCGS(w, len) != 0) return -1;
    memcpy(w.data + w.size, data, len);
    w.size += len;
    return 0;
}

int MemReadBlockCGS(char *base, size_t size, size_t &pos, bam_block *block) {
    if (pos >= size) return -1;
    if (pos + BLOCK_HEADER_LENGTH > size) return -1;
    uint16_t bsize = *(uint16_t *)(base + pos + 16);
    bsize += 1;
    if (pos + bsize > size) return -1;
    memcpy(block->data, base + pos, bsize);
    block->length = bsize;
    block->pos = 0;
    block->errcode = 0;
    pos += bsize;
    return bsize;
}

bool MemGetlineCGS(MemReader &r, kstring_t *ks) {
    if (r.pos >= r.size) return false;
    size_t start = r.pos;
    while (r.pos < r.size && r.base[r.pos] != '\n') r.pos++;
    size_t len = r.pos - start;
    if (len > 0 && r.base[r.pos - 1] == '\r') len--;
    if (len + 1 > ks->m) {
        char *new_s = (char *)realloc(ks->s, len + 1);
        if (!new_s) return false;
        ks->s = new_s;
        ks->m = len + 1;
    }
    memcpy(ks->s, r.base + start, len);
    ks->s[len] = '\0';
    ks->l = len;
    if (r.pos < r.size && r.base[r.pos] == '\n') r.pos++;
    return true;
}

void InitEmptyCompParaCGS(Comp_Para *para, int block_id) {
    para->block_id = block_id;
    para->input_records = nullptr;
    para->n_records = 0;
    para->un_comp_block = nullptr;
    para->un_comp_size = 0;
    para->output_block = nullptr;
    para->output_size = 0;
    para->status = -1;
}

int CheckCgsResources() {
#ifdef PLATFORM_SUNWAY
    return 0;
#else
    return -1;
#endif
}

CgsCrossAllocator::CgsCrossAllocator() {}

CgsCrossAllocator::~CgsCrossAllocator() {
#ifdef PLATFORM_SUNWAY
    for (size_t i = allocations_.size(); i > 0; --i) {
        _sw_xfree(allocations_[i - 1]);
    }
#else
    for (size_t i = allocations_.size(); i > 0; --i) {
        free(allocations_[i - 1]);
    }
#endif
}

void *CgsCrossAllocator::allocate(size_t alignment, size_t size) {
    if (alignment == 0) alignment = 64;
    size_t total_size = size + alignment + sizeof(void *);
#ifdef PLATFORM_SUNWAY
    unsigned char *raw = (unsigned char *)_sw_xmalloc(total_size);
#else
    unsigned char *raw = (unsigned char *)malloc(total_size);
#endif
    if (!raw) return nullptr;
    uintptr_t raw_addr = (uintptr_t)raw;
    uintptr_t aligned_addr = (raw_addr + sizeof(void *) + alignment - 1) & ~(uintptr_t)(alignment - 1);
    unsigned char *aligned = (unsigned char *)aligned_addr;
    void **storage = (void **)(aligned - sizeof(void *));
    *storage = raw;
    allocations_.push_back(raw);
    return aligned;
}

int AllocBlockSet(CgsCrossAllocator &alloc, CgsBlockSet *set, int n) {
    set->n = n;
    set->blocks = (bam_block *)alloc.allocate(64, (size_t)n * sizeof(bam_block));
    set->data = (unsigned char *)alloc.allocate(64, (size_t)n * BGZF_MAX_BLOCK_SIZE);
    if (!set->blocks || !set->data) return -1;
    memset(set->blocks, 0, (size_t)n * sizeof(bam_block));
    for (int i = 0; i < n; ++i) {
        set->blocks[i].data = set->data + (size_t)i * BGZF_MAX_BLOCK_SIZE;
        set->blocks[i].block_id = i;
    }
    return 0;
}

int AllocRecordSet(CgsCrossAllocator &alloc, CgsRecordSet *set, int total_records) {
    set->capacity = total_records;
    set->records = (bam1_t *)alloc.allocate(64, (size_t)total_records * sizeof(bam1_t));
    set->data = (unsigned char *)alloc.allocate(64, (size_t)total_records * INIT_DATA_SIZE);
    set->ptrs = (bam1_t **)alloc.allocate(64, (size_t)total_records * sizeof(bam1_t *));
    set->bam_lens = (uint32_t *)alloc.allocate(64, (size_t)total_records * sizeof(uint32_t));
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

void ResetRecordPointersForBlock(CgsRecordSet *set, int block_id, int records_per_block) {
    int base = block_id * records_per_block;
    int end = std::min(base + records_per_block, set->capacity);
    for (int r = base; r < end; ++r) {
        set->records[r].data = set->data + (size_t)r * INIT_DATA_SIZE;
        set->records[r].m_data = INIT_DATA_SIZE;
        set->records[r].l_data = 0;
        set->records[r].mempolicy = BAM_USER_OWNS_DATA;
        set->ptrs[r] = &set->records[r];
    }
}

void ResetRecordRange(CgsRecordSet *set, int begin, int count) {
    int end = std::min(begin + count, set->capacity);
    for (int i = begin; i < end; ++i) {
        set->records[i].data = set->data + (size_t)i * INIT_DATA_SIZE;
        set->records[i].m_data = INIT_DATA_SIZE;
        set->records[i].l_data = 0;
        set->records[i].mempolicy = BAM_USER_OWNS_DATA;
        set->ptrs[i] = &set->records[i];
    }
}

int AllocPackWorkspace(CgsCrossAllocator &alloc, CgsPackWorkspace *workspace, int capacity, int n_plans) {
    workspace->capacity = capacity;
    workspace->max_plans = n_plans;
    workspace->records = (bam1_t **)alloc.allocate(64, (size_t)capacity * sizeof(bam1_t *));
    workspace->plans = (CgsPackedPlan *)alloc.allocate(64, (size_t)n_plans * sizeof(CgsPackedPlan));
    if (!workspace->records || !workspace->plans) return -1;
    ResetPackWorkspace(workspace);
    return 0;
}

void ResetPackWorkspace(CgsPackWorkspace *workspace) {
    workspace->active_blocks = 0;
    workspace->total_records = 0;
    workspace->current_begin = workspace->records;
    workspace->current_records = 0;
    workspace->current_len = 0;
}

int SealCurrentPackBlock(CgsPackWorkspace *workspace) {
    if (workspace->current_records == 0) return 0;
    if (workspace->active_blocks >= workspace->max_plans) return -1;
    CgsPackedPlan &plan = workspace->plans[workspace->active_blocks++];
    plan.records = workspace->current_begin;
    plan.n_records = workspace->current_records;
    plan.total_len = workspace->current_len;
    workspace->current_begin = workspace->records + workspace->total_records;
    workspace->current_records = 0;
    workspace->current_len = 0;
    return 0;
}

int AppendRecordToPackWorkspace(CgsPackWorkspace *workspace, bam1_t *record, uint32_t bam_len) {
    const uint32_t packed_len = bam_len + 4;
    if (packed_len > BGZF_BLOCK_SIZE) return -1;
    if (workspace->current_records > 0 && workspace->current_len + packed_len > BGZF_BLOCK_SIZE) return 1;
    if (workspace->total_records >= workspace->capacity) return -2;
    if (workspace->current_records == 0) {
        workspace->current_begin = workspace->records + workspace->total_records;
    }
    workspace->records[workspace->total_records++] = record;
    workspace->current_records++;
    workspace->current_len += packed_len;
    return 0;
}

namespace {

int NormalizeFormatCGS(int format) {
    if (CgsIsBamLikeFormat(format)) return bam;
    if (CgsIsSamLikeFormat(format)) return sam;
    return format;
}

int LoadFileToMemoryCGS(const std::string &path, const char *mode, char **data, size_t *size) {
    FILE *f = fopen(path.c_str(), mode);
    if (!f) {
        perror("fopen");
        return -1;
    }
    if (fseek(f, 0, SEEK_END) != 0) {
        fclose(f);
        return -1;
    }
    long file_size = ftell(f);
    if (file_size < 0) {
        fclose(f);
        return -1;
    }
    if (fseek(f, 0, SEEK_SET) != 0) {
        fclose(f);
        return -1;
    }
    *size = (size_t)file_size;
    *data = (char *)malloc(*size ? *size : 1);
    if (!*data) {
        fclose(f);
        return -1;
    }
    if (*size > 0 && fread(*data, 1, *size, f) != *size) {
        fclose(f);
        free(*data);
        *data = nullptr;
        return -1;
    }
    fclose(f);
    return 0;
}

int SkipSamHeaderInMemoryCGS(MemReader &reader) {
    kstring_t tmp;
    tmp.l = 0;
    tmp.m = MAX_SAM_LINE_SIZE;
    tmp.s = (char *)malloc(MAX_SAM_LINE_SIZE);
    if (!tmp.s) return -1;

    while (true) {
        size_t old_pos = reader.pos;
        if (!MemGetlineCGS(reader, &tmp)) break;
        if (tmp.l == 0) continue;
        if (tmp.s[0] != '@') {
            reader.pos = old_pos;
            break;
        }
    }
    free(tmp.s);
    return 0;
}

} // namespace

int ProcessSwBamCGS(CmdInfo *cmd_info) {
    double t0 = GetTime();
    double t_header = 0;
    double t_total = 0;
    int exit_code = 1;
    samFile *sin = nullptr;
    samFile *sout = nullptr;
    sam_hdr_t *hdr = nullptr;
    char *input_mem = nullptr;
    size_t input_size = 0;
    MemReader reader = {};
    MemWriter mem_writer = {};
    bool ran_body = false;
    int input_format = -1;
    int output_format = -1;
    bool bam_to_bam = false;
    bool bam_to_sam = false;
    bool sam_to_bam = false;

    if (cmd_info->validate_bounds_) {
        fprintf(stderr, "ERROR: --validate-bounds is not supported by RabbitBAM-CGS. Use RabbitBAM-X for checked paths.\n");
        return 1;
    }

    bool filter_requested = CgsHasBamFilterRequest(cmd_info);
    BamFilterOptions bam_filter = CgsBuildBamFilterOptions(cmd_info);

    sin = sam_open(cmd_info->in_file_name_.c_str(), "r");
    if (!sin) {
        fprintf(stderr, "Error opening input file %s\n", cmd_info->in_file_name_.c_str());
        goto cleanup;
    }
    sout = sam_open(cmd_info->out_file_name_.c_str(), CgsOutputOpenMode(cmd_info->out_file_name_));
    if (!sout) {
        fprintf(stderr, "Error opening output file %s\n", cmd_info->out_file_name_.c_str());
        goto cleanup;
    }
    printf("open the files cost %lf---\n", GetTime() - t0);

    t0 = GetTime();
    hdr = sam_hdr_read(sin);
    if (!hdr) {
        fprintf(stderr, "Error reading header from input file %s\n", cmd_info->in_file_name_.c_str());
        goto cleanup;
    }

    input_format = NormalizeFormatCGS(sin->format.format);
    output_format = NormalizeFormatCGS(sout->format.format);
    bam_to_bam = input_format == bam && output_format == bam;
    bam_to_sam = input_format == bam && output_format == sam;
    sam_to_bam = input_format == sam && output_format == bam;

    if (!bam_to_bam && !bam_to_sam && !sam_to_bam) {
        fprintf(stderr, "Unsupported CGS conversion: input format %d -> output format %d\n",
                input_format, output_format);
        goto cleanup;
    }

    if (filter_requested && !bam_to_bam) {
        fprintf(stderr, "ERROR: BAM filtering options are only supported for CGS BAM -> BAM.\n");
        goto cleanup;
    }
    if (cmd_info->min_mapq_ > cmd_info->max_mapq_ && cmd_info->max_mapq_ >= 0) {
        fprintf(stderr, "ERROR: --min-mapq cannot be greater than --max-mapq.\n");
        goto cleanup;
    }
    if (cmd_info->min_read_len_ > cmd_info->max_read_len_ && cmd_info->max_read_len_ >= 0) {
        fprintf(stderr, "ERROR: --min-read-len cannot be greater than --max-read-len.\n");
        goto cleanup;
    }
    if (bam_to_bam && !cmd_info->ref_name_.empty()) {
        int ref_tid = sam_hdr_name2tid(hdr, cmd_info->ref_name_.c_str());
        if (ref_tid < 0) {
            fprintf(stderr, "ERROR: reference name '%s' does not exist in the BAM header.\n",
                    cmd_info->ref_name_.c_str());
            goto cleanup;
        }
        bam_filter.ref_tid = ref_tid;
    }

    if (sam_hdr_write(sout, hdr) != 0) {
        fprintf(stderr, "Error writing header to output file %s\n", cmd_info->out_file_name_.c_str());
        goto cleanup;
    }
    t_header = GetTime() - t0;
    t_total += t_header;
    printf("Complete the head cost %lf\n", t_header);

    t0 = GetTime();
    if (input_format == bam) {
        if (LoadFileToMemoryCGS(cmd_info->in_file_name_, "rb", &input_mem, &input_size) != 0) {
            fprintf(stderr, "Failed to load BAM file into memory\n");
            goto cleanup;
        }
        printf("Loaded BAM into memory: %.2f MB\n", input_size / 1024.0 / 1024.0);
        reader.base = input_mem;
        reader.size = input_size;
        reader.pos = sin->fp.bgzf->block_address;
        size_t writer_estimate = bam_to_sam ? input_size * 5 : input_size;
        InitMemWriterCGS(mem_writer, writer_estimate ? writer_estimate : 64 * 1024 * 1024);
        // InitMemWriterCGS(mem_writer, 8ULL * 1024 * 1024 * 1024);
        printf("Initialized %s memory writer: %.2f MB\n",
               output_format == bam ? "BAM" : "SAM",
               mem_writer.capacity / 1024.0 / 1024.0);
    } else {
        if (LoadFileToMemoryCGS(cmd_info->in_file_name_, "r", &input_mem, &input_size) != 0) {
            fprintf(stderr, "Failed to load SAM file into memory\n");
            goto cleanup;
        }
        printf("Loaded SAM into memory: %.2f MB\n", input_size / 1024.0 / 1024.0);
        reader.base = input_mem;
        reader.size = input_size;
        reader.pos = 0;
        if (SkipSamHeaderInMemoryCGS(reader) != 0) {
            fprintf(stderr, "Failed to skip SAM header in memory\n");
            goto cleanup;
        }
        size_t writer_estimate = input_size / 2;
        InitMemWriterCGS(mem_writer, writer_estimate ? writer_estimate : 64 * 1024 * 1024);
        printf("Initialized BAM memory writer: %.2f MB\n", mem_writer.capacity / 1024.0 / 1024.0);
    }
    if (!mem_writer.data) {
        fprintf(stderr, "Failed to allocate CGS memory writer\n");
        goto cleanup;
    }
    printf("Complete the memory cost %lf\n", GetTime() - t0);

    {
        double tbody = GetTime();
        if (bam_to_bam) {
            printf("Enable CGS BAM2BAM mode (1 MPE + %d CPEs)!!!\n", CGS_NB);
            if (filter_requested) printf("Enable CGS BAM filtering options.\n");
            if (FusedBamToBamCGS(reader, mem_writer, bam_filter) != 0) goto cleanup;
        } else if (bam_to_sam) {
            printf("Enable CGS BAM2SAM mode (1 MPE + %d CPEs)!!!\n", CGS_NB);
            if (FusedBamToSamCGS(reader, mem_writer, hdr) != 0) goto cleanup;
        } else {
            printf("Enable CGS SAM2BAM mode (1 MPE + %d CPEs)!!!\n", CGS_NB);
            if (FusedSamToBamCGS(reader, mem_writer, hdr) != 0) goto cleanup;
        }
        ran_body = true;
        exit_code = 0;
        t_total += GetTime() - tbody;
        printf("Complete the body cost %lf\n", GetTime() - tbody);
        printf("Complete the total cost %lf---\n", t_total);
    }

    if (ran_body && exit_code == 0) {
        double dump_t0 = GetTime();
        if (output_format == bam) {
            if (bgzf_flush(sout->fp.bgzf) != 0) {
                fprintf(stderr, "Error flushing BAM output header stream\n");
                exit_code = 1;
                goto cleanup;
            }
            if (hwrite(sout->fp.bgzf->fp, mem_writer.data, mem_writer.size) != (ssize_t)mem_writer.size) {
                fprintf(stderr, "Error dumping CGS memory output to BAM file\n");
                exit_code = 1;
                goto cleanup;
            }
        } else {
            if (hflush(sout->fp.hfile) != 0) {
                fprintf(stderr, "Error flushing SAM output header stream\n");
                exit_code = 1;
                goto cleanup;
            }
            if (hwrite(sout->fp.hfile, mem_writer.data, mem_writer.size) != (ssize_t)mem_writer.size) {
                fprintf(stderr, "Error dumping CGS memory output to SAM file\n");
                exit_code = 1;
                goto cleanup;
            }
        }
        printf("Dump memory to output file cost %lf\n", GetTime() - dump_t0);
    }

cleanup:
    if (mem_writer.data) free(mem_writer.data);
    if (input_mem) free(input_mem);

    t0 = GetTime();
    if (hdr) sam_hdr_destroy(hdr);
    if (sout) {
        int ret = hts_close(sout);
        if (ret < 0) fprintf(stderr, "Error closing output.\n");
    }
    if (sin) {
        int ret = hts_close(sin);
        if (ret < 0) fprintf(stderr, "Error closing input.\n");
    }
    printf("close the files cost %lf---\n", GetTime() - t0);
    return exit_code;
}

int ProcessBamToBamCGS(CmdInfo *cmd_info) {
    return ProcessSwBamCGS(cmd_info);
}
