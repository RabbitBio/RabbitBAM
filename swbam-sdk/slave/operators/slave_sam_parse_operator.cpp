#include "swbam/bam_types.h"
#include <climits>
#include <cstring>
#include <htslib/hts_endian.h>
#include "sam_parse.h"
#include "swbam/cpe_codec.h"
#define MPI_SLAVE_RPCC() swbam_cpe_cycle_now()

extern "C" void slave_mpi_copy_and_count(void *arg) {
    MpiSamParseBatch *batch = (MpiSamParseBatch *)arg;
    int tid = _PEN;
    MpiSamParseChunk *chunk = &batch->chunks[tid];

    if (chunk->src_len == 0) {
        chunk->text_len = 0;
        chunk->count    = 0;
        chunk->status = 0;
        chunk->max_line_len = 0;
        chunk->estimated_bam_data = 0;
        chunk->record_index = 0;
        chunk->actual_value = 0;
        chunk->limit_value = 0;
        chunk->limit_id = BOUNDS_LIMIT_NONE;
        chunk->parse_fast_records = 0;
        chunk->parse_fallback_records = 0;
        chunk->parse_total_cycles = 0;
        chunk->parse_core_cycles = 0;
        chunk->parse_aux_cycles = 0;
        chunk->parse_cg_cycles = 0;
        chunk->parse_fallback_cycles = 0;
        return;
    }

    memcpy(chunk->text_buf, chunk->src_ptr, chunk->src_len);
    chunk->text_buf[chunk->src_len] = '\0';
    chunk->text_len = chunk->src_len;
    chunk->status = 0;
    chunk->max_line_len = 0;
    chunk->estimated_bam_data = 0;
    chunk->record_index = 0;
    chunk->actual_value = 0;
    chunk->limit_value = 0;
    chunk->limit_id = BOUNDS_LIMIT_NONE;
    chunk->parse_fast_records = 0;
    chunk->parse_fallback_records = 0;
    chunk->parse_total_cycles = 0;
    chunk->parse_core_cycles = 0;
    chunk->parse_aux_cycles = 0;
    chunk->parse_cg_cycles = 0;
    chunk->parse_fallback_cycles = 0;

    int count = 0;
    char *p = chunk->text_buf;
    char *end = chunk->text_buf + chunk->text_len;
    while (p < end) {
        char *line_end = p;
        while (line_end < end && *line_end != '\n' && *line_end != '\0') line_end++;

        int line_len = (int)(line_end - p);
        if (line_len > 0 && p[line_len - 1] == '\r') line_len--;
        if (line_len > 0) {
            count++;
            if ((size_t)line_len > chunk->max_line_len) {
                chunk->max_line_len = (size_t)line_len;
            }
            chunk->estimated_bam_data += (size_t)line_len + 64u;
        }

        if (line_end < end && *line_end == '\n') line_end++;
        p = line_end;
    }
    chunk->count = count;
}

extern "C" void slave_mpi_sam_parse_chunk(void *arg) {
    MpiSamParseBatch *batch = (MpiSamParseBatch *)arg;
    int tid = _PEN;
    MpiSamParseChunk *chunk = &batch->chunks[tid];

    chunk->parse_fast_records = 0;
    chunk->parse_fallback_records = 0;
    chunk->parse_total_cycles = 0;
    chunk->parse_core_cycles = 0;
    chunk->parse_aux_cycles = 0;
    chunk->parse_cg_cycles = 0;
    chunk->parse_fallback_cycles = 0;
    chunk->status = 0;
    chunk->bam_data_used = 0;
    chunk->record_index = 0;
    chunk->actual_value = 0;
    chunk->limit_value = 0;
    chunk->limit_id = BOUNDS_LIMIT_NONE;
    if (chunk->text_len == 0 || chunk->count == 0) return;
    if (chunk->bam_data == NULL || chunk->bam_data_capacity == 0) {
        chunk->status = -3;
        chunk->limit_id = BOUNDS_LIMIT_INIT_DATA_SIZE;
        chunk->limit_value = 0;
        chunk->actual_value = chunk->estimated_bam_data;
        return;
    }

    int valid_count = 0;
    size_t arena_used = 0;
    char *ptr = chunk->text_buf;
    char *end = chunk->text_buf + chunk->text_len;
    MpiSamParseFastCache fast_cache;
    memset(&fast_cache, 0, sizeof(fast_cache));

    while (ptr < end) {
        char *eol = ptr;
        while (eol < end && *eol != '\n' && *eol != '\0') eol++;

        int line_len = (int)(eol - ptr);
        kstring_t ks;
        ks.s = ptr;
        ks.l = line_len;
        ks.m = line_len + 1;

        if (ks.l > 0 && ks.s[ks.l - 1] == '\r') {
            ks.l--;
        }

        if (ks.l > 0) {
            char saved_char = ks.s[ks.l];
            ks.s[ks.l] = '\0';

	            if (valid_count >= chunk->count || valid_count >= MAX_BAMS_PER_CHUNK) {
	                chunk->status = -3;
	                chunk->limit_id = BOUNDS_LIMIT_MAX_BAMS_PER_CHUNK;
	                chunk->limit_value = MAX_BAMS_PER_CHUNK;
	                chunk->actual_value = valid_count + 1;
	                chunk->record_index = valid_count;
	                ks.s[ks.l] = saved_char;
	                break;
	            }
	            if (arena_used >= chunk->bam_data_capacity) {
	                chunk->status = -3;
	                chunk->limit_id = BOUNDS_LIMIT_INIT_DATA_SIZE;
	                chunk->limit_value = chunk->bam_data_capacity;
	                chunk->actual_value = arena_used + (size_t)ks.l + 64u;
	                chunk->record_index = valid_count;
	                ks.s[ks.l] = saved_char;
	                break;
	            }
	            bam1_t *b = chunk->bams + valid_count;
	            size_t arena_remain = chunk->bam_data_capacity - arena_used;
	            b->data = chunk->bam_data + arena_used;
	            b->m_data = arena_remain > UINT32_MAX ? UINT32_MAX : (uint32_t)arena_remain;
	            b->l_data = 0;
	            b->mempolicy = BAM_USER_OWNS_DATA;
            MpiSamParseFastTiming *timing_ptr = nullptr;
#if MPI_SAM_PARSE_DETAIL
            MpiSamParseFastTiming timing;
            timing.core_cycles = 0;
            timing.aux_cycles = 0;
            timing.cg_cycles = 0;
            timing_ptr = &timing;
#endif
            int ret = sam_parse1_mpi_fast(&ks, (sam_hdr_t *)batch->hdr, b, &fast_cache, timing_ptr);
            int parsed_by_fast = (ret == 0);
#if MPI_SAM_PARSE_DETAIL
            uint64_t fallback_cycles = 0;
#endif
            if (ret != 0) {
#if MPI_SAM_PARSE_DETAIL
                uint64_t fallback_t0 = MPI_SLAVE_RPCC();
#endif
	                b->data = chunk->bam_data + arena_used;
	                b->m_data = arena_remain > UINT32_MAX ? UINT32_MAX : (uint32_t)arena_remain;
	                b->l_data = 0;
	                b->mempolicy = BAM_USER_OWNS_DATA;
                ret = sam_parse1(&ks, (sam_hdr_t *)batch->hdr, b);
#if MPI_SAM_PARSE_DETAIL
                uint64_t fallback_t1 = MPI_SLAVE_RPCC();
                fallback_cycles = fallback_t1 - fallback_t0;
                chunk->parse_fallback_cycles += fallback_cycles;
#endif
            }
#if MPI_SAM_PARSE_DETAIL
            chunk->parse_core_cycles += timing.core_cycles;
            chunk->parse_aux_cycles += timing.aux_cycles;
            chunk->parse_cg_cycles += timing.cg_cycles;
            chunk->parse_total_cycles += timing.core_cycles + timing.aux_cycles +
                                         timing.cg_cycles + fallback_cycles;
#endif

            ks.s[ks.l] = saved_char;

	            if (ret >= 0) {
	                if (parsed_by_fast) {
	                    chunk->parse_fast_records++;
	                } else {
	                    chunk->parse_fallback_records++;
	                }
	                if ((size_t)b->l_data > arena_remain) {
	                    chunk->status = -3;
	                    chunk->limit_id = BOUNDS_LIMIT_INIT_DATA_SIZE;
	                    chunk->limit_value = arena_remain;
	                    chunk->actual_value = b->l_data;
	                    chunk->record_index = valid_count;
	                    break;
	                }
	                if (chunk->bam_offsets) {
	                    chunk->bam_offsets[valid_count] = (uint32_t)arena_used;
	                }
	                chunk->bam_lens[valid_count] = (uint32_t)(b->l_data - b->core.l_extranul + 32);
	                arena_used += ((size_t)b->l_data + 7u) & ~(size_t)7u;
	                chunk->bam_data_used = arena_used;
	                valid_count++;
	            } else {
	                chunk->status = -2;
	                chunk->record_index = valid_count;
	                chunk->actual_value = ret;
	                chunk->limit_id = BOUNDS_LIMIT_GENERIC_RUNTIME_ERROR;
	                break;
	            }
        }

        if (eol < end && *eol == '\n') eol++;
        ptr = eol;
    }
    chunk->count = valid_count;
}
