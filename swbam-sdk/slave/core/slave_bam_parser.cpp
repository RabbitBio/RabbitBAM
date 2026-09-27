#include "swbam/cpe_bam_parser.h"

#include <climits>
#include <cstdio>
#include <cstring>

#include <htslib/sam.h>

#if defined(__GNUC__)
#define SWBAM_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define SWBAM_UNLIKELY(x) (x)
#endif

int swbam_cpe_read_bam_record(bam_block *decoded, bam1_t *record_out) {
    if (SWBAM_UNLIKELY(decoded->pos >= decoded->length)) return -1;
    bam1_core_t *core = &record_out->core;
    const unsigned int start = decoded->pos;
    const unsigned int remaining = decoded->length - start;
    const uint8_t *record = decoded->data + start;
    const uint8_t *payload;
    uint32_t raw_qname_len;
    uint32_t rest_len;
    int32_t block_len;
    uint32_t fields[8], new_data_len;
    int cigar_changed = 0;

    record_out->l_data = 0;
    if (SWBAM_UNLIKELY(remaining < 4)) return -2;
    memcpy(&block_len, record, 4);
    if (SWBAM_UNLIKELY(block_len < 32)) return -4;
    if (SWBAM_UNLIKELY(remaining < 36)) return -3;
    if (SWBAM_UNLIKELY((uint64_t)block_len + 4 > remaining)) return -4;

    memcpy(fields, record + 4, 32);
    core->tid = fields[0];
    core->pos = (int32_t)fields[1];
    core->bin = fields[2] >> 16;
    core->qual = fields[2] >> 8 & 0xff;
    core->l_qname = fields[2] & 0xff;
    core->l_extranul = (-core->l_qname) & 3;
    core->flag = fields[3] >> 16;
    core->n_cigar = fields[3] & 0xffff;
    core->l_qseq = fields[4];
    core->mtid = fields[5];
    core->mpos = (int32_t)fields[6];
    core->isize = (int32_t)fields[7];

    raw_qname_len = core->l_qname;
    new_data_len = block_len - 32 + core->l_extranul;
    if (SWBAM_UNLIKELY(new_data_len > INT_MAX || core->l_qseq < 0 ||
                       raw_qname_len < 1)) {
        return -4;
    }
    if (SWBAM_UNLIKELY(((uint64_t)core->n_cigar << 2) + raw_qname_len +
                       core->l_extranul +
                       (((uint64_t)core->l_qseq + 1) >> 1) +
                       core->l_qseq > new_data_len)) {
        return -4;
    }
    if (SWBAM_UNLIKELY(new_data_len > record_out->m_data)) {
        if (record_out->m_data == INIT_DATA_SIZE &&
            realloc_bam_data(record_out, new_data_len) == 0) {
            record_out->l_data = new_data_len;
            goto data_ready;
        }
        record_out->l_data = new_data_len;
        return -5;
    }
    record_out->l_data = new_data_len;

data_ready:
    payload = record + 36;
    if (core->l_extranul == 0 && payload[raw_qname_len - 1] == '\0') {
        memcpy(record_out->data, payload, new_data_len);
    } else {
        const int qname_has_nul = payload[raw_qname_len - 1] == '\0';
        memcpy(record_out->data, payload, raw_qname_len);
        if (SWBAM_UNLIKELY(!qname_has_nul) &&
            SWBAM_UNLIKELY(fixup_missing_qname_nul(record_out) < 0)) {
            return -4;
        }
        switch (core->l_extranul) {
        case 3:
            record_out->data[core->l_qname + 2] = '\0';
            /* fall through */
        case 2:
            record_out->data[core->l_qname + 1] = '\0';
            /* fall through */
        case 1:
            record_out->data[core->l_qname] = '\0';
            break;
        default:
            break;
        }
        core->l_qname += core->l_extranul;
        rest_len = (uint32_t)block_len - 32 - raw_qname_len;
        if (SWBAM_UNLIKELY(record_out->l_data < core->l_qname ||
                           (uint64_t)core->l_qname + rest_len >
                               (uint64_t)record_out->l_data)) {
            return -4;
        }
        memcpy(record_out->data + core->l_qname,
               payload + raw_qname_len, rest_len);
    }
    decoded->pos = start + 4 + (unsigned int)block_len;

    if (core->n_cigar != 0 && core->tid >= 0 && core->pos >= 0) {
        uint32_t *cigar = bam_get_cigar(record_out);
        if (SWBAM_UNLIKELY(bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP &&
                            bam_cigar_oplen(cigar[0]) == core->l_qseq)) {
            cigar_changed = bam_tag2cigar(record_out, 0, 0);
            if (SWBAM_UNLIKELY(cigar_changed < 0)) return -4;
        }
    }
    if (cigar_changed > 0 && core->n_cigar > 0) {
        hts_pos_t ref_len, query_len;
        bam_cigar2rqlens(core->n_cigar, bam_get_cigar(record_out),
                         &ref_len, &query_len);
        if ((core->flag & BAM_FUNMAP) || ref_len == 0) ref_len = 1;
        core->bin = hts_reg2bin(core->pos, core->pos + ref_len, 14, 5);
        if (SWBAM_UNLIKELY(core->l_qseq > 0 &&
                           !(core->flag & BAM_FUNMAP) &&
                           query_len != core->l_qseq)) {
            fprintf(stderr, "CIGAR/query length mismatch for %s\n",
                    bam_get_qname(record_out));
            return -4;
        }
    }
    return 4 + block_len;
}
