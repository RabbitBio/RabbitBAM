#include <cstring>
#include <stdint.h>

#include "BamTools.h"
#include "libdeflate.h"

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

#if defined(__GNUC__)
#define RB_COLLATE_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define RB_COLLATE_UNLIKELY(x) (x)
#endif

static struct libdeflate_decompressor *g_collate_decompressors[64] = {0};

static inline unsigned long collate_cycle_now() {
#ifdef PLATFORM_SUNWAY
    unsigned long counter = 0;
    asm volatile("rcsr %0, 4" : "=r"(counter));
    return counter;
#else
    return 0;
#endif
}

static inline uint16_t collate_read_le16(const unsigned char *p) {
    return (uint16_t)p[0] | ((uint16_t)p[1] << 8);
}

static inline uint32_t collate_read_le32(const unsigned char *p) {
    return (uint32_t)p[0] |
           ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) |
           ((uint32_t)p[3] << 24);
}

static inline uint32_t collate_hash_wang(uint32_t key) {
    key += ~(key << 15);
    key ^= key >> 10;
    key += key << 3;
    key ^= key >> 6;
    key += ~(key << 11);
    key ^= key >> 16;
    return key;
}

static inline uint32_t collate_hash_qname(
        const unsigned char *qname, size_t max_len,
        uint16_t *qname_len) {
    size_t len = 0;
    uint32_t h = max_len > 0 ? qname[0] : 0;
    if (h == 0) {
        *qname_len = 0;
        return 0;
    }
    len = 1;
    while (len < max_len && qname[len] != 0) {
        h = (h << 5) - h + qname[len];
        ++len;
    }
    *qname_len = (uint16_t)len;
    return collate_hash_wang(h);
}

static inline struct libdeflate_decompressor *
collate_get_decompressor(int id, uint64_t *alloc_cycles) {
    struct libdeflate_decompressor *z = g_collate_decompressors[id];
    if (z) return z;
    unsigned long t0 = collate_cycle_now();
    z = libdeflate_alloc_decompressor();
    if (alloc_cycles) {
        *alloc_cycles +=
            (uint64_t)(collate_cycle_now() - t0);
    }
    g_collate_decompressors[id] = z;
    return z;
}

static int collate_decode_block(
        bam_block *input, bam_block *output,
        struct libdeflate_decompressor *z,
        uint64_t *inflate_cycles,
        uint64_t *crc_cycles) {
    if (!input || !output || !z ||
        input->length < BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH) {
        return -1;
    }
    size_t capacity =
        output->length > 0 &&
        output->length <= BGZF_MAX_BLOCK_SIZE
            ? (size_t)output->length
            : (size_t)BGZF_MAX_BLOCK_SIZE;
    size_t output_len = capacity;
    unsigned long t0 = collate_cycle_now();
    int ret = libdeflate_deflate_decompress(
        z, input->data + BLOCK_HEADER_LENGTH,
        input->length - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH,
        output->data, output_len, &output_len);
    if (inflate_cycles) {
        *inflate_cycles +=
            (uint64_t)(collate_cycle_now() - t0);
    }
    if (ret != 0) return -1;

    uint32_t expected_crc = collate_read_le32(
        input->data + input->length - 8);
    t0 = collate_cycle_now();
    uint32_t actual_crc =
        libdeflate_crc32(0, output->data, output_len);
    if (crc_cycles) {
        *crc_cycles +=
            (uint64_t)(collate_cycle_now() - t0);
    }
    if (actual_crc != expected_crc) return -1;
    output->length = (unsigned int)output_len;
    output->pos = 0;
    output->errcode = 0;
    return 0;
}

extern "C" void slave_mpi_collate_extract(
        MpiCollateExtractPara paras[64]) {
    const int id = _PEN;
    MpiCollateExtractPara *para = &paras[id];
    para->raw_used = 0;
    para->n_records = 0;
    para->decomp_alloc_cycles = 0;
    para->decomp_inflate_cycles = 0;
    para->decomp_crc_cycles = 0;
    para->decomp_parse_cycles = 0;
    para->decomp_total_cycles = 0;
    if (para->status != 0 || !para->input_block ||
        !para->un_comp_block || !para->raw_arena ||
        !para->records || para->n_bins <= 0) {
        return;
    }

    unsigned long total_t0 = collate_cycle_now();
    struct libdeflate_decompressor *z =
        collate_get_decompressor(id, &para->decomp_alloc_cycles);
    if (!z || collate_decode_block(
                  para->input_block, para->un_comp_block, z,
                  &para->decomp_inflate_cycles,
                  &para->decomp_crc_cycles) != 0) {
        para->status = -2;
        para->decomp_total_cycles =
            (uint64_t)(collate_cycle_now() - total_t0);
        return;
    }

    unsigned long parse_t0 = collate_cycle_now();
    size_t pos = 0;
    int count = 0;
    while (pos < para->un_comp_block->length) {
        if (RB_COLLATE_UNLIKELY(
                para->un_comp_block->length - pos < 4)) {
            para->status = -4;
            break;
        }
        uint32_t block_size =
            collate_read_le32(para->un_comp_block->data + pos);
        size_t raw_len = (size_t)block_size + 4;
        if (RB_COLLATE_UNLIKELY(
                block_size < 32 ||
                raw_len > para->un_comp_block->length - pos)) {
            para->status = -4;
            break;
        }
        if (RB_COLLATE_UNLIKELY(
                count >= para->record_capacity ||
                raw_len > para->raw_capacity - para->raw_used)) {
            para->status = -3;
            break;
        }

        const unsigned char *record =
            para->un_comp_block->data + pos;
        const unsigned char *body = record + 4;
        const uint8_t l_qname = body[8];
        const uint16_t flag = collate_read_le16(body + 14);
        if (RB_COLLATE_UNLIKELY(
                l_qname == 0 || raw_len < 36u + l_qname)) {
            para->status = -2;
            break;
        }
        const unsigned char *qname = record + 36;
        uint16_t qname_len = 0;
        uint32_t hash =
            collate_hash_qname(qname, l_qname, &qname_len);
        if (RB_COLLATE_UNLIKELY(
                qname_len >= l_qname || qname[qname_len] != 0)) {
            para->status = -2;
            break;
        }

        MpiCollateRecordMetaShared *meta =
            &para->records[count];
        meta->hash = hash;
        meta->bin = hash % (uint32_t)para->n_bins;
        meta->qname_len = qname_len;
        meta->flag_order =
            (uint8_t)((flag >> 6) & 3u);
        meta->pad = 0;
        meta->raw_len = (uint32_t)raw_len;
        meta->raw_offset =
            para->raw_base_offset + para->raw_used;
        meta->global_order =
            ((uint64_t)para->global_block_index << 32) |
            (uint32_t)count;
        if (para->raw_arena + para->raw_used != record) {
            memcpy(para->raw_arena + para->raw_used,
                   record, raw_len);
        }
        para->raw_used += raw_len;
        pos += raw_len;
        ++count;
    }
    para->record_index = count;
    para->n_records = count;
    para->decomp_parse_cycles =
        (uint64_t)(collate_cycle_now() - parse_t0);
    para->decomp_total_cycles =
        (uint64_t)(collate_cycle_now() - total_t0);
    if (para->status == 0 && pos != para->un_comp_block->length) {
        para->status = -4;
    }
}
