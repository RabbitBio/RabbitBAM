#include <cstring>
#include <stdint.h>

#include "BamTools.h"
#include "libdeflate.h"

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

#if defined(__GNUC__)
#define RB_SORT_LIKELY(x) __builtin_expect(!!(x), 1)
#define RB_SORT_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define RB_SORT_LIKELY(x) (x)
#define RB_SORT_UNLIKELY(x) (x)
#endif

#ifndef RABBITBAM_MPI_LEVEL1_NICE_LEN
#define RABBITBAM_MPI_LEVEL1_NICE_LEN 4
#endif

#ifndef RABBITBAM_MPI_LEVEL1_SINGLE_PROBE
#define RABBITBAM_MPI_LEVEL1_SINGLE_PROBE 1
#endif

#ifndef RABBITBAM_MPI_LEVEL1_STRIDE2_PROBE
#define RABBITBAM_MPI_LEVEL1_STRIDE2_PROBE 1
#endif

#if defined(PLATFORM_SUNWAY) && defined(RABBITBAM_ENABLE_SUNWAY_CRC16_LDM)
extern "C" void libdeflate_crc32_sunway_set_ldm_enabled(int enabled);
#endif

static inline unsigned long sort_slave_cycle_now() {
#ifdef PLATFORM_SUNWAY
    unsigned long counter = 0;
    asm volatile("rcsr %0, 4" : "=r"(counter));
    return counter;
#else
    return 0;
#endif
}

static inline uint16_t sort_read_le16(const unsigned char *p) {
    return (uint16_t)p[0] | ((uint16_t)p[1] << 8);
}

static inline uint32_t sort_read_le32(const unsigned char *p) {
    return (uint32_t)p[0] |
           ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) |
           ((uint32_t)p[3] << 24);
}

static inline int32_t sort_read_le32s(const unsigned char *p) {
    return (int32_t)sort_read_le32(p);
}

static struct libdeflate_compressor *g_sort_mpi_compressors[64] = {0};
static int g_sort_mpi_compressor_levels[64] = {0};
static struct libdeflate_decompressor *g_sort_decompressors[64] = {0};

static inline struct libdeflate_decompressor *sort_get_reused_decompressor(
        int id, uint64_t *alloc_cycles) {
    struct libdeflate_decompressor *z = nullptr;
    unsigned long t0 = 0;
    if (id >= 0 && id < 64) {
        z = g_sort_decompressors[id];
        if (z) return z;
        t0 = sort_slave_cycle_now();
        z = libdeflate_alloc_decompressor();
        if (alloc_cycles) *alloc_cycles += (uint64_t)(sort_slave_cycle_now() - t0);
        g_sort_decompressors[id] = z;
        return z;
    }
    t0 = sort_slave_cycle_now();
    z = libdeflate_alloc_decompressor();
    if (alloc_cycles) *alloc_cycles += (uint64_t)(sort_slave_cycle_now() - t0);
    return z;
}

static inline struct libdeflate_compressor *sort_get_reused_mpi_compressor(
        int id, int level, uint64_t *alloc_cycles) {
    struct libdeflate_compressor *z = nullptr;
    unsigned long t0 = 0;
    if (id >= 0 && id < 64) {
        z = g_sort_mpi_compressors[id];
        if (z && g_sort_mpi_compressor_levels[id] == level) return z;
        t0 = sort_slave_cycle_now();
        if (z) libdeflate_free_compressor(z);
        z = libdeflate_alloc_compressor(level);
        if (level == 1 && z) {
            libdeflate_set_level1_nice_match_length(z, RABBITBAM_MPI_LEVEL1_NICE_LEN);
#if RABBITBAM_MPI_LEVEL1_SINGLE_PROBE
            libdeflate_set_level1_single_probe(z, 1);
#endif
#if RABBITBAM_MPI_LEVEL1_STRIDE2_PROBE
            libdeflate_set_level1_stride2_probe(z, 1);
#endif
        }
        if (alloc_cycles) *alloc_cycles += (uint64_t)(sort_slave_cycle_now() - t0);
        if (z) {
            g_sort_mpi_compressors[id] = z;
            g_sort_mpi_compressor_levels[id] = level;
        } else {
            g_sort_mpi_compressors[id] = nullptr;
            g_sort_mpi_compressor_levels[id] = 0;
        }
        return z;
    }

    t0 = sort_slave_cycle_now();
    z = libdeflate_alloc_compressor(level);
    if (level == 1 && z) {
        libdeflate_set_level1_nice_match_length(z, RABBITBAM_MPI_LEVEL1_NICE_LEN);
#if RABBITBAM_MPI_LEVEL1_SINGLE_PROBE
        libdeflate_set_level1_single_probe(z, 1);
#endif
#if RABBITBAM_MPI_LEVEL1_STRIDE2_PROBE
        libdeflate_set_level1_stride2_probe(z, 1);
#endif
    }
    if (alloc_cycles) *alloc_cycles += (uint64_t)(sort_slave_cycle_now() - t0);
    return z;
}

static int sort_bgzf_uncompress_reuse(uint8_t *dst, size_t *dlen,
                                      const uint8_t *src, size_t slen,
                                      uint32_t expected_crc,
                                      struct libdeflate_decompressor *z,
                                      uint64_t *inflate_cycles,
                                      uint64_t *crc_cycles) {
    if (!z) return -1;
    unsigned long inflate_t0 = sort_slave_cycle_now();
    int ret = libdeflate_deflate_decompress(z, src, slen, dst, *dlen, dlen);
    if (inflate_cycles) *inflate_cycles += (uint64_t)(sort_slave_cycle_now() - inflate_t0);
    if (ret != 0) return -1;

    unsigned long crc_t0 = sort_slave_cycle_now();
    uint32_t crc = libdeflate_crc32(0, (unsigned char *)dst, *dlen);
    if (crc_cycles) *crc_cycles += (uint64_t)(sort_slave_cycle_now() - crc_t0);
    return crc == expected_crc ? 0 : -2;
}

static int sort_block_decode_func_reuse(struct bam_block *comp, struct bam_block *un_comp,
                                        struct libdeflate_decompressor *z,
                                        uint64_t *inflate_cycles,
                                        uint64_t *crc_cycles) {
    size_t out_capacity = (un_comp->length > 0 && un_comp->length <= BGZF_MAX_BLOCK_SIZE)
                          ? (size_t)un_comp->length
                          : (size_t)BGZF_MAX_BLOCK_SIZE;
    un_comp->pos = 0;
    uint32_t crc = sort_read_le32((unsigned char *)comp->data + comp->length - 8);
    size_t un_comp_len = out_capacity;
    int ret = sort_bgzf_uncompress_reuse(un_comp->data, &un_comp_len,
                                         comp->data + BLOCK_HEADER_LENGTH,
                                         comp->length - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH,
                                         crc, z, inflate_cycles, crc_cycles);
    un_comp->length = (unsigned int)un_comp_len;
    if (ret != 0) un_comp->errcode |= BGZF_ERR_ZLIB;
    return ret;
}

static int sort_bgzf_compress_reuse(void *_dst, size_t *dlen,
                                    const void *src, size_t slen,
                                    int level,
                                    struct libdeflate_compressor *z,
                                    uint64_t *deflate_cycles,
                                    uint64_t *footer_cycles) {
    if (slen == 0) return 0;
    if (level != 0 && !z) return -1;

    uint8_t *dst = (uint8_t *)_dst;
    if (level == 0) {
        if (*dlen < slen + 5 + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH) return -1;
        dst[BLOCK_HEADER_LENGTH] = 1;
        packInt16(&dst[BLOCK_HEADER_LENGTH + 1], (uint16_t)slen);
        packInt16(&dst[BLOCK_HEADER_LENGTH + 3], (uint16_t)~(uint16_t)slen);
        memcpy(dst + BLOCK_HEADER_LENGTH + 5, src, slen);
        *dlen = slen + 5 + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    } else {
        unsigned long deflate_t0 = sort_slave_cycle_now();
        size_t clen = libdeflate_deflate_compress(
            z, src, slen, dst + BLOCK_HEADER_LENGTH,
            *dlen - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH);
        if (deflate_cycles) *deflate_cycles += (uint64_t)(sort_slave_cycle_now() - deflate_t0);
        if (clen <= 0) return -1;
        *dlen = clen + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    }

    unsigned long footer_t0 = sort_slave_cycle_now();
    memcpy(dst, g_magic, BLOCK_HEADER_LENGTH);
    packInt16(&dst[16], *dlen - 1);
#if defined(PLATFORM_SUNWAY) && defined(RABBITBAM_ENABLE_SUNWAY_CRC16_LDM)
    libdeflate_crc32_sunway_set_ldm_enabled(level != 0 && z != nullptr);
#endif
    uint32_t crc = libdeflate_crc32(0, src, slen);
#if defined(PLATFORM_SUNWAY) && defined(RABBITBAM_ENABLE_SUNWAY_CRC16_LDM)
    libdeflate_crc32_sunway_set_ldm_enabled(0);
#endif
    packInt32((uint8_t *)&dst[*dlen - 8], crc);
    packInt32((uint8_t *)&dst[*dlen - 4], slen);
    if (footer_cycles) *footer_cycles += (uint64_t)(sort_slave_cycle_now() - footer_t0);
    return 0;
}

extern "C" void slave_mpi_sort_extract_raw(MpiSortExtractPara paras[64]) {
    int id = _PEN;
    MpiSortExtractPara *para = &paras[id];

    para->decomp_alloc_cycles = 0;
    para->decomp_inflate_cycles = 0;
    para->decomp_crc_cycles = 0;
    para->decomp_parse_cycles = 0;
    para->decomp_total_cycles = 0;
    para->raw_used = 0;
    para->n_records = 0;

    bam_block *comp = para->input_block;
    bam_block *un_comp = para->un_comp_block;
    if (comp == NULL || para->status != 0) return;

    unsigned long total_t0 = sort_slave_cycle_now();
    struct libdeflate_decompressor *z =
        sort_get_reused_decompressor(id, &para->decomp_alloc_cycles);
    if (!z) {
        para->status = -2;
        para->decomp_total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
        return;
    }

    if (sort_block_decode_func_reuse(comp, un_comp, z,
                                     &para->decomp_inflate_cycles,
                                     &para->decomp_crc_cycles) != 0) {
        para->status = -2;
        para->decomp_total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
        return;
    }

    unsigned long parse_t0 = sort_slave_cycle_now();
    size_t pos = 0;
    int count = 0;
    while (pos < un_comp->length) {
        if (RB_SORT_UNLIKELY(un_comp->length - pos < 4)) {
            para->status = -4;
            para->record_index = count;
            para->actual_value = (long long)(un_comp->length - pos);
            para->limit_value = 4;
            para->limit_id = BOUNDS_LIMIT_BGZF_RECORD_SIZE;
            para->decomp_parse_cycles = (uint64_t)(sort_slave_cycle_now() - parse_t0);
            para->decomp_total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
            return;
        }
        uint32_t block_size = sort_read_le32(un_comp->data + pos);
        if (RB_SORT_UNLIKELY(block_size < 32)) {
            para->status = -2;
            para->record_index = count;
            para->decomp_parse_cycles = (uint64_t)(sort_slave_cycle_now() - parse_t0);
            para->decomp_total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
            return;
        }
        size_t raw_len = (size_t)block_size + 4u;
        if (RB_SORT_UNLIKELY(raw_len > un_comp->length - pos)) {
            para->status = -4;
            para->record_index = count;
            para->actual_value = (long long)raw_len;
            para->limit_value = (long long)(un_comp->length - pos);
            para->limit_id = BOUNDS_LIMIT_BGZF_RECORD_SIZE;
            para->decomp_parse_cycles = (uint64_t)(sort_slave_cycle_now() - parse_t0);
            para->decomp_total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
            return;
        }
        if (RB_SORT_UNLIKELY(raw_len > para->raw_capacity - para->raw_used)) {
            para->status = -3;
            para->record_index = count;
            para->actual_value = (long long)raw_len;
            para->limit_value = (long long)(para->raw_capacity - para->raw_used);
            para->limit_id = BOUNDS_LIMIT_INIT_DATA_SIZE;
            para->decomp_parse_cycles = (uint64_t)(sort_slave_cycle_now() - parse_t0);
            para->decomp_total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
            return;
        }
        if (RB_SORT_UNLIKELY(count >= para->record_capacity)) {
            para->status = -3;
            para->record_index = count;
            para->actual_value = count + 1;
            para->limit_value = para->record_capacity;
            para->limit_id = BOUNDS_LIMIT_MAX_RECORDS_PER_BLOCK;
            para->decomp_parse_cycles = (uint64_t)(sort_slave_cycle_now() - parse_t0);
            para->decomp_total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
            return;
        }

        const unsigned char *body = un_comp->data + pos + 4;
        MpiSortRecordMetaShared *meta = &para->records[count];
        meta->tid = sort_read_le32s(body);
        meta->pos = sort_read_le32s(body + 4);
        meta->flag = sort_read_le16(body + 14);
        meta->pad = 0;
        meta->raw_len = (uint32_t)raw_len;
        meta->pad2 = 0;
        meta->raw_offset = para->raw_base_offset + (uint64_t)para->raw_used;
        meta->global_order = ((uint64_t)para->global_block_index << 32) | (uint32_t)count;
        if (para->raw_arena + para->raw_used != un_comp->data + pos) {
            memcpy(para->raw_arena + para->raw_used, un_comp->data + pos, raw_len);
        }
        para->raw_used += raw_len;
        pos += raw_len;
        count++;
    }

    para->n_records = count;
    para->status = 0;
    para->decomp_parse_cycles = (uint64_t)(sort_slave_cycle_now() - parse_t0);
    para->decomp_total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
}

extern "C" void slave_mpi_sort_compress_payload(MpiSortRawCompressPara paras[64]) {
    int id = _PEN;
    MpiSortRawCompressPara *para = &paras[id];

    if (para->status != 0 || para->un_comp_block == nullptr || para->un_comp_size <= 0) return;
    int compress_level = para->compress_level;
    if (compress_level != 0 && compress_level != 1 && compress_level != 6) {
        para->status = -2;
        return;
    }

    para->compress_pack_cycles = 0;
    para->compress_alloc_cycles = 0;
    para->compress_deflate_cycles = 0;
    para->compress_footer_cycles = 0;
    para->compress_total_cycles = 0;

    unsigned long total_t0 = sort_slave_cycle_now();
    struct libdeflate_compressor *z = nullptr;
    if (compress_level != 0) {
        z = sort_get_reused_mpi_compressor(id, compress_level, &para->compress_alloc_cycles);
    }
    if (compress_level != 0 && !z) {
        para->status = -2;
        para->compress_total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
        return;
    }

    bam_block *uncompressed = para->un_comp_block;
    bam_block *compressed = para->output_block;
    size_t comp_size = BGZF_MAX_BLOCK_SIZE;
    int ret = sort_bgzf_compress_reuse(compressed->data, &comp_size,
                                       uncompressed->data, (size_t)para->un_comp_size,
                                       compress_level, z,
                                       &para->compress_deflate_cycles,
                                       &para->compress_footer_cycles);
    compressed->length = ret == 0 ? (int)comp_size : -1;
    compressed->pos = 0;
    compressed->errcode = 0;
    compressed->block_id = para->block_id;

    if (compressed->length <= 0) {
        para->status = -2;
        para->compress_total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
        return;
    }

    para->output_size = compressed->length;
    para->compress_total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
    para->status = 0;
}

extern "C" void slave_mpi_sort_bucket_pack(MpiSortBucketPackPara paras[64]) {
    int id = _PEN;
    MpiSortBucketPackPara *para = &paras[id];

    para->pack_cycles = 0;
    para->total_cycles = 0;
    if (para->status != 0) return;
    if (para->record_begin >= para->record_end) {
        para->status = 0;
        return;
    }
    if (para->local_records == nullptr || para->local_raw == nullptr ||
        para->bucket_ids == nullptr || para->record_meta_in_bucket == nullptr ||
        para->record_raw_in_bucket == nullptr || para->send_meta_displs == nullptr ||
        para->send_raw_displs == nullptr || para->send_meta == nullptr ||
        para->send_raw == nullptr) {
        para->status = -2;
        return;
    }

    unsigned long total_t0 = sort_slave_cycle_now();
    unsigned long pack_t0 = sort_slave_cycle_now();
    for (size_t i = para->record_begin; i < para->record_end; ++i) {
        const MpiSortRecordMetaShared *src = &para->local_records[i];
        int bucket = para->bucket_ids[i];
        uint64_t meta_pos = (uint64_t)para->send_meta_displs[bucket] +
                            para->record_meta_in_bucket[i];
        uint64_t raw_in_bucket = para->record_raw_in_bucket[i];
        uint64_t raw_pos = (uint64_t)para->send_raw_displs[bucket] + raw_in_bucket;
        uint64_t src_end = src->raw_offset + (uint64_t)src->raw_len;
        uint64_t dst_end = raw_pos + (uint64_t)src->raw_len;

        if (RB_SORT_UNLIKELY(src_end > (uint64_t)para->local_raw_size ||
                             meta_pos >= para->send_meta_capacity ||
                             dst_end > para->send_raw_capacity)) {
            para->status = -3;
            para->record_index = (int)i;
            para->actual_value = (long long)(dst_end > (uint64_t)para->local_raw_size ? dst_end : src_end);
            para->limit_value = (long long)(dst_end > para->send_raw_capacity
                                            ? para->send_raw_capacity
                                            : para->local_raw_size);
            para->limit_id = BOUNDS_LIMIT_GENERIC_RUNTIME_ERROR;
            para->pack_cycles = (uint64_t)(sort_slave_cycle_now() - pack_t0);
            para->total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
            return;
        }

        MpiSortRecordMetaShared dst = *src;
        dst.raw_offset = raw_in_bucket;
        para->send_meta[meta_pos] = dst;
        memcpy(para->send_raw + raw_pos, para->local_raw + src->raw_offset, src->raw_len);
    }
    para->pack_cycles = (uint64_t)(sort_slave_cycle_now() - pack_t0);
    para->total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
    para->status = 0;
}

extern "C" void slave_mpi_sort_range_pack(MpiSortRangePackPara paras[64]) {
    int id = _PEN;
    MpiSortRangePackPara *para = &paras[id];

    para->pack_cycles = 0;
    para->total_cycles = 0;
    if (para->status != 0) return;
    if (para->record_begin >= para->record_end) {
        para->status = 0;
        return;
    }
    if (para->local_records == nullptr || para->local_raw == nullptr ||
        para->record_raw_offsets == nullptr || para->send_meta == nullptr ||
        para->send_raw == nullptr) {
        para->status = -2;
        return;
    }

    unsigned long total_t0 = sort_slave_cycle_now();
    unsigned long pack_t0 = sort_slave_cycle_now();
    for (size_t i = para->record_begin; i < para->record_end; ++i) {
        const MpiSortRecordMetaShared *src = &para->local_records[i];
        size_t out_index = para->output_record_begin + (i - para->record_begin);
        uint64_t raw_pos = para->record_raw_offsets[out_index];
        uint64_t src_end = src->raw_offset + (uint64_t)src->raw_len;
        uint64_t dst_end = raw_pos + (uint64_t)src->raw_len;
        if (RB_SORT_UNLIKELY(src_end > (uint64_t)para->local_raw_size ||
                             out_index >= para->send_meta_capacity ||
                             dst_end > para->send_raw_capacity)) {
            para->status = -3;
            para->record_index = (int)i;
            para->actual_value = (long long)(src_end > dst_end ? src_end : dst_end);
            para->limit_value = (long long)(src_end > (uint64_t)para->local_raw_size
                                            ? para->local_raw_size
                                            : para->send_raw_capacity);
            para->limit_id = BOUNDS_LIMIT_GENERIC_RUNTIME_ERROR;
            para->pack_cycles = (uint64_t)(sort_slave_cycle_now() - pack_t0);
            para->total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
            return;
        }
        MpiSortRecordMetaShared dst = *src;
        dst.raw_offset = raw_pos;
        para->send_meta[out_index] = dst;
        memcpy(para->send_raw + raw_pos, para->local_raw + src->raw_offset, src->raw_len);
    }
    para->pack_cycles = (uint64_t)(sort_slave_cycle_now() - pack_t0);
    para->total_cycles = (uint64_t)(sort_slave_cycle_now() - total_t0);
    para->status = 0;
}
