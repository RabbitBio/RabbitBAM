#include "BamTools.h"


//读取一个压缩块，放到j中
int read_block(BGZF *fp, struct bam_block *j) {
    uint8_t header[BLOCK_HEADER_LENGTH], *compressed_block;
    int count, size = 0, block_length, remaining;

    // NOTE: Guaranteed to be compressed as we block multi-threading in
    // uncompressed mode.  However it may be gzip compression instead
    // of bgzf.

    // Reading compressed file
    int64_t block_address;
    block_address = htell(fp->fp);

    if (fp->cache_size && load_block_from_cache(fp, block_address)) return 0;
    count = hpeek(fp->fp, header, sizeof(header));
    if (count == 0) // no data read
        return -1;
    int ret;
    if (count != sizeof(header) || (ret = check_header(header)) == -2) {
        j->errcode |= BGZF_ERR_HEADER;
        return -1;
    }
    if (ret == -1) {
        j->errcode |= BGZF_ERR_MT;
        return -1;
    }

    count = hread(fp->fp, header, sizeof(header));
    if (count != sizeof(header)) // no data read
        return -1;

    size = count;
    block_length = unpackInt16((uint8_t * ) & header[16]) + 1; // +1 because when writing this number, we used "-1"
    if (block_length < BLOCK_HEADER_LENGTH) {
        j->errcode |= BGZF_ERR_HEADER;
        return -1;
    }

    compressed_block = (uint8_t *) j->data;
    memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
    remaining = block_length - BLOCK_HEADER_LENGTH;
    count = hread(fp->fp, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
    if (count != remaining) {
        j->errcode |= BGZF_ERR_IO;
        return -1;
    }
    size += count;
    j->length = block_length;
    j->block_address = block_address;
    j->errcode = 0;
    return 0;
}


int sam_write1_sw(samFile *fp, const sam_hdr_t *h, const bam1_t *b){

    //调用 sam_format1() 将 BAM 结构转换为一行 SAM 文本，写入 kstring_t fp->line
    if (sam_format1(h, b, &fp->line) < 0) return -1;
    kputc('\n', &fp->line);
    //写入 SAM 文本数据
    if ( hwrite(fp->fp.hfile, fp->line.s, fp->line.l) != fp->line.l ) return -1;      

    return fp->line.l;

}


//--------------------------------------------------------------------------------------------------------------


int sam_read1_sw(samFile *fp, const sam_hdr_t *h, const bam1_t *b){
    int ret = sam_read1_sam(fp, h, b);
    return ret;
}
