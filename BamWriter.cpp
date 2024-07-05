//
// Created by 赵展 on 2022/8/17.
//

#include "BamWriter.h"
int rabbit_write_deflate_block(BGZF *fp, bam_write_block* write_block){
    size_t comp_size = BGZF_MAX_BLOCK_SIZE;
    int ret;
    if ( !fp->is_gzip )
        ret = rabbit_bgzf_compress(write_block->compressed_data, &comp_size, write_block->uncompressed_data, write_block->block_offset, fp->compress_level);
    else
        ret = rabbit_bgzf_gzip_compress(fp, write_block->compressed_data, &comp_size, write_block->uncompressed_data, write_block->block_offset, fp->compress_level);

    if ( ret != 0 )
    {
        hts_log_debug("Compression error %d", ret);
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    return comp_size;
}
int rabbit_bgzf_flush(BGZF *fp,bam_write_block* write_block)
{
    //TODO 此处可能会出现问题
    while (write_block->block_offset > 0) {
        int block_length;
        printf("Write Block Offset : %d\n",write_block->block_offset);
        block_length = rabbit_write_deflate_block(fp, write_block);
        if (block_length < 0) {
            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block_length, NULL));
            return -1;
        }
        if (write_block== nullptr){
            printf("Write Game Over!!!\n");
        }
        if (hwrite(fp->fp, write_block->compressed_data, block_length) != block_length) {
            printf("Write Failed\n");
            hts_log_error("File write failed (wrong size)");
            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
            return -1;
        }


        write_block->block_offset=0;
        fp->block_address += block_length;
    }
    write_block->block_offset=0;
    return 0;
}
int rabbit_bgzf_mul_flush(BGZF *fp,BamWriteCompress *bam_write_compress,bam_write_block* &write_block)
{
//    printf("Try to input One Uncompressed data\n");
//    InputBlockNum++;
//    printf("Write Block Offset : %d\n",write_block->block_offset);
    bam_write_compress->inputUnCompressData(write_block);
    write_block=bam_write_compress->getEmpty();
//    printf("Get Another Empty Block Block Num : %d\n",write_block->block_num);
    return 0;
}
int rabbit_bgzf_write(BGZF *fp,bam_write_block* &write_block,const void *data, size_t length)
{
    const uint8_t *input = (const uint8_t*)data;
    ssize_t remaining = length;
//    assert(fp->is_write);
    while (remaining > 0) {
        uint8_t* buffer = (uint8_t*)write_block->uncompressed_data;
        int copy_length = BGZF_BLOCK_SIZE - write_block->block_offset;
        if (copy_length > remaining) copy_length = remaining;
        memcpy(buffer + write_block->block_offset, input, copy_length);
        write_block->block_offset += copy_length;
        input += copy_length;
        remaining -= copy_length;
        if (write_block->block_offset == BGZF_BLOCK_SIZE) {
            //BUG  write block 不是 引用，所以没有改变
            if (rabbit_bgzf_flush(fp,write_block) != 0) return -1;
        }
    }
    return length - remaining;
}
int rabbit_bgzf_mul_write(BGZF *fp, BamWriteCompress *bam_write_compress,bam_write_block* &write_block,const void *data, size_t length)
{
    const uint8_t *input = (const uint8_t*)data;
    ssize_t remaining = length;
//    assert(fp->is_write);
    while (remaining > 0) {
        uint8_t* buffer = (uint8_t*)write_block->uncompressed_data;
//        printf("In rabbit bgzf mul write block Block Num : %d\n",write_block->block_num);
        int copy_length = BGZF_BLOCK_SIZE - write_block->block_offset;
        if (copy_length > remaining) copy_length = remaining;
//        printf("In rabbit bgzf mul write block Block Num : %d\n",write_block->block_num);
        memcpy(buffer + write_block->block_offset, input, copy_length);
        write_block->block_offset += copy_length;
//        printf("In rabbit bgzf mul write block Block Num : %d\n",write_block->block_num);
        input += copy_length;
        remaining -= copy_length;
        if (write_block->block_offset == BGZF_BLOCK_SIZE) {
            if (rabbit_bgzf_mul_flush(fp,bam_write_compress,write_block) != 0) return -1;
        }
    }
    return length - remaining;
}
int rabbit_bgzf_flush_try(BGZF *fp, bam_write_block* write_block,ssize_t size)
{
    if (write_block->block_offset + size > BGZF_BLOCK_SIZE) {
        return rabbit_bgzf_flush(fp,write_block);
    }
    return 0;
}
int rabbit_bgzf_mul_flush_try(BGZF *fp,BamWriteCompress* bam_write_compress,bam_write_block* &write_block,ssize_t size)
{
    if (write_block->block_offset + size > BGZF_BLOCK_SIZE) {
        return rabbit_bgzf_mul_flush(fp,bam_write_compress,write_block);
    }
    return 0;
}
int bam_write_pack(BGZF *fp,BamWriteCompress *bam_write_compress){
    bam_write_block* block;
    while(1){
//        printf("Try to Get Compress Data\n");
//
        block=bam_write_compress->getCompressData();
//        printf("Has Get One Compress Data %d\n",block->block_num);
        if (block == nullptr){
            return 0;
        }
        int block_offset = block->block_offset;
//        printf("BlockLength = %d\n",block_offset);
//        bam_write_compress->wait_num++;
//        printf("Block Num = %d\n",bam_write_compress->wait_num);
//        std::this_thread::sleep_for(std::chrono::nanoseconds(5));
        if (block->block_length < 0) {
            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block->block_length, NULL));
            return -1;
        }
        if (hwrite(fp->fp, block->compressed_data, block->block_length) != block->block_length) {
//            printf("Write Failed\n");
            hts_log_error("File write failed (wrong size)");
            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
            return -1;
        }
//        printf("Has write One Block\n");
        block->block_offset=0;
        fp->block_address += block->block_length;
        bam_write_compress->backEmpty(block);
    }

}
void bam_write_compress_pack(BGZF *fp,BamWriteCompress *bam_write_compress){
//    printf("Start Compress\n");
    bam_write_block * block;
    while (1){
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
//        printf("Has Start Try to Get One Uncompress\n");
        block=bam_write_compress->getUnCompressData();
//        printf("Has get One Uncompress data\n");

        //printf("%d is get One compressed data\n",id);
        if (block == nullptr) {
            //printf("%d is Over\n",id);
            break;
        }
//        printf("This Uncompress data block num : %d\n",block->block_num);
        /*
         * 压缩
         */
//        int block_num=block->block_num;
        block->block_length = rabbit_write_deflate_block(fp, block);
//        printf("Has Compress One Block\n");
        bam_write_compress->inputCompressData(block);
//        printf("Can,t input Compress Data\n");
//        while (!compress->tryinputUnCompressData(un_comp,comp.second)){
//            std::this_thread::sleep_for(std::chrono::milliseconds(1));
//        }
    }

//    if (block->block_length>0){
//        bam_write_compress->inputUnCompressData(block);
//    }
//    printf("One Compress Thread has been over!\n");
    bam_write_compress->CompressThreadComplete();
}

int rabbit_bam_write_test(BGZF *fp,bam_write_block* write_block,bam1_t *b){
    const bam1_core_t *c = &b->core;
    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
    int i, ok;
    if (c->l_qname - c->l_extranul > 255) {
        hts_log_error("QNAME \"%s\" is longer than 254 characters", bam_get_qname(b));
        errno = EOVERFLOW;
        return -1;
    }
    if (c->n_cigar > 0xffff) block_len += 16; // "16" for "CGBI", 4-byte tag length and 8-byte fake CIGAR
    if (c->pos > INT_MAX ||
        c->mpos > INT_MAX ||
        c->isize < INT_MIN || c->isize > INT_MAX) {
        hts_log_error("Positional data is too large for BAM format");
        return -1;
    }
    x[0] = c->tid;
    x[1] = c->pos;
    x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | (c->l_qname - c->l_extranul);
    if (c->n_cigar > 0xffff) x[3] = (uint32_t)c->flag << 16 | 2;
    else x[3] = (uint32_t)c->flag << 16 | (c->n_cigar & 0xffff);
    x[4] = c->l_qseq;
    x[5] = c->mtid;
    x[6] = c->mpos;
    x[7] = c->isize;
    ok = (rabbit_bgzf_flush_try(fp, write_block, 4 + block_len) >= 0);
    if (fp->is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
        y = block_len;
        if (ok) ok = (rabbit_bgzf_write(fp, write_block,ed_swap_4p(&y), 4) >= 0);
        swap_data(c, b->l_data, b->data, 1);
    } else {
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, &block_len, 4) >= 0);
        }
    }
    if (ok) {
        ok = (rabbit_bgzf_write(fp, write_block, x, 32) >= 0);
    }
    if (ok) ok = (rabbit_bgzf_write(fp, write_block, b->data, c->l_qname - c->l_extranul) >= 0);
    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
        if (ok) ok = (rabbit_bgzf_write(fp, write_block, b->data + c->l_qname, b->l_data - c->l_qname) >= 0);
    } else { // with long CIGAR, insert a fake CIGAR record and move the real CIGAR to the CG:B,I tag
        uint8_t buf[8];
        uint32_t cigar_st, cigar_en, cigar[2];
        hts_pos_t cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
        if (cigreflen >= (1<<28)) {
            // Length of reference covered is greater than the biggest
            // CIGAR operation currently allowed.
            hts_log_error("Record %s with %d CIGAR ops and ref length %" PRIhts_pos
                                  " cannot be written in BAM.  Try writing SAM or CRAM instead.\n",
                          bam_get_qname(b), c->n_cigar, cigreflen);
            return -1;
        }
        cigar_st = (uint8_t*)bam_get_cigar(b) - b->data;
        cigar_en = cigar_st + c->n_cigar * 4;
        cigar[0] = (uint32_t)c->l_qseq << 4 | BAM_CSOFT_CLIP;
        cigar[1] = (uint32_t)cigreflen << 4 | BAM_CREF_SKIP;
        u32_to_le(cigar[0], buf);
        u32_to_le(cigar[1], buf + 4);
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, buf, 8) >= 0); // write cigar: <read_length>S<ref_length>N
        }
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, &b->data[cigar_en], b->l_data - cigar_en) >= 0); // write data after CIGAR
        }
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, "CGBI", 4) >= 0); // write CG:B,I
        }
        u32_to_le(c->n_cigar, buf);
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, buf, 4) >= 0); // write the true CIGAR length
        }
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, &b->data[cigar_st], c->n_cigar * 4) >= 0); // write the real CIGAR
        }
    }
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    return ok? 4 + block_len : -1;
}
int rabbit_bam_write_mul_test(BGZF *fp,BamWriteCompress *bam_write_compress,bam_write_block* &write_block,bam1_t *b){
    const bam1_core_t *c = &b->core;
    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
    int i, ok;
    if (c->l_qname - c->l_extranul > 255) {
        hts_log_error("QNAME \"%s\" is longer than 254 characters", bam_get_qname(b));
        errno = EOVERFLOW;
        return -1;
    }
    if (c->n_cigar > 0xffff) block_len += 16; // "16" for "CGBI", 4-byte tag length and 8-byte fake CIGAR


    if (c->pos > INT_MAX ||
        c->mpos > INT_MAX ||
        c->isize < INT_MIN || c->isize > INT_MAX) {
        hts_log_error("Positional data is too large for BAM format");
        return -1;
    }
    x[0] = c->tid;
    x[1] = c->pos;
    x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | (c->l_qname - c->l_extranul);
    if (c->n_cigar > 0xffff) x[3] = (uint32_t)c->flag << 16 | 2;
    else x[3] = (uint32_t)c->flag << 16 | (c->n_cigar & 0xffff);
    x[4] = c->l_qseq;
    x[5] = c->mtid;
    x[6] = c->mpos;
    x[7] = c->isize;
    ok = (rabbit_bgzf_mul_flush_try(fp,bam_write_compress, write_block, 4 + block_len) >= 0);
    if (fp->is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
        y = block_len;
        if (ok) ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block,ed_swap_4p(&y), 4) >= 0);
        swap_data(c, b->l_data, b->data, 1);
    } else {
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp,bam_write_compress, write_block, &block_len, 4) >= 0);
        }
    }
    if (ok) {
        ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, x, 32) >= 0);
    }
    if (ok) ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, b->data, c->l_qname - c->l_extranul) >= 0);
    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
        if (ok) ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, b->data + c->l_qname, b->l_data - c->l_qname) >= 0);
    } else { // with long CIGAR, insert a fake CIGAR record and move the real CIGAR to the CG:B,I tag
        uint8_t buf[8];
        uint32_t cigar_st, cigar_en, cigar[2];
        hts_pos_t cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
        if (cigreflen >= (1<<28)) {
            // Length of reference covered is greater than the biggest
            // CIGAR operation currently allowed.
            hts_log_error("Record %s with %d CIGAR ops and ref length %" PRIhts_pos
                                  " cannot be written in BAM.  Try writing SAM or CRAM instead.\n",
                          bam_get_qname(b), c->n_cigar, cigreflen);
            return -1;
        }
        cigar_st = (uint8_t*)bam_get_cigar(b) - b->data;
        cigar_en = cigar_st + c->n_cigar * 4;
        cigar[0] = (uint32_t)c->l_qseq << 4 | BAM_CSOFT_CLIP;
        cigar[1] = (uint32_t)cigreflen << 4 | BAM_CREF_SKIP;
        u32_to_le(cigar[0], buf);
        u32_to_le(cigar[1], buf + 4);
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp,bam_write_compress, write_block, buf, 8) >= 0); // write cigar: <read_length>S<ref_length>N
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp,bam_write_compress,write_block, &b->data[cigar_en], b->l_data - cigar_en) >= 0); // write data after CIGAR
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp,bam_write_compress,write_block, "CGBI", 4) >= 0); // write CG:B,I
        }
        u32_to_le(c->n_cigar, buf);
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, buf, 4) >= 0); // write the true CIGAR length
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress,write_block, &b->data[cigar_st], c->n_cigar * 4) >= 0); // write the real CIGAR
        }
    }
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    return ok? 4 + block_len : -1;
}




void benchmark_write_pack(BamCompleteBlock* completeBlock,samFile *output,sam_hdr_t *hdr,int level){

    uint8_t* compress_block_test = new uint8_t[BGZF_BLOCK_SIZE];
    uint8_t* uncompress_block_test = new uint8_t[BGZF_BLOCK_SIZE];
    output->fp.bgzf->block_offset=0;
    output->fp.bgzf->uncompressed_block=uncompress_block_test;
    output->fp.bgzf->compressed_block=compress_block_test;
    output->fp.bgzf->compress_level=level;
    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    if (sam_hdr_write(output, hdr) != 0){
        printf("HDR Write False!\n");
        return ;
    }
    bam_complete_block* un_comp;
    long long ans = 0;
    long long res = 0;
    bam_write_block *write_block=new bam_write_block();
    write_block->block_offset=0;
    write_block->uncompressed_data=new uint8_t[BGZF_BLOCK_SIZE];
    write_block->compressed_data=new uint8_t[BGZF_BLOCK_SIZE];
    write_block->status=0;
    while (1){
        un_comp = completeBlock->getCompleteBlock();
        if (un_comp == nullptr){
            break;
        }
//        printf("assign over block length is %d\n",un_comp->length);
        int ret;
        while ((ret=(read_bam(un_comp,b,0)))>=0) {
//            printf("One Bam1_t Size is %d\n",ret);
//            printf("This Bam1_t Char Number is %d\n",b->core.l_qseq);
/*
 *  尝试单线程输出
 */
//            sam_write1(output,hdr,b);
            rabbit_bam_write_test(output->fp.bgzf,write_block,b);
            ans++;
        }
        res++;
        completeBlock->backEmpty(un_comp);
    }
    rabbit_bgzf_flush(output->fp.bgzf,write_block);
    printf("Bam1_t Number is %lld\n",ans);
    printf("Block  Number is %lld\n",res);
}

void benchmark_write_mul_pack(BamCompleteBlock* completeBlock,BamWriteCompress* bam_write_compress,samFile *output,sam_hdr_t *hdr,int level){

    output->fp.bgzf->block_offset=0;
    output->fp.bgzf->compress_level=level;
    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    if (sam_hdr_write(output, hdr) != 0){
        printf("HDR Write False!\n");
        return ;
    }
    bam_complete_block* un_comp;
    long long ans = 0;
    long long res = 0;
    bam_write_block *write_block=bam_write_compress->getEmpty();
//    printf("Write Pack Num : %d\n",write_block->block_num);
//    printf("Write Pack Block Offset : %d\n",write_block->block_offset);
//    printf("Has get One Empty!\n");
    int bam_num=1;
    while (1){
        un_comp = completeBlock->getCompleteBlock();
        if (un_comp == nullptr){
            break;
        }
//        printf("assign over block length is %d\n",un_comp->length);
        int ret;
        while ((ret=(read_bam(un_comp,b,0)))>=0) {
//            printf("One Bam1_t Size is %d\n",ret);
//            printf("This Bam1_t Char Number is %d\n",b->core.l_qseq);
/*
 *  尝试单线程输出
 */
//            sam_write1(output,hdr,b);
//            printf("Try to output one Bam1_t : %d\n",bam_num++);
            rabbit_bam_write_mul_test(output->fp.bgzf,bam_write_compress,write_block,b);
            ans++;
        }
        res++;
        completeBlock->backEmpty(un_comp);
    }
    if (write_block->block_offset>0) {
//        rabbit_bgzf_mul_flush(output->fp.bgzf,bam_write_compress,write_block);
        bam_write_compress->inputUnCompressData(write_block);
    }
//    printf("The Input Block Num : %d\n",InputBlockNum);
    bam_write_compress->WriteComplete();
//    printf("Bam1_t Number is %lld\n",ans);
//    printf("Block  Number is %lld\n",res);
}
void BamWriter::write(bam1_t* b) {
    rabbit_bam_write_mul_test(output->fp.bgzf, bam_write_compress, write_block, b);

}

void BamWriter::bam_write(bam1_t* b){
    rabbit_bam_write_mul_test(output->fp.bgzf, bam_write_compress, write_block, b);
}
//BamWriter::BamWriter(std::string filename,sam_hdr_t *hdr,int level){
//
//    if ((output=sam_open(filename.c_str(),"w+"))==NULL){
//        printf("Can`t open this file!\n");
//        //TODO 处理一下无法打开的清空
//    }
//    if (sam_hdr_write(output, hdr) != 0) {
//        printf("HDR Write False!\n");
//        //TODO 处理一下无法输出的情况
//    }
//
//    n_thread_write = 1;
//    bam_write_compress = new BamWriteCompress(4000,n_thread_write);
//
//
//    write_compress_thread = new std::thread*[n_thread_write];
//    for (int i=0;i<n_thread_write;i++) write_compress_thread[i] = new std::thread(&bam_write_compress_pack,output->fp.bgzf,bam_write_compress);
//
//    write_output_thread = new std::thread(&bam_write_pack,output->fp.bgzf,bam_write_compress);
//
//
//
//    output->fp.bgzf->block_offset=0;
//    output->fp.bgzf->compress_level=level;
//
//    write_block=bam_write_compress->getEmpty();
//}
//BamWriter::BamWriter(std::string filename,sam_hdr_t *hdr, int threadNumber,int level){
//
//    if ((output=sam_open(filename.c_str(),"wb"))==NULL){
//        printf("Can`t open this file!\n");
//        //TODO 处理一下无法打开的清空
//    }
//    if (sam_hdr_write(output, hdr) != 0) {
//        printf("HDR Write False!\n");
//        //TODO 处理一下无法输出的情况
//    }
//
//    n_thread_write = threadNumber;
//    bam_write_compress = new BamWriteCompress(4000,n_thread_write);
//
//
//    write_compress_thread = new std::thread*[n_thread_write];
//    for (int i=0;i<n_thread_write;i++) write_compress_thread[i] = new std::thread(&bam_write_compress_pack,output->fp.bgzf,bam_write_compress);
//
//    write_output_thread = new std::thread(&bam_write_pack,output->fp.bgzf,bam_write_compress);
//
//
//
//    output->fp.bgzf->block_offset=0;
//    output->fp.bgzf->compress_level=level;
//
//    write_block=bam_write_compress->getEmpty();
//}

BamWriter::BamWriter(int threadNumber,  int level,int BufferSize){

    n_thread_write = threadNumber;
    bam_write_compress = new BamWriteCompress(BufferSize,n_thread_write);



}


BamWriter::BamWriter(std::string filename, int threadNumber,  int level,int BufferSize){

    if ((output=sam_open(filename.c_str(),"wb"))==NULL){
        printf("Can`t open this file!\n");
        //TODO 处理一下无法打开的清空
    }

    n_thread_write = threadNumber;
    bam_write_compress = new BamWriteCompress(BufferSize,n_thread_write);


    write_compress_thread = new std::thread*[n_thread_write];
    for (int i=0;i<n_thread_write;i++) write_compress_thread[i] = new std::thread(&bam_write_compress_pack,output->fp.bgzf,bam_write_compress);

    write_output_thread = new std::thread(&bam_write_pack,output->fp.bgzf,bam_write_compress);



    output->fp.bgzf->block_offset=0;
    output->fp.bgzf->compress_level=level;

    write_block=bam_write_compress->getEmpty();
}

BamWriter::BamWriter(std::string filename,sam_hdr_t *hdr, int threadNumber,  int level,int BufferSize){

    if ((output=sam_open(filename.c_str(),"wb"))==NULL){
        printf("Can`t open this file!\n");
        //TODO 处理一下无法打开的清空
    }
    if (sam_hdr_write(output, hdr) != 0) {
        printf("HDR Write False!\n");
        //TODO 处理一下无法输出的情况
    }

    n_thread_write = threadNumber;
    bam_write_compress = new BamWriteCompress(BufferSize,n_thread_write);


    write_compress_thread = new std::thread*[n_thread_write];
    for (int i=0;i<n_thread_write;i++) write_compress_thread[i] = new std::thread(&bam_write_compress_pack,output->fp.bgzf,bam_write_compress);

    write_output_thread = new std::thread(&bam_write_pack,output->fp.bgzf,bam_write_compress);



    output->fp.bgzf->block_offset=0;
    output->fp.bgzf->compress_level=level;

    write_block=bam_write_compress->getEmpty();
}

void BamWriter::hdr_write(sam_hdr_t* hdr){
    if (sam_hdr_write(output, hdr) != 0) {
        printf("HDR Write False!\n");
        //TODO 处理一下无法输出的情况
    }
}


void BamWriter::set_output(samFile *output){
    this->output=output;
    write_compress_thread = new std::thread*[n_thread_write];
    for (int i=0;i<n_thread_write;i++) write_compress_thread[i] = new std::thread(&bam_write_compress_pack,output->fp.bgzf,bam_write_compress);

    write_output_thread = new std::thread(&bam_write_pack,output->fp.bgzf,bam_write_compress);



    output->fp.bgzf->block_offset=0;
    output->fp.bgzf->compress_level=6;

    write_block=bam_write_compress->getEmpty();
}

void BamWriter::over(){

    if (write_block->block_offset>0) {
        bam_write_compress->inputUnCompressData(write_block);
        write_block=bam_write_compress->getEmpty();
    }
    bam_write_compress->WriteComplete();
//    sleep(100);
    for (int i=0;i<n_thread_write;i++) write_compress_thread[i]->join();
    write_output_thread->join();
//    sam_close(output);
}

// TODO
/*
 * 这部分需要放在析构函数里面
 *
 * if (write_block->block_offset>0) {
        InputBlockNum++;
        bam_write_compress->inputUnCompressData(write_block);
    }
    bam_write_compress->WriteComplete();
 */