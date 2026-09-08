#include "BamReader.h"

// RabbitBAM's sources rely on `using namespace std;` (previously pulled in
// transitively via Buffer.h, which the vendored reader no longer includes).
using namespace std;

// sam_hdr_read() may leave records in its already-decompressed BGZF block.
// Preserve those bytes before the raw reader starts at the next physical block.
// In particular, the parallel sort writer packs headers and records together.
static bam_block *buffered_record_block(BGZF *fp, BamCompress *compress) {
    if (fp->block_offset >= fp->block_length) return nullptr;
    bam_block *block = compress->getEmpty();
    block->errcode = 0;
    block->pos = 0;
    block->block_address = fp->block_address;
    block->length = fp->block_length - fp->block_offset;
    memcpy(block->data,
           static_cast<unsigned char *>(fp->uncompressed_block) + fp->block_offset,
           block->length);
    const auto split = find_divide_pos_and_get_read_number(block);
    block->split_pos = split.first;
    block->bam_number = split.second;
    return block;
}


void read_pack(BGZF *fp, BamRead *read) {
    bam_block *b;
    b = read->getEmpty();
    int count = 0;
    while (read_block(fp, b) == 0) {
        read->inputBlock(b);
        b = read->getEmpty();
    }
    read->ReadComplete();
}

void compress_pack(BamRead *read, BamCompress *compress) {
    pair < bam_block * , int > comp;
    bam_block *un_comp = compress->getEmpty();
    while (1) {
        comp = read->getReadBlock();
        if (comp.second < 0) {
            break;
        }
        block_decode_func(comp.first, un_comp);
        read->backBlock(comp.first);

        std::pair<int, int> tmp_pair = find_divide_pos_and_get_read_number(un_comp);
        un_comp->split_pos = tmp_pair.first, un_comp->bam_number = tmp_pair.second;
        compress->inputUnCompressData(un_comp, comp.second);
        un_comp = compress->getEmpty();
    }
    compress->CompressThreadComplete();
}

void assign_pack(BamCompress *compress, BamCompleteBlock *completeBlock,
                 bam_block *initial_block) {
    bam_block *un_comp = nullptr;
    bam_complete_block *assign_block = completeBlock->getEmpty();
    int need_block_len = 0, input_length = 0;
    int last_use_block_length = 0;
    bool isclean = true;
    int ret = -1;
    while (1) {
        if (isclean && un_comp != nullptr) {
            compress->backEmpty(un_comp);
        }
        if (initial_block != nullptr) {
            un_comp = initial_block;
            initial_block = nullptr;
        } else {
            un_comp = compress->getUnCompressData();
        }
        if (un_comp == nullptr) {
            break;
        }
        ret = un_comp->split_pos;
        need_block_len = ret;
        if (assign_block->length + need_block_len > assign_block->data_size) {
            completeBlock->inputCompleteBlock(assign_block);
            assign_block = completeBlock->getEmpty();
        }
        if (ret != un_comp->length) {

            memcpy(assign_block->data + assign_block->length, un_comp->data, ret * sizeof(char));
            assign_block->length += ret;
            completeBlock->inputCompleteBlock(assign_block);
            assign_block = completeBlock->getEmpty();

            memcpy(assign_block->data + assign_block->length, un_comp->data + ret,
                   (un_comp->length - ret) * sizeof(char));
            assign_block->length += (un_comp->length - ret);
            compress->backEmpty(un_comp);
            un_comp = compress->getUnCompressData();
            memcpy(assign_block->data + assign_block->length, un_comp->data, un_comp->length * sizeof(char));
            assign_block->length += un_comp->length;


            int divide_pos = 0;
            int ret = 0;
            uint32_t x[8], new_l_data;
            while (divide_pos < assign_block->length) {
                Rabbit_memcpy(&ret, assign_block->data + divide_pos, 4);
                if (ret >= 32) {
                    if (divide_pos + 4 + 32 > assign_block->length) {
                        compress->backEmpty(un_comp);
                        un_comp = compress->getUnCompressData();
                        if (assign_block->length + un_comp->length > assign_block->data_size) {
                            change_data_size(assign_block);
                        }
                        memcpy(assign_block->data + assign_block->length, un_comp->data,
                               un_comp->length * sizeof(char));
                        assign_block->length += un_comp->length;
                    }
                    Rabbit_memcpy(x, assign_block->data + divide_pos + 4, 32);
                    int pos = (int32_t) x[1];
                    int l_qname = x[2] & 0xff;
                    int l_extranul = (l_qname % 4 != 0) ? (4 - l_qname % 4) : 0;
                    int n_cigar = x[3] & 0xffff;
                    int l_qseq = x[4];
                    new_l_data = ret - 32 + l_extranul;//block_len + c->l_extranul
                    if (new_l_data > INT_MAX || l_qseq < 0 || l_qname < 1) {
                        divide_pos += 4 + 32;
                        continue;
                    }
                    if (((uint64_t) n_cigar << 2) + l_qname + l_extranul
                        + (((uint64_t) l_qseq + 1) >> 1) + l_qseq > (uint64_t) new_l_data) {
                        divide_pos += 4 + 32;
                        continue;
                    }
                    while (divide_pos + 4 + 32 + l_qname > assign_block->length) {
                        compress->backEmpty(un_comp);
                        un_comp = compress->getUnCompressData();
                        if (assign_block->length + un_comp->length > assign_block->data_size) {
                            change_data_size(assign_block);
                        }
                        memcpy(assign_block->data + assign_block->length, un_comp->data,
                               un_comp->length * sizeof(char));
                        assign_block->length += un_comp->length;
                    }
                    char fg_char;
                    Rabbit_memcpy(&fg_char, assign_block->data + divide_pos + 4 + 32 + l_qname - 1, 1);
                    if (fg_char != '\0') {
                    }
                    if (fg_char != '\0' && l_extranul <= 0 && new_l_data > INT_MAX - 4) {

                        while (divide_pos + 4 + 32 + l_qname > assign_block->length) {
                            compress->backEmpty(un_comp);
                            un_comp = compress->getUnCompressData();
                            if (assign_block->length + un_comp->length > assign_block->data_size) {
                                change_data_size(assign_block);
                            }
                            memcpy(assign_block->data + assign_block->length, un_comp->data,
                                   un_comp->length * sizeof(char));
                            assign_block->length += un_comp->length;
                        }
                        divide_pos += 4 + 32 + l_qname;
                        continue;
                    }

                    while (divide_pos + 4 + ret > assign_block->length) {
                        compress->backEmpty(un_comp);
                        un_comp = compress->getUnCompressData();
                        if (assign_block->length + un_comp->length > assign_block->data_size) {
                            change_data_size(assign_block);
                        }
                        memcpy(assign_block->data + assign_block->length, un_comp->data,
                               un_comp->length * sizeof(char));
                        assign_block->length += un_comp->length;
                    }
                    divide_pos += 4 + ret;
                } else {
                    if (divide_pos + 4 > assign_block->length) {
                        compress->backEmpty(un_comp);
                        un_comp = compress->getUnCompressData();
                        if (assign_block->length + un_comp->length > assign_block->data_size) {
                            change_data_size(assign_block);
                        }
                        memcpy(assign_block->data + assign_block->length, un_comp->data,
                               un_comp->length * sizeof(char));
                        assign_block->length += un_comp->length;
                    }
                    divide_pos += 4;
                }

            }

            completeBlock->inputCompleteBlock(assign_block);
            assign_block = completeBlock->getEmpty();
        } else {
            if (ret != un_comp->length) {
                break;
            }
            memcpy(assign_block->data + assign_block->length, un_comp->data, ret * sizeof(char));
            assign_block->length += ret;
            last_use_block_length = 0;
            isclean = true;
        }
    }
    if (assign_block->length != 0) {
        completeBlock->inputCompleteBlock(assign_block);
    }
    completeBlock->is_over();

}

void compress_test_pack(BamCompress *compress) {
    bam_block *un_comp = nullptr;
    while (1) {
        un_comp = compress->getUnCompressData();
        if (un_comp == nullptr) {
            break;
        }
        compress->backEmpty(un_comp);

    }
}


BamReader::BamReader(std::string file_name, int n_thread, bool is_tgs) {


    if ((sin = sam_open(file_name.c_str(), "r")) == NULL) {
        printf("Can`t open this file!\n");

    }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
    }


    read = new BamRead(1024);
    this->n_thread = n_thread;
    compress = new BamCompress(1024, n_thread);
    completeBlock = new BamCompleteBlock(256);

    bam_block *initial_block = buffered_record_block(sin->fp.bgzf, compress);

    read_thread = new thread(&read_pack, sin->fp.bgzf, read);
    compress_thread = new thread *[n_thread];

    for (int i = 0; i < n_thread; i++) {

        compress_thread[i] = new thread(&compress_pack, read, compress);

    }
    // NOTE: assign to the MEMBER assign_thread (was shadowed by a local
    // declaration here, leaving the member uninitialised so ~BamReader()'s
    // assign_thread->join() dereferenced garbage and crashed).
    assign_thread = new thread(&assign_pack, compress, completeBlock, initial_block);


    if(is_tgs) {
        un_comp = completeBlock->getCompleteBlock();
    }
//
//#ifdef use_parallel_read
//#else
//    un_comp = completeBlock->getCompleteBlock();
//#endif


}


BamReader::BamReader(std::string file_name, int read_block, int compress_block, int compress_complete_block,
                     int n_thread, bool is_tgs) {

    if ((sin = sam_open(file_name.c_str(), "r")) == NULL) {
        printf("Can`t open this file!\n");

    }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
    }

    read = new BamRead(read_block);
    this->n_thread = n_thread;
    compress = new BamCompress(compress_block, n_thread);
    completeBlock = new BamCompleteBlock(compress_block);

    bam_block *initial_block = buffered_record_block(sin->fp.bgzf, compress);

    read_thread = new thread(&read_pack, sin->fp.bgzf, read);
    compress_thread = new thread *[n_thread];

    for (int i = 0; i < n_thread; i++) {
        compress_thread[i] = new thread(&compress_pack, read, compress);
    }
    // Assign to the MEMBER assign_thread (see note in the other constructor).
    assign_thread = new thread(&assign_pack, compress, completeBlock, initial_block);

    if(is_tgs) {
        un_comp = completeBlock->getCompleteBlock();
    }

//#ifdef use_parallel_read
//#else
//    un_comp = completeBlock->getCompleteBlock();
//#endif
}


sam_hdr_t *BamReader::getHeader() {
    return hdr;
}

bool BamReader::getBam1_t(bam1_t *b) {
    int ret;
    while (un_comp != nullptr) {
        if ((ret = (read_bam(un_comp, b, 0))) >= 0) {
            return true;
        } else {
            completeBlock->backEmpty(un_comp);
            un_comp = completeBlock->getCompleteBlock();
        }
    }
    return false;
}
