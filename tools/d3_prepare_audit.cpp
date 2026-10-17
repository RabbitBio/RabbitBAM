#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/kseq.h>
#include <htslib/sam.h>
#include <zlib.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>

KSEQ_INIT(gzFile, gzread)

namespace {
void Require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

std::string PairName(const kseq_t *read) {
    std::string name(read->name.s, read->name.l);
    if (name.size() >= 2 && name[name.size() - 2] == '/' &&
        (name.back() == '1' || name.back() == '2')) name.resize(name.size() - 2);
    return name;
}

void CheckFastq(const char *r1, const char *r2, uint64_t expected) {
    gzFile files[2] = {gzopen(r1, "rb"), gzopen(r2, "rb")};
    Require(files[0] && files[1], "cannot open FASTQ pair");
    gzbuffer(files[0], 1u << 20);
    gzbuffer(files[1], 1u << 20);
    kseq_t *reads[2] = {kseq_init(files[0]), kseq_init(files[1])};
    uint64_t pairs = 0, bases = 0;
    size_t minimum = SIZE_MAX, maximum = 0;
    for (;;) {
        const int a = kseq_read(reads[0]), b = kseq_read(reads[1]);
        if (a == -1 && b == -1) break;
        Require(a > 0 && b > 0, "truncated, unequal, or empty FASTQ pair");
        Require(PairName(reads[0]) == PairName(reads[1]), "FASTQ mate QNAME mismatch");
        for (int i = 0; i < 2; ++i) {
            Require(reads[i]->seq.l == reads[i]->qual.l, "FASTQ quality length mismatch");
            bases += reads[i]->seq.l;
            minimum = std::min(minimum, reads[i]->seq.l);
            maximum = std::max(maximum, reads[i]->seq.l);
        }
        ++pairs;
        if (pairs % 10000000 == 0)
            fprintf(stderr, "fastq_audit pairs=%llu\n", (unsigned long long)pairs);
    }
    for (int i = 0; i < 2; ++i) {
        int error = Z_OK;
        gzerror(files[i], &error);
        Require(error == Z_OK || error == Z_STREAM_END, "FASTQ gzip/CRC error");
        kseq_destroy(reads[i]);
        Require(gzclose(files[i]) == Z_OK, "FASTQ close error");
    }
    Require(pairs == expected, "FASTQ pair count differs from accession metadata");
    printf("{\"pairs\":%llu,\"reads\":%llu,\"bases\":%llu,"
           "\"min_length\":%zu,\"max_length\":%zu,\"paired_qnames\":true,\"gzip_crc\":true}\n",
           (unsigned long long)pairs, (unsigned long long)(pairs * 2),
           (unsigned long long)bases, pairs ? minimum : 0, maximum);
}

uint32_t Le32(const unsigned char *p) {
    return uint32_t(p[0]) | uint32_t(p[1]) << 8 |
           uint32_t(p[2]) << 16 | uint32_t(p[3]) << 24;
}

void CheckBam(const char *path, uint64_t expected_primary) {
    BGZF *input = bgzf_open(path, "r");
    Require(input != nullptr, "cannot open BAM");
    Require(bgzf_check_EOF(input) == 1, "missing BGZF EOF");
    bam_hdr_t *header = bam_hdr_read(input);
    bam1_t *record = bam_init1();
    hts_md5_context *digest = hts_md5_init();
    Require(header && record && digest, "BAM header or allocation failed");
    uint64_t records = 0, primary = 0, paired = 0, proper = 0, mapped = 0;
    uint64_t secondary = 0, supplementary = 0, duplicates = 0, raw_bytes = 0;
    uint32_t min_length = UINT32_MAX, max_length = 0, max_cigar = 0, max_record = 0;
    int32_t previous_tid = -1;
    hts_pos_t previous_pos = -1;
    bool unplaced = false;
    for (;;) {
        const int peek = bgzf_peek(input);
        if (peek == -1) break;
        Require(peek >= 0, "BGZF inflate/CRC failure");
        const size_t available = size_t(input->block_length - input->block_offset);
        Require(available >= 4, "BAM record length crosses BGZF block");
        const auto *raw = static_cast<const unsigned char *>(input->uncompressed_block)
                          + input->block_offset;
        const uint32_t length = Le32(raw);
        Require(length >= 32 && uint64_t(length) + 4 <= available,
                "BAM record crosses BGZF block; prepare common input before benchmarking");
        hts_md5_update(digest, raw, size_t(length) + 4);
        raw_bytes += uint64_t(length) + 4;
        max_record = std::max(max_record, length + 4);
        Require(bam_read1(input, record) >= 0, "invalid BAM record");
        const auto &core = record->core;
        Require(core.tid >= -1 && core.tid < header->n_targets, "invalid BAM reference ID");
        Require(core.mtid >= -1 && core.mtid < header->n_targets, "invalid mate reference ID");
        if (core.tid < 0) {
            unplaced = true;
        } else {
            Require(!unplaced && (core.tid > previous_tid ||
                    (core.tid == previous_tid && core.pos >= previous_pos)),
                    "BAM is not coordinate sorted");
            previous_tid = core.tid;
            previous_pos = core.pos;
        }
        ++records;
        primary += !(core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY));
        paired += !!(core.flag & BAM_FPAIRED);
        proper += !!(core.flag & BAM_FPROPER_PAIR);
        mapped += !(core.flag & BAM_FUNMAP);
        secondary += !!(core.flag & BAM_FSECONDARY);
        supplementary += !!(core.flag & BAM_FSUPPLEMENTARY);
        duplicates += !!(core.flag & BAM_FDUP);
        min_length = std::min(min_length, uint32_t(core.l_qseq));
        max_length = std::max(max_length, uint32_t(core.l_qseq));
        max_cigar = std::max(max_cigar, core.n_cigar);
        if (records % 10000000 == 0)
            fprintf(stderr, "bam_audit records=%llu\n", (unsigned long long)records);
    }
    Require(!input->errcode, "BGZF read error");
    Require(primary == expected_primary, "BAM primary count differs from input reads");
    Require(duplicates == 0, "initial BAM unexpectedly contains duplicate flags");
    unsigned char md5[16];
    char hex[33];
    hts_md5_final(md5, digest);
    hts_md5_hex(hex, md5);
    printf("{\"records\":%llu,\"primary\":%llu,\"paired\":%llu,\"proper_pair\":%llu,"
           "\"mapped\":%llu,\"secondary\":%llu,\"supplementary\":%llu,\"duplicates\":%llu,"
           "\"min_read_length\":%u,\"max_read_length\":%u,\"max_cigar_ops\":%u,"
           "\"max_record_bytes\":%u,\"references\":%d,\"raw_record_bytes\":%llu,"
           "\"record_md5\":\"%s\",\"coordinate_sorted\":true,"
           "\"record_aligned_bgzf\":true,\"crc\":true,\"eof\":true}\n",
           (unsigned long long)records, (unsigned long long)primary,
           (unsigned long long)paired, (unsigned long long)proper,
           (unsigned long long)mapped, (unsigned long long)secondary,
           (unsigned long long)supplementary, (unsigned long long)duplicates,
           records ? min_length : 0, max_length, max_cigar, max_record,
           header->n_targets, (unsigned long long)raw_bytes, hex);
    hts_md5_destroy(digest);
    bam_destroy1(record);
    bam_hdr_destroy(header);
    Require(bgzf_close(input) == 0, "BAM close error");
}
} // namespace

int main(int argc, char **argv) {
    try {
        if (argc == 5 && std::string(argv[1]) == "fastq")
            CheckFastq(argv[2], argv[3], std::stoull(argv[4]));
        else if (argc == 4 && std::string(argv[1]) == "bam")
            CheckBam(argv[2], std::stoull(argv[3]));
        else throw std::runtime_error("usage: d3_prepare_audit fastq R1 R2 pairs | bam input.bam primary_reads");
        return 0;
    } catch (const std::exception &error) {
        fprintf(stderr, "ERROR: %s\n", error.what());
        return 1;
    }
}
