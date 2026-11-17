// Exercise the production halo probe and CPE extraction with tiny BAM files.
#include "../swbam-app/algorithms/markdup/swbam_markdup_mpi.cpp"
#include <unistd.h>
#include <zlib.h>

static void Check(bool ok, const char *message) {
    if (ok) return;
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    fprintf(stderr, "[rank %d] markdup boundary regression FAILED: %s\n",
            rank, message);
    MPI_Abort(MPI_COMM_WORLD, 2);
    abort();
}

static std::vector<unsigned char> RunAndDecode(
        const swbam::BamInputBackend &input,
        const std::vector<swbam::BgzfBlockSpan> &spans,
        int rank, int ranks, const MdOnePassBoundaryPlan *plan) {
    MemWriter storage = {};
    swbam::MemoryRankBodySink prefix(&storage);
    swbam_usage::OrderedBody body(&prefix, 16u << 20);
    MpiMarkdupStats stats = {};
    uint64_t records = 0;
    int has = 0, first_tid = -1, first_pos = -1, last_tid = -1, last_pos = -1;
    Check(MdMarkdupOnePassToSink(input, spans.data(), spans.size(),
              (long long)rank * spans.size(), rank, ranks, 0, 300,
              256ull << 20, 1, &body, plan, &records, &has,
              &first_tid, &first_pos, &last_tid, &last_pos, &stats) == 0 &&
          body.Flush() == 0 && records == spans.size(), "one-pass fixture output");
    std::vector<unsigned char> compressed((size_t)body.size());
    Check(body.ReadAt(0, compressed.data(), compressed.size()) == 0,
          "ordered fixture readback");
    std::vector<unsigned char> result;
    unsigned char decoded[BGZF_MAX_BLOCK_SIZE];
    size_t offset = 0;
    while (offset < compressed.size()) {
        Check(compressed.size() - offset >= 18, "fixture BGZF header");
        const size_t bytes = (size_t)compressed[offset + 16] +
            ((size_t)compressed[offset + 17] << 8) + 1;
        Check(bytes <= compressed.size() - offset, "fixture BGZF length");
        z_stream stream = {};
        Check(inflateInit2(&stream, 31) == Z_OK, "fixture inflate init");
        stream.next_in = compressed.data() + offset;
        stream.avail_in = (uInt)bytes;
        stream.next_out = decoded;
        stream.avail_out = sizeof(decoded);
        Check(inflate(&stream, Z_FINISH) == Z_STREAM_END &&
              stream.total_in == bytes, "fixture BGZF inflate and CRC");
        result.insert(result.end(), decoded, decoded + stream.total_out);
        Check(inflateEnd(&stream) == Z_OK, "fixture inflate cleanup");
        offset += bytes;
    }
    free(storage.data);
    return result;
}

static void TestBoundary(int test, int rank, int ranks, const char *directory) {
    fprintf(stderr, "[rank %d] boundary fixture case=%d begin\n", rank, test);
    fflush(stderr);
    const std::string path = std::string(directory) + "/markdup-halo." +
        std::to_string((long long)getpid()) + "." + std::to_string(rank) + ".bam";
    const char *text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000000\n";
    sam_hdr_t *header = sam_hdr_parse(strlen(text), text);
    samFile *file = sam_open(path.c_str(), "wb1");
    Check(header && file, "fixture open");
    fprintf(stderr, "[rank %d] boundary fixture case=%d open\n", rank, test);
    fflush(stderr);
    Check(sam_hdr_write(file, header) == 0 && bgzf_flush(file->fp.bgzf) == 0,
          "fixture header");
    bam1_t *record = bam_init1();
    Check(record != nullptr, "fixture record allocation");
    const std::string seq(100, 'A');
    const uint32_t cigar = bam_cigar_gen(100, BAM_CMATCH);
    const int blocks = test == 2 || test == 3 ? 768 :
        test == 5 || test == 8 ? 128 : 64;
    for (int block = 0; block < blocks; ++block) {
        const bool unmapped =
            (test == 1 && rank == ranks - 1 && block >= 32) ||
            (test == 6 && rank == ranks - 1) || test == 7;
        int position = rank * 10000 + block * (blocks == 64 ? 100 : 10);
        if (test == 2 && rank == 1 && block < 384) position = rank * 10000 + block;
        if (test == 3 && rank == 0 && block >= 384) position = 9990;
        if (test == 4)
            position = rank == 0 ? 9400 + block * 9 : rank * 10000 + block * 10;
        if (test == 5) {
            position = rank == 0 ? 9000 + block * 7 : rank * 10000 + block * 10;
            if (rank == 1 && block >= 30 && block < 64) position = 10300;
            if (rank == 0 && block == 126) position = 9889;
        }
        if (test == 8) {
            position = rank * 10000 + (block == 0 ? 999 :
                block == 1 ? 4999 : block <= 64 ? 5299 : 6000 + block * 10);
        }
        char name[64];
        snprintf(name, sizeof(name), "r%d-b%d", rank, block);
        const int flag = unmapped ? BAM_FUNMAP :
            test == 5 ? BAM_FPAIRED | BAM_FREAD1 | BAM_FMREVERSE :
            test == 8 && block != 1 && block != 64 ? BAM_FSECONDARY : 0;
        const uint32_t limit_cigar[] = {
            bam_cigar_gen(300, BAM_CHARD_CLIP), bam_cigar_gen(1, BAM_CMATCH)};
        const bool limit_read = test == 8 && (block == 1 || block == 64);
        const uint32_t *record_cigar = limit_read
            ? limit_cigar + (block == 1 ? 1 : 0) : &cigar;
        const size_t record_cigars = limit_read && block == 64 ? 2 : 1;
        const size_t record_length = limit_read ? 1 : seq.size();
        Check(bam_set1(record, strlen(name), name, flag,
                  unmapped ? -1 : 0, unmapped ? -1 : position, 60,
                  unmapped ? 0 : record_cigars, unmapped ? nullptr : record_cigar,
                  test == 5 ? 0 : -1, test == 5 ? 50000 : -1, 0,
                  record_length, seq.c_str(), nullptr, 0) >= 0,
              "fixture record fields");
        if (limit_read)
            memset(bam_get_qual(record), block == 1 ? 20 : 40, record_length);
        if (test == 5) {
            const int32_t mate_score = 3000;
            Check(bam_aux_append(record, "MC", 'Z', 5, (const uint8_t *)"100M") == 0 &&
                  bam_aux_append(record, "ms", 'i', sizeof(mate_score),
                                 (const uint8_t *)&mate_score) == 0,
                  "paired fixture mate tags");
        }
        Check(sam_write1(file, header, record) >= 0 && bgzf_flush(file->fp.bgzf) == 0,
              "fixture block write");
    }
    bam_destroy1(record);
    sam_hdr_destroy(header);
    Check(sam_close(file) == 0, "fixture close");
    fprintf(stderr, "[rank %d] boundary fixture case=%d written\n", rank, test);
    fflush(stderr);

    swbam::MemoryBamInput input;
    Check(input.Open(path) == 0, "fixture input");
    std::vector<swbam::BgzfBlockSpan> spans;
    Check(input.ScanBlocks(&spans) == 0 && spans.size() == (size_t)blocks, "fixture spans");
    fprintf(stderr, "[rank %d] boundary fixture case=%d loaded\n", rank, test);
    fflush(stderr);
    MdOnePassBoundaryPlan plan;
    MpiMarkdupStats stats = {};
    Check(MdBuildOnePassBoundaryPlan(input, spans.data(), spans.size(),
              rank * blocks, rank, ranks, 0, 300, 256ull * 1024 * 1024,
              &plan, &stats) == 0, "boundary probe");
    Check(plan.enabled == (test == 6 || test == 7 ? 0 : 1),
          "unplaced-only ranges fall back without rejecting the input");
    if (plan.enabled)
        Check(plan.partial == (test == 2 || test == 3 || test == 5),
              "fixed probe distinguishes complete coverage from partial coverage");
    Check(plan.first_probe_end_block == (uint64_t)(rank * blocks + 32) &&
          plan.last_probe_begin_block == (uint64_t)(rank * blocks + blocks - 32),
          "probe stays fixed at first/last 32 blocks");
    if (plan.partial) {
        Check(plan.conservative && plan.closed_owner_keys.empty(),
              "partial coverage retains boundary guards and certifies no closed key");
        for (const auto &decision : plan.decisions)
            Check(decision.second == kMdOnePassDuplicate, "partial prepass is DUP-only");
        for (const auto &id : plan.candidates)
            Check(id.kind != 2, "partial prepass retains all paired markers");
    }
    if (test == 5 && rank == 0) {
        Check(plan.decisions.count((uint64_t)126 << 32) == 0 &&
              plan.decisions.at((uint64_t)127 << 32) == kMdOnePassDuplicate,
              "partial paired group retains its winner and finalizes the QNAME loser");
    }
    if (test == 4 && rank == 1) {
        const uint64_t fringe = (uint64_t)(rank * blocks + 25) << 32;
        Check(plan.decisions.count(fringe) == 0, "fringe winner must not receive final KEEP");
        const MdBoundaryCandidateId id = {fringe, 1};
        Check(plan.candidates.count(id) == 0, "fringe winner must remain in the main window");
    }
    if (test == 3 || test == 5) {
        MdOnePassBoundaryPlan fallback;
        const auto actual = RunAndDecode(input, spans, rank, ranks, &plan);
        const auto reference = RunAndDecode(input, spans, rank, ranks, &fallback);
        Check(actual == reference, "partial prepass preserves every output record byte and order");
    }
    if (test >= 6) {
        const auto actual = RunAndDecode(input, spans, rank, ranks, &plan);
        size_t offset = 0, records = 0, duplicates = 0;
        while (offset < actual.size()) {
            Check(actual.size() - offset >= 36, "fixture record header");
            const size_t bytes = (size_t)le_to_u32(actual.data() + offset) + 4;
            Check(bytes >= 36 && bytes <= actual.size() - offset,
                  "fixture record length");
            const uint16_t flag = (uint16_t)(
                le_to_u32(actual.data() + offset + 16) >> 16);
            const bool duplicate = (flag & BAM_FDUP) != 0;
            Check(duplicate == (test == 8 && records == 1),
                  "exact duplicate member across the key+L batch boundary");
            if ((test == 6 && rank == ranks - 1) || test == 7)
                Check((flag & BAM_FUNMAP) != 0, "unmapped tail flags are preserved");
            duplicates += duplicate;
            ++records;
            offset += bytes;
        }
        Check(records == (size_t)blocks && duplicates == (test == 8 ? 1u : 0u),
              "all mapped/unmapped fixture records are preserved");
    }
    input.Close();
    Check(unlink(path.c_str()) == 0, "fixture cleanup");
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        printf("PASS markdup boundary ranks=%d case=%d enabled=%d\n",
               ranks, test, plan.enabled);
        fflush(stdout);
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    Check(argc == 2 && ranks >= 2, "usage: >=2 ranks, existing temp directory");
    fprintf(stderr, "[rank %d] boundary regression athread init\n", rank);
    fflush(stderr);
    athread_init();
    fprintf(stderr, "[rank %d] boundary regression athread ready\n", rank);
    fflush(stderr);
    for (int test = 0; test < 9; ++test) TestBoundary(test, rank, ranks, argv[1]);
    athread_halt();
    MPI_Finalize();
    return 0;
}
