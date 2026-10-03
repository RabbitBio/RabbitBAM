// Exercise the production halo probe and CPE extraction with tiny BAM files.
#include "../swbam-app/algorithms/markdup/swbam_markdup_mpi.cpp"
#include <unistd.h>

static void Check(bool ok, const char *message) {
    if (ok) return;
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    fprintf(stderr, "[rank %d] markdup boundary regression FAILED: %s\n",
            rank, message);
    MPI_Abort(MPI_COMM_WORLD, 2);
    abort();
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
    for (int block = 0; block < 64; ++block) {
        const bool unmapped = test == 1 && rank == ranks - 1 && block >= 32;
        int position = rank * 10000 + block * 100;
        if (test == 2 && rank == 1 && block < 32) position = rank * 10000 + block;
        if (test == 3 && rank == 0 && block >= 32) position = 9990;
        char name[64];
        snprintf(name, sizeof(name), "r%d-b%d", rank, block);
        Check(bam_set1(record, strlen(name), name, unmapped ? BAM_FUNMAP : 0,
                  unmapped ? -1 : 0, unmapped ? -1 : position, 60,
                  unmapped ? 0 : 1, unmapped ? nullptr : &cigar,
                  -1, -1, 0, seq.size(), seq.c_str(), nullptr, 0) >= 0,
              "fixture record fields");
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
    Check(input.ScanBlocks(&spans) == 0 && spans.size() == 64, "fixture spans");
    fprintf(stderr, "[rank %d] boundary fixture case=%d loaded\n", rank, test);
    fflush(stderr);
    MdOnePassBoundaryPlan plan;
    MpiMarkdupStats stats = {};
    Check(MdBuildOnePassBoundaryPlan(input, spans.data(), spans.size(),
              rank * 64, rank, ranks, 0, 300, 256ull * 1024 * 1024,
              &plan, &stats) == 0, "boundary probe");
    Check(plan.enabled == (test < 2), "enabled for valid halos only");
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
    for (int test = 0; test < 4; ++test) TestBoundary(test, rank, ranks, argv[1]);
    athread_halt();
    MPI_Finalize();
    return 0;
}
