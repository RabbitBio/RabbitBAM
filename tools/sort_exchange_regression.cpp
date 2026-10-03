// Exercise the production internal byte exchange without running BAM kernels.
#include "../swbam-app/algorithms/sort/swbam_sort_optimized_mpi.cpp"

static void Check(bool ok, const char *message) {
    if (ok) return;
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    fprintf(stderr, "[rank %d] sort exchange regression FAILED: %s\n",
            rank, message);
    MPI_Abort(MPI_COMM_WORLD, 2);
    abort();
}

static size_t Bytes(int test, int src, int dst, int ranks) {
    const size_t chunk = kSortExchangeChunk;
    if (test == 0) return 0;
    if (test == 1) {
        if (src == 0 && dst == 1) return 1;
        if (src == 1 && dst == 2) return 7;
        return 0;
    }
    if (test == 2) {
        if (dst != (src + 1) % ranks) return 0;
        if (src == 0) return 2 * chunk + 17;
        if (src == 1) return 1;
        if (src == 2) return chunk + 3;
        return 13 + src;
    }
    if (src == 0 && dst == ranks - 1) return chunk + 1;
    return (src * 7 + dst * 3) % 5;
}

static unsigned char Pattern(int src, int dst, size_t index, int pass) {
    return (unsigned char)(src * 17 + dst * 11 + index * 13 + pass);
}

static void TestExchange(int test, int rank, int ranks) {
    for (int step = 1; step < ranks; ++step) {
        const int send_to = (rank + step) % ranks;
        const int recv_from = (rank - step + ranks) % ranks;
        const size_t send_size = Bytes(test, rank, send_to, ranks);
        const size_t recv_size = Bytes(test, recv_from, rank, ranks);
        std::vector<unsigned char> send(send_size);
        std::vector<unsigned char> recv(recv_size);
        // Reuse the production tags across metadata/raw calls and rounds.
        for (int pass = 0; pass < 2; ++pass) {
            for (size_t i = 0; i < send_size; ++i)
                send[i] = Pattern(rank, send_to, i, pass);
            std::fill(recv.begin(), recv.end(), 0xff);
            Check(MpiSortSendRecvBytes(send_to,
                      send_size ? (const char *)send.data() : nullptr, send_size,
                      recv_from, recv_size ? (char *)recv.data() : nullptr, recv_size,
                      pass ? kSortTagRaw : kSortTagMeta) == 0, "byte exchange");
            for (size_t i = 0; i < recv_size; ++i)
                Check(recv[i] == Pattern(recv_from, rank, i, pass), "payload bytes");
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        printf("PASS sort byte exchange ranks=%d case=%d\n", ranks, test);
        fflush(stdout);
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    Check(ranks >= 3, "requires at least three MPI ranks");
    for (int test = 0; test < 4; ++test) TestExchange(test, rank, ranks);
    MPI_Finalize();
    return 0;
}
