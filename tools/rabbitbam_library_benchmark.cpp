#include "BamReader.h"
#include "BamWriter.h"

#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

static uint64_t Number(const char *text) {
    char *end = nullptr;
    errno = 0;
    const unsigned long long value = strtoull(text, &end, 10);
    if (errno || !text[0] || text[0] == '-' || *end) {
        throw std::runtime_error("invalid integer argument");
    }
    return value;
}

int main(int argc, char **argv) {
    try {
        if (argc != 7) {
            std::cerr << "usage: rabbitbam-library count|copy INPUT OUTPUT|- NR NW EXPECTED\n";
            return 2;
        }
        const std::string mode(argv[1]), input(argv[2]), output(argv[3]);
        const int nr = (int)Number(argv[4]), nw = (int)Number(argv[5]);
        const uint64_t expected = Number(argv[6]);
        struct stat info;
        if ((mode != "count" && mode != "copy") || nr < 1 || nr > 120 ||
            nw < 0 || nr + nw > 120 || !expected || stat(input.c_str(), &info) ||
            !S_ISREG(info.st_mode) || (mode == "copy" && (!nw || output == input || output == "-")) ||
            (mode == "count" && (nw || output != "-"))) {
            std::cerr << "invalid mode, input, thread allocation, or expected count\n";
            return 2;
        }
        const auto begin = std::chrono::steady_clock::now();
        omp_set_dynamic(0);
        // Match upstream parallel_test: same pools, parser and serializer APIs.
        BamReader *reader = new BamReader(input, nr, false);
        BamWriter *writer = mode == "copy"
            ? new BamWriter(output, reader->getHeader(), nw, 1, 4 << 10, false)
            : nullptr;
        const int capacity = 40 << 10;
        std::vector<bam1_t *> records[THREAD_NUM_P];
        int allocation_failed = 0;
#pragma omp parallel for num_threads(THREAD_NUM_P) reduction(|:allocation_failed) schedule(static)
        for (int lane = 0; lane < THREAD_NUM_P; ++lane) {
            records[lane].resize(capacity);
            for (int i = 0; i < capacity; ++i) {
                records[lane][i] = bam_init1();
                if (!records[lane][i]) allocation_failed = 1;
            }
        }
        if (allocation_failed) return 1;
        uint64_t count = 0, batches = 0;
        while (true) {
            std::vector<bam1_t *> batch = reader->getBam1_t_parallel(records);
            if (batch.empty()) break;
            count += batch.size();
            ++batches;
            if (writer) writer->write_parallel(batch);
        }
        if (writer) writer->over_parallel();
        for (int lane = 0; lane < THREAD_NUM_P; ++lane) {
            for (bam1_t *record : records[lane]) bam_destroy1(record);
        }
        // Upstream's example leaves reader/writer teardown to process exit.
        // Its reader destructor is not a supported safe close path in this revision.
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - begin).count();
        std::cout << "BENCH records=" << count << " batches=" << batches
                  << " nr=" << nr << " nw=" << nw
                  << " parser_threads=4 serializer_threads=" << (writer ? 4 : 0)
                  << " level=1 internal_s=" << seconds << std::endl;
        if (count != expected) {
            std::cerr << "record count mismatch: " << count << " != " << expected << '\n';
            return 1;
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
