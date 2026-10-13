#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <mpi.h>

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <string>
#include <unistd.h>
#include <vector>

namespace {
const size_t kChunk = 8u * 1024u * 1024u;
volatile uint64_t guard = 0;

bool AllOk(bool ok) {
    int local = ok ? 1 : 0, global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return global != 0;
}

void Require(bool ok, const char *stage, int rank) {
    if (!ok) {
        fprintf(stderr, "[rank %d] IO_ERROR stage=%s errno=%d message=%s\n",
                rank, stage, errno, strerror(errno));
        fflush(stderr);
    }
    if (!AllOk(ok)) MPI_Abort(MPI_COMM_WORLD, 1);
}

size_t Number(const char *text, size_t maximum) {
    char *end = nullptr;
    errno = 0;
    const unsigned long long value = strtoull(text, &end, 10);
    return errno || end == text || *end || !value || value > maximum
        ? 0 : static_cast<size_t>(value);
}

char *Align(void *raw) {
    return reinterpret_cast<char *>((reinterpret_cast<uintptr_t>(raw) + 4095u) &
                                     ~uintptr_t(4095u));
}

void Fill(char *data, size_t bytes, int rank) {
    uint64_t x = 0x9e3779b97f4a7c15ull ^ static_cast<uint64_t>(rank + 1);
    uint64_t *words = reinterpret_cast<uint64_t *>(data);
    for (size_t i = 0; i < bytes / sizeof(uint64_t); ++i) {
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        words[i] = x * 2685821657736338717ull;
    }
}

// Keep repeated copies observable without including validation in the timer.
__attribute__((noinline)) void Copy(char *dst, const char *src, char *scratch,
                                   size_t bytes, size_t chunk, bool assemble) {
    for (size_t pos = 0; pos < bytes; pos += chunk) {
        const size_t count = std::min(chunk, bytes - pos);
        if (assemble) {
            memcpy(scratch, src + pos, count);
            asm volatile("" : : "r"(scratch) : "memory");
            memcpy(dst + pos, scratch, count);
        } else {
            memcpy(dst + pos, src + pos, count);
        }
        asm volatile("" : : "r"(dst + pos) : "memory");
    }
}

void Report(const char *name, int repeat, uint64_t bytes, double seconds,
            int rank, int ranks) {
    std::vector<double> times(static_cast<size_t>(ranks));
    MPI_Gather(&seconds, 1, MPI_DOUBLE, times.data(), 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    if (rank != 0) return;
    const double slowest = *std::max_element(times.begin(), times.end());
    const double fastest = *std::min_element(times.begin(), times.end());
    for (int r = 0; r < ranks; ++r) {
        printf("IO_RANK name=%s repeat=%d rank=%d bytes=%llu seconds=%.6f GBps=%.6f\n",
               name, repeat, r, static_cast<unsigned long long>(bytes), times[r],
               static_cast<double>(bytes) / times[r] / 1e9);
    }
    printf("IO_SUMMARY name=%s repeat=%d bytes_per_rank=%llu max_seconds=%.6f "
           "aggregate_GBps=%.6f rank_min_GBps=%.6f rank_max_GBps=%.6f\n",
           name, repeat, static_cast<unsigned long long>(bytes), slowest,
           static_cast<double>(bytes) * ranks / slowest / 1e9,
           static_cast<double>(bytes) / slowest / 1e9,
           static_cast<double>(bytes) / fastest / 1e9);
    fflush(stdout);
}

bool Transfer(int fd, char *data, size_t bytes, bool write, off_t base = 0) {
    for (size_t pos = 0; pos < bytes;) {
        const size_t count = std::min(kChunk, bytes - pos);
        const ssize_t n = write ? pwrite(fd, data + pos, count, base + pos)
                                : pread(fd, data + pos, count, base + pos);
        if (n < 0 && errno == EINTR) continue;
        if (n <= 0) { if (n == 0) errno = EIO; return false; }
        pos += static_cast<size_t>(n);
    }
    return true;
}

void ReadCase(const char *name, int repeat, int fd, char *dst, const char *src,
              size_t bytes, int rank, int ranks) {
    MPI_Barrier(MPI_COMM_WORLD);
    const double t0 = MPI_Wtime();
    const bool ok = Transfer(fd, dst, bytes, false);
    const double seconds = MPI_Wtime() - t0;
    Require(ok, name, rank);
    Require(memcmp(dst, src, bytes) == 0, "read-content", rank);
    Report(name, repeat, bytes, seconds, rank, ranks);
}
} // namespace

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    int rank = 0, ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    if (argc != 7) {
        if (rank == 0) fprintf(stderr,
            "Usage: probe DIRECTORY MEMORY_MIB FILE_MIB REPEATS COPY_ITERATIONS EXPECTED_RANKS\n");
        MPI_Finalize();
        return 2;
    }
    const size_t mib = Number(argv[2], 4096);
    const size_t file_mib = Number(argv[3], 4096);
    const size_t repeats = Number(argv[4], 10);
    const size_t iterations = Number(argv[5], 1024);
    const size_t expected = Number(argv[6], 1024);
    Require(mib && file_mib && file_mib <= mib && repeats && iterations &&
            expected == static_cast<size_t>(ranks), "arguments", rank);
    const size_t bytes = mib * 1024u * 1024u;
    const size_t file_bytes = file_mib * 1024u * 1024u;
    void *src_raw = malloc(bytes + 4095);
    void *dst_raw = malloc(bytes + 4095);
    void *scratch_raw = malloc(kChunk + 4095);
    Require(src_raw && dst_raw && scratch_raw, "malloc", rank);
    char *src = Align(src_raw), *dst = Align(dst_raw), *scratch = Align(scratch_raw);
    Fill(src, bytes, rank);
    memset(dst, 0, bytes);
    memset(scratch, 0, kChunk);
    char host[MPI_MAX_PROCESSOR_NAME] = {};
    int host_length = 0;
    MPI_Get_processor_name(host, &host_length);
    std::vector<char> hosts(static_cast<size_t>(ranks) * MPI_MAX_PROCESSOR_NAME);
    MPI_Gather(host, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, hosts.data(),
               MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("IO_CONFIG ranks=%d memory_MiB=%zu file_MiB=%zu repeats=%zu "
               "iterations=%zu io_chunk_bytes=%zu units=decimal_GB_per_second\n",
               ranks, mib, file_mib, repeats, iterations, kChunk);
        printf("IO_NOTES bandwidth_counts_payload_not_double_memory_traffic; "
               "allocation_pretouch_open_validation_cleanup_outside_timers; "
               "writes_include_sync; fadvise_is_advisory; O_DIRECT_bypasses_client_cache_only\n");
        for (int r = 0; r < ranks; ++r)
            printf("IO_HOST rank=%d name=%s\n", r,
                   hosts.data() + static_cast<size_t>(r) * MPI_MAX_PROCESSOR_NAME);
        fflush(stdout);
    }

    for (int repeat = 1; repeat <= static_cast<int>(repeats); ++repeat) {
        const char *names[] = {"memcpy_64KiB", "memcpy_8MiB", "assemble_8MiB"};
        for (int mode = 0; mode < 3; ++mode) {
            const size_t chunk = mode == 0 ? 64u * 1024u : kChunk;
            Copy(dst, src, scratch, bytes, chunk, mode == 2);
            MPI_Barrier(MPI_COMM_WORLD);
            const double t0 = MPI_Wtime();
            for (size_t i = 0; i < iterations; ++i)
                Copy(dst, src, scratch, bytes, chunk, mode == 2);
            const double seconds = MPI_Wtime() - t0;
            Require(memcmp(dst, src, bytes) == 0, "copy-content", rank);
            guard ^= static_cast<unsigned char>(dst[bytes - 1]);
            Report(names[mode], repeat, static_cast<uint64_t>(bytes) * iterations,
                   seconds, rank, ranks);
        }

        const std::string path = std::string(argv[1]) + "/rank-" +
            std::to_string(rank) + "-r" + std::to_string(repeat) + ".bin";
        int fd = open(path.c_str(), O_CREAT | O_EXCL | O_RDWR, 0600);
        Require(fd >= 0, "file-open-exclusive", rank);
        MPI_Barrier(MPI_COMM_WORLD);
        double t0 = MPI_Wtime();
        bool ok = Transfer(fd, src, file_bytes, true);
        if (ok) ok = fsync(fd) == 0;
        double seconds = MPI_Wtime() - t0;
        Require(ok, "posix-write-fsync", rank);
        Report("posix_write_fsync", repeat, file_bytes, seconds, rank, ranks);
        Require(close(fd) == 0, "file-close", rank);
        fd = open(path.c_str(), O_RDONLY);
        Require(fd >= 0, "file-reopen", rank);
        const int advise = posix_fadvise(fd, 0, file_bytes, POSIX_FADV_DONTNEED);
        if (advise) fprintf(stderr, "[rank %d] IO_HINT fadvise_error=%d\n", rank, advise);
        ReadCase("posix_read_after_fadvise", repeat, fd, dst, src,
                 file_bytes, rank, ranks);
        ReadCase("posix_read_repeat", repeat, fd, dst, src, file_bytes, rank, ranks);
        Require(close(fd) == 0, "read-close", rank);
#ifdef O_DIRECT
        const int direct = open(path.c_str(), O_RDONLY | O_DIRECT);
        const int direct_error = direct < 0 ? errno : 0;
        if (AllOk(direct >= 0)) {
            MPI_Barrier(MPI_COMM_WORLD);
            t0 = MPI_Wtime();
            ok = Transfer(direct, dst, file_bytes, false);
            seconds = MPI_Wtime() - t0;
            if (AllOk(ok)) {
                Require(memcmp(dst, src, file_bytes) == 0, "direct-read-content", rank);
                Report("posix_read_O_DIRECT", repeat, file_bytes, seconds, rank, ranks);
            } else if (rank == 0) {
                printf("IO_SKIP name=posix_read_O_DIRECT repeat=%d reason=read-unsupported\n", repeat);
            }
        } else {
            fprintf(stderr, "[rank %d] IO_SKIP name=posix_read_O_DIRECT repeat=%d open_errno=%d\n",
                    rank, repeat, direct_error);
        }
        Require(direct < 0 || close(direct) == 0, "direct-close", rank);
#else
        if (rank == 0) printf("IO_SKIP name=posix_read_O_DIRECT reason=not-compiled\n");
#endif
        Require(unlink(path.c_str()) == 0, "rank-file-cleanup", rank);

        const std::string shared = std::string(argv[1]) + "/shared-r" +
            std::to_string(repeat) + ".bin";
        MPI_File file = MPI_FILE_NULL;
        int status = MPI_File_open(MPI_COMM_WORLD, const_cast<char *>(shared.c_str()),
            MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
        Require(status == MPI_SUCCESS, "mpi-file-open", rank);
        MPI_File_set_errhandler(file, MPI_ERRORS_RETURN);
        Require(MPI_File_set_size(file, static_cast<MPI_Offset>(file_bytes) * ranks)
                    == MPI_SUCCESS, "mpi-file-size", rank);
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = MPI_Wtime();
        ok = true;
        for (size_t pos = 0; pos < file_bytes;) {
            const int count = static_cast<int>(std::min(kChunk, file_bytes - pos));
            MPI_Status written;
            status = MPI_File_write_at(file,
                static_cast<MPI_Offset>(file_bytes) * rank + pos,
                src + pos, count, MPI_BYTE, &written);
            int actual = 0;
            if (status != MPI_SUCCESS || MPI_Get_count(&written, MPI_BYTE, &actual)
                    != MPI_SUCCESS || actual != count) { ok = false; break; }
            pos += static_cast<size_t>(count);
        }
        Require(ok, "mpi-file-write", rank);
        ok = MPI_File_sync(file) == MPI_SUCCESS;
        seconds = MPI_Wtime() - t0;
        Require(ok, "mpi-file-sync", rank);
        Report("mpi_shared_write_sync", repeat, file_bytes, seconds, rank, ranks);
        Require(MPI_File_close(&file) == MPI_SUCCESS, "mpi-file-close", rank);
        fd = open(shared.c_str(), O_RDONLY);
        Require(fd >= 0, "shared-verify-open", rank);
        Require(Transfer(fd, dst, file_bytes, false,
                         static_cast<off_t>(file_bytes) * rank), "shared-verify-read", rank);
        Require(memcmp(dst, src, file_bytes) == 0, "shared-content", rank);
        Require(close(fd) == 0, "shared-verify-close", rank);
        MPI_Barrier(MPI_COMM_WORLD);
        Require(rank != 0 || unlink(shared.c_str()) == 0, "shared-cleanup", rank);
    }
    free(src_raw);
    free(dst_raw);
    free(scratch_raw);
    if (rank == 0) {
        printf("IO_PROBE_COMPLETE ranks=%d status=0\n", ranks);
        fflush(stdout);
    }
    MPI_Finalize();
    return 0;
}
