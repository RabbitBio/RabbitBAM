#include "cgs_memory_probe_common.h"

#include <athread.h>

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <string>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

extern "C" void slave_cgs_memory_probe(CgsMemoryProbeArgs *args);

namespace {

struct Options {
    uint64_t allocate_bytes;
    uint64_t touch_stride;
    uint64_t read_buffer_bytes;
    std::string file_path;

    Options()
        : allocate_bytes(1ULL << 30),
          touch_stride(4096),
          read_buffer_bytes(8ULL << 20) {}
};

uint64_t ParseSize(const char *text, bool *ok) {
    if (ok) *ok = false;
    if (!text || !*text) return 0;
    errno = 0;
    char *end = nullptr;
    const unsigned long long value = std::strtoull(text, &end, 10);
    if (errno != 0 || end == text) return 0;

    uint64_t scale = 1;
    if (*end != '\0') {
        if (end[1] != '\0') return 0;
        switch (*end) {
        case 'k':
        case 'K': scale = 1ULL << 10; break;
        case 'm':
        case 'M': scale = 1ULL << 20; break;
        case 'g':
        case 'G': scale = 1ULL << 30; break;
        case 't':
        case 'T': scale = 1ULL << 40; break;
        default: return 0;
        }
    }
    if (value > ULLONG_MAX / scale) return 0;
    if (ok) *ok = true;
    return (uint64_t)value * scale;
}

void PrintUsage(const char *program) {
    std::fprintf(
        stderr,
        "Usage: %s [--allocate 1G] [--touch-stride 4K] "
        "[--file PATH] [--read-buffer 8M]\n",
        program ? program : "swbam-cgs-memory-probe");
}

int ParseOptions(int argc, char **argv, Options *options) {
    if (!options) return -1;
    for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        if (std::strcmp(arg, "--help") == 0 ||
            std::strcmp(arg, "-h") == 0) {
            PrintUsage(argv[0]);
            return 1;
        }
        if (i + 1 >= argc) return -1;
        const char *value = argv[++i];
        bool ok = false;
        if (std::strcmp(arg, "--allocate") == 0) {
            options->allocate_bytes = ParseSize(value, &ok);
        } else if (std::strcmp(arg, "--touch-stride") == 0) {
            options->touch_stride = ParseSize(value, &ok);
        } else if (std::strcmp(arg, "--read-buffer") == 0) {
            options->read_buffer_bytes = ParseSize(value, &ok);
        } else if (std::strcmp(arg, "--file") == 0) {
            options->file_path = value;
            ok = !options->file_path.empty();
        } else {
            return -1;
        }
        if (!ok) return -1;
    }
    if (options->allocate_bytes == 0 || options->touch_stride == 0 ||
        options->read_buffer_bytes == 0) {
        return -1;
    }
    return 0;
}

double GiB(uint64_t bytes) {
    return (double)bytes / (double)(1ULL << 30);
}

double Now() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

void ResetResults(CgsMemoryProbeArgs *args) {
    std::memset(args->counts, 0,
                SWBAM_CGS_PROBE_THREADS * sizeof(uint64_t));
    std::memset(args->checksums, 0,
                SWBAM_CGS_PROBE_THREADS * sizeof(uint64_t));
    std::memset(args->seen, 0,
                SWBAM_CGS_PROBE_THREADS * sizeof(int));
}

int SpawnProbe(CgsMemoryProbeArgs *args, const char *name,
               bool cgs, double *seconds) {
    ResetResults(args);
    const double begin = Now();
    const int spawn_status = cgs
        ? __real_athread_spawn_cgs(
              (void *)slave_cgs_memory_probe, args, 1)
        : __real_athread_spawn(
              (void *)slave_cgs_memory_probe, args, 1);
    if (spawn_status != 0) {
        std::fprintf(stderr, "ERROR: %s spawn_cgs status=%d\n",
                     name, spawn_status);
        return -1;
    }
    const int join_status = cgs ? athread_join_cgs() : athread_join();
    const double end = Now();
    if (join_status != 0) {
        std::fprintf(stderr, "ERROR: %s join_cgs status=%d\n",
                     name, join_status);
        return -1;
    }
    if (seconds) *seconds = end - begin;
    return 0;
}

int VerifyWorkers(const CgsMemoryProbeArgs *args,
                  uint64_t expected_count) {
    uint64_t actual_count = 0;
    for (int i = 0; i < args->worker_count; ++i) {
        if (args->seen[i] != i + 1) {
            std::fprintf(stderr,
                         "ERROR: CPE %d did not report, seen=%d\n",
                         i, args->seen[i]);
            return -1;
        }
        actual_count += args->counts[i];
    }
    if (actual_count != expected_count) {
        std::fprintf(stderr,
                     "ERROR: CPE work count mismatch expected=%llu "
                     "actual=%llu\n",
                     (unsigned long long)expected_count,
                     (unsigned long long)actual_count);
        return -1;
    }
    return 0;
}

int TouchAndVerify(CgsMemoryProbeArgs *args, double *cpe_seconds,
                   double *verify_seconds, bool cgs) {
    args->mode = SWBAM_CGS_PROBE_TOUCH;
    args->use_cgs_tid = cgs ? 1 : 0;
    args->worker_count = cgs
        ? SWBAM_CGS_PROBE_THREADS : SWBAM_CGS_PROBE_PES_PER_GROUP;
    if (SpawnProbe(args, cgs ? "cgs touch" : "single-CG touch",
                   cgs, cpe_seconds) != 0) {
        return -1;
    }
    const uint64_t pages =
        (args->bytes + args->stride - 1ULL) / args->stride;
    if (VerifyWorkers(args, pages) != 0) return -1;

    const double begin = Now();
    for (uint64_t page = 0; page < pages; ++page) {
        const uint64_t offset = page * args->stride;
        if (offset < args->bytes &&
            args->data[offset] != CgsMemoryProbePattern(page)) {
            std::fprintf(stderr,
                         "ERROR: cross memory mismatch page=%llu "
                         "offset=%llu expected=%u actual=%u\n",
                         (unsigned long long)page,
                         (unsigned long long)offset,
                         (unsigned)CgsMemoryProbePattern(page),
                         (unsigned)args->data[offset]);
            return -1;
        }
    }
    const double end = Now();
    if (verify_seconds) *verify_seconds = end - begin;
    return 0;
}

int LoadFile(const Options &options, unsigned char *data,
             uint64_t capacity, uint64_t *loaded, double *seconds) {
    if (loaded) *loaded = 0;
    if (options.file_path.empty()) return 0;
    const int fd = open(options.file_path.c_str(), O_RDONLY);
    if (fd < 0) {
        std::fprintf(stderr, "ERROR: open %s: %s\n",
                     options.file_path.c_str(), std::strerror(errno));
        return -1;
    }
    struct stat st;
    if (fstat(fd, &st) != 0 || st.st_size < 0) {
        std::fprintf(stderr, "ERROR: stat %s: %s\n",
                     options.file_path.c_str(), std::strerror(errno));
        close(fd);
        return -1;
    }
    const uint64_t file_size = (uint64_t)st.st_size;
    if (file_size > capacity) {
        std::fprintf(stderr,
                     "ERROR: file %.3f GiB exceeds cross allocation "
                     "%.3f GiB\n",
                     GiB(file_size), GiB(capacity));
        close(fd);
        return -1;
    }

    uint64_t offset = 0;
    const double begin = Now();
    while (offset < file_size) {
        const size_t request = (size_t)std::min(
            options.read_buffer_bytes, file_size - offset);
        const ssize_t got = pread(
            fd, data + offset, request, (off_t)offset);
        if (got <= 0) {
            std::fprintf(stderr, "ERROR: pread %s at %llu: %s\n",
                         options.file_path.c_str(),
                         (unsigned long long)offset,
                         got == 0 ? "unexpected EOF" :
                                    std::strerror(errno));
            close(fd);
            return -1;
        }
        offset += (uint64_t)got;
    }
    const double end = Now();
    close(fd);
    if (loaded) *loaded = offset;
    if (seconds) *seconds = end - begin;
    return 0;
}

uint64_t HostChecksum(const unsigned char *data, uint64_t bytes) {
    uint64_t checksum = 0;
    for (uint64_t offset = 0; offset < bytes; ++offset) {
        checksum += CgsMemoryProbeContribution(data[offset], offset);
    }
    return checksum;
}

int ChecksumFile(CgsMemoryProbeArgs *args, uint64_t bytes,
                 double *cpe_seconds, double *host_seconds,
                 uint64_t *checksum) {
    args->bytes = bytes;
    args->mode = SWBAM_CGS_PROBE_CHECKSUM;
    args->use_cgs_tid = 1;
    args->worker_count = SWBAM_CGS_PROBE_THREADS;
    if (SpawnProbe(args, "checksum", true, cpe_seconds) != 0) return -1;
    if (VerifyWorkers(args, bytes) != 0) return -1;
    uint64_t cpe_checksum = 0;
    for (int i = 0; i < SWBAM_CGS_PROBE_THREADS; ++i) {
        cpe_checksum += args->checksums[i];
    }

    const double begin = Now();
    const uint64_t host_checksum = HostChecksum(args->data, bytes);
    const double end = Now();
    if (host_seconds) *host_seconds = end - begin;
    if (checksum) *checksum = cpe_checksum;
    if (cpe_checksum != host_checksum) {
        std::fprintf(stderr,
                     "ERROR: file checksum mismatch cpe=%llu host=%llu\n",
                     (unsigned long long)cpe_checksum,
                     (unsigned long long)host_checksum);
        return -1;
    }
    return 0;
}

}  // namespace

int main(int argc, char **argv) {
    Options options;
    const int parse_status = ParseOptions(argc, argv, &options);
    if (parse_status != 0) {
        if (parse_status < 0) PrintUsage(argv[0]);
        return parse_status > 0 ? 0 : 2;
    }

    athread_init_cgs();
    std::printf(
        "cgs_resources threads=%d max_threads=%d pemask=0x%lx "
        "share_size=%d cross_size=%d priv_size=%d cache_size=%d\n",
        athread_get_num_threads(), athread_get_max_threads(),
        athread_get_pemask(), athread_get_share_size(),
        athread_get_cross_size(), athread_get_priv_size(),
        athread_get_cache_size());

    unsigned char *data = static_cast<unsigned char *>(
        _sw_xmalloc((size_t)options.allocate_bytes));
    CgsMemoryProbeArgs *args = static_cast<CgsMemoryProbeArgs *>(
        _sw_xmalloc(sizeof(CgsMemoryProbeArgs)));
    uint64_t *counts = static_cast<uint64_t *>(
        _sw_xmalloc(SWBAM_CGS_PROBE_THREADS * sizeof(uint64_t)));
    uint64_t *checksums = static_cast<uint64_t *>(
        _sw_xmalloc(SWBAM_CGS_PROBE_THREADS * sizeof(uint64_t)));
    int *seen = static_cast<int *>(
        _sw_xmalloc(SWBAM_CGS_PROBE_THREADS * sizeof(int)));
    if (!data || !args || !counts || !checksums || !seen) {
        std::fprintf(stderr,
                     "ERROR: _sw_xmalloc failed for %.3f GiB; "
                     "increase -cross_size or lower --allocate\n",
                     GiB(options.allocate_bytes));
        if (seen) _sw_xfree(seen);
        if (checksums) _sw_xfree(checksums);
        if (counts) _sw_xfree(counts);
        if (args) _sw_xfree(args);
        if (data) _sw_xfree(data);
        athread_halt();
        return 1;
    }

    args->data = data;
    args->bytes = options.allocate_bytes;
    args->stride = options.touch_stride;
    args->counts = counts;
    args->checksums = checksums;
    args->seen = seen;
    args->mode = 0;
    args->worker_count = SWBAM_CGS_PROBE_THREADS;
    args->use_cgs_tid = 1;

    std::printf(
        "cross_allocation address=%p bytes=%llu size=%.3f GiB "
        "touch_stride=%llu\n",
        (void *)data, (unsigned long long)options.allocate_bytes,
        GiB(options.allocate_bytes),
        (unsigned long long)options.touch_stride);

    int status = 0;
    double cpe_touch = 0.0;
    double host_verify = 0.0;
    if (TouchAndVerify(args, &cpe_touch, &host_verify, true) != 0) {
        status = 1;
    } else {
        const double touched_gib = GiB(options.allocate_bytes);
        std::printf(
            "cross_touch status=ok cpe_threads=%d cpe_time=%.6f s "
            "host_verify=%.6f s address_span_rate=%.3f GiB/s\n",
            SWBAM_CGS_PROBE_THREADS, cpe_touch, host_verify,
            touched_gib / std::max(cpe_touch, 1e-12));
    }

    if (status == 0) {
        const uint64_t original_bytes = args->bytes;
        args->bytes = std::min<uint64_t>(original_bytes, 1ULL << 30);
        double single_touch = 0.0;
        double single_verify = 0.0;
        if (TouchAndVerify(args, &single_touch, &single_verify,
                           false) != 0) {
            status = 1;
        } else {
            std::printf(
                "single_cg_cross_access status=ok cpe_threads=%d "
                "span=%.3f GiB cpe_time=%.6f s "
                "host_verify=%.6f s\n",
                SWBAM_CGS_PROBE_PES_PER_GROUP, GiB(args->bytes),
                single_touch, single_verify);
        }
        args->bytes = original_bytes;
    }

    if (status == 0 && !options.file_path.empty()) {
        uint64_t loaded = 0;
        double read_seconds = 0.0;
        if (LoadFile(options, data, options.allocate_bytes,
                     &loaded, &read_seconds) != 0) {
            status = 1;
        } else {
            double cpe_checksum_time = 0.0;
            double host_checksum_time = 0.0;
            uint64_t checksum = 0;
            if (ChecksumFile(args, loaded, &cpe_checksum_time,
                             &host_checksum_time, &checksum) != 0) {
                status = 1;
            } else {
                std::printf(
                    "cross_file_load status=ok file=%s bytes=%llu "
                    "read_time=%.6f s read_rate=%.3f GiB/s "
                    "cpe_checksum_time=%.6f s "
                    "host_checksum_time=%.6f s checksum=%llu\n",
                    options.file_path.c_str(),
                    (unsigned long long)loaded, read_seconds,
                    GiB(loaded) / std::max(read_seconds, 1e-12),
                    cpe_checksum_time, host_checksum_time,
                    (unsigned long long)checksum);
            }
        }
    }

    _sw_xfree(seen);
    _sw_xfree(checksums);
    _sw_xfree(counts);
    _sw_xfree(args);
    _sw_xfree(data);
    athread_halt();
    return status;
}
