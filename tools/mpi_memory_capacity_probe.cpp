#include <mpi.h>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <map>
#include <string>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <vector>

namespace {

struct Options {
    uint64_t allocate_bytes;
    uint64_t chunk_bytes;
    uint64_t touch_stride;
    std::string file_path;
    std::string file_mode;
    uint64_t file_limit;
    uint64_t read_buffer_bytes;

    Options()
        : allocate_bytes(512ULL << 20),
          chunk_bytes(256ULL << 20),
          touch_stride(4096),
          file_mode("replicated"),
          file_limit(0),
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
        "Usage: %s [--allocate 512M] [--chunk 256M] "
        "[--touch-stride 4K] [--file PATH] "
        "[--file-mode replicated|partitioned] "
        "[--file-limit SIZE] [--read-buffer 8M]\n",
        program ? program : "swbam-mpi-memory-probe");
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
        } else if (std::strcmp(arg, "--chunk") == 0) {
            options->chunk_bytes = ParseSize(value, &ok);
        } else if (std::strcmp(arg, "--touch-stride") == 0) {
            options->touch_stride = ParseSize(value, &ok);
        } else if (std::strcmp(arg, "--file") == 0) {
            options->file_path = value;
            ok = !options->file_path.empty();
        } else if (std::strcmp(arg, "--file-mode") == 0) {
            options->file_mode = value;
            ok = options->file_mode == "replicated" ||
                 options->file_mode == "partitioned";
        } else if (std::strcmp(arg, "--file-limit") == 0) {
            options->file_limit = ParseSize(value, &ok);
        } else if (std::strcmp(arg, "--read-buffer") == 0) {
            options->read_buffer_bytes = ParseSize(value, &ok);
        } else {
            return -1;
        }
        if (!ok) return -1;
    }
    if (options->chunk_bytes == 0 || options->touch_stride == 0 ||
        options->read_buffer_bytes == 0) {
        return -1;
    }
    return 0;
}

double BytesToGiB(uint64_t bytes) {
    return (double)bytes / (double)(1ULL << 30);
}

uint64_t SaturatingMul(uint64_t a, uint64_t b) {
    return b != 0 && a > UINT64_MAX / b ? UINT64_MAX : a * b;
}

uint64_t ReadMemAvailableBytes() {
    FILE *fp = std::fopen("/proc/meminfo", "r");
    if (!fp) return 0;
    char key[64];
    unsigned long long kib = 0;
    char unit[32];
    uint64_t result = 0;
    while (std::fscanf(fp, "%63s %llu %31s", key, &kib, unit) == 3) {
        if (std::strcmp(key, "MemAvailable:") == 0) {
            result = (uint64_t)kib << 10;
            break;
        }
    }
    std::fclose(fp);
    return result;
}

uint64_t SystemTotalBytes() {
    struct sysinfo info;
    if (sysinfo(&info) != 0) return 0;
    return SaturatingMul((uint64_t)info.totalram, (uint64_t)info.mem_unit);
}

uint64_t AddressSpaceLimitBytes() {
    struct rlimit limit;
    if (getrlimit(RLIMIT_AS, &limit) != 0 || limit.rlim_cur == RLIM_INFINITY) {
        return 0;
    }
    return (uint64_t)limit.rlim_cur;
}

void PrintTopology(int rank, int world_size, const char *processor_name) {
    const int width = MPI_MAX_PROCESSOR_NAME + 1;
    std::vector<char> all_names;
    if (rank == 0) all_names.resize((size_t)world_size * width, 0);
    char local[MPI_MAX_PROCESSOR_NAME + 1];
    std::memset(local, 0, sizeof(local));
    std::strncpy(local, processor_name ? processor_name : "unknown",
                 sizeof(local) - 1);
    MPI_Gather(local, width, MPI_CHAR,
               rank == 0 ? all_names.data() : nullptr,
               width, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (rank != 0) return;

    std::map<std::string, std::vector<int> > nodes;
    for (int i = 0; i < world_size; ++i) {
        const char *name = all_names.data() + (size_t)i * width;
        nodes[name].push_back(i);
    }
    std::printf("topology world_ranks=%d unique_nodes=%zu\n",
                world_size, nodes.size());
    for (std::map<std::string, std::vector<int> >::const_iterator it =
             nodes.begin(); it != nodes.end(); ++it) {
        std::printf("  node=%s ranks=", it->first.c_str());
        for (size_t i = 0; i < it->second.size(); ++i) {
            std::printf("%s%d", i == 0 ? "" : ",", it->second[i]);
        }
        std::printf(" count=%zu\n", it->second.size());
    }
    std::fflush(stdout);
}

void PrintRankEnvironment(int rank, int world_size,
                          const char *processor_name) {
    const uint64_t total = SystemTotalBytes();
    const uint64_t available = ReadMemAvailableBytes();
    const uint64_t as_limit = AddressSpaceLimitBytes();
#ifdef PLATFORM_SUNWAY
    const int share_size = athread_get_share_size();
    const int cross_size = athread_get_cross_size();
    const int priv_size = athread_get_priv_size();
    const int cache_size = athread_get_cache_size();
#else
    const int share_size = -1;
    const int cross_size = -1;
    const int priv_size = -1;
    const int cache_size = -1;
#endif

    for (int r = 0; r < world_size; ++r) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == r) {
            std::printf(
                "[rank %d] node=%s pid=%ld linux_private_total=%.3f GiB "
                "linux_private_available=%.3f GiB rlimit_as=%s "
                "share_size=%d cross_size=%d priv_size=%d "
                "cache_size=%d\n",
                rank, processor_name ? processor_name : "unknown",
                (long)getpid(), BytesToGiB(total), BytesToGiB(available),
                as_limit == 0 ? "unlimited" : "limited",
                share_size, cross_size, priv_size, cache_size);
            if (as_limit != 0) {
                std::printf("[rank %d] rlimit_as_bytes=%llu (%.3f GiB)\n",
                            rank, (unsigned long long)as_limit,
                            BytesToGiB(as_limit));
            }
            std::fflush(stdout);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

struct AllocationResult {
    std::vector<void *> buffers;
    std::vector<size_t> sizes;
    uint64_t allocated;
    double seconds;
    int ok;

    AllocationResult() : allocated(0), seconds(0.0), ok(1) {}
};

AllocationResult AllocateAndTouch(const Options &options, int rank) {
    AllocationResult result;
    const double begin = MPI_Wtime();
    while (result.allocated < options.allocate_bytes) {
        const uint64_t remaining = options.allocate_bytes - result.allocated;
        const size_t bytes = (size_t)std::min(remaining, options.chunk_bytes);
        void *memory = nullptr;
        if (posix_memalign(&memory, 4096, bytes) != 0 || !memory) {
            result.ok = 0;
            break;
        }
        unsigned char *data = static_cast<unsigned char *>(memory);
        for (size_t offset = 0; offset < bytes;
             offset += (size_t)options.touch_stride) {
            data[offset] = (unsigned char)(rank + offset / options.touch_stride);
        }
        if (bytes != 0) data[bytes - 1] = (unsigned char)rank;
        result.buffers.push_back(memory);
        result.sizes.push_back(bytes);
        result.allocated += bytes;
    }
    result.seconds = MPI_Wtime() - begin;
    return result;
}

void FreeAllocation(AllocationResult *result) {
    if (!result) return;
    for (size_t i = 0; i < result->buffers.size(); ++i) {
        std::free(result->buffers[i]);
    }
    result->buffers.clear();
    result->sizes.clear();
}

int RunFileRead(const Options &options, int rank, int world_size) {
    if (options.file_path.empty()) return 0;
    int fd = open(options.file_path.c_str(), O_RDONLY);
    struct stat st;
    int local_ok = fd >= 0 && fstat(fd, &st) == 0 && st.st_size >= 0;
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    if (!global_ok) {
        if (rank == 0) {
            std::fprintf(stderr, "ERROR: cannot open probe file %s\n",
                         options.file_path.c_str());
        }
        if (fd >= 0) close(fd);
        return -1;
    }

    uint64_t logical_size = (uint64_t)st.st_size;
    if (options.file_limit != 0) {
        logical_size = std::min(logical_size, options.file_limit);
    }
    uint64_t begin = 0;
    uint64_t end = logical_size;
    if (options.file_mode == "partitioned") {
        begin = logical_size * (uint64_t)rank / (uint64_t)world_size;
        end = logical_size * (uint64_t)(rank + 1) / (uint64_t)world_size;
    }

    std::vector<unsigned char> buffer((size_t)options.read_buffer_bytes);
    uint64_t offset = begin;
    uint64_t bytes_read = 0;
    uint64_t checksum = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    const double t0 = MPI_Wtime();
    while (offset < end) {
        const size_t request = (size_t)std::min(
            (uint64_t)buffer.size(), end - offset);
        const ssize_t got = pread(fd, buffer.data(), request, (off_t)offset);
        if (got <= 0) {
            local_ok = 0;
            break;
        }
        bytes_read += (uint64_t)got;
        offset += (uint64_t)got;
        checksum ^= (uint64_t)buffer[0] +
                    ((uint64_t)buffer[(size_t)got - 1] << 8);
    }
    const double seconds = MPI_Wtime() - t0;
    close(fd);

    uint64_t total_read = 0;
    double max_seconds = 0.0;
    MPI_Reduce(&bytes_read, &total_read, 1, MPI_UNSIGNED_LONG_LONG,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&seconds, &max_seconds, 1, MPI_DOUBLE,
               MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    if (rank == 0) {
        const double gbps = (double)total_read /
            std::max(max_seconds, 1e-12) / 1e9;
        std::printf(
            "file_read mode=%s logical_size=%.3f GiB "
            "physical_total=%.3f GiB max_time=%.6f s "
            "aggregate=%.3f GB/s status=%s "
            "note=page_cache_may_affect_result\n",
            options.file_mode.c_str(), BytesToGiB(logical_size),
            BytesToGiB(total_read), max_seconds, gbps,
            global_ok ? "ok" : "failed");
        std::fflush(stdout);
    }
    (void)checksum;
    return global_ok ? 0 : -1;
}

}  // namespace

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0;
    int world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    Options options;
    const int parse_result = ParseOptions(argc, argv, &options);
    if (parse_result != 0) {
        if (rank == 0 && parse_result < 0) PrintUsage(argv[0]);
        MPI_Finalize();
        return parse_result > 0 ? 0 : 2;
    }

    char processor[MPI_MAX_PROCESSOR_NAME + 1];
    int processor_length = 0;
    std::memset(processor, 0, sizeof(processor));
    MPI_Get_processor_name(processor, &processor_length);
    processor[std::min(processor_length, MPI_MAX_PROCESSOR_NAME)] = '\0';

    PrintTopology(rank, world_size, processor);
    PrintRankEnvironment(rank, world_size, processor);

    MPI_Barrier(MPI_COMM_WORLD);
    AllocationResult allocation = AllocateAndTouch(options, rank);
    double max_seconds = 0.0;
    uint64_t min_allocated = 0;
    uint64_t total_allocated = 0;
    int global_ok = 0;
    MPI_Reduce(&allocation.seconds, &max_seconds, 1, MPI_DOUBLE,
               MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&allocation.allocated, &min_allocated, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&allocation.allocated, &total_allocated, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Allreduce(&allocation.ok, &global_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    if (rank == 0) {
        std::printf(
            "allocation requested_per_rank=%.3f GiB "
            "minimum_per_rank=%.3f GiB total=%.3f GiB "
            "max_touch_time=%.6f s status=%s\n",
            BytesToGiB(options.allocate_bytes), BytesToGiB(min_allocated),
            BytesToGiB(total_allocated), max_seconds,
            global_ok ? "ok" : "failed");
        std::fflush(stdout);
    }

    const int file_status = RunFileRead(options, rank, world_size);
    FreeAllocation(&allocation);
    MPI_Finalize();
    return global_ok && file_status == 0 ? 0 : 1;
}
