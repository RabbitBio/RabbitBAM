#ifndef SWBAM_USAGE_NATIVE_STORAGE_H
#define SWBAM_USAGE_NATIVE_STORAGE_H
#include "swbam/mpi_runtime.h"
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
namespace swbam_usage {
typedef swbam::AdaptiveRankBodySink BodySink;
typedef swbam::SegmentedRankBodySource BodySequence;
class Session {
public:
    Session() : rank_(0) {}
    int Prepare(int *argc, char ***argv, bool example = false) {
        int n = 1;
        for (int i = 1; i < *argc; ++i) {
            const char *arg = (*argv)[i];
            if (!strcmp(arg, "--storage")) {
                if (++i >= *argc || strcmp((*argv)[i], "posix")) return Unsupported();
            } else if (!strcmp(arg, "--storage=posix")) {
            } else if (!strncmp(arg, "--ram-", 6) || !strcmp(arg, "--no-dump") || !strncmp(arg, "--storage=", 10)) {
                return Unsupported();
            } else if (example && (!strcmp(arg, "-o") || !strcmp(arg, "--output"))) {
                if (++i >= *argc) return -1;
                output_ = (*argv)[i];
            } else {
                if (example && arg[0] != '-' && input_.empty()) input_ = arg;
                (*argv)[n++] = (*argv)[i];
            }
        }
        *argc = n; (*argv)[n] = nullptr; return 0;
    }
    int EnterMpi() { return MPI_Comm_rank(MPI_COMM_WORLD, &rank_) == MPI_SUCCESS ? 0 : -1; }
    void LeaveMpi() {}
    int Finish(int status, double seconds) {
        printf("[swbam-e2e rank=%d] storage=posix seconds=%.6f status=%d\n", rank_, seconds, status);
        return status;
    }
    bool ram() const { return false; }
    const std::string &input_path() const { return input_; }
    const std::string &output_path() const { return output_; }
    swbam::SamInputBackend *sam_input() { return nullptr; }
private:
    int Unsupported() {
        fprintf(stderr, "ERROR: RAM mode requires SWBAM_BENCHMARK_SUPPORT_DIR at build time\n");
        return -1;
    }
    int rank_;
    std::string input_, output_;
};
inline bool IsRam() { return false; }
inline void BeginMpiInterval() {}
inline void EndMpiInterval() {}
inline double WallTime() {
    timespec time;
    clock_gettime(CLOCK_MONOTONIC, &time);
    return time.tv_sec + time.tv_nsec * 1e-9;
}
inline Session *CurrentSession() { return nullptr; }
inline int OpenInput(swbam::mpi::MpiBamInput &input, const std::string &path,
                      const std::string &backend, const std::string &limit) {
    return input.Open(path, backend, limit);
}
inline int WriteOutput(const std::string &path, const swbam::RankBodySource &body,
                        size_t prefix_size, const void *prefix, size_t suffix_size, const void *suffix) {
    swbam::mpi::DistributedBamOutput output;
    return output.Prepare(body, prefix_size, suffix_size) ||
        output.WriteMpiIo(path, body, prefix, suffix) ? -1 : 0;
}
inline int BuildBamPrefix(sam_hdr_t *header, int level, swbam::MemoryBamOutput *prefix) {
    if (!header || !prefix) return -1;
    const char *text = sam_hdr_str(header);
    const size_t length = sam_hdr_length(header);
    if (length > UINT32_MAX || (length && !text)) return -1;
    std::vector<unsigned char> raw = {'B', 'A', 'M', 1};
    auto put32 = [&](uint32_t value) {
        for (int i = 0; i < 4; ++i) raw.push_back((value >> (i * 8)) & 255);
    };
    put32(length);
    if (length) raw.insert(raw.end(), text, text + length);
    const int refs = sam_hdr_nref(header);
    if (refs < 0) return -1;
    put32(refs);
    for (int i = 0; i < refs; ++i) {
        const char *name = sam_hdr_tid2name(header, i);
        const hts_pos_t size = sam_hdr_tid2len(header, i);
        if (!name || size < 0 || size > UINT32_MAX || strlen(name) >= UINT32_MAX) return -1;
        put32(strlen(name) + 1);
        raw.insert(raw.end(), name, name + strlen(name) + 1);
        put32(size);
    }
    prefix->Clear();
    for (size_t pos = 0; pos < raw.size();) {
        unsigned char block[BGZF_MAX_BLOCK_SIZE];
        size_t bytes = sizeof(block);
        const size_t count = std::min<size_t>(BGZF_BLOCK_SIZE, raw.size() - pos);
        if (bgzf_compress(block, &bytes, raw.data() + pos, count, level) ||
            prefix->Append(block, bytes)) return -1;
        pos += count;
    }
    return 0;
}
inline int OpenTemporary(const std::string &path) {
    int fd = open(path.c_str(), O_RDWR | O_CREAT | O_EXCL, 0600);
    if (fd >= 0 && unlink(path.c_str())) { close(fd); return -1; }
    return fd;
}
inline int CloseTemporary(int fd) { return close(fd); }
inline ssize_t ReadTemporary(int fd, void *data, size_t bytes, off_t offset) { return pread(fd, data, bytes, offset); }
inline ssize_t WriteTemporary(int fd, const void *data, size_t bytes, off_t offset) { return pwrite(fd, data, bytes, offset); }
}
#endif
