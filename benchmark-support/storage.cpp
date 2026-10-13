#include "storage.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <map>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

namespace swbam_usage {
namespace {
Session *current = nullptr;
std::map<int, std::shared_ptr<RamFile> > temporary;
int next_handle = 0x40000000;
uint64_t temp_read_bytes = 0, temp_write_bytes = 0;
size_t work_memory = 0;

double Now() {
    timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}
int Transfer(int fd, uint64_t offset, void *data, size_t bytes, bool write) {
    char *p = static_cast<char *>(data);
    while (bytes) {
        const size_t count = std::min<size_t>(bytes, 8u * 1024u * 1024u);
        const ssize_t n = write ? pwrite(fd, p, count, offset) : pread(fd, p, count, offset);
        if (n < 0 && errno == EINTR) continue;
        if (n <= 0) return -1;
        p += n; offset += n; bytes -= n;
    }
    return 0;
}
int Mark(const std::string &path) {
    int fd = open(path.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0600);
    return fd < 0 ? -1 : close(fd);
}
bool Exists(const std::string &path) { return access(path.c_str(), F_OK) == 0; }
struct Manifest {
    uint64_t magic, size, body, prefix, count;
    int kind, ranks;
};
struct BlockHeader {
    uint64_t offset;
    unsigned char bytes[18];
};
struct Boundary {
    uint64_t offset, bytes;
};
const uint64_t kManifestMagic = 0x535742414d52414dULL;

// Only raw bytes are staged. Headers and partition plans are rebuilt in MPI.
struct Staged {
    Manifest manifest = {};
    std::shared_ptr<RamFile> prefix, catalog, body;
    uint64_t local_offset = 0;
    std::vector<Boundary> boundaries;
    char *Range(uint64_t offset, size_t bytes) const {
        if (offset < local_offset || offset - local_offset > body->size() ||
            bytes > body->size() - (offset - local_offset)) return nullptr;
        return body->Range(offset - local_offset, bytes);
    }
};

class RamBamInput : public swbam::BamInputBackend {
public:
    explicit RamBamInput(Staged *staged) : staged_(staged) {}
    int Open(const std::string &) {
        static const unsigned char eof[28] = {
            31,139,8,4,0,0,0,0,0,255,6,0,66,67,2,0,27,0,3,0,0,0,0,0,0,0,0,0
        };
        std::vector<char> bytes(staged_->prefix->size() + sizeof(eof));
        if (staged_->prefix->ReadAt(0, bytes.data(), staged_->prefix->size())) return -1;
        memcpy(bytes.data() + staged_->prefix->size(), eof, sizeof(eof));
        return parsed_.OpenMemoryCopy(bytes.data(), bytes.size());
    }
    void Close() { parsed_.Close(); }
    const sam_hdr_t *header() const { return parsed_.header(); }
    int format() const { return parsed_.format(); }
    size_t size() const { return staged_->manifest.size; }
    uint64_t body_offset() const { return parsed_.body_offset(); }
    char *ResidentRange(uint64_t offset, size_t bytes) const {
        return staged_->Range(offset, bytes);
    }
    int ScanBlocks(std::vector<swbam::BgzfBlockSpan> *spans) const {
        spans->clear();
        uint64_t expected = body_offset();
        for (uint64_t i = 0; i < staged_->manifest.count; ++i) {
            BlockHeader h;
            if (staged_->catalog->ReadAt(i * sizeof(h), &h, sizeof(h))) return -1;
            const unsigned char *p = h.bytes;
            const uint32_t bytes = uint32_t(p[16]) + (uint32_t(p[17]) << 8) + 1;
            if (h.offset != expected || p[0] != 31 || p[1] != 139 || p[2] != 8 ||
                p[3] != 4 || p[10] != 6 || p[11] != 0 || p[12] != 'B' ||
                p[13] != 'C' || p[14] != 2 || p[15] != 0 || bytes < 26 ||
                bytes > BGZF_MAX_BLOCK_SIZE || expected > size() || bytes > size() - expected)
                return -1;
            spans->push_back({expected, bytes});
            expected += bytes;
        }
        return expected <= size() && size() - expected <= 28 ? 0 : -1;
    }
    int ReadBatch(const swbam::BgzfBlockSpan *spans, size_t count,
                  swbam::BgzfBlockBatch *batch) const {
        if (count > batch->capacity()) return -1;
        for (size_t i = 0; i < count; ++i) {
            const size_t bytes = spans[i].compressed_size;
            char *source = ResidentRange(spans[i].offset, bytes);
            if (!source || bytes > BGZF_MAX_BLOCK_SIZE) { errno = ENODATA; return -1; }
            bam_block &b = batch->blocks()[i];
            memcpy(b.data, source, bytes);
            b.length = bytes; b.pos = 0; b.errcode = 0; b.block_id = i;
        }
        return 0;
    }
private:
    Staged *staged_;
    swbam::MemoryBamInput parsed_;
};

class RamSamInput : public swbam::SamInputBackend {
public:
    explicit RamSamInput(Staged *staged) : staged_(staged), header_(nullptr), body_(0) {}
    ~RamSamInput() { Close(); }
    int Open(const std::string &) {
        Close();
        const char *p = staged_->prefix->Range(0, staged_->prefix->size());
        if (!p) return -1;
        size_t pos = 0;
        while (pos < staged_->prefix->size() && p[pos] == '@') {
            const char *end = static_cast<const char *>(memchr(p + pos, '\n', staged_->prefix->size() - pos));
            if (!end) return -1;
            pos = end - p + 1;
        }
        body_ = pos;
        header_ = sam_hdr_parse(pos, p);
        return header_ && body_ == staged_->manifest.body ? 0 : -1;
    }
    void Close() { if (header_) sam_hdr_destroy(header_); header_ = nullptr; }
    const sam_hdr_t *header() const { return header_; }
    size_t size() const { return staged_->manifest.size; }
    uint64_t body_offset() const { return body_; }
    char *ResidentRange(uint64_t offset, size_t bytes) const { return staged_->Range(offset, bytes); }
    int ReadRange(uint64_t offset, size_t bytes, std::vector<char> *data) const {
        char *p = ResidentRange(offset, bytes);
        if (!p && bytes) return -1;
        data->resize(bytes);
        if (bytes) memcpy(data->data(), p, bytes);
        return 0;
    }
    int SplitRanges(int ranks, std::vector<long long> *offsets,
                    std::vector<long long> *lengths) const {
        if (ranks != staged_->manifest.ranks) return -1;
        std::vector<uint64_t> cuts(ranks + 1, body_);
        cuts.back() = size();
        uint64_t cursor = 0;
        for (int r = 1; r < ranks; ++r) {
            const uint64_t bytes = size() - body_;
            const uint64_t target = body_ + bytes / ranks * r + bytes % ranks * r / ranks;
            const Boundary &window = staged_->boundaries[r - 1];
            if (target == 0 || target == size()) { cuts[r] = target; continue; }
            std::vector<char> raw(window.bytes);
            if (window.offset != target - 1 ||
                staged_->catalog->ReadAt(cursor, raw.data(), raw.size())) return -1;
            cursor += window.bytes;
            if (raw[0] == '\n') cuts[r] = target;
            else {
                const char *end = static_cast<const char *>(memchr(raw.data() + 1, '\n', raw.size() - 1));
                cuts[r] = end ? window.offset + (end - raw.data()) + 1 : size();
            }
        }
        offsets->resize(ranks); lengths->resize(ranks);
        for (int r = 0; r < ranks; ++r) {
            if (cuts[r + 1] < cuts[r]) return -1;
            (*offsets)[r] = cuts[r]; (*lengths)[r] = cuts[r + 1] - cuts[r];
        }
        return 0;
    }
private:
    Staged *staged_;
    sam_hdr_t *header_;
    uint64_t body_;
};
}

struct Session::Impl {
    bool ram = false, no_dump = false, invalid_layout = false;
    int rank = 0, ranks = 1;
    size_t capacity = 0;
    uint64_t requested_work = 0, node_allowance = 0;
    double timeout = 7200;
    std::string directory, input, output;
    RamPool pool;
    Staged staged;
    std::unique_ptr<RamBamInput> bam;
    std::unique_ptr<RamSamInput> sam;
    std::vector<std::shared_ptr<RamFile> > final_parts;
    std::shared_ptr<RamFile> envelope;
    uint64_t output_size = 0, body_offset = 0, suffix_offset = 0;
    size_t prefix_bytes = 0, suffix_bytes = 0;

    std::string Marker(const std::string &name, int r) const {
        return directory + "/" + name + "." + std::to_string(r);
    }
    int Wait(const std::string &path, bool errors = true) const {
        const double start = Now();
        while (!Exists(path)) {
            if (errors)
                for (int r = 0; r < ranks; ++r)
                    if (Exists(Marker("error", r))) { errno = EIO; return -1; }
            if (Now() - start > timeout) {
                fprintf(stderr, "ERROR: RAM rendezvous timed out waiting for %s (expected ranks=%d)\n",
                        path.c_str(), ranks);
                errno = ETIMEDOUT; return -1;
            }
            usleep(100000);
        }
        return 0;
    }
    int StagePart(int fd, uint64_t source, uint64_t bytes, RamFile *file) {
        if (bytes > SIZE_MAX || file->ReserveContiguous(bytes)) return -1;
        std::vector<char> buffer(std::min<uint64_t>(bytes, 8u * 1024u * 1024u));
        for (uint64_t pos = 0; pos < bytes;) {
            const size_t count = std::min<uint64_t>(bytes - pos, buffer.size());
            if (Transfer(fd, source + pos, buffer.data(), count, false) ||
                file->WriteAt(pos, buffer.data(), count)) return -1;
            pos += count;
        }
        return 0;
    }
    int MakeManifest() {
        Manifest m = {};
        m.magic = kManifestMagic; m.ranks = ranks;
        samFile *probe = sam_open(input.c_str(), "r");
        if (!probe) return -1;
        const int format = probe->format.format;
        const int close_result = sam_close(probe);
        if (close_result) return -1;
        if (format != ::bam && format != ::sam) { errno = EINVAL; return -1; }
        const bool is_bam = format == ::bam;
        swbam::PosixBamInput bi;
        swbam::PosixSamInput si;
        if (is_bam ? bi.Open(input) : si.Open(input)) return -1;
        m.kind = is_bam ? 1 : 2;
        m.size = is_bam ? bi.size() : si.size();
        m.body = is_bam ? bi.body_offset() : si.body_offset();
        m.prefix = std::min<uint64_t>(m.size, m.body + BGZF_MAX_BLOCK_SIZE);
        int fd = open(input.c_str(), O_RDONLY);
        if (fd < 0) return -1;
        std::vector<BlockHeader> headers;
        std::vector<Boundary> windows;
        std::vector<char> raw_windows;
        bool ok = true;
        if (is_bam) {
            std::vector<char> cache(8u * 1024u * 1024u);
            uint64_t cached = UINT64_MAX;
            size_t available = 0;
            uint64_t offset = m.body;
            while (ok && offset < m.size) {
                if (cached == UINT64_MAX || offset < cached || offset + 28 > cached + available) {
                    cached = offset; available = std::min<uint64_t>(cache.size(), m.size - offset);
                    ok = Transfer(fd, cached, cache.data(), available, false) == 0;
                }
                if (!ok || m.size - offset < 18) { ok = false; break; }
                const unsigned char *p = reinterpret_cast<const unsigned char *>(cache.data() + offset - cached);
                const uint32_t bytes = uint32_t(p[16]) + (uint32_t(p[17]) << 8) + 1;
                if (p[0] != 31 || p[1] != 139 || bytes < 26 || bytes > BGZF_MAX_BLOCK_SIZE || bytes > m.size - offset) {
                    ok = false; break;
                }
                if (bytes == 28 && offset + bytes == m.size && p[18] == 3 && p[19] == 0) break;
                BlockHeader h = {}; h.offset = offset; memcpy(h.bytes, p, 18);
                headers.push_back(h); offset += bytes;
            }
            m.count = headers.size();
            m.prefix = m.body;
            if (!headers.empty())
                m.prefix += uint32_t(headers[0].bytes[16]) +
                            (uint32_t(headers[0].bytes[17]) << 8) + 1;
        } else {
            for (int r = 1; ok && r < ranks; ++r) {
                const uint64_t body_bytes = m.size - m.body;
                const uint64_t target = m.body + body_bytes / ranks * r + body_bytes % ranks * r / ranks;
                Boundary w = {target ? target - 1 : 0, 0};
                if (target && target != m.size) {
                    char c = 0;
                    do {
                        if (Transfer(fd, w.offset + w.bytes, &c, 1, false)) { ok = false; break; }
                        raw_windows.push_back(c); ++w.bytes;
                    } while (c != '\n' && w.offset + w.bytes < m.size);
                }
                windows.push_back(w);
            }
            m.count = windows.size();
        }
        int out = -1;
        if (ok) out = open((directory + "/manifest").c_str(), O_CREAT | O_EXCL | O_WRONLY, 0600);
        if (out < 0) ok = false;
        uint64_t pos = sizeof(m);
        if (ok) ok = Transfer(out, 0, &m, sizeof(m), true) == 0;
        std::vector<char> prefix(m.prefix);
        if (ok) ok = Transfer(fd, 0, prefix.data(), prefix.size(), false) == 0 &&
                      Transfer(out, pos, prefix.data(), prefix.size(), true) == 0;
        pos += prefix.size();
        if (ok && is_bam) ok = Transfer(out, pos, headers.data(), headers.size() * sizeof(BlockHeader), true) == 0;
        if (ok && !is_bam) {
            ok = Transfer(out, pos, windows.data(), windows.size() * sizeof(Boundary), true) == 0;
            pos += windows.size() * sizeof(Boundary);
            if (ok) ok = Transfer(out, pos, raw_windows.data(), raw_windows.size(), true) == 0;
        }
        if (out >= 0 && close(out)) ok = false;
        close(fd);
        return ok ? Mark(directory + "/manifest.ready") : -1;
    }
    int Load() {
        int fd = open((directory + "/manifest").c_str(), O_RDONLY);
        if (fd < 0) return -1;
        Manifest &m = staged.manifest;
        bool ok = Transfer(fd, 0, &m, sizeof(m), false) == 0 && m.magic == kManifestMagic && m.ranks == ranks;
        staged.prefix.reset(new RamFile(&pool));
        staged.catalog.reset(new RamFile(&pool));
        staged.body.reset(new RamFile(&pool));
        if (ok) ok = StagePart(fd, sizeof(m), m.prefix, staged.prefix.get()) == 0;
        uint64_t catalog_start = sizeof(m) + m.prefix;
        uint64_t begin = m.body, end = m.body;
        if (ok && m.kind == 1) {
            ok = m.count <= SIZE_MAX / sizeof(BlockHeader) &&
                 StagePart(fd, catalog_start, m.count * sizeof(BlockHeader), staged.catalog.get()) == 0;
            const uint64_t first = m.count * rank / ranks, last = m.count * (rank + 1) / ranks;
            if (ok && first < last) {
                BlockHeader a, b;
                ok = staged.catalog->ReadAt(first * sizeof(a), &a, sizeof(a)) == 0 &&
                     staged.catalog->ReadAt((last - 1) * sizeof(b), &b, sizeof(b)) == 0;
                begin = a.offset; end = b.offset + b.bytes[16] + (uint32_t(b.bytes[17]) << 8) + 1;
            }
        } else if (ok && m.kind == 2) {
            staged.boundaries.resize(m.count);
            ok = Transfer(fd, catalog_start, staged.boundaries.data(), m.count * sizeof(Boundary), false) == 0;
            catalog_start += m.count * sizeof(Boundary);
            uint64_t bytes = 0;
            for (const Boundary &w : staged.boundaries) bytes += w.bytes;
            if (ok) ok = StagePart(fd, catalog_start, bytes, staged.catalog.get()) == 0;
            if (ok) {
                RamSamInput split(&staged);
                std::vector<long long> offsets, lengths;
                ok = split.Open(input) == 0 && split.SplitRanges(ranks, &offsets, &lengths) == 0;
                if (ok) { begin = offsets[rank]; end = begin + lengths[rank]; }
            }
        } else ok = false;
        close(fd);
        staged.local_offset = begin;
        fd = ok ? open(input.c_str(), O_RDONLY) : -1;
        if (fd < 0) return -1;
        ok = StagePart(fd, begin, end - begin, staged.body.get()) == 0;
        close(fd);
        return ok ? 0 : -1;
    }
};

Session::Session() : impl_(new Impl) {}
Session::~Session() {
    temporary.clear();
    if (current == this) current = nullptr;
}
bool Session::ram() const { return impl_->ram; }
const std::string &Session::input_path() const { return impl_->input; }
const std::string &Session::output_path() const { return impl_->output; }
RamPool *Session::pool() { return &impl_->pool; }
Session *CurrentSession() { return current; }
bool IsRam() { return current && current->ram(); }
double WallTime() { return Now(); }
size_t WorkMemoryLimit() { return work_memory; }
size_t CacheMemoryLimit() {
    return std::min<size_t>(256u << 20, std::max<size_t>(8u << 20, work_memory / 16));
}

int Session::Prepare(int *argc, char ***argv, bool example) {
    current = this;
    std::vector<char *> keep;
    keep.push_back((*argv)[0]);
    std::string storage = "posix", capacity;
    bool legacy = false, ram_option = false;
    for (int i = 1; i < *argc; ++i) {
        std::string arg = (*argv)[i];
        const size_t eq = arg.find('=');
        const std::string key = arg.substr(0, eq);
        const bool common = key == "--storage" || key == "--work-memory" || key == "--ram-capacity" ||
            key == "--ram-ranks" || key == "--ram-session" || key == "--ram-ready-timeout";
        if (common || (example && (key == "-o" || key == "--output"))) {
            if (key.find("--ram-") == 0) ram_option = true;
            std::string value;
            if (eq != std::string::npos) value = arg.substr(eq + 1);
            else if (++i < *argc) value = (*argv)[i];
            else { fprintf(stderr, "ERROR: missing value for %s\n", key.c_str()); return -1; }
            if (key == "--storage") storage = value;
            else if (key == "--work-memory") {
                if (swbam::ParseByteSize(value, &impl_->requested_work) || !impl_->requested_work) return -1;
            }
            else if (key == "--ram-capacity") capacity = value;
            else if (key == "--ram-session") impl_->directory = value;
            else if (key == "--ram-ranks") {
                char *end = nullptr; long n = strtol(value.c_str(), &end, 10);
                if (!end || *end || n < 1 || n > INT_MAX) return -1;
                impl_->ranks = n;
            } else if (key == "--ram-ready-timeout") {
                char *end = nullptr; impl_->timeout = strtod(value.c_str(), &end);
                if (!end || *end || !(impl_->timeout > 0)) return -1;
            } else impl_->output = value;
            continue;
        }
        if (key == "--no-dump") {
            if (eq != std::string::npos) { fprintf(stderr, "ERROR: --no-dump takes no value\n"); return -1; }
            impl_->no_dump = true; continue;
        }
        if (key.find("--io-") == 0 || key.find("--rank-body-") == 0) legacy = true;
        if (example && arg.size() && arg[0] != '-' && impl_->input.empty()) impl_->input = arg;
        if (!example && (key == "-i" || key == "--inFile"))
            impl_->input = eq != std::string::npos ? arg.substr(eq + 1) : i + 1 < *argc ? (*argv)[i + 1] : "";
        if (!example && (key == "-o" || key == "--outFile"))
            impl_->output = eq != std::string::npos ? arg.substr(eq + 1) : i + 1 < *argc ? (*argv)[i + 1] : "";
        keep.push_back((*argv)[i]);
    }
    if (storage != "posix" && storage != "ram") { fprintf(stderr, "ERROR: --storage must be posix or ram\n"); return -1; }
    impl_->ram = storage == "ram";
    const char *node_memory = getenv("SWBAM_NODE_MEMORY");
    if (node_memory && (swbam::ParseByteSize(node_memory, &impl_->node_allowance) ||
                        !impl_->node_allowance)) return -1;
    for (size_t i = 0; i < keep.size(); ++i) (*argv)[i] = keep[i];
    *argc = keep.size(); (*argv)[*argc] = nullptr;
    if (!ram()) {
        if (impl_->no_dump || ram_option) {
            fprintf(stderr, "ERROR: RAM options require --storage ram\n"); return -1;
        }
        return 0;
    }
    uint64_t parsed = 0;
    const char *rank = getenv("RMS_RANK");
    char *end = nullptr;
    const long rank_number = rank ? strtol(rank, &end, 10) : 0;
    if ((rank && (!end || *end)) || (!rank && impl_->ranks > 1) || rank_number < 0 || rank_number >= impl_->ranks ||
        legacy || impl_->input.empty() || impl_->directory.empty() || impl_->directory[0] != '/' ||
        capacity.empty() || swbam::ParseByteSize(capacity, &parsed) || parsed > SIZE_MAX || parsed < RamPool::kPage ||
        (!impl_->output.empty() && impl_->input == impl_->output)) {
        fprintf(stderr, "ERROR: RAM mode requires explicit capacity, ranks, a unique absolute session directory and input; no legacy I/O overrides\n"); return -1;
    }
    for (int i = 1; i < *argc; ++i)
        if (!strcmp((*argv)[i], "dedup-workflow") || !strcmp((*argv)[i], "dedup-pipeline")) {
            fprintf(stderr, "ERROR: RAM mode supports separate commands, not dedup wrappers\n"); return -1;
        }
    impl_->rank = rank_number; impl_->capacity = parsed;
    if (impl_->rank == 0) {
        if (mkdir(impl_->directory.c_str(), 0700)) { perror("RAM session mkdir (must not exist)"); return -1; }
    } else if (impl_->Wait(impl_->directory, false)) return -1;
    const double start = Now();
    int result = 0;
    try {
        if (impl_->pool.Initialize(parsed)) result = -1;
        if (!result && impl_->rank == 0) result = impl_->MakeManifest();
        if (!result) result = impl_->Wait(impl_->directory + "/manifest.ready");
        if (!result) result = impl_->Load();
    } catch (...) { errno = ENOMEM; result = -1; }
    if (result) {
        fprintf(stderr, "[rank %d] ERROR: RAM preparation failed: %s\n", impl_->rank, strerror(errno));
        Mark(impl_->Marker("error", impl_->rank)); return -1;
    }
    if (Mark(impl_->Marker("ready", impl_->rank))) return -1;
    for (int r = 0; r < impl_->ranks; ++r)
        if (impl_->Wait(impl_->Marker("ready", r))) return -1;
    if (impl_->rank == 0) { unlink((impl_->directory + "/manifest").c_str()); unlink((impl_->directory + "/manifest.ready").c_str()); }
    fprintf(stderr, "[ram-prepare rank=%d] input_bytes=%llu capacity=%zu seconds=%.6f (outside MPI interval)\n",
            impl_->rank, (unsigned long long)impl_->staged.body->size(), pool()->capacity(), Now() - start);
    return 0;
}
int Session::EnterMpi() {
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    const int matches = !ram() || (rank == impl_->rank && ranks == impl_->ranks);
    if (!swbam::mpi::AllRanksOk(matches)) {
        impl_->invalid_layout = true;
        fprintf(stderr, "ERROR: preloaded rank layout does not match MPI communicator\n"); return -1;
    }
    impl_->rank = rank; impl_->ranks = ranks;
    return swbam::mpi::ResolveWorkMemory(impl_->capacity, impl_->requested_work,
                                         impl_->node_allowance, &work_memory);
}
void Session::LeaveMpi() { impl_->bam.reset(); impl_->sam.reset(); temporary.clear(); }
swbam::BamInputBackend *Session::bam_input() {
    if (!ram() || impl_->staged.manifest.kind != 1) return nullptr;
    if (!impl_->bam) impl_->bam.reset(new RamBamInput(&impl_->staged));
    return impl_->bam.get();
}
swbam::SamInputBackend *Session::sam_input() {
    if (!ram() || impl_->staged.manifest.kind != 2) return nullptr;
    if (!impl_->sam) impl_->sam.reset(new RamSamInput(&impl_->staged));
    return impl_->sam.get();
}
int Session::input_format() const { return impl_->staged.manifest.kind == 1 ? bam : sam; }
int InputFormat(const std::string &path) {
    if (IsRam()) return current->input_format();
    samFile *input = sam_open(path.c_str(), "r");
    if (!input) return -1;
    const int format = input->format.format;
    return sam_close(input) ? -1 : format;
}
int OpenInput(swbam::mpi::MpiBamInput &input, const std::string &path,
              const std::string &backend, const std::string &limit) {
    if (!IsRam()) return input.Open(path, backend, limit);
    swbam::BamInputBackend *b = current->bam_input();
    const int ok = b && b->Open(path) == 0 && input.Attach(b, true) == 0;
    return swbam::mpi::AllRanksOk(ok) ? 0 : -1;
}
int BodySink::Open(const std::string &backend, uint64_t limit, size_t reserve,
                   const std::string &prefix, const std::string &directory) {
    Close();
    if (IsRam()) {
        try { file_.reset(new RamFile(current->pool())); }
        catch (...) { errno = ENOMEM; return -1; }
        return 0;
    }
    const char *tmp = getenv("TMPDIR");
    return native_.Open(backend, limit, reserve, prefix,
                        directory.empty() ? (tmp ? tmp : "/tmp") : directory);
}
int BodySink::Close() { file_.reset(); return native_.Close(); }
int BodySink::Append(const void *data, size_t bytes) {
    const int result = file_ ? file_->WriteAt(file_->size(), data, bytes) : native_.Append(data, bytes);
    if (file_ && result)
        fprintf(stderr, "ERROR: RAM output append failed: %s; capacity=%zu used=%zu output=%llu\n",
                strerror(errno), current->pool()->capacity(), current->pool()->used(),
                (unsigned long long)file_->size());
    return result;
}
int BodySink::Flush() { return file_ ? 0 : native_.Flush(); }
int BodySink::ReadAt(uint64_t offset, void *data, size_t bytes) const {
    return file_ ? file_->ReadAt(offset, data, bytes) : native_.ReadAt(offset, data, bytes);
}
uint64_t BodySink::size() const { return file_ ? file_->size() : native_.size(); }
bool BodySink::is_memory() const { return file_ || native_.is_memory(); }
int BodySink::ReleaseBefore(uint64_t offset) const {
    return file_ ? file_->ReleaseBefore(offset) : native_.ReleaseBefore(offset);
}
int BodySequence::Add(BodySink *body) {
    if (!body || native_.Add(body)) return -1;
    if (IsRam()) {
        if (!body->file()) return -1;
        try { files_.push_back(body->file()); }
        catch (...) { errno = ENOMEM; return -1; }
    }
    return 0;
}

int Session::RetainOutput(const std::string &path, const swbam::RankBodySource &body,
                          uint64_t prefix_size, const void *prefix,
                          size_t suffix_size, const void *suffix) {
    swbam::mpi::DistributedBamOutput layout;
    if (layout.Prepare(body, prefix_size, suffix_size)) return -1;
    int ok = impl_->final_parts.empty();
    try { if (ok) {
        impl_->final_parts.push_back(std::make_shared<RamFile>(pool()));
        impl_->envelope.reset(new RamFile(pool()));
        impl_->output = path;
        impl_->prefix_bytes = layout.layout().prefix_size;
        impl_->suffix_bytes = layout.layout().suffix_size;
        impl_->body_offset = impl_->prefix_bytes + layout.local_body_offset();
        impl_->suffix_offset = impl_->prefix_bytes + layout.layout().total_body_size;
        impl_->output_size = layout.layout().total_size;
    } } catch (...) { errno = ENOMEM; ok = 0; }
    if (!swbam::mpi::AllRanksOk(ok)) return -1;
    class RamOutput : public swbam::mpi::PositionedOutput {
    public:
        explicit RamOutput(Impl *impl, uint64_t bytes) : p_(impl), bytes_(bytes) {}
        int WriteAt(uint64_t offset, const void *data, uint64_t bytes) {
            if (bytes > SIZE_MAX) return -1;
            if (p_->rank == 0 && offset == 0 && bytes == p_->prefix_bytes && bytes)
                return p_->envelope->WriteAt(0, data, bytes);
            if (p_->rank == 0 && offset == p_->suffix_offset && bytes == p_->suffix_bytes)
                return p_->envelope->WriteAt(p_->prefix_bytes, data, bytes);
            if (offset < p_->body_offset || offset - p_->body_offset > bytes_ ||
                bytes > bytes_ - (offset - p_->body_offset)) return -1;
            return p_->final_parts[0]->WriteAt(offset - p_->body_offset, data, bytes);
        }
        int Sync() { return 0; }
    private:
        Impl *p_;
        uint64_t bytes_;
    } output(impl_.get(), body.size());
    const double start = Now();
    const int result = layout.Write(&output, body, prefix, suffix, true);
    fprintf(stdout, "[final-output rank=%d] body_read_bytes=%llu body_write_bytes=%llu seconds=%.6f\n",
            impl_->rank, (unsigned long long)body.size(),
            (unsigned long long)impl_->final_parts[0]->size(), Now() - start);
    return result;
}
int WriteOutput(const std::string &path, const swbam::RankBodySource &body,
                size_t prefix_size, const void *prefix,
                size_t suffix_size, const void *suffix) {
    if (IsRam()) return current->RetainOutput(path, body, prefix_size, prefix, suffix_size, suffix);
    swbam::mpi::DistributedBamOutput output;
    return output.Prepare(body, prefix_size, suffix_size) ||
           output.WriteMpiIo(path, body, prefix, suffix, true) ? -1 : 0;
}
int BuildBamPrefix(sam_hdr_t *header, int level, swbam::MemoryBamOutput *prefix) {
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
int Session::Finish(int status, double seconds, double processing_seconds) {
    if (PhysicalDataIoCalls()) status = 1;
#ifdef SWBAM_RAM_IO_AUDIT
    fprintf(stdout, "[ram-io-audit rank=%d] physical_data_io_calls=%d\n", impl_->rank, PhysicalDataIoCalls());
#endif
    // MPI is finalized: the cached rank and control-file rendezvous are intentional.
    if (processing_seconds >= 0) {
        fprintf(stdout, "[swbam-processing rank=%d] storage=%s seconds=%.6f status=%d\n",
                impl_->rank, ram() ? "ram" : "posix", processing_seconds, status);
    }
    fprintf(stdout, "[swbam-e2e rank=%d] storage=%s seconds=%.6f status=%d\n",
            impl_->rank, ram() ? "ram" : "posix", seconds, status);
    fflush(stdout);
    if (!ram() || impl_->invalid_layout) return status;
    fprintf(stdout, "[ram-storage rank=%d] capacity=%zu peak=%zu temp_read_bytes=%llu temp_write_bytes=%llu output_bytes=%llu\n",
            impl_->rank, pool()->capacity(), pool()->peak(),
            (unsigned long long)temp_read_bytes, (unsigned long long)temp_write_bytes,
            (unsigned long long)impl_->output_size);
    if (status) Mark(impl_->Marker("error", impl_->rank));
    if (Mark(impl_->Marker("finished", impl_->rank))) return 1;
    for (int r = 0; r < impl_->ranks; ++r)
        if (impl_->Wait(impl_->Marker("finished", r))) return 1;
    for (int r = 0; r < impl_->ranks; ++r)
        if (Exists(impl_->Marker("error", r))) return 1;
    if (impl_->no_dump || impl_->final_parts.empty()) return status;
    int fd = -1;
    if (impl_->rank == 0) {
        fd = open(impl_->output.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0666);
        if (fd < 0 || ftruncate(fd, impl_->output_size)) status = 1;
        if (fd >= 0) close(fd);
        if (!status && Mark(impl_->directory + "/output.ready")) status = 1;
        if (status) Mark(impl_->Marker("error", impl_->rank));
    }
    if (impl_->Wait(impl_->directory + "/output.ready")) return 1;
    fd = open(impl_->output.c_str(), O_WRONLY);
    std::vector<char> buffer(8u * 1024u * 1024u);
    if (fd < 0) status = 1;
    auto dump = [&](RamFile &file, uint64_t from, uint64_t to, uint64_t bytes) {
        for (uint64_t pos = 0; !status && pos < bytes;) {
            const size_t count = std::min<uint64_t>(buffer.size(), bytes - pos);
            if (file.ReadAt(from + pos, buffer.data(), count) || Transfer(fd, to + pos, buffer.data(), count, true)) status = 1;
            pos += count;
        }
    };
    uint64_t offset = impl_->body_offset;
    for (const std::shared_ptr<RamFile> &part : impl_->final_parts) {
        dump(*part, 0, offset, part->size());
        offset += part->size();
    }
    if (impl_->rank == 0) {
        dump(*impl_->envelope, 0, 0, impl_->prefix_bytes);
        dump(*impl_->envelope, impl_->prefix_bytes, impl_->suffix_offset, impl_->suffix_bytes);
    }
    if (fd >= 0) { if (fsync(fd)) status = 1; if (close(fd)) status = 1; }
    if (status) fprintf(stderr, "[rank %d] ERROR: post-MPI output save failed for %s: %s\n",
                        impl_->rank, impl_->output.c_str(), strerror(errno));
    if (status) Mark(impl_->Marker("error", impl_->rank));
    if (Mark(impl_->Marker("dumped", impl_->rank))) return 1;
    for (int r = 0; r < impl_->ranks; ++r)
        if (impl_->Wait(impl_->Marker("dumped", r))) return 1;
    for (int r = 0; r < impl_->ranks; ++r)
        if (Exists(impl_->Marker("error", r))) return 1;
    if (impl_->rank == 0) printf("RAM output saved after MPI_Finalize: %s bytes=%llu\n", impl_->output.c_str(), (unsigned long long)impl_->output_size);
    return status;
}
int OpenTemporary(const std::string &path) {
    if (IsRam()) {
        const int handle = next_handle++;
        try { temporary[handle].reset(new RamFile(current->pool())); }
        catch (...) { temporary.erase(handle); errno = ENOMEM; return -1; }
        return handle;
    }
    int fd = open(path.c_str(), O_RDWR | O_CREAT | O_EXCL, 0600);
    if (fd >= 0 && unlink(path.c_str())) { close(fd); return -1; }
    return fd;
}
int CloseTemporary(int handle) {
    auto it = temporary.find(handle);
    if (it == temporary.end()) {
        if (IsRam()) { errno = EBADF; return -1; }
        return close(handle);
    }
    temporary.erase(it); return 0;
}
ssize_t ReadTemporary(int handle, void *data, size_t bytes, off_t offset) {
    auto it = temporary.find(handle);
    if (it == temporary.end()) {
        if (IsRam()) { errno = EBADF; return -1; }
        return pread(handle, data, bytes, offset);
    }
    if (offset < 0) { errno = EINVAL; return -1; }
    if ((uint64_t)offset >= it->second->size()) return 0;
    bytes = std::min<uint64_t>(bytes, it->second->size() - offset);
    if (it->second->ReadAt(offset, data, bytes)) return -1;
    temp_read_bytes += bytes; return bytes;
}
ssize_t WriteTemporary(int handle, const void *data, size_t bytes, off_t offset) {
    auto it = temporary.find(handle);
    if (it == temporary.end()) {
        if (IsRam()) { errno = EBADF; return -1; }
        return pwrite(handle, data, bytes, offset);
    }
    if (offset < 0) { errno = EINVAL; return -1; }
    if (it->second->WriteAt(offset, data, bytes)) {
        const int error = errno;
        fprintf(stderr, "ERROR: RAM temporary write failed: %s; capacity=%zu used=%zu\n",
                strerror(error), current->pool()->capacity(), current->pool()->used());
        fflush(stderr);
        // Peers may already be waiting in the next exchange; never leave them stuck.
        if (error == ENOSPC || error == ENOMEM) MPI_Abort(MPI_COMM_WORLD, error);
        errno = error; return -1;
    }
    temp_write_bytes += bytes; return bytes;
}
}
