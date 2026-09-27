#include "swbam/io.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace swbam {

namespace {

int EnsureMemoryCapacity(MemWriter *writer, size_t add) {
    if (!writer || add > SIZE_MAX - writer->size) return -1;
    const size_t required = writer->size + add;
    if (required <= writer->capacity) return 0;
    size_t next = writer->capacity
        ? writer->capacity : 64u * 1024u * 1024u;
    while (required > next) {
        if (next > SIZE_MAX / 2) return -1;
        next *= 2;
    }
    char *data = static_cast<char *>(realloc(writer->data, next));
    if (!data) return -1;
    writer->data = data;
    writer->capacity = next;
    return 0;
}

int WriteAll(int fd, const unsigned char *data, size_t size) {
    size_t done = 0;
    while (done < size) {
        const size_t request = std::min(
            size - done, static_cast<size_t>(INT_MAX));
        const ssize_t written = write(fd, data + done, request);
        if (written < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        if (written == 0) return -1;
        done += static_cast<size_t>(written);
    }
    return 0;
}

int ReadAllAt(int fd, uint64_t offset,
              unsigned char *data, size_t size) {
    size_t done = 0;
    while (done < size) {
        const size_t request = std::min(
            size - done, static_cast<size_t>(INT_MAX));
        if (offset + done > static_cast<uint64_t>(
                std::numeric_limits<off_t>::max())) {
            return -1;
        }
        const ssize_t count = pread(
            fd, data + done, request,
            static_cast<off_t>(offset + done));
        if (count < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        if (count == 0) return -1;
        done += static_cast<size_t>(count);
    }
    return 0;
}

} // namespace

MemoryRankBodySink::MemoryRankBodySink(MemWriter *writer)
    : writer_(writer) {}

int MemoryRankBodySink::Append(const void *data, size_t size) {
    if (!writer_ || (size > 0 && !data) ||
        EnsureMemoryCapacity(writer_, size) != 0) {
        return -1;
    }
    if (size > 0) {
        memcpy(writer_->data + writer_->size, data, size);
        writer_->size += size;
    }
    return 0;
}

uint64_t MemoryRankBodySink::size() const {
    return writer_ ? static_cast<uint64_t>(writer_->size) : 0;
}

int MemoryRankBodySink::ReadAt(
        uint64_t offset, void *data, size_t size) const {
    if (!writer_ || (size > 0 && !data) ||
        offset > writer_->size || size > writer_->size -
            static_cast<size_t>(offset)) {
        return -1;
    }
    if (size > 0) {
        memcpy(data, writer_->data + static_cast<size_t>(offset), size);
    }
    return 0;
}

SpoolRankBodySink::SpoolRankBodySink()
    : fd_(-1), size_(0), remove_on_close_(true) {}

SpoolRankBodySink::~SpoolRankBodySink() {
    Close();
}

int SpoolRankBodySink::Open(
        const std::string &path, bool remove_on_close) {
    if (path.empty()) return -1;
    Close();
    const int fd = open(
        path.c_str(), O_CREAT | O_TRUNC | O_RDWR, 0600);
    if (fd < 0) return -1;
    fd_ = fd;
    size_ = 0;
    path_ = path;
    remove_on_close_ = remove_on_close;
    return 0;
}

int SpoolRankBodySink::OpenTemporary(
        const std::string &directory,
        const std::string &prefix) {
    if (directory.empty() || prefix.empty()) return -1;
    Close();

    std::string pattern = directory;
    if (pattern[pattern.size() - 1] != '/') pattern += '/';
    pattern += prefix;
    pattern += ".XXXXXX";
    std::vector<char> writable(pattern.begin(), pattern.end());
    writable.push_back('\0');
    const int fd = mkstemp(writable.data());
    if (fd < 0) return -1;
    if (unlink(writable.data()) != 0) {
        const int saved_errno = errno;
        close(fd);
        errno = saved_errno;
        return -1;
    }

    fd_ = fd;
    size_ = 0;
    path_.assign(writable.data());
    remove_on_close_ = false;
    return 0;
}

int SpoolRankBodySink::Close() {
    int ret = 0;
    if (fd_ >= 0 && close(fd_) != 0) ret = -1;
    fd_ = -1;
    if (remove_on_close_ && !path_.empty() &&
        unlink(path_.c_str()) != 0 && errno != ENOENT) {
        ret = -1;
    }
    path_.clear();
    size_ = 0;
    remove_on_close_ = true;
    return ret;
}

int SpoolRankBodySink::Append(const void *data, size_t size) {
    if (fd_ < 0 || (size > 0 && !data) ||
        size > UINT64_MAX - size_) {
        return -1;
    }
    if (size > 0 && WriteAll(
            fd_, static_cast<const unsigned char *>(data), size) != 0) {
        return -1;
    }
    size_ += static_cast<uint64_t>(size);
    return 0;
}

int SpoolRankBodySink::Flush() {
    return fd_ >= 0 ? 0 : -1;
}

uint64_t SpoolRankBodySink::size() const {
    return size_;
}

int SpoolRankBodySink::ReadAt(
        uint64_t offset, void *data, size_t size) const {
    if (fd_ < 0 || (size > 0 && !data) ||
        offset > size_ || size > size_ - offset) {
        return -1;
    }
    return size > 0
        ? ReadAllAt(fd_, offset,
                    static_cast<unsigned char *>(data), size) : 0;
}

AdaptiveRankBodySink::AdaptiveRankBodySink()
    : memory_sink_(&memory_), memory_limit_(0), opened_(false) {
    memset(&memory_, 0, sizeof(memory_));
}

AdaptiveRankBodySink::~AdaptiveRankBodySink() {
    Close();
}

int AdaptiveRankBodySink::Open(
        const std::string &requested_backend,
        uint64_t memory_limit, size_t initial_capacity,
        const std::string &spool_prefix,
        const std::string &spool_directory) {
    Close();
    if (requested_backend != "memory" && spool_directory.empty()) {
        fprintf(stderr,
                "ERROR: --rank-body-temp-dir is required when rank "
                "body backend is '%s'.\n",
                requested_backend.c_str());
        return -1;
    }
    if ((requested_backend != "memory" &&
         requested_backend != "spool" &&
         requested_backend != "auto") ||
        spool_prefix.empty() ||
        (requested_backend == "auto" && memory_limit == 0)) {
        return -1;
    }
    requested_backend_ = requested_backend;
    selected_backend_ = requested_backend == "spool"
        ? "spool" : "memory";
    spool_prefix_ = spool_prefix;
    spool_directory_ = spool_directory;
    memory_limit_ = memory_limit;
    opened_ = true;

    if (selected_backend_ == "spool") {
        if (OpenSpool() != 0) {
            Close();
            return -1;
        }
        return 0;
    }

    size_t capacity = initial_capacity
        ? initial_capacity : 64u * 1024u * 1024u;
    if (requested_backend_ == "auto" &&
        capacity > memory_limit_) {
        capacity = memory_limit_ > SIZE_MAX
            ? SIZE_MAX : static_cast<size_t>(memory_limit_);
    }
    if (capacity == 0) capacity = 1;
    memory_.data = static_cast<char *>(malloc(capacity));
    if (!memory_.data) {
        Close();
        return -1;
    }
    memory_.capacity = capacity;
    memory_.size = 0;
    return 0;
}

int AdaptiveRankBodySink::OpenSpool() {
    if (spool_sink_.OpenTemporary(
            spool_directory_, spool_prefix_) != 0) {
        fprintf(stderr,
                "ERROR: cannot create rank body spool in '%s': %s.\n",
                spool_directory_.c_str(), strerror(errno));
        return -1;
    }
    return 0;
}

int AdaptiveRankBodySink::SpillToSpool() {
    if (selected_backend_ == "spool") return 0;
    if (OpenSpool() != 0) return -1;
    if (memory_.size > 0 &&
        AppendToSpool(memory_.data, memory_.size) != 0) {
        spool_sink_.Close();
        return -1;
    }
    free(memory_.data);
    memset(&memory_, 0, sizeof(memory_));
    selected_backend_ = "spool";
    return 0;
}

int AdaptiveRankBodySink::AppendToSpool(
        const void *data, size_t size) {
    if (spool_sink_.Append(data, size) == 0) return 0;
    fprintf(stderr,
            "ERROR: rank body spool write failed in '%s': %s; "
            "check --rank-body-temp-dir capacity.\n",
            spool_directory_.c_str(), strerror(errno));
    return -1;
}

int AdaptiveRankBodySink::Close() {
    int ret = spool_sink_.Close();
    free(memory_.data);
    memset(&memory_, 0, sizeof(memory_));
    requested_backend_.clear();
    selected_backend_.clear();
    spool_prefix_.clear();
    spool_directory_.clear();
    memory_limit_ = 0;
    opened_ = false;
    return ret;
}

int AdaptiveRankBodySink::Append(const void *data, size_t size) {
    if (!opened_ || (size > 0 && !data)) return -1;
    if (requested_backend_ == "auto" &&
        selected_backend_ == "memory") {
        const uint64_t current = memory_sink_.size();
        if (size > UINT64_MAX - current ||
            current + static_cast<uint64_t>(size) > memory_limit_) {
            if (SpillToSpool() != 0) return -1;
        }
    }
    return selected_backend_ == "memory"
        ? memory_sink_.Append(data, size)
        : AppendToSpool(data, size);
}

int AdaptiveRankBodySink::Flush() {
    if (!opened_) return -1;
    return selected_backend_ == "memory"
        ? memory_sink_.Flush() : spool_sink_.Flush();
}

uint64_t AdaptiveRankBodySink::size() const {
    if (!opened_) return 0;
    return selected_backend_ == "memory"
        ? memory_sink_.size() : spool_sink_.size();
}

int AdaptiveRankBodySink::ReadAt(
        uint64_t offset, void *data, size_t size) const {
    if (!opened_) return -1;
    return selected_backend_ == "memory"
        ? memory_sink_.ReadAt(offset, data, size)
        : spool_sink_.ReadAt(offset, data, size);
}

bool AdaptiveRankBodySink::is_memory() const {
    return opened_ && selected_backend_ == "memory";
}

MemWriter *AdaptiveRankBodySink::memory_writer() {
    return is_memory() ? &memory_ : nullptr;
}

const MemWriter *AdaptiveRankBodySink::memory_writer() const {
    return is_memory() ? &memory_ : nullptr;
}

SegmentedRankBodySource::SegmentedRankBodySource()
    : total_size_(0) {}

void SegmentedRankBodySource::Clear() {
    segments_.clear();
    offsets_.clear();
    total_size_ = 0;
}

int SegmentedRankBodySource::Add(const RankBodySource *source) {
    if (!source || source->size() > UINT64_MAX - total_size_) {
        return -1;
    }
    offsets_.push_back(total_size_);
    segments_.push_back(source);
    total_size_ += source->size();
    return 0;
}

int SegmentedRankBodySource::ReadAt(
        uint64_t offset, void *data, size_t size) const {
    if ((size > 0 && !data) || offset > total_size_ ||
        static_cast<uint64_t>(size) > total_size_ - offset) {
        return -1;
    }
    unsigned char *output = static_cast<unsigned char *>(data);
    size_t remaining = size;
    uint64_t position = offset;
    for (size_t i = 0; remaining > 0 && i < segments_.size(); ++i) {
        const uint64_t begin = offsets_[i];
        const uint64_t end = begin + segments_[i]->size();
        if (position >= end) continue;
        if (position < begin) return -1;
        const uint64_t available = end - position;
        const size_t count = static_cast<size_t>(std::min<uint64_t>(
            available, static_cast<uint64_t>(remaining)));
        if (segments_[i]->ReadAt(
                position - begin, output, count) != 0) {
            return -1;
        }
        output += count;
        remaining -= count;
        position += count;
    }
    return remaining == 0 ? 0 : -1;
}

bool SegmentedRankBodySource::is_memory() const {
    for (size_t i = 0; i < segments_.size(); ++i) {
        if (!segments_[i]->is_memory()) return false;
    }
    return true;
}

int AppendBgzfBlock(RankBodySink *sink, const bam_block *block) {
    if (!sink || !block || !block->data ||
        block->length == 0 || block->length > BGZF_MAX_BLOCK_SIZE) {
        return -1;
    }
    return sink->Append(block->data, block->length);
}

int ParseByteSize(const std::string &text, uint64_t *value) {
    if (!value || text.empty()) return -1;
    errno = 0;
    char *end = nullptr;
    const unsigned long long base =
        strtoull(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str()) return -1;
    while (*end == ' ' || *end == '\t') ++end;
    uint64_t multiplier = 1;
    if (*end != '\0') {
        const char suffix = *end++;
        if (suffix == 'K' || suffix == 'k') multiplier = 1024ULL;
        else if (suffix == 'M' || suffix == 'm') {
            multiplier = 1024ULL * 1024ULL;
        } else if (suffix == 'G' || suffix == 'g') {
            multiplier = 1024ULL * 1024ULL * 1024ULL;
        } else if (suffix == 'T' || suffix == 't') {
            multiplier = 1024ULL * 1024ULL * 1024ULL * 1024ULL;
        } else {
            return -1;
        }
        if (*end == 'B' || *end == 'b') ++end;
        while (*end == ' ' || *end == '\t') ++end;
        if (*end != '\0') return -1;
    }
    if (base > UINT64_MAX / multiplier) return -1;
    *value = static_cast<uint64_t>(base) * multiplier;
    return 0;
}

} // namespace swbam
