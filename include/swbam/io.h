#ifndef SWBAM_IO_H
#define SWBAM_IO_H

#include "BamTools.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace swbam {

struct BgzfBlockSpan {
    uint64_t offset;
    uint32_t compressed_size;
};

class BgzfBlockBatch {
public:
    BgzfBlockBatch();
    ~BgzfBlockBatch();

    int Allocate(size_t capacity);
    int Attach(bam_block *blocks, size_t capacity);
    void Release();

    bam_block *blocks() { return blocks_; }
    const bam_block *blocks() const { return blocks_; }
    size_t capacity() const { return capacity_; }

private:
    BgzfBlockBatch(const BgzfBlockBatch &);
    BgzfBlockBatch &operator=(const BgzfBlockBatch &);

    bam_block *blocks_;
    unsigned char *data_;
    size_t capacity_;
    bool owns_storage_;
};

class BamInputBackend {
public:
    virtual ~BamInputBackend() {}

    virtual int Open(const std::string &path) = 0;
    virtual void Close() = 0;
    virtual const sam_hdr_t *header() const = 0;
    virtual int format() const = 0;
    virtual size_t size() const = 0;
    virtual uint64_t body_offset() const = 0;
    virtual int ScanBlocks(std::vector<BgzfBlockSpan> *spans) const = 0;
    virtual int ReadBatch(const BgzfBlockSpan *spans, size_t count,
                          BgzfBlockBatch *batch) const = 0;
};

class MemoryBamInput : public BamInputBackend {
public:
    MemoryBamInput();
    ~MemoryBamInput();

    int Open(const std::string &path);
    int OpenMemoryCopy(const void *data, size_t size);
    int Load(const std::string &path);
    int LoadMemoryCopy(const void *data, size_t size);
    int ParseHeader();
    void Close();
    const sam_hdr_t *header() const { return header_; }
    int format() const { return format_; }
    size_t size() const { return size_; }
    uint64_t body_offset() const { return body_offset_; }
    int ScanBlocks(std::vector<BgzfBlockSpan> *spans) const;
    int ReadBatch(const BgzfBlockSpan *spans, size_t count,
                  BgzfBlockBatch *batch) const;

    char *data() { return data_; }
    const char *data() const { return data_; }

private:
    MemoryBamInput(const MemoryBamInput &);
    MemoryBamInput &operator=(const MemoryBamInput &);

    char *data_;
    size_t size_;
    samFile *input_;
    hFILE *pending_hfile_;
    sam_hdr_t *header_;
    int format_;
    uint64_t body_offset_;
};

class PosixBamInput : public BamInputBackend {
public:
    PosixBamInput();
    ~PosixBamInput();

    int Open(const std::string &path);
    void Close();
    const sam_hdr_t *header() const;
    int format() const;
    size_t size() const;
    uint64_t body_offset() const;
    int ScanBlocks(std::vector<BgzfBlockSpan> *spans) const;
    int ReadBatch(const BgzfBlockSpan *spans, size_t count,
                  BgzfBlockBatch *batch) const;

private:
    class Impl;
    Impl *impl_;

    PosixBamInput(const PosixBamInput &);
    PosixBamInput &operator=(const PosixBamInput &);
};

class SamInputBackend {
public:
    virtual ~SamInputBackend() {}
    virtual int Open(const std::string &path) = 0;
    virtual void Close() = 0;
    virtual const sam_hdr_t *header() const = 0;
    virtual size_t size() const = 0;
    virtual uint64_t body_offset() const = 0;
    virtual int SplitRanges(int range_count,
                            std::vector<long long> *offsets,
                            std::vector<long long> *lengths) const = 0;
    virtual int ReadRange(uint64_t offset, size_t length,
                          std::vector<char> *data) const = 0;
};

class PosixSamInput : public SamInputBackend {
public:
    PosixSamInput();
    ~PosixSamInput();

    int Open(const std::string &path);
    void Close();
    const sam_hdr_t *header() const;
    size_t size() const;
    uint64_t body_offset() const;
    int SplitRanges(int range_count,
                    std::vector<long long> *offsets,
                    std::vector<long long> *lengths) const;
    int ReadRange(uint64_t offset, size_t length,
                  std::vector<char> *data) const;

private:
    class Impl;
    Impl *impl_;

    PosixSamInput(const PosixSamInput &);
    PosixSamInput &operator=(const PosixSamInput &);
};

class BgzfSpanBatchReader {
public:
    BgzfSpanBatchReader(const BamInputBackend *input,
                        const BgzfBlockSpan *spans, size_t count);

    int ReadNext(BgzfBlockBatch *batch, size_t *count);
    size_t remaining() const { return span_count_ - position_; }

private:
    const BamInputBackend *input_;
    const BgzfBlockSpan *spans_;
    size_t span_count_;
    size_t position_;
};

class BamOutputBackend {
public:
    virtual ~BamOutputBackend() {}

    virtual int Write(const void *data, size_t size) = 0;
    virtual int Flush() = 0;
    virtual uint64_t bytes_written() const = 0;
};

int WriteBgzfBlocks(BamOutputBackend *output,
                    const bam_block *blocks, size_t count);
int WriteBgzfEof(BamOutputBackend *output);

class MemoryBamOutput : public BamOutputBackend {
public:
    MemoryBamOutput();

    void Clear();
    int Reserve(size_t capacity);
    int Write(const void *data, size_t size);
    int Flush() { return 0; }
    uint64_t bytes_written() const { return data_.size(); }

    // Compatibility helpers for existing memory-based stages.
    int Append(const void *data, size_t size);
    int AppendBlocks(const bam_block *blocks, size_t count);
    int AppendBgzfEof();

    const unsigned char *data() const {
        return data_.empty() ? nullptr : data_.data();
    }
    size_t size() const { return data_.size(); }

private:
    std::vector<unsigned char> data_;
};

class PosixBamOutput : public BamOutputBackend {
public:
    PosixBamOutput();
    ~PosixBamOutput();

    int Open(const std::string &path);
    int Close();
    int Write(const void *data, size_t size);
    int Flush();
    uint64_t bytes_written() const;

private:
    class Impl;
    Impl *impl_;

    PosixBamOutput(const PosixBamOutput &);
    PosixBamOutput &operator=(const PosixBamOutput &);
};

class RankBodySource {
public:
    virtual ~RankBodySource() {}

    virtual uint64_t size() const = 0;
    virtual int ReadAt(uint64_t offset, void *data,
                       size_t size) const = 0;
    virtual bool is_memory() const = 0;
};

class RankBodySink : public RankBodySource {
public:
    virtual ~RankBodySink() {}

    virtual int Append(const void *data, size_t size) = 0;
    virtual int Flush() = 0;
};

class MemoryRankBodySink : public RankBodySink {
public:
    explicit MemoryRankBodySink(MemWriter *writer);

    int Append(const void *data, size_t size);
    int Flush() { return 0; }
    uint64_t size() const;
    int ReadAt(uint64_t offset, void *data, size_t size) const;
    bool is_memory() const { return true; }

    MemWriter *writer() { return writer_; }
    const MemWriter *writer() const { return writer_; }

private:
    MemWriter *writer_;
};

class SpoolRankBodySink : public RankBodySink {
public:
    SpoolRankBodySink();
    ~SpoolRankBodySink();

    int Open(const std::string &path, bool remove_on_close = true);
    int OpenTemporary(const std::string &directory,
                      const std::string &prefix);
    int Close();
    int Append(const void *data, size_t size);
    int Flush();
    uint64_t size() const;
    int ReadAt(uint64_t offset, void *data, size_t size) const;
    bool is_memory() const { return false; }
    const std::string &path() const { return path_; }

private:
    SpoolRankBodySink(const SpoolRankBodySink &);
    SpoolRankBodySink &operator=(const SpoolRankBodySink &);

    int fd_;
    uint64_t size_;
    std::string path_;
    bool remove_on_close_;
};

class AdaptiveRankBodySink : public RankBodySink {
public:
    AdaptiveRankBodySink();
    ~AdaptiveRankBodySink();

    int Open(const std::string &requested_backend,
             uint64_t memory_limit, size_t initial_capacity,
             const std::string &spool_prefix,
             const std::string &spool_directory);
    int Close();
    int Append(const void *data, size_t size);
    int Flush();
    uint64_t size() const;
    int ReadAt(uint64_t offset, void *data, size_t size) const;
    bool is_memory() const;

    const std::string &requested_backend() const {
        return requested_backend_;
    }
    const std::string &selected_backend() const {
        return selected_backend_;
    }
    MemWriter *memory_writer();
    const MemWriter *memory_writer() const;

private:
    AdaptiveRankBodySink(const AdaptiveRankBodySink &);
    AdaptiveRankBodySink &operator=(const AdaptiveRankBodySink &);

    int OpenSpool();
    int SpillToSpool();
    int AppendToSpool(const void *data, size_t size);

    MemWriter memory_;
    MemoryRankBodySink memory_sink_;
    SpoolRankBodySink spool_sink_;
    std::string requested_backend_;
    std::string selected_backend_;
    std::string spool_prefix_;
    std::string spool_directory_;
    uint64_t memory_limit_;
    bool opened_;
};

class SegmentedRankBodySource : public RankBodySource {
public:
    SegmentedRankBodySource();

    void Clear();
    int Add(const RankBodySource *source);
    uint64_t size() const { return total_size_; }
    int ReadAt(uint64_t offset, void *data, size_t size) const;
    bool is_memory() const;
    size_t segment_count() const { return segments_.size(); }

private:
    std::vector<const RankBodySource *> segments_;
    std::vector<uint64_t> offsets_;
    uint64_t total_size_;
};

int AppendBgzfBlock(RankBodySink *sink, const bam_block *block);
int ParseByteSize(const std::string &text, uint64_t *value);

int LoadFileToMemory(const std::string &path, char **data, size_t *size);
int ScanBgzfBlocksInMemory(const char *base, size_t size, long long body_start,
                           std::vector<long long> *offsets,
                           std::vector<long long> *lengths);
int SelectBlockRangeFromMemory(char *base, size_t input_size,
                               const std::vector<long long> &offsets,
                               const std::vector<long long> &lengths,
                               long long begin, long long end,
                               char **data, size_t *size);
int FindSamBodyStartInMemory(const char *base, size_t size,
                             long long *body_start);
int SplitSamRangesInMemory(const char *base, size_t size,
                           long long body_start, int range_count,
                           std::vector<long long> *offsets,
                           std::vector<long long> *lengths);

} // namespace swbam

#endif
