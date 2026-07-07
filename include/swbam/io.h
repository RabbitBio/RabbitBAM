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
