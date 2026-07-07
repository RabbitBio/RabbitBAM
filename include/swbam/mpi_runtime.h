#ifndef SWBAM_MPI_RUNTIME_H
#define SWBAM_MPI_RUNTIME_H

#include "swbam/io.h"

#include <mpi.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace swbam {
namespace mpi {

struct MpiBamInputPlan {
    std::vector<BgzfBlockSpan> blocks;
    size_t rank_begin;
    size_t rank_end;
    uint64_t body_offset;

    MpiBamInputPlan()
        : rank_begin(0), rank_end(0), body_offset(0) {}

    const BgzfBlockSpan *rank_spans() const {
        return rank_begin < rank_end ? blocks.data() + rank_begin : nullptr;
    }
    size_t rank_block_count() const { return rank_end - rank_begin; }
};

class MpiIoBamInput : public BamInputBackend {
public:
    explicit MpiIoBamInput(MPI_Comm comm = MPI_COMM_WORLD);
    ~MpiIoBamInput();

    int Open(const std::string &path);
    void Close();
    const sam_hdr_t *header() const { return metadata_.header(); }
    int format() const { return metadata_.format(); }
    size_t size() const { return metadata_.size(); }
    uint64_t body_offset() const { return metadata_.body_offset(); }
    int ScanBlocks(std::vector<BgzfBlockSpan> *spans) const;
    int ReadBatch(const BgzfBlockSpan *spans, size_t count,
                  BgzfBlockBatch *batch) const;

private:
    MpiIoBamInput(const MpiIoBamInput &);
    MpiIoBamInput &operator=(const MpiIoBamInput &);

    MPI_Comm comm_;
    MPI_File file_;
    PosixBamInput metadata_;
};

class MpiIoSamInput : public SamInputBackend {
public:
    explicit MpiIoSamInput(MPI_Comm comm = MPI_COMM_WORLD);
    ~MpiIoSamInput();

    int Open(const std::string &path);
    void Close();
    const sam_hdr_t *header() const { return metadata_.header(); }
    size_t size() const { return metadata_.size(); }
    uint64_t body_offset() const { return metadata_.body_offset(); }
    int SplitRanges(int range_count,
                    std::vector<long long> *offsets,
                    std::vector<long long> *lengths) const;
    int ReadRange(uint64_t offset, size_t length,
                  std::vector<char> *data) const;

private:
    MpiIoSamInput(const MpiIoSamInput &);
    MpiIoSamInput &operator=(const MpiIoSamInput &);

    MPI_Comm comm_;
    MPI_File file_;
    PosixSamInput metadata_;
};

class MpiBamInput {
public:
    MpiBamInput();
    ~MpiBamInput();

    int Open(const std::string &path,
             const std::string &requested_backend,
             const std::string &memory_limit,
             MPI_Comm comm = MPI_COMM_WORLD);
    void Close();

    BamInputBackend *backend() { return backend_; }
    const BamInputBackend *backend() const { return backend_; }
    const std::string &selected_backend() const { return selected_backend_; }
    double data_open_cost() const { return data_open_cost_; }
    double header_open_cost() const { return header_open_cost_; }
    uint64_t auto_memory_budget() const { return auto_memory_budget_; }

private:
    MpiBamInput(const MpiBamInput &);
    MpiBamInput &operator=(const MpiBamInput &);

    BamInputBackend *backend_;
    std::string selected_backend_;
    double data_open_cost_;
    double header_open_cost_;
    uint64_t auto_memory_budget_;
};

class MpiFileOutput {
public:
    MpiFileOutput();
    ~MpiFileOutput();

    int Open(const std::string &path, uint64_t final_size,
             MPI_Comm comm = MPI_COMM_WORLD);
    int WriteAt(uint64_t offset, const void *data, uint64_t size);
    int Sync();
    int Close();
    bool is_open() const { return file_ != MPI_FILE_NULL; }

private:
    MpiFileOutput(const MpiFileOutput &);
    MpiFileOutput &operator=(const MpiFileOutput &);

    MPI_File file_;
};

int AllRanksOk(int local_ok, MPI_Comm comm = MPI_COMM_WORLD);
double ReduceMaxCost(double local_cost, int root = 0,
                     MPI_Comm comm = MPI_COMM_WORLD);
int PrepareMpiBamInputPlan(const BamInputBackend &input,
                           MpiBamInputPlan *plan,
                           int root = 0,
                           MPI_Comm comm = MPI_COMM_WORLD);

} // namespace mpi
} // namespace swbam

#endif
