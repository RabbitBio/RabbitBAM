#ifndef SWBAM_USAGE_STORAGE_H
#define SWBAM_USAGE_STORAGE_H

#include "ram_file.h"
#include "swbam/mpi_runtime.h"

#include <memory>
#include <sys/types.h>

namespace swbam_usage {

// One session owns staged storage until after MPI_Finalize. No SDK dependency
// points back to this usage-side module.
class Session {
public:
    Session();
    ~Session();
    int Prepare(int *argc, char ***argv, bool sdk_example = false);
    int EnterMpi();
    void LeaveMpi();
    int Finish(int status, double seconds);
    bool ram() const;
    const std::string &input_path() const;
    const std::string &output_path() const;
    swbam::BamInputBackend *bam_input();
    swbam::SamInputBackend *sam_input();
    int input_format() const;
    RamPool *pool();
    int RetainOutput(const std::string &path,
                     const swbam::RankBodySource &body,
                     uint64_t prefix_size, const void *prefix,
                     size_t suffix_size, const void *suffix);
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

Session *CurrentSession();
bool IsRam();
double WallTime();
int InputFormat(const std::string &path);
#ifdef SWBAM_RAM_IO_AUDIT
void BeginMpiInterval();
void EndMpiInterval();
int PhysicalDataIoCalls();
#else
inline void BeginMpiInterval() {}
inline void EndMpiInterval() {}
inline int PhysicalDataIoCalls() { return 0; }
#endif
int OpenInput(swbam::mpi::MpiBamInput &input, const std::string &path,
              const std::string &backend, const std::string &limit);
int WriteOutput(const std::string &path, const swbam::RankBodySource &body,
                size_t prefix_size, const void *prefix,
                size_t suffix_size, const void *suffix);
int BuildBamPrefix(sam_hdr_t *header, int level, swbam::MemoryBamOutput *prefix);

class BodySink : public swbam::RankBodySink {
public:
    int Open(const std::string &backend, uint64_t limit, size_t reserve,
             const std::string &prefix, const std::string &directory);
    int Close();
    int Append(const void *data, size_t size);
    int Flush();
    int ReadAt(uint64_t offset, void *data, size_t size) const;
    uint64_t size() const;
    bool is_memory() const;
    MemWriter *memory_writer() { return file_ ? nullptr : native_.memory_writer(); }
    const MemWriter *memory_writer() const { return file_ ? nullptr : native_.memory_writer(); }
    std::shared_ptr<RamFile> file() const { return file_; }
private:
    swbam::AdaptiveRankBodySink native_;
    std::shared_ptr<RamFile> file_;
};

class BodySequence : public swbam::RankBodySource {
public:
    int Add(BodySink *body);
    int ReadAt(uint64_t offset, void *data, size_t size) const { return native_.ReadAt(offset, data, size); }
    uint64_t size() const { return native_.size(); }
    bool is_memory() const { return native_.is_memory(); }
    const std::vector<std::shared_ptr<RamFile> > &files() const { return files_; }
private:
    swbam::SegmentedRankBodySource native_;
    std::vector<std::shared_ptr<RamFile> > files_;
};

int OpenTemporary(const std::string &path);
int CloseTemporary(int handle);
ssize_t ReadTemporary(int handle, void *data, size_t size, off_t offset);
ssize_t WriteTemporary(int handle, const void *data, size_t size, off_t offset);
}
#endif
