#include "swbam/mpi_runtime.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace swbam {
namespace mpi {
namespace {

int ParseSize(const std::string &text, uint64_t *value) {
    if (!value || text.empty()) return -1;
    errno = 0;
    char *end = nullptr;
    const unsigned long long base = strtoull(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str()) return -1;
    uint64_t multiplier = 1;
    if (*end != '\0') {
        const char suffix = *end++;
        if (*end != '\0') return -1;
        if (suffix == 'K' || suffix == 'k') multiplier = 1024ULL;
        else if (suffix == 'M' || suffix == 'm') multiplier = 1024ULL * 1024ULL;
        else if (suffix == 'G' || suffix == 'g') multiplier = 1024ULL * 1024ULL * 1024ULL;
        else if (suffix == 'T' || suffix == 't') multiplier = 1024ULL * 1024ULL * 1024ULL * 1024ULL;
        else return -1;
    }
    if (base > UINT64_MAX / multiplier) return -1;
    *value = (uint64_t)base * multiplier;
    return 0;
}

int SelectBackend(const std::string &path,
                  const std::string &requested,
                  const std::string &memory_limit,
                  int comm_size,
                  int *selected,
                  uint64_t *budget) {
    if (!selected || !budget) return -1;
    *budget = 0;
    if (requested == "memory") {
        *selected = 0;
        return 0;
    }
    if (requested == "posix") {
        *selected = 1;
        return 0;
    }
    if (requested == "mpiio") {
        *selected = 2;
        return 0;
    }
    if (requested != "auto") return -1;

    uint64_t configured = 0;
    if (ParseSize(memory_limit, &configured) != 0 || configured == 0) {
        return -1;
    }
    uint64_t available_per_rank = UINT64_MAX;
    const long pages = sysconf(_SC_AVPHYS_PAGES);
    const long page_size = sysconf(_SC_PAGESIZE);
    if (pages > 0 && page_size > 0 && comm_size > 0) {
        const unsigned long long available =
            (unsigned long long)pages * (unsigned long long)page_size;
        available_per_rank = (uint64_t)(available / 2 /
            (unsigned long long)comm_size);
    }
    *budget = configured < available_per_rank
        ? configured : available_per_rank;

    struct stat file_stat;
    if (stat(path.c_str(), &file_stat) != 0 || file_stat.st_size < 0) {
        return -1;
    }
    *selected = (uint64_t)file_stat.st_size <= *budget ? 0 : 1;
    return 0;
}

} // namespace

MpiIoBamInput::MpiIoBamInput(MPI_Comm comm)
    : comm_(comm), file_(MPI_FILE_NULL) {}

MpiIoBamInput::~MpiIoBamInput() {
    Close();
}

int MpiIoBamInput::Open(const std::string &path) {
    Close();
    int local_ok = metadata_.Open(path) == 0 ? 1 : 0;
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, comm_);
    if (!global_ok) {
        metadata_.Close();
        return -1;
    }

    int ret = MPI_File_open(
        comm_, const_cast<char *>(path.c_str()), MPI_MODE_RDONLY,
        MPI_INFO_NULL, &file_);
    if (ret != MPI_SUCCESS) {
        file_ = MPI_FILE_NULL;
        metadata_.Close();
        return -1;
    }
    MPI_File_set_errhandler(file_, MPI_ERRORS_RETURN);
    return 0;
}

void MpiIoBamInput::Close() {
    if (file_ != MPI_FILE_NULL) {
        MPI_File_close(&file_);
        file_ = MPI_FILE_NULL;
    }
    metadata_.Close();
}

int MpiIoBamInput::ScanBlocks(
        std::vector<BgzfBlockSpan> *spans) const {
    return metadata_.ScanBlocks(spans);
}

int MpiIoBamInput::ReadBatch(
        const BgzfBlockSpan *spans, size_t count,
        BgzfBlockBatch *batch) const {
    if (file_ == MPI_FILE_NULL || !batch || count > batch->capacity() ||
        (count > 0 && !spans) || count > static_cast<size_t>(INT_MAX)) {
        return -1;
    }
    if (count == 0) return 0;
    bam_block *blocks = batch->blocks();
    if (!blocks) return -1;

    bool contiguous = true;
    uint64_t next_offset = spans[0].offset;
    std::vector<int> lengths(count);
    std::vector<MPI_Aint> displacements(count);
    MPI_Aint base_address = 0;
    MPI_Get_address(blocks[0].data, &base_address);
    for (size_t i = 0; i < count; ++i) {
        const uint64_t offset = spans[i].offset;
        const size_t length = spans[i].compressed_size;
        if (length == 0 || length > BGZF_MAX_BLOCK_SIZE ||
            offset > static_cast<uint64_t>(
                std::numeric_limits<MPI_Offset>::max()) ||
            offset > static_cast<uint64_t>(metadata_.size()) ||
            length > metadata_.size() - static_cast<size_t>(offset)) {
            return -1;
        }
        if (i > 0 && offset != next_offset) contiguous = false;
        next_offset = offset + static_cast<uint64_t>(length);
        lengths[i] = static_cast<int>(length);
        MPI_Aint address = 0;
        MPI_Get_address(blocks[i].data, &address);
        displacements[i] = address - base_address;
    }

    if (contiguous) {
        MPI_Datatype memory_type = MPI_DATATYPE_NULL;
        if (MPI_Type_create_hindexed(
                static_cast<int>(count), lengths.data(),
                displacements.data(), MPI_BYTE,
                &memory_type) != MPI_SUCCESS) {
            return -1;
        }
        int ret = MPI_Type_commit(&memory_type);
        MPI_Status status;
        if (ret == MPI_SUCCESS) {
            ret = MPI_File_read_at(
                file_, static_cast<MPI_Offset>(spans[0].offset),
                blocks[0].data, 1, memory_type, &status);
        }
        int actual = 0;
        if (ret == MPI_SUCCESS &&
            (MPI_Get_count(&status, memory_type, &actual) != MPI_SUCCESS ||
             actual != 1)) {
            ret = MPI_ERR_IO;
        }
        MPI_Type_free(&memory_type);
        if (ret != MPI_SUCCESS) return -1;
    } else {
        for (size_t i = 0; i < count; ++i) {
            MPI_Status status;
            int ret = MPI_File_read_at(
                file_, static_cast<MPI_Offset>(spans[i].offset),
                blocks[i].data, lengths[i], MPI_BYTE, &status);
            int actual = 0;
            if (ret != MPI_SUCCESS ||
                MPI_Get_count(&status, MPI_BYTE, &actual) != MPI_SUCCESS ||
                actual != lengths[i]) {
                return -1;
            }
        }
    }

    for (size_t i = 0; i < count; ++i) {
        blocks[i].length = static_cast<unsigned int>(lengths[i]);
        blocks[i].pos = 0;
        blocks[i].errcode = 0;
        blocks[i].block_address = static_cast<int64_t>(spans[i].offset);
        blocks[i].block_id = static_cast<int>(i);
    }
    return 0;
}

MpiIoSamInput::MpiIoSamInput(MPI_Comm comm)
    : comm_(comm), file_(MPI_FILE_NULL) {}

MpiIoSamInput::~MpiIoSamInput() {
    Close();
}

int MpiIoSamInput::Open(const std::string &path) {
    Close();
    int local_ok = metadata_.Open(path) == 0 ? 1 : 0;
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, comm_);
    if (!global_ok) {
        metadata_.Close();
        return -1;
    }
    int ret = MPI_File_open(
        comm_, const_cast<char *>(path.c_str()), MPI_MODE_RDONLY,
        MPI_INFO_NULL, &file_);
    if (ret != MPI_SUCCESS) {
        file_ = MPI_FILE_NULL;
        metadata_.Close();
        return -1;
    }
    MPI_File_set_errhandler(file_, MPI_ERRORS_RETURN);
    return 0;
}

void MpiIoSamInput::Close() {
    if (file_ != MPI_FILE_NULL) {
        MPI_File_close(&file_);
        file_ = MPI_FILE_NULL;
    }
    metadata_.Close();
}

int MpiIoSamInput::SplitRanges(
        int range_count, std::vector<long long> *offsets,
        std::vector<long long> *lengths) const {
    return metadata_.SplitRanges(range_count, offsets, lengths);
}

int MpiIoSamInput::ReadRange(uint64_t offset, size_t length,
                             std::vector<char> *data) const {
    if (file_ == MPI_FILE_NULL || !data || offset > metadata_.size() ||
        length > metadata_.size() - static_cast<size_t>(offset) ||
        offset > static_cast<uint64_t>(
            std::numeric_limits<MPI_Offset>::max()) ||
        static_cast<uint64_t>(length) > UINT64_MAX - offset ||
        offset + static_cast<uint64_t>(length) > static_cast<uint64_t>(
            std::numeric_limits<MPI_Offset>::max())) {
        return -1;
    }
    try {
        data->resize(length);
    } catch (...) {
        return -1;
    }
    uint64_t done = 0;
    while (done < length) {
        const int chunk = static_cast<int>(std::min<uint64_t>(
            static_cast<uint64_t>(length) - done, INT_MAX));
        MPI_Status status;
        int ret = MPI_File_read_at(
            file_, static_cast<MPI_Offset>(offset + done),
            data->data() + static_cast<size_t>(done),
            chunk, MPI_BYTE, &status);
        int actual = 0;
        if (ret != MPI_SUCCESS ||
            MPI_Get_count(&status, MPI_BYTE, &actual) != MPI_SUCCESS ||
            actual != chunk) {
            return -1;
        }
        done += static_cast<uint64_t>(chunk);
    }
    return 0;
}

MpiBamInput::MpiBamInput()
    : backend_(nullptr), data_open_cost_(0.0), header_open_cost_(0.0),
      auto_memory_budget_(0) {}

MpiBamInput::~MpiBamInput() {
    Close();
}

int MpiBamInput::Open(const std::string &path,
                      const std::string &requested_backend,
                      const std::string &memory_limit,
                      MPI_Comm comm) {
    Close();
    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    int selected = -1;
    int selection_ok = 1;
    unsigned long long budget = 0;
    if (rank == 0) {
        uint64_t selected_budget = 0;
        if (SelectBackend(
                path, requested_backend, memory_limit, comm_size,
                &selected, &selected_budget) != 0) {
            selection_ok = 0;
        }
        budget = (unsigned long long)selected_budget;
    }
    MPI_Bcast(&selection_ok, 1, MPI_INT, 0, comm);
    if (!selection_ok) return -1;
    MPI_Bcast(&selected, 1, MPI_INT, 0, comm);
    MPI_Bcast(&budget, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
    auto_memory_budget_ = (uint64_t)budget;

    int local_ok = 1;
    if (selected == 0) {
        MemoryBamInput *memory = new MemoryBamInput();
        const double data_t0 = GetTime();
        if (memory->Load(path) != 0) local_ok = 0;
        data_open_cost_ = GetTime() - data_t0;
        const double header_t0 = GetTime();
        if (local_ok && memory->ParseHeader() != 0) local_ok = 0;
        header_open_cost_ = GetTime() - header_t0;
        backend_ = memory;
        selected_backend_ = "memory";
    } else if (selected == 1) {
        PosixBamInput *posix = new PosixBamInput();
        const double data_t0 = GetTime();
        if (posix->Open(path) != 0) local_ok = 0;
        data_open_cost_ = GetTime() - data_t0;
        header_open_cost_ = 0.0;
        backend_ = posix;
        selected_backend_ = "posix";
    } else if (selected == 2) {
        MpiIoBamInput *mpiio = new MpiIoBamInput(comm);
        const double data_t0 = GetTime();
        if (mpiio->Open(path) != 0) local_ok = 0;
        data_open_cost_ = GetTime() - data_t0;
        header_open_cost_ = 0.0;
        backend_ = mpiio;
        selected_backend_ = "mpiio";
    } else {
        local_ok = 0;
    }

    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, comm);
    if (!global_ok && requested_backend == "auto" && selected == 0) {
        Close();
        auto_memory_budget_ = (uint64_t)budget;
        PosixBamInput *posix = new PosixBamInput();
        const double data_t0 = GetTime();
        local_ok = posix->Open(path) == 0 ? 1 : 0;
        data_open_cost_ = GetTime() - data_t0;
        header_open_cost_ = 0.0;
        backend_ = posix;
        selected_backend_ = "posix";
        MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, comm);
    }
    if (!global_ok) {
        Close();
        return -1;
    }
    return 0;
}

void MpiBamInput::Close() {
    if (backend_) {
        backend_->Close();
        delete backend_;
        backend_ = nullptr;
    }
    selected_backend_.clear();
    data_open_cost_ = 0.0;
    header_open_cost_ = 0.0;
    auto_memory_budget_ = 0;
}

MpiFileOutput::MpiFileOutput() : file_(MPI_FILE_NULL) {}

MpiFileOutput::~MpiFileOutput() {
    Close();
}

int MpiFileOutput::Open(const std::string &path, uint64_t final_size,
                        MPI_Comm comm) {
    if (file_ != MPI_FILE_NULL || path.empty() ||
        final_size > static_cast<uint64_t>(
            std::numeric_limits<MPI_Offset>::max())) {
        return -1;
    }
    int ret = MPI_File_open(
        comm, const_cast<char *>(path.c_str()),
        MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_);
    if (ret != MPI_SUCCESS) {
        file_ = MPI_FILE_NULL;
        return -1;
    }
    ret = MPI_File_set_size(file_, static_cast<MPI_Offset>(final_size));
    if (ret != MPI_SUCCESS) {
        MPI_File_close(&file_);
        file_ = MPI_FILE_NULL;
        return -1;
    }
    return 0;
}

int MpiFileOutput::WriteAt(uint64_t offset, const void *data,
                           uint64_t size) {
    if (file_ == MPI_FILE_NULL || (size > 0 && !data) ||
        offset > static_cast<uint64_t>(
            std::numeric_limits<MPI_Offset>::max()) ||
        size > UINT64_MAX - offset ||
        offset + size > static_cast<uint64_t>(
            std::numeric_limits<MPI_Offset>::max())) {
        return -1;
    }
    const unsigned char *bytes =
        static_cast<const unsigned char *>(data);
    uint64_t written = 0;
    while (written < size) {
        const int chunk = static_cast<int>(
            std::min<uint64_t>(size - written, INT_MAX));
        MPI_Status status;
        int ret = MPI_File_write_at(
            file_, static_cast<MPI_Offset>(offset + written),
            const_cast<unsigned char *>(bytes + written),
            chunk, MPI_BYTE, &status);
        if (ret != MPI_SUCCESS) return -1;
        int actual = 0;
        if (MPI_Get_count(&status, MPI_BYTE, &actual) != MPI_SUCCESS ||
            actual != chunk) {
            return -1;
        }
        written += static_cast<uint64_t>(chunk);
    }
    return 0;
}

int MpiFileOutput::Sync() {
    return file_ != MPI_FILE_NULL && MPI_File_sync(file_) == MPI_SUCCESS
        ? 0 : -1;
}

int MpiFileOutput::Close() {
    if (file_ == MPI_FILE_NULL) return 0;
    int ret = MPI_File_close(&file_);
    file_ = MPI_FILE_NULL;
    return ret == MPI_SUCCESS ? 0 : -1;
}

DistributedBamOutput::DistributedBamOutput()
    : comm_(MPI_COMM_NULL), rank_(0), comm_size_(0), root_(0),
      prepared_(false) {}

int DistributedBamOutput::Prepare(
        const RankBodySource &local_body,
        uint64_t prefix_size, uint64_t suffix_size,
        int root, MPI_Comm comm) {
    if (PrepareBody(local_body, root, comm) != 0) return -1;
    return SetEnvelope(prefix_size, suffix_size);
}

int DistributedBamOutput::PrepareBody(
        const RankBodySource &local_body,
        int root, MPI_Comm comm) {
    prepared_ = false;
    layout_ = DistributedBamLayout();
    comm_ = comm;
    root_ = root;
    if (comm_ == MPI_COMM_NULL ||
        MPI_Comm_rank(comm_, &rank_) != MPI_SUCCESS ||
        MPI_Comm_size(comm_, &comm_size_) != MPI_SUCCESS ||
        root_ < 0 || root_ >= comm_size_) {
        return -1;
    }

    unsigned long long local_size =
        static_cast<unsigned long long>(local_body.size());
    std::vector<unsigned long long> gathered(
        static_cast<size_t>(comm_size_), 0);
    if (MPI_Allgather(
            &local_size, 1,
            MPI_UNSIGNED_LONG_LONG, gathered.data(), 1,
            MPI_UNSIGNED_LONG_LONG, comm_) != MPI_SUCCESS) {
        return -1;
    }

    layout_.body_sizes.resize(static_cast<size_t>(comm_size_));
    layout_.body_offsets.resize(static_cast<size_t>(comm_size_));
    uint64_t total = 0;
    int local_ok = 1;
    for (int i = 0; i < comm_size_; ++i) {
        const uint64_t size = static_cast<uint64_t>(
            gathered[static_cast<size_t>(i)]);
        layout_.body_sizes[static_cast<size_t>(i)] = size;
        layout_.body_offsets[static_cast<size_t>(i)] = total;
        if (size > UINT64_MAX - total) {
            local_ok = 0;
            break;
        }
        total += size;
    }
    layout_.total_body_size = total;
    if (!AllRanksOk(local_ok, comm_)) return -1;
    return 0;
}

int DistributedBamOutput::SetEnvelope(
        uint64_t prefix_size, uint64_t suffix_size) {
    if (comm_ == MPI_COMM_NULL || comm_size_ <= 0) return -1;
    unsigned long long root_sizes[2] = {0, 0};
    if (rank_ == root_) {
        root_sizes[0] = static_cast<unsigned long long>(prefix_size);
        root_sizes[1] = static_cast<unsigned long long>(suffix_size);
    }
    if (MPI_Bcast(root_sizes, 2, MPI_UNSIGNED_LONG_LONG,
                  root_, comm_) != MPI_SUCCESS) {
        return -1;
    }
    layout_.prefix_size = static_cast<uint64_t>(root_sizes[0]);
    layout_.suffix_size = static_cast<uint64_t>(root_sizes[1]);
    int local_ok = 1;
    if (layout_.prefix_size >
            UINT64_MAX - layout_.total_body_size ||
        layout_.prefix_size + layout_.total_body_size >
            UINT64_MAX - layout_.suffix_size) {
        local_ok = 0;
    } else {
        layout_.total_size = layout_.prefix_size +
            layout_.total_body_size + layout_.suffix_size;
    }
    if (!AllRanksOk(local_ok, comm_)) return -1;
    prepared_ = true;
    return 0;
}

uint64_t DistributedBamOutput::local_body_offset() const {
    if (!prepared_ || rank_ < 0 || rank_ >= comm_size_) return 0;
    return layout_.body_offsets[static_cast<size_t>(rank_)];
}

int DistributedBamOutput::GatherToRoot(
        const RankBodySource &local_body,
        const void *prefix, const void *suffix,
        char **root_data, size_t *root_size, int mpi_tag) const {
    if (!prepared_ || !root_data || !root_size || mpi_tag < 0 ||
        local_body.size() !=
            layout_.body_sizes[static_cast<size_t>(rank_)]) {
        return -1;
    }
    *root_data = nullptr;
    *root_size = 0;
    int local_ok = 1;
    if (rank_ == root_) {
        if (layout_.total_size > static_cast<uint64_t>(SIZE_MAX) ||
            (layout_.prefix_size > 0 && !prefix) ||
            (layout_.suffix_size > 0 && !suffix)) {
            local_ok = 0;
        } else {
            *root_size = static_cast<size_t>(layout_.total_size);
            *root_data = *root_size
                ? static_cast<char *>(malloc(*root_size)) : nullptr;
            if (*root_size && !*root_data) local_ok = 0;
        }
    }

    const size_t chunk_capacity = 8u * 1024u * 1024u;
    std::vector<unsigned char> chunk;
    if (rank_ != root_ && local_body.size() > 0) {
        try {
            chunk.resize(static_cast<size_t>(std::min<uint64_t>(
                local_body.size(), chunk_capacity)));
        } catch (...) {
            local_ok = 0;
        }
    }
    if (!AllRanksOk(local_ok, comm_)) {
        if (rank_ == root_) {
            free(*root_data);
            *root_data = nullptr;
            *root_size = 0;
        }
        return -1;
    }

    if (rank_ == root_) {
        if (layout_.prefix_size > 0) {
            memcpy(*root_data, prefix,
                   static_cast<size_t>(layout_.prefix_size));
        }
        const uint64_t local_offset = layout_.prefix_size +
            layout_.body_offsets[static_cast<size_t>(root_)];
        if (local_body.size() > 0 &&
            local_body.ReadAt(
                0, *root_data + static_cast<size_t>(local_offset),
                static_cast<size_t>(local_body.size())) != 0) {
            local_ok = 0;
        }
        for (int src = 0; local_ok && src < comm_size_; ++src) {
            if (src == root_) continue;
            uint64_t received = 0;
            const uint64_t size =
                layout_.body_sizes[static_cast<size_t>(src)];
            while (received < size) {
                const int count = static_cast<int>(std::min<uint64_t>(
                    size - received, chunk_capacity));
                MPI_Status status;
                const uint64_t output_offset = layout_.prefix_size +
                    layout_.body_offsets[static_cast<size_t>(src)] +
                    received;
                if (MPI_Recv(
                        *root_data + static_cast<size_t>(output_offset),
                        count, MPI_BYTE, src, mpi_tag, comm_,
                        &status) != MPI_SUCCESS) {
                    local_ok = 0;
                    break;
                }
                int actual = 0;
                if (MPI_Get_count(&status, MPI_BYTE, &actual) !=
                        MPI_SUCCESS || actual != count) {
                    local_ok = 0;
                    break;
                }
                received += static_cast<uint64_t>(count);
            }
        }
        if (local_ok && layout_.suffix_size > 0) {
            memcpy(*root_data + static_cast<size_t>(
                       layout_.prefix_size + layout_.total_body_size),
                   suffix, static_cast<size_t>(layout_.suffix_size));
        }
    } else {
        uint64_t sent = 0;
        while (local_ok && sent < local_body.size()) {
            const size_t count = static_cast<size_t>(std::min<uint64_t>(
                local_body.size() - sent, chunk.size()));
            if (local_body.ReadAt(sent, chunk.data(), count) != 0 ||
                MPI_Send(chunk.data(), static_cast<int>(count), MPI_BYTE,
                         root_, mpi_tag, comm_) != MPI_SUCCESS) {
                local_ok = 0;
                break;
            }
            sent += count;
        }
    }
    if (!AllRanksOk(local_ok, comm_)) {
        if (rank_ == root_) {
            free(*root_data);
            *root_data = nullptr;
            *root_size = 0;
        }
        return -1;
    }
    return 0;
}

int DistributedBamOutput::WriteMpiIo(
        const std::string &path,
        const RankBodySource &local_body,
        const void *prefix, const void *suffix) const {
    if (!prepared_ || path.empty() ||
        local_body.size() !=
            layout_.body_sizes[static_cast<size_t>(rank_)] ||
        (rank_ == root_ && layout_.prefix_size > 0 && !prefix) ||
        (rank_ == root_ && layout_.suffix_size > 0 && !suffix)) {
        return -1;
    }
    int local_ok = 1;
    MpiFileOutput output;
    if (output.Open(path, layout_.total_size, comm_) != 0) local_ok = 0;
    if (!AllRanksOk(local_ok, comm_)) {
        output.Close();
        return -1;
    }

    if (rank_ == root_ && layout_.prefix_size > 0 &&
        output.WriteAt(0, prefix, layout_.prefix_size) != 0) {
        local_ok = 0;
    }
    const size_t chunk_capacity = 8u * 1024u * 1024u;
    std::vector<unsigned char> chunk;
    if (local_body.size() > 0) {
        try {
            chunk.resize(static_cast<size_t>(std::min<uint64_t>(
                local_body.size(), chunk_capacity)));
        } catch (...) {
            local_ok = 0;
        }
    }
    uint64_t copied = 0;
    while (local_ok && copied < local_body.size()) {
        const size_t count = static_cast<size_t>(std::min<uint64_t>(
            local_body.size() - copied, chunk.size()));
        if (local_body.ReadAt(copied, chunk.data(), count) != 0 ||
            output.WriteAt(
                layout_.prefix_size + local_body_offset() + copied,
                chunk.data(), count) != 0) {
            local_ok = 0;
            break;
        }
        copied += count;
    }
    if (rank_ == root_ && layout_.suffix_size > 0 &&
        output.WriteAt(
            layout_.prefix_size + layout_.total_body_size,
            suffix, layout_.suffix_size) != 0) {
        local_ok = 0;
    }
    const int writes_ok = AllRanksOk(local_ok, comm_);
    if (writes_ok && output.Sync() != 0) local_ok = 0;
    if (output.Close() != 0) local_ok = 0;
    return AllRanksOk(local_ok, comm_) ? 0 : -1;
}

int AllRanksOk(int local_ok, MPI_Comm comm) {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, comm);
    return global_ok;
}

double ReduceMaxCost(double local_cost, int root, MPI_Comm comm) {
    double max_cost = 0.0;
    MPI_Reduce(&local_cost, &max_cost, 1, MPI_DOUBLE, MPI_MAX, root, comm);
    return max_cost;
}

int PrepareMpiBamInputPlan(const BamInputBackend &input,
                           MpiBamInputPlan *plan,
                           int root, MPI_Comm comm) {
    if (!plan) return -1;

    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);
    if (root < 0 || root >= comm_size) return -1;

    int local_ok = 1;
    unsigned long long body_offset = 0;
    long long block_count = 0;
    if (rank == root) {
        plan->blocks.clear();
        if (input.ScanBlocks(&plan->blocks) != 0 ||
            plan->blocks.size() > static_cast<size_t>(INT_MAX)) {
            local_ok = 0;
        } else {
            body_offset = static_cast<unsigned long long>(input.body_offset());
            block_count = static_cast<long long>(plan->blocks.size());
        }
    }

    MPI_Bcast(&local_ok, 1, MPI_INT, root, comm);
    if (!local_ok) return -1;
    MPI_Bcast(&body_offset, 1, MPI_UNSIGNED_LONG_LONG, root, comm);
    MPI_Bcast(&block_count, 1, MPI_LONG_LONG, root, comm);
    if (block_count < 0 || block_count > INT_MAX) return -1;

    if (rank != root) plan->blocks.resize(static_cast<size_t>(block_count));
    std::vector<unsigned long long> offsets(static_cast<size_t>(block_count));
    std::vector<unsigned int> lengths(static_cast<size_t>(block_count));
    if (rank == root) {
        for (size_t i = 0; i < plan->blocks.size(); ++i) {
            offsets[i] = static_cast<unsigned long long>(plan->blocks[i].offset);
            lengths[i] = static_cast<unsigned int>(
                plan->blocks[i].compressed_size);
        }
    }
    if (block_count > 0) {
        MPI_Bcast(offsets.data(), static_cast<int>(block_count),
                  MPI_UNSIGNED_LONG_LONG, root, comm);
        MPI_Bcast(lengths.data(), static_cast<int>(block_count),
                  MPI_UNSIGNED, root, comm);
    }
    if (rank != root) {
        for (size_t i = 0; i < plan->blocks.size(); ++i) {
            plan->blocks[i].offset = static_cast<uint64_t>(offsets[i]);
            plan->blocks[i].compressed_size = static_cast<uint32_t>(lengths[i]);
        }
    }

    plan->body_offset = static_cast<uint64_t>(body_offset);
    size_t n = plan->blocks.size();
    plan->rank_begin = n * static_cast<size_t>(rank) /
                       static_cast<size_t>(comm_size);
    plan->rank_end = n * static_cast<size_t>(rank + 1) /
                     static_cast<size_t>(comm_size);
    return 0;
}

} // namespace mpi
} // namespace swbam
