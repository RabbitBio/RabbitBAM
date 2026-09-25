#include "swbam/cpe_sam_write_pipeline.h"

#include <cstdio>
#include <cstring>
#include <stdint.h>
#include <utility>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

extern "C" void slave_sam_format();

namespace swbam {
namespace cpe {
namespace {

int AllocateBatch(SamFormatBatch *batch, const sam_hdr_t *header) {
    memset(batch, 0, sizeof(*batch));
    batch->hdr = header;
    for (int i = 0; i < 64; ++i) {
        kstring_t &line = batch->core_out_lines[i];
        line.m = MAX_SAM_FORMAT_CORE_BUFFER_SIZE;
        line.s = reinterpret_cast<char *>(
            aligned_alloc_custom(64, line.m));
        if (!line.s) return -1;
        batch->required_capacity[i] = line.m;
    }
    return 0;
}

void ReleaseBatch(SamFormatBatch *batch) {
    if (!batch) return;
    for (int i = 0; i < 64; ++i) {
        if (batch->core_out_lines[i].s) {
            aligned_free_custom(reinterpret_cast<unsigned char *>(
                batch->core_out_lines[i].s));
            batch->core_out_lines[i].s = nullptr;
        }
    }
    aligned_free_custom(reinterpret_cast<unsigned char *>(batch));
}

int EnsureCapacity(SamFormatBatch *batch, int core, size_t need) {
    kstring_t &line = batch->core_out_lines[core];
    if (need <= line.m) return 0;
    size_t capacity = line.m ? line.m : MAX_SAM_FORMAT_CORE_BUFFER_SIZE;
    while (capacity < need) {
        if (capacity > SIZE_MAX / 2) {
            capacity = need;
            break;
        }
        capacity *= 2;
    }
    char *data = reinterpret_cast<char *>(aligned_alloc_custom(64, capacity));
    if (!data) return -1;
    if (line.l) memcpy(data, line.s, line.l);
    aligned_free_custom(reinterpret_cast<unsigned char *>(line.s));
    line.s = data;
    line.m = capacity;
    batch->required_capacity[core] = capacity;
    return 0;
}

void ResetBatch(SamFormatBatch *batch, int count) {
    batch->count = count;
    for (int i = 0; i < 64; ++i) {
        batch->core_out_lines[i].l = 0;
        batch->status[i] = 0;
        batch->required_capacity[i] = batch->core_out_lines[i].m;
    }
}

int ReserveForRecords(SamFormatBatch *batch, int count) {
    for (int core = 0; core < 64; ++core) {
        const int base = count >> 6;
        const int remainder = count & 63;
        const int start = core < remainder ? core * (base + 1)
                                           : remainder + core * base;
        const int end = start + base + (core < remainder);
        size_t need = 1;
        for (int i = start; i < end; ++i) {
            const bam1_t *record = batch->bams[i];
            if (record) {
                need += static_cast<size_t>(record->l_data) * 8u +
                        static_cast<size_t>(record->core.l_qseq) * 2u + 1024u;
            }
        }
        if (EnsureCapacity(batch, core, need) != 0) return -1;
    }
    return 0;
}

} // namespace

CpeSamWriteSession::CpeSamWriteSession()
    : sink_(nullptr), current_(nullptr), pending_(nullptr),
      has_pending_(false), kernel_(DefaultSamWriteKernelSpec()) {}

SamWriteKernelSpec DefaultSamWriteKernelSpec() {
    return {reinterpret_cast<void *>(slave_sam_format), BATCH_SIZE};
}

CpeSamWriteSession::~CpeSamWriteSession() {
    Close();
}

void CpeSamWriteSession::Close() {
    ReleaseBatch(current_);
    ReleaseBatch(pending_);
    current_ = nullptr;
    pending_ = nullptr;
    sink_ = nullptr;
    has_pending_ = false;
}

int CpeSamWriteSession::Initialize(RankBodySink *sink,
                                   const sam_hdr_t *header) {
    return Initialize(sink, header, DefaultSamWriteKernelSpec());
}

int CpeSamWriteSession::Initialize(RankBodySink *sink,
                                   const sam_hdr_t *header,
                                   const SamWriteKernelSpec &kernel) {
    if (!sink || !header || current_ || pending_ || !kernel.format_entry ||
        !kernel.batch_capacity || kernel.batch_capacity > BATCH_SIZE) return -1;
    kernel_ = kernel;
    sink_ = sink;
    current_ = reinterpret_cast<SamFormatBatch *>(
        aligned_alloc_custom(64, sizeof(SamFormatBatch)));
    pending_ = reinterpret_cast<SamFormatBatch *>(
        aligned_alloc_custom(64, sizeof(SamFormatBatch)));
    if (current_) memset(current_, 0, sizeof(*current_));
    if (pending_) memset(pending_, 0, sizeof(*pending_));
    if (!current_ || !pending_) return -1;
    if (AllocateBatch(current_, header) != 0 ||
        AllocateBatch(pending_, header) != 0) return -1;
    return 0;
}

bam1_t **CpeSamWriteSession::RecordSlots() {
    return current_ ? current_->bams : nullptr;
}

int CpeSamWriteSession::SubmitRecords(bam1_t *const *records, size_t count) {
    if (!current_ || (count && !records) || count > capacity()) return -1;
    if (count && records != current_->bams)
        memcpy(current_->bams, records, count * sizeof(*records));
    return Submit(count);
}

int CpeSamWriteSession::FlushPending() {
    if (!has_pending_) return 0;
    const double t0 = GetTime();
    for (int i = 0; i < 64; ++i) {
        const kstring_t &line = pending_->core_out_lines[i];
        if (line.l && sink_->Append(line.s, line.l) != 0) {
            fprintf(stderr, "ERROR: failed to append SAM text to rank body sink.\n");
            return -1;
        }
    }
    timing_.write += GetTime() - t0;
    has_pending_ = false;
    return 0;
}

int CpeSamWriteSession::FlushPrevious() { return FlushPending(); }

int CpeSamWriteSession::Submit(size_t count) {
    return SubmitWithPrefetch(count, nullptr, nullptr);
}

int CpeSamWriteSession::SubmitWithPrefetch(
        size_t count, CpeBatchPrefetch prefetch, void *context) {
    if (!current_ || count > capacity()) return -1;
    if (!count) return prefetch ? prefetch(context) : 0;
    ResetBatch(current_, static_cast<int>(count));
    if (ReserveForRecords(current_, static_cast<int>(count)) != 0) return -1;
    const double t0 = GetTime();
#ifdef PLATFORM_SUNWAY
    __real_athread_spawn(kernel_.format_entry, current_, 1);
    const int prefetch_status = prefetch ? prefetch(context) : 0;
    athread_join();
    timing_.format += GetTime() - t0;
    if (prefetch_status != 0) return -1;
    if (FlushPending() != 0) return -1;
    bool retry = false;
    for (int i = 0; i < 64; ++i) {
        if (current_->status[i] == BOUNDS_LIMIT_MAX_SAM_FORMAT_CORE_BUFFER_SIZE) {
            if (EnsureCapacity(current_, i, current_->required_capacity[i]) != 0)
                return -1;
            retry = true;
        } else if (current_->status[i] != 0) {
            fprintf(stderr, "ERROR: SAM formatter core %d failed: %d.\n",
                    i, current_->status[i]);
            return -1;
        }
    }
    if (retry) {
        ResetBatch(current_, static_cast<int>(count));
        const double retry_t0 = GetTime();
        __real_athread_spawn(kernel_.format_entry, current_, 1);
        athread_join();
        timing_.format += GetTime() - retry_t0;
        for (int i = 0; i < 64; ++i) {
            if (current_->status[i] != 0) return -1;
        }
    }
#else
    return -1;
#endif
    std::swap(current_, pending_);
    has_pending_ = true;
    timing_.batches++;
    return 0;
}

int CpeSamWriteSession::Finish() { return FlushPending(); }

} // namespace cpe
} // namespace swbam
