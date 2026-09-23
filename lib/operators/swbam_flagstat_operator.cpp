#include "swbam/operators/flagstat.h"

#include "BamTools.h"

#include <cstdio>
#include <cstring>
#include <stdint.h>

extern "C" {
    void slave_mpi_flagstat_count();
}

namespace swbam {
namespace operators {

namespace {

const int kFlagstatBatchBlocks = 64;

#define SWBAM_FLAGSTAT_ID_MATCH(public_id, cpe_id)                         \
    static_assert(static_cast<int>(public_id) == static_cast<int>(cpe_id), \
                  "flagstat ID mismatch")

SWBAM_FLAGSTAT_ID_MATCH(kFlagstatTotal, RB_FLAGSTAT_TOTAL);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatPrimary, RB_FLAGSTAT_PRIMARY);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatSecondary, RB_FLAGSTAT_SECONDARY);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatSupplementary, RB_FLAGSTAT_SUPPLEMENTARY);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatDuplicates, RB_FLAGSTAT_DUPLICATES);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatPrimaryDuplicates,
                        RB_FLAGSTAT_PRIMARY_DUPLICATES);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatMapped, RB_FLAGSTAT_MAPPED);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatPrimaryMapped, RB_FLAGSTAT_PRIMARY_MAPPED);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatPaired, RB_FLAGSTAT_PAIRED);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatRead1, RB_FLAGSTAT_READ1);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatRead2, RB_FLAGSTAT_READ2);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatProperlyPaired,
                        RB_FLAGSTAT_PROPERLY_PAIRED);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatPairMapped, RB_FLAGSTAT_PAIR_MAPPED);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatSingletons, RB_FLAGSTAT_SINGLETONS);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatDiffChr, RB_FLAGSTAT_DIFF_CHR);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatDiffChrMapq5, RB_FLAGSTAT_DIFF_CHR_MAPQ5);
SWBAM_FLAGSTAT_ID_MATCH(kFlagstatCounterCount, RB_FLAGSTAT_COUNTER_COUNT);

#undef SWBAM_FLAGSTAT_ID_MATCH

void AccumulateDecompDetail(const MpiFlagstatCountPara *paras,
                            int active_blocks,
                            double kernel_wall,
                            FlagstatMetrics *metrics) {
    if (!metrics || active_blocks <= 0 || kernel_wall <= 0.0) return;

    const MpiFlagstatCountPara *critical = nullptr;
    uint64_t critical_total = 0;
    for (int b = 0; b < active_blocks; ++b) {
        if (paras[b].decomp_total_cycles >= critical_total) {
            critical_total = paras[b].decomp_total_cycles;
            critical = &paras[b];
        }
    }
    if (!critical || critical_total == 0) {
        metrics->decomp_other += kernel_wall;
        return;
    }

    const double scale = kernel_wall / static_cast<double>(critical_total);
    const double alloc = scale * critical->decomp_alloc_cycles;
    const double inflate = scale * critical->decomp_inflate_cycles;
    const double crc = scale * critical->decomp_crc_cycles;
    const double parse = scale * critical->decomp_parse_cycles;
    double other = kernel_wall - alloc - inflate - crc - parse;
    if (other < 0.0) other = 0.0;

    metrics->decomp_alloc += alloc;
    metrics->decomp_inflate += inflate;
    metrics->decomp_crc += crc;
    metrics->decomp_parse += parse;
    metrics->decomp_other += other;
}

void MergeSlice(const MpiFlagstatCountSlice &slice,
                FlagstatCounts *counts) {
    for (int i = 0; i < kFlagstatCounterCount; ++i) {
        counts->values[i][0] += slice.values[i][0];
        counts->values[i][1] += slice.values[i][1];
    }
}

} // namespace

FlagstatCounts::FlagstatCounts() {
    std::memset(values, 0, sizeof(values));
}

FlagstatMetrics::FlagstatMetrics()
    : input_blocks(0), batch_count(0), total_records(0),
      read(0.0), kernel(0.0), decomp_alloc(0.0),
      decomp_inflate(0.0), decomp_crc(0.0), decomp_parse(0.0),
      decomp_other(0.0), post_process(0.0), total(0.0) {}

class FlagstatOperator::Impl {
public:
    Impl(FlagstatCounts *counts_in, FlagstatMetrics *metrics_in)
        : counts(counts_in), metrics(metrics_in), count_slices(nullptr),
          scratch_data(nullptr) {
        std::memset(paras, 0, sizeof(paras));
    }

    FlagstatCounts *counts;
    FlagstatMetrics *metrics;
    MpiFlagstatCountPara paras[kFlagstatBatchBlocks];
    MpiFlagstatCountSlice *count_slices;
    unsigned char *scratch_data;
};

FlagstatOperator::FlagstatOperator(FlagstatCounts *counts,
                                   FlagstatMetrics *metrics)
    : impl_(new Impl(counts, metrics)) {}

FlagstatOperator::~FlagstatOperator() {
    Shutdown();
    delete impl_;
}

const char *FlagstatOperator::name() const {
    return "flagstat";
}

size_t FlagstatOperator::batch_capacity() const {
    return kFlagstatBatchBlocks;
}

int FlagstatOperator::Initialize() {
    if (!impl_ || !impl_->counts) return -1;
    Shutdown();
    *impl_->counts = FlagstatCounts();
    if (impl_->metrics) *impl_->metrics = FlagstatMetrics();
    impl_->count_slices =
        reinterpret_cast<MpiFlagstatCountSlice *>(aligned_alloc_custom(
            64, static_cast<size_t>(kFlagstatBatchBlocks) *
                    sizeof(MpiFlagstatCountSlice)));
    impl_->scratch_data = aligned_alloc_custom(
        64, static_cast<size_t>(kFlagstatBatchBlocks) *
                MPI_BAM_BLOCK_ARENA_SIZE);
    return impl_->count_slices && impl_->scratch_data ? 0 : -1;
}

void FlagstatOperator::Shutdown() {
    if (!impl_) return;
    if (impl_->count_slices) {
        aligned_free_custom(
            reinterpret_cast<unsigned char *>(impl_->count_slices));
        impl_->count_slices = nullptr;
    }
    if (impl_->scratch_data) {
        aligned_free_custom(impl_->scratch_data);
        impl_->scratch_data = nullptr;
    }
}

int FlagstatOperator::Prepare(const BgzfBlockBatch &compressed,
                              BgzfBlockBatch *decoded,
                              size_t active_blocks) {
    if (!impl_ || !impl_->count_slices || !impl_->scratch_data ||
        !decoded || active_blocks > kFlagstatBatchBlocks) {
        return -1;
    }
    std::memset(impl_->count_slices, 0,
                static_cast<size_t>(kFlagstatBatchBlocks) *
                    sizeof(MpiFlagstatCountSlice));
    for (int b = 0; b < kFlagstatBatchBlocks; ++b) {
        MpiFlagstatCountPara &para = impl_->paras[b];
        para.block_id = b;
        para.scratch_data = impl_->scratch_data +
            static_cast<size_t>(b) * MPI_BAM_BLOCK_ARENA_SIZE;
        para.scratch_capacity = MPI_BAM_BLOCK_ARENA_SIZE;
        para.counts = &impl_->count_slices[b];
        para.n_total_records = 0;
        para.record_index = 0;
        para.actual_value = 0;
        para.limit_value = 0;
        para.limit_id = BOUNDS_LIMIT_NONE;
        para.decomp_alloc_cycles = 0;
        para.decomp_inflate_cycles = 0;
        para.decomp_crc_cycles = 0;
        para.decomp_parse_cycles = 0;
        para.decomp_total_cycles = 0;
        if (static_cast<size_t>(b) < active_blocks) {
            para.input_block = const_cast<bam_block *>(
                &compressed.blocks()[b]);
            para.un_comp_block = &decoded->blocks()[b];
            para.status = 0;
        } else {
            para.input_block = nullptr;
            para.un_comp_block = nullptr;
            para.status = -1;
        }
    }
    return 0;
}

void *FlagstatOperator::kernel_entry() const {
    return reinterpret_cast<void *>(slave_mpi_flagstat_count);
}

void *FlagstatOperator::kernel_arguments() {
    return impl_ ? impl_->paras : nullptr;
}

void FlagstatOperator::ObserveKernel(double wall_seconds,
                                     size_t active_blocks) {
    if (!impl_) return;
    AccumulateDecompDetail(impl_->paras, static_cast<int>(active_blocks),
                           wall_seconds, impl_->metrics);
}

int FlagstatOperator::Validate(size_t active_blocks) const {
    if (!impl_ || active_blocks > kFlagstatBatchBlocks) return -1;
    for (size_t b = 0; b < active_blocks; ++b) {
        const MpiFlagstatCountPara &para = impl_->paras[b];
        if (para.status == 0) continue;
        if (para.status == -3) {
            std::fprintf(
                stderr,
                "ERROR: MPI flagstat capacity exceeded on input block %zu. "
                "limit_id=%d limit=%lld actual=%lld record=%d.\n",
                b, para.limit_id, para.limit_value,
                para.actual_value, para.record_index);
        } else {
            std::fprintf(
                stderr,
                "ERROR: MPI flagstat count failed on input block %zu "
                "with status %d.\n",
                b, para.status);
        }
        return -1;
    }
    return 0;
}

int FlagstatOperator::PostProcessBatch(size_t active_blocks,
                                       long long *records_processed) {
    if (!impl_ || !impl_->counts || active_blocks > kFlagstatBatchBlocks) {
        return -1;
    }
    long long records = 0;
    for (size_t b = 0; b < active_blocks; ++b) {
        records += impl_->paras[b].n_total_records;
        MergeSlice(impl_->count_slices[b], impl_->counts);
    }
    if (records_processed) *records_processed = records;
    return 0;
}

int FlagstatOperator::Finish() {
    return 0;
}

int RunFlagstatPipeline(const BamInputBackend &input,
                        const BgzfBlockSpan *spans,
                        size_t span_count,
                        FlagstatCounts *counts,
                        FlagstatMetrics *metrics,
                        const cpe::CpeReadPipelineOptions &options) {
    if (!counts) return -1;

    FlagstatOperator op(counts, metrics);
    cpe::CpeReadPipelineTiming timing;
    const int ret = cpe::RunCpeReadPipeline(
        input, spans, span_count, &op, &timing, options);
    if (metrics) {
        metrics->input_blocks = timing.input_blocks;
        metrics->batch_count = timing.batch_count;
        metrics->total_records = timing.total_records;
        metrics->read = timing.read;
        metrics->kernel = timing.kernel;
        metrics->post_process = timing.post_process;
        metrics->total = timing.total;
    }
    return ret;
}

} // namespace operators
} // namespace swbam
