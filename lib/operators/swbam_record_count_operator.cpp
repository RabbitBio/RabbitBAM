#include "swbam/operators/record_count.h"
#include "swbam/operators/bam_read_batch.h"

#include "BamTools.h"

#include <cstdio>

extern "C" void slave_swbam_record_count();

namespace swbam {
namespace operators {
namespace {

class RecordCountOperator : public BamReadBatchOperator {
public:
    RecordCountOperator()
        : BamReadBatchOperator({"record-count",
              reinterpret_cast<void *>(slave_swbam_record_count), 64}),
          scratch_(nullptr) {}
    ~RecordCountOperator() { Shutdown(); }
    int Initialize() {
        scratch_ = aligned_alloc_custom(64, 64 * MPI_BAM_BLOCK_ARENA_SIZE);
        return scratch_ ? 0 : -1;
    }
    void Shutdown() {
        if (scratch_) aligned_free_custom(scratch_);
        scratch_ = nullptr;
    }

    int Prepare(const BgzfBlockBatch &compressed,
                BgzfBlockBatch *decoded, size_t count) {
        if (!decoded || !scratch_ || count > 64) return -1;
        BindBamReadBatch(paras_, 64, compressed, decoded, count);
        for (size_t i = 0; i < 64; ++i) {
            paras_[i].scratch_data = scratch_ +
                i * MPI_BAM_BLOCK_ARENA_SIZE;
            paras_[i].scratch_capacity = MPI_BAM_BLOCK_ARENA_SIZE;
            paras_[i].records = 0;
        }
        return 0;
    }

    void *kernel_arguments() { return paras_; }
    void ObserveKernel(double, size_t) {}

    int Validate(size_t count) const {
        for (size_t i = 0; i < count; ++i) {
            if (paras_[i].status != 0) {
                fprintf(stderr, "ERROR: record count failed on BGZF block %zu (status %d).\n",
                        i, paras_[i].status);
                return -1;
            }
        }
        return 0;
    }

    int PostProcessBatch(size_t count, long long *records_processed) {
        long long records = 0;
        for (size_t i = 0; i < count; ++i) records += paras_[i].records;
        *records_processed = records;
        return 0;
    }
    int Finish() { return 0; }

private:
    RecordCountPara paras_[64];
    unsigned char *scratch_;
};

} // namespace

int RunRecordCountPipeline(const BamInputBackend &input,
                           const BgzfBlockSpan *spans, size_t span_count,
                           RecordCountMetrics *metrics) {
    if (!metrics) return -1;
    *metrics = RecordCountMetrics();
    RecordCountOperator op;
    cpe::CpeReadPipelineTiming timing;
    const int status = cpe::RunCpeReadPipeline(
        input, spans, span_count, &op, &timing);
    metrics->records = timing.total_records;
    metrics->blocks = timing.input_blocks;
    metrics->read = timing.read;
    metrics->kernel = timing.kernel;
    metrics->total = timing.total;
    return status;
}

} // namespace operators
} // namespace swbam
