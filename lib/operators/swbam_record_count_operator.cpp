#include "swbam/operators/record_count.h"

#include <cstdio>

extern "C" void slave_swbam_record_count();

namespace swbam {
namespace operators {
namespace {

class RecordCountOperator : public cpe::CpeBatchOperator {
public:
    const char *name() const { return "record-count"; }
    size_t batch_capacity() const { return 64; }
    int Initialize() { return 0; }
    void Shutdown() {}

    int Prepare(const BgzfBlockBatch &compressed,
                BgzfBlockBatch *decoded, size_t count) {
        if (!decoded || count > 64) return -1;
        for (size_t i = 0; i < 64; ++i) {
            paras_[i].compressed = i < count
                ? const_cast<bam_block *>(&compressed.blocks()[i]) : nullptr;
            paras_[i].decoded = i < count ? &decoded->blocks()[i] : nullptr;
            paras_[i].records = 0;
            paras_[i].status = i < count ? 0 : -1;
        }
        return 0;
    }

    void *kernel_entry() const {
        return reinterpret_cast<void *>(slave_swbam_record_count);
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
