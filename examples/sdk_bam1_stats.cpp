#include "swbam/swbam.h"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <memory>
#include <vector>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

namespace {

double NowSeconds() {
    typedef std::chrono::steady_clock Clock;
    return std::chrono::duration<double>(
        Clock::now().time_since_epoch()).count();
}

class Bam1StatsBatchPostProcessor : public swbam::cpe::Bam1BatchPostProcessor {
public:
    Bam1StatsBatchPostProcessor()
        : records_(0), mapped_(0), duplicates_(0), mapq_ge_30_(0),
          nm_records_(0), nm_sum_(0) {}

    int PostProcessBam1Batch(const bam1_t *const *records, size_t count) {
        if (!records && count != 0) return -1;
        for (size_t i = 0; i < count; ++i) {
            const bam1_t *record = records[i];
            if (!record || !record->data || !bam_get_qname(record)) {
                return -1;
            }
            if (record->core.n_cigar > 0 && !bam_get_cigar(record)) {
                return -1;
            }
            if (record->core.l_qseq > 0 &&
                (!bam_get_seq(record) || !bam_get_qual(record))) {
                return -1;
            }

            ++records_;
            if ((record->core.flag & BAM_FUNMAP) == 0) ++mapped_;
            if ((record->core.flag & BAM_FDUP) != 0) ++duplicates_;
            if (record->core.qual >= 30) ++mapq_ge_30_;

            uint8_t *nm = bam_aux_get(record, "NM");
            if (nm) {
                ++nm_records_;
                nm_sum_ += bam_aux2i(nm);
            }
        }
        return 0;
    }

    long long records() const { return records_; }
    long long mapped() const { return mapped_; }
    long long duplicates() const { return duplicates_; }
    long long mapq_ge_30() const { return mapq_ge_30_; }
    long long nm_records() const { return nm_records_; }
    long long nm_sum() const { return nm_sum_; }

private:
    long long records_;
    long long mapped_;
    long long duplicates_;
    long long mapq_ge_30_;
    long long nm_records_;
    long long nm_sum_;
};

} // namespace

int main(int argc, char **argv) {
    if (argc != 2 && argc != 3) {
        fprintf(stderr, "Usage: %s input.bam [--memory-io]\n", argv[0]);
        return 2;
    }
    const bool memory_io = argc == 3 &&
        strcmp(argv[2], "--memory-io") == 0;
    if (argc == 3 && !memory_io) {
        fprintf(stderr, "Unknown option: %s\n", argv[2]);
        return 2;
    }

    int exit_code = 1;
    std::unique_ptr<swbam::BamInputBackend> input(
        memory_io
            ? static_cast<swbam::BamInputBackend *>(
                  new swbam::MemoryBamInput())
            : static_cast<swbam::BamInputBackend *>(
                  new swbam::PosixBamInput()));
    std::vector<swbam::BgzfBlockSpan> spans;
    Bam1StatsBatchPostProcessor post_processor;
    swbam::cpe::CpeReadPipelineTiming timing;
    swbam::cpe::ComposableDecodeMetrics decode_metrics;
    swbam::cpe::ComposableBam1Metrics bam1_metrics;
    double open_seconds = 0.0;
    double scan_seconds = 0.0;
    double phase_t0 = 0.0;

#ifdef PLATFORM_SUNWAY
    athread_init();
#endif
    phase_t0 = NowSeconds();
    if (input->Open(argv[1]) != 0) {
        fprintf(stderr, "Failed to open BAM: %s\n", argv[1]);
        goto cleanup;
    }
    open_seconds = NowSeconds() - phase_t0;
    phase_t0 = NowSeconds();
    if (input->ScanBlocks(&spans) != 0) {
        fprintf(stderr, "Failed to scan BAM: %s\n", argv[1]);
        goto cleanup;
    }
    scan_seconds = NowSeconds() - phase_t0;
    if (swbam::cpe::RunComposableBam1Pipeline(
            *input, spans.empty() ? nullptr : spans.data(), spans.size(),
            &post_processor, &timing, &decode_metrics, &bam1_metrics) != 0) {
        fprintf(stderr, "Composable bam1_t pipeline failed\n");
        goto cleanup;
    }

    printf("records %lld\n", post_processor.records());
    printf("mapped %lld\n", post_processor.mapped());
    printf("duplicates %lld\n", post_processor.duplicates());
    printf("mapq_ge_30 %lld\n", post_processor.mapq_ge_30());
    printf("nm_records %lld\n", post_processor.nm_records());
    printf("nm_sum %lld\n", post_processor.nm_sum());
    printf("bgzf_blocks %zu\n", spans.size());
    printf("decoded_bytes %lld\n", decode_metrics.decoded_bytes);
    printf("encoded_record_bytes %lld\n", bam1_metrics.encoded_bytes);
    printf("bam1_data_bytes %lld\n", bam1_metrics.data_bytes);
    printf("peak_batch_records %zu\n", bam1_metrics.peak_batch_records);
    printf("io_backend %s\n", memory_io ? "memory" : "posix");
    printf("timing_open %.6f\n", open_seconds);
    printf("timing_scan %.6f\n", scan_seconds);
    printf("timing_core %.6f\n", timing.total);
    printf("timing_read %.6f\n", timing.read);
    printf("timing_cpe %.6f\n", timing.kernel);
    printf("timing_pipeline_post_process %.6f\n", timing.post_process);
    printf("timing_materialize %.6f\n", bam1_metrics.materialize);
    printf("timing_bam1_post_process %.6f\n", bam1_metrics.post_process);
    exit_code = 0;

cleanup:
    input->Close();
#ifdef PLATFORM_SUNWAY
    athread_halt();
#endif
    return exit_code;
}
