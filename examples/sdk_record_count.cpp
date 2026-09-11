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

uint32_t ReadLe32(const unsigned char *data) {
    return (uint32_t)data[0] |
           ((uint32_t)data[1] << 8) |
           ((uint32_t)data[2] << 16) |
           ((uint32_t)data[3] << 24);
}

class CountConsumer : public swbam::cpe::RawBamRecordConsumer {
public:
    CountConsumer() : records_(0), mapped_(0), duplicates_(0) {}

    int ConsumeRaw(const swbam::cpe::RawBamRecordView *records,
                   size_t count) {
        for (size_t i = 0; i < count; ++i) {
            if (!records[i].encoded || records[i].encoded_size < 20) {
                return -1;
            }
            const uint16_t flag = (uint16_t)(
                ReadLe32(records[i].encoded + 16) >> 16);
            ++records_;
            if ((flag & BAM_FUNMAP) == 0) ++mapped_;
            if ((flag & BAM_FDUP) != 0) ++duplicates_;
        }
        return 0;
    }

    long long records() const { return records_; }
    long long mapped() const { return mapped_; }
    long long duplicates() const { return duplicates_; }

private:
    long long records_;
    long long mapped_;
    long long duplicates_;
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
    CountConsumer consumer;
    swbam::cpe::CpeReadPipelineTiming timing;
    swbam::cpe::GenericDecodeMetrics decode_metrics;
    swbam::cpe::GenericRawBamMetrics raw_metrics;
    double open_seconds = 0.0;
    double scan_seconds = 0.0;
    double phase_t0 = 0.0;

#ifdef PLATFORM_SUNWAY
    athread_init();
#endif
    phase_t0 = NowSeconds();
    if (input->Open(argv[1]) != 0) {
        fprintf(stderr, "Failed to open or scan BAM: %s\n", argv[1]);
        goto cleanup;
    }
    open_seconds = NowSeconds() - phase_t0;
    phase_t0 = NowSeconds();
    if (input->ScanBlocks(&spans) != 0) {
        fprintf(stderr, "Failed to scan BAM: %s\n", argv[1]);
        goto cleanup;
    }
    scan_seconds = NowSeconds() - phase_t0;
    if (swbam::cpe::RunGenericRawBamPipeline(
            *input, spans.empty() ? nullptr : spans.data(), spans.size(),
            &consumer, &timing, &decode_metrics, &raw_metrics) != 0) {
        fprintf(stderr, "Generic raw BAM pipeline failed\n");
        goto cleanup;
    }
    printf("records %lld\n", consumer.records());
    printf("mapped %lld\n", consumer.mapped());
    printf("duplicates %lld\n", consumer.duplicates());
    printf("bgzf_blocks %zu\n", spans.size());
    printf("decoded_bytes %lld\n", decode_metrics.decoded_bytes);
    printf("io_backend %s\n", memory_io ? "memory" : "posix");
    printf("timing_open %.6f\n", open_seconds);
    printf("timing_scan %.6f\n", scan_seconds);
    printf("timing_core %.6f\n", timing.total);
    printf("timing_read %.6f\n", timing.read);
    printf("timing_cpe %.6f\n", timing.kernel);
    printf("timing_consume %.6f\n", timing.consume);
    exit_code = 0;

cleanup:
    input->Close();
#ifdef PLATFORM_SUNWAY
    athread_halt();
#endif
    return exit_code;
}
