#include "swbam/swbam.h"

#include <cstdio>
#include <cstring>
#include <vector>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

namespace {

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
    if (argc != 2) {
        fprintf(stderr, "Usage: %s input.bam\n", argv[0]);
        return 2;
    }

    int exit_code = 1;
    swbam::PosixBamInput input;
    std::vector<swbam::BgzfBlockSpan> spans;
    CountConsumer consumer;
    swbam::cpe::CpeReadPipelineTiming timing;
    swbam::cpe::GenericDecodeMetrics decode_metrics;
    swbam::cpe::GenericRawBamMetrics raw_metrics;

#ifdef PLATFORM_SUNWAY
    athread_init();
#endif
    if (input.Open(argv[1]) != 0 || input.ScanBlocks(&spans) != 0) {
        fprintf(stderr, "Failed to open or scan BAM: %s\n", argv[1]);
        goto cleanup;
    }
    if (swbam::cpe::RunGenericRawBamPipeline(
            input, spans.empty() ? nullptr : spans.data(), spans.size(),
            &consumer, &timing, &decode_metrics, &raw_metrics) != 0) {
        fprintf(stderr, "Generic raw BAM pipeline failed\n");
        goto cleanup;
    }
    printf("records %lld\n", consumer.records());
    printf("mapped %lld\n", consumer.mapped());
    printf("duplicates %lld\n", consumer.duplicates());
    printf("bgzf_blocks %zu\n", spans.size());
    exit_code = 0;

cleanup:
    input.Close();
#ifdef PLATFORM_SUNWAY
    athread_halt();
#endif
    return exit_code;
}
