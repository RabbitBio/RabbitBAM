#include "swbam/swbam.h"

#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <vector>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

namespace {

int ParseMapq(const char *text, int *value) {
    if (!text || !value) return -1;
    errno = 0;
    char *end = nullptr;
    const long parsed = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' ||
        parsed < 0 || parsed > 255) {
        return -1;
    }
    *value = (int)parsed;
    return 0;
}

} // namespace

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s input.bam output.bam min_mapq\n",
                argv[0]);
        return 2;
    }

    swbam::BamFilterOptions options =
        swbam::DefaultBamFilterOptions();
    if (ParseMapq(argv[3], &options.min_mapq) != 0) {
        fprintf(stderr, "Invalid min_mapq: %s\n", argv[3]);
        return 2;
    }

    int exit_code = 1;
    swbam::PosixBamInput input;
    swbam::PosixBamOutput output;
    std::vector<swbam::BgzfBlockSpan> spans;
    swbam::cpe::RawBamWriter writer;
    swbam::cpe::RawBamFilterConsumer filter(options, &writer);
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
    if (output.Open(argv[2]) != 0 ||
        writer.InitializeBam(&output, input.header(), 1) != 0) {
        fprintf(stderr, "Failed to initialize BAM output: %s\n", argv[2]);
        goto cleanup;
    }
    if (swbam::cpe::RunGenericRawBamPipeline(
            input, spans.empty() ? nullptr : spans.data(), spans.size(),
            &filter, &timing, &decode_metrics, &raw_metrics) != 0 ||
        writer.Finish() != 0 || output.Close() != 0) {
        fprintf(stderr, "Generic BAM filter pipeline failed\n");
        goto cleanup;
    }
    printf("total_records %lld\n", filter.metrics().total_records);
    printf("kept_records %lld\n", filter.metrics().kept_records);
    printf("dropped_records %lld\n", filter.metrics().dropped_records);
    exit_code = 0;

cleanup:
    output.Close();
    input.Close();
#ifdef PLATFORM_SUNWAY
    athread_halt();
#endif
    return exit_code;
}
