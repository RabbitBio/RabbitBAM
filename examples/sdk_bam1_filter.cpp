#include "swbam/swbam.h"

#include <chrono>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
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

int DumpMemoryOutput(const char *path,
                     const swbam::MemoryBamOutput &output) {
    FILE *file = fopen(path, "wb");
    if (!file) return -1;
    size_t offset = 0;
    while (offset < output.size()) {
        const size_t written = fwrite(
            output.data() + offset, 1, output.size() - offset, file);
        if (written == 0) {
            fclose(file);
            return -1;
        }
        offset += written;
    }
    return fclose(file) == 0 ? 0 : -1;
}

struct FilterMetrics {
    long long total;
    long long kept;
    long long dropped;
    double filter;

    FilterMetrics() : total(0), kept(0), dropped(0), filter(0.0) {}
};

class MapqFilterBatchPostProcessor : public swbam::cpe::Bam1BatchPostProcessor {
public:
    MapqFilterBatchPostProcessor(int min_mapq,
                       swbam::cpe::Bam1BatchPostProcessor *downstream)
        : min_mapq_(min_mapq), downstream_(downstream) {
        kept_.reserve(64 * 1024);
    }

    int PostProcessBam1Batch(const bam1_t *const *records, size_t count) {
        if (!downstream_ || (!records && count != 0)) return -1;
        const double filter_t0 = NowSeconds();
        kept_.clear();
        for (size_t i = 0; i < count; ++i) {
            if (!records[i]) return -1;
            if (records[i]->core.qual >= min_mapq_) {
                kept_.push_back(records[i]);
            }
        }
        metrics_.total += (long long)count;
        metrics_.kept += (long long)kept_.size();
        metrics_.dropped +=
            (long long)count - (long long)kept_.size();
        metrics_.filter += NowSeconds() - filter_t0;
        return downstream_->PostProcessBam1Batch(
            kept_.empty() ? nullptr : kept_.data(), kept_.size());
    }

    const FilterMetrics &metrics() const { return metrics_; }

private:
    int min_mapq_;
    swbam::cpe::Bam1BatchPostProcessor *downstream_;
    FilterMetrics metrics_;
    std::vector<const bam1_t *> kept_;
};

} // namespace

int main(int argc, char **argv) {
    if (argc != 4 && argc != 5) {
        fprintf(stderr,
                "Usage: %s input.bam output.bam min_mapq [--memory-io]\n",
                argv[0]);
        return 2;
    }
    const bool memory_io = argc == 5 &&
        strcmp(argv[4], "--memory-io") == 0;
    if (argc == 5 && !memory_io) {
        fprintf(stderr, "Unknown option: %s\n", argv[4]);
        return 2;
    }
    int min_mapq = 0;
    if (ParseMapq(argv[3], &min_mapq) != 0) {
        fprintf(stderr, "Invalid min_mapq: %s\n", argv[3]);
        return 2;
    }

    int exit_code = 1;
    std::unique_ptr<swbam::BamInputBackend> input(
        memory_io
            ? static_cast<swbam::BamInputBackend *>(
                  new swbam::MemoryBamInput())
            : static_cast<swbam::BamInputBackend *>(
                  new swbam::PosixBamInput()));
    swbam::MemoryBamOutput memory_output;
    swbam::PosixBamOutput posix_output;
    swbam::BamOutputBackend *output = memory_io
        ? static_cast<swbam::BamOutputBackend *>(&memory_output)
        : static_cast<swbam::BamOutputBackend *>(&posix_output);
    std::vector<swbam::BgzfBlockSpan> spans;
    swbam::cpe::Bam1Writer writer;
    MapqFilterBatchPostProcessor filter(min_mapq, &writer);
    swbam::cpe::CpeReadPipelineTiming timing;
    swbam::cpe::ComposableDecodeMetrics decode_metrics;
    swbam::cpe::ComposableBam1Metrics bam1_metrics;
    double open_seconds = 0.0;
    double scan_seconds = 0.0;
    double core_seconds = 0.0;
    double dump_seconds = 0.0;
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
    if (memory_io) {
        if (memory_output.Reserve(input->size()) != 0) {
            fprintf(stderr, "Failed to reserve memory BAM output\n");
            goto cleanup;
        }
    } else if (posix_output.Open(argv[2]) != 0) {
        fprintf(stderr, "Failed to initialize BAM output: %s\n", argv[2]);
        goto cleanup;
    }

    phase_t0 = NowSeconds();
    if (writer.InitializeBam(output, input->header(), 1) != 0) {
        fprintf(stderr, "Failed to initialize bam1_t writer\n");
        goto cleanup;
    }
    if (swbam::cpe::RunComposableBam1Pipeline(
            *input, spans.empty() ? nullptr : spans.data(), spans.size(),
            &filter, &timing, &decode_metrics, &bam1_metrics) != 0 ||
        writer.Finish() != 0 || (!memory_io && posix_output.Close() != 0)) {
        fprintf(stderr, "Composable bam1_t filter pipeline failed\n");
        goto cleanup;
    }
    core_seconds = NowSeconds() - phase_t0;

    if (memory_io) {
        phase_t0 = NowSeconds();
        if (DumpMemoryOutput(argv[2], memory_output) != 0) {
            fprintf(stderr, "Failed to dump BAM output: %s\n", argv[2]);
            goto cleanup;
        }
        dump_seconds = NowSeconds() - phase_t0;
    }

    printf("total_records %lld\n", filter.metrics().total);
    printf("kept_records %lld\n", filter.metrics().kept);
    printf("dropped_records %lld\n", filter.metrics().dropped);
    printf("decoded_bytes %lld\n", decode_metrics.decoded_bytes);
    printf("encoded_record_bytes %lld\n",
           writer.metrics().encoded_bytes);
    printf("output_bytes %llu\n",
           (unsigned long long)output->bytes_written());
    printf("io_backend %s\n", memory_io ? "memory" : "posix");
    printf("timing_open %.6f\n", open_seconds);
    printf("timing_scan %.6f\n", scan_seconds);
    printf("timing_core %.6f\n", core_seconds);
    printf("timing_read_pipeline %.6f\n", timing.total);
    printf("timing_materialize %.6f\n", bam1_metrics.materialize);
    printf("timing_filter %.6f\n", filter.metrics().filter);
    printf("timing_encode %.6f\n", writer.metrics().encode);
    printf("timing_pack %.6f\n", writer.raw_metrics().pack);
    printf("timing_compress %.6f\n", writer.raw_metrics().kernel);
    printf("timing_writer %.6f\n", writer.raw_metrics().total);
    printf("timing_dump %.6f\n", dump_seconds);
    exit_code = 0;

cleanup:
    posix_output.Close();
    input->Close();
#ifdef PLATFORM_SUNWAY
    athread_halt();
#endif
    return exit_code;
}
