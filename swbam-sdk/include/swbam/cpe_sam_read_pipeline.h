#ifndef SWBAM_CPE_SAM_READ_PIPELINE_H
#define SWBAM_CPE_SAM_READ_PIPELINE_H

#include "swbam/bam_types.h"

#include <cstddef>

namespace swbam {
namespace cpe {

struct CpeSamReadTiming {
    long long input_chunks;
    long long chunk_groups;
    long long parse_fast_records;
    long long parse_fallback_records;
    double split;
    double copy_count;
    double parse;
    double parse_core;
    double parse_aux;
    double parse_cg;
    double parse_fallback;
    double parse_other;

    CpeSamReadTiming();
};

class SamParsedBatchPostProcessor {
public:
    virtual ~SamParsedBatchPostProcessor() {}

    // Parsed records and their data arenas are borrowed until this returns.
    virtual int PostProcessParsedBatch(const MpiSamParseChunk *chunks,
                                       size_t count) = 0;

    // Called on the MPE while the copy/count CPE kernel runs. An output
    // adapter may use this opportunity to flush its prior compressed batch.
    virtual int FlushPendingOutput() { return 0; }
};

// Copy/count and parse run as two existing CPE phases. Custom entries must
// accept MpiSamParseBatch and produce the same chunk/status/arena contract.
// Applications usually customize only the MPE post-processor.
struct SamReadKernelSpec {
    void *copy_count_entry;
    void *parse_entry;
};

SamReadKernelSpec DefaultSamReadKernelSpec();

int RunCpeSamReadPipeline(MemReader *reader, const sam_hdr_t *header,
                          SamParsedBatchPostProcessor *post_processor,
                          CpeSamReadTiming *timing,
                          const SamReadKernelSpec &kernels);

int RunCpeSamReadPipeline(MemReader *reader, const sam_hdr_t *header,
                          SamParsedBatchPostProcessor *post_processor,
                          CpeSamReadTiming *timing);

} // namespace cpe
} // namespace swbam

#endif
