#include "swbam_mpi.h"
#include "swbam/operators/bam_transform.h"

int OptimizedBamToBamMPI(MemReader &reader, MemWriter &mem_writer,
                          const BamFilterOptions &filter, int compress_level,
                          MpiBamToBamStats *stats) {
    swbam::MemoryRankBodySink sink(&mem_writer);
    return swbam::operators::RunBamTransformPipeline(
        &reader, &sink, filter, compress_level, stats);
}

int OptimizedBamToBamMPI(MemReader &reader, swbam::RankBodySink &sink,
                          const BamFilterOptions &filter, int compress_level,
                          MpiBamToBamStats *stats) {
    return swbam::operators::RunBamTransformPipeline(
        &reader, &sink, filter, compress_level, stats);
}

int OptimizedBamToBamMPI(const swbam::BamInputBackend &input,
                          const swbam::BgzfBlockSpan *spans, size_t count,
                          MemWriter &mem_writer,
                          const BamFilterOptions &filter, int compress_level,
                          MpiBamToBamStats *stats) {
    swbam::MemoryRankBodySink sink(&mem_writer);
    return swbam::operators::RunBamTransformPipeline(
        input, spans, count, &sink, filter, compress_level, stats);
}

int OptimizedBamToBamMPI(const swbam::BamInputBackend &input,
                          const swbam::BgzfBlockSpan *spans, size_t count,
                          swbam::RankBodySink &sink,
                          const BamFilterOptions &filter, int compress_level,
                          MpiBamToBamStats *stats) {
    return swbam::operators::RunBamTransformPipeline(
        input, spans, count, &sink, filter, compress_level, stats);
}
