#include "swbam_mpi.h"

int OptimizedBamToSamMPI(MemReader &reader, MemWriter &mem_writer,
                         sam_hdr_t *header, MpiBamToSamStats *stats) {
    swbam::MemoryRankBodySink sink(&mem_writer);
    return swbam::operators::RunBamToSamPipeline(
        &reader, &sink, header, stats);
}

int OptimizedBamToSamMPI(MemReader &reader,
                         swbam::RankBodySink &sink,
                         sam_hdr_t *header, MpiBamToSamStats *stats) {
    return swbam::operators::RunBamToSamPipeline(
        &reader, &sink, header, stats);
}

int OptimizedBamToSamMPI(const swbam::BamInputBackend &input,
                         const swbam::BgzfBlockSpan *spans, size_t span_count,
                         MemWriter &mem_writer, sam_hdr_t *header,
                         MpiBamToSamStats *stats) {
    swbam::MemoryRankBodySink sink(&mem_writer);
    return swbam::operators::RunBamToSamPipeline(
        input, spans, span_count, &sink, header, stats);
}

int OptimizedBamToSamMPI(const swbam::BamInputBackend &input,
                         const swbam::BgzfBlockSpan *spans, size_t span_count,
                         swbam::RankBodySink &sink, sam_hdr_t *header,
                         MpiBamToSamStats *stats) {
    return swbam::operators::RunBamToSamPipeline(
        input, spans, span_count, &sink, header, stats);
}
