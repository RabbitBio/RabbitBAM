#include "swbam_mpi.h"
#include "swbam/operators/sam_to_bam.h"

int OptimizedSamToBamMPI(MemReader &reader, MemWriter &mem_writer,
                         sam_hdr_t *header, int compress_level,
                         MpiSamToBamStats *stats) {
    swbam::MemoryRankBodySink sink(&mem_writer);
    return swbam::operators::RunSamToBamPipeline(
        &reader, &sink, header, compress_level, stats);
}

int OptimizedSamToBamMPI(MemReader &reader, swbam::RankBodySink &sink,
                         sam_hdr_t *header, int compress_level,
                         MpiSamToBamStats *stats) {
    return swbam::operators::RunSamToBamPipeline(
        &reader, &sink, header, compress_level, stats);
}
