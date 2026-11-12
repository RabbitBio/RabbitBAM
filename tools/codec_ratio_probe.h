#ifndef SWBAM_CODEC_RATIO_PROBE_H
#define SWBAM_CODEC_RATIO_PROBE_H

#include "swbam/bam_types.h"

struct CodecRatioPara {
    bam_block *input;
    bam_block *decoded;
    bam_block *original;
    bam_block *optimized;
    int status;
};

extern "C" void slave_codec_ratio_probe(CodecRatioPara *parameters);

#endif
