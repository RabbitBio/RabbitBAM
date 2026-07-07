#ifndef SWBAM_CPE_BAM_PARSER_H
#define SWBAM_CPE_BAM_PARSER_H

#include "BamTools.h"

// Parses one BAM record from a decoded BGZF block. The caller owns record data.
// Returns bytes consumed, -1 at block end, -5 for insufficient record capacity,
// or another negative value for malformed input.
int swbam_cpe_read_bam_record(bam_block *decoded, bam1_t *record);

#endif
