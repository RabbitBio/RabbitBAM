#ifndef SWBAM_SWBAM_H
#define SWBAM_SWBAM_H

#include "swbam/bam1.h"
#include "swbam/bam1_writer.h"
#include "swbam/cpe_pipeline.h"
#include "swbam/cpe_sam_read_pipeline.h"
#include "swbam/cpe_sam_write_pipeline.h"
#include "swbam/cpe_write_pipeline.h"
#include "swbam/cpe_record_write_adapter.h"
#include "swbam/composable_compress.h"
#include "swbam/composable_decode.h"
#include "swbam/io.h"
#include "swbam/mpi_runtime.h"
#include "swbam/operators/flagstat.h"
#include "swbam/operators/stats_basic.h"
#include "swbam/operators/bam_transform.h"
#include "swbam/operators/bam_to_sam.h"
#include "swbam/operators/sam_to_bam.h"
#include "swbam/raw_bam.h"
#include "swbam/raw_bam_filter.h"
#include "swbam/raw_bam_writer.h"

#endif
