#ifndef SWBAM_SWBAM_H
#define SWBAM_SWBAM_H

#include "swbam/cpe_pipeline.h"
#include "swbam/cpe_sam_read_pipeline.h"
#include "swbam/cpe_sam_write_pipeline.h"
#include "swbam/cpe_write_pipeline.h"
#include "swbam/cpe_record_write_adapter.h"
#include "swbam/io.h"
#include "swbam/operators/bgzf_compress.h"
#include "swbam/mpi_runtime.h"
#include "swbam/operators/flagstat.h"
#include "swbam/operators/bam_read_batch.h"
#include "swbam/operators/record_count.h"
#include "swbam/operators/stats_basic.h"
#include "swbam/operators/bam_transform.h"
#include "swbam/operators/bam_to_sam.h"
#include "swbam/operators/sam_to_bam.h"

#endif
