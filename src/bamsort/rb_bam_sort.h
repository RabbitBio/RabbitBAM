#ifndef RB_BAM_SORT_H_
#define RB_BAM_SORT_H_

// Shared coordinate-sort + BAM writer extracted from RabbitBin.
// Used by `rabbitbam sortbam` and available to library callers.

#include <string>
#include <vector>

#include <htslib/sam.h>

// Stable coordinate sort (samtools-equivalent key + stable tie-break).
void rb_stable_coord_sort(std::vector<bam1_t *> &recs);

// Write already-sorted records as BAM to out_path ("-"/empty = stdout).
// Sets @HD SO:coordinate, optionally builds <out>.bai. Returns 0 on success.
int rb_write_sorted_bam(sam_hdr_t *hdr, std::vector<bam1_t *> &recs,
                        const std::string &out_path, int threads, int level,
                        bool write_index);

// Build <bam_path>.bai for a finished, coordinate-sorted BAM (htslib indexer).
// Run this on the MAIN thread: htslib's index thread pool does not engage when
// called from a background worker thread.
int rb_index_bam(const std::string &bam_path, int threads);

// Destroy all records in parallel and clear the vector.
void rb_free_records(std::vector<bam1_t *> &recs, int threads);

#endif  // RB_BAM_SORT_H_
