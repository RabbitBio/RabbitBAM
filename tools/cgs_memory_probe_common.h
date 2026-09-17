#ifndef SWBAM_TOOLS_CGS_MEMORY_PROBE_COMMON_H
#define SWBAM_TOOLS_CGS_MEMORY_PROBE_COMMON_H

#include <stdint.h>

enum {
    SWBAM_CGS_PROBE_GROUPS = 6,
    SWBAM_CGS_PROBE_PES_PER_GROUP = 64,
    SWBAM_CGS_PROBE_THREADS =
        SWBAM_CGS_PROBE_GROUPS * SWBAM_CGS_PROBE_PES_PER_GROUP
};

enum CgsMemoryProbeMode {
    SWBAM_CGS_PROBE_TOUCH = 1,
    SWBAM_CGS_PROBE_CHECKSUM = 2
};

struct CgsMemoryProbeArgs {
    unsigned char *data;
    uint64_t bytes;
    uint64_t stride;
    uint64_t *counts;
    uint64_t *checksums;
    int *seen;
    int mode;
    int worker_count;
    int use_cgs_tid;
};

static inline unsigned char CgsMemoryProbePattern(uint64_t page) {
    return (unsigned char)((page * 131ULL + 17ULL) & 0xffULL);
}

static inline uint64_t CgsMemoryProbeContribution(
        unsigned char value, uint64_t offset) {
    return ((uint64_t)value + 1ULL) *
        (offset + 0x9e3779b97f4a7c15ULL);
}

#endif
