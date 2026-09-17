#include "cgs_memory_probe_common.h"

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

extern "C" void slave_cgs_memory_probe(CgsMemoryProbeArgs *args) {
#ifdef PLATFORM_SUNWAY
    const int tid = args && args->use_cgs_tid
        ? (int)_CGN * SWBAM_CGS_PROBE_PES_PER_GROUP + (int)_PEN
        : (int)_PEN;
#else
    const int tid = 0;
#endif
    if (!args || tid < 0 || args->worker_count <= 0 ||
        tid >= args->worker_count ||
        args->worker_count > SWBAM_CGS_PROBE_THREADS) {
        return;
    }

    uint64_t count = 0;
    uint64_t checksum = 0;
    if (args->mode == SWBAM_CGS_PROBE_TOUCH) {
        const uint64_t stride = args->stride ? args->stride : 4096ULL;
        const uint64_t pages =
            (args->bytes + stride - 1ULL) / stride;
        const uint64_t begin =
            pages * (uint64_t)tid / (uint64_t)args->worker_count;
        const uint64_t end =
            pages * (uint64_t)(tid + 1) /
            (uint64_t)args->worker_count;
        for (uint64_t page = begin; page < end; ++page) {
            const uint64_t offset = page * stride;
            if (offset < args->bytes) {
                args->data[offset] = CgsMemoryProbePattern(page);
                ++count;
            }
        }
    } else if (args->mode == SWBAM_CGS_PROBE_CHECKSUM) {
        const uint64_t begin =
            args->bytes * (uint64_t)tid /
            (uint64_t)args->worker_count;
        const uint64_t end =
            args->bytes * (uint64_t)(tid + 1) /
            (uint64_t)args->worker_count;
        for (uint64_t offset = begin; offset < end; ++offset) {
            checksum += CgsMemoryProbeContribution(
                args->data[offset], offset);
        }
        count = end - begin;
    }

    args->counts[tid] = count;
    args->checksums[tid] = checksum;
    args->seen[tid] = tid + 1;
}
