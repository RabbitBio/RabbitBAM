#ifndef SWBAM_OPERATORS_BAM_READ_BATCH_H
#define SWBAM_OPERATORS_BAM_READ_BATCH_H

#include "swbam/cpe_pipeline.h"

namespace swbam {
namespace operators {

struct BamReadKernelSpec {
    const char *name;
    void *entry;
    size_t batch_capacity;
};

// The spec describes one CPE launch per batch. Derived classes only own
// operator-specific scratch, result slices and batch post-processing.
class BamReadBatchOperator : public cpe::CpeBatchOperator {
public:
    explicit BamReadBatchOperator(const BamReadKernelSpec &spec)
        : spec_(spec) {}

    const char *name() const { return spec_.name; }
    size_t batch_capacity() const { return spec_.batch_capacity; }
    void *kernel_entry() const { return spec_.entry; }

private:
    BamReadKernelSpec spec_;
};

// Common batch boundary for CPE BAM-read operators. Each operator still owns
// its output slices and CPE parameter layout.
template <typename Para>
inline void BindBamReadBatch(Para *paras, size_t slots,
                             const BgzfBlockBatch &compressed,
                             BgzfBlockBatch *decoded, size_t active) {
    for (size_t i = 0; i < slots; ++i) {
        paras[i].input_block = i < active
            ? const_cast<bam_block *>(&compressed.blocks()[i]) : nullptr;
        paras[i].un_comp_block = i < active
            ? &decoded->blocks()[i] : nullptr;
        paras[i].status = i < active ? 0 : -1;
    }
}

} // namespace operators
} // namespace swbam

#endif
