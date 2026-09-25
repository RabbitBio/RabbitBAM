#include "swbam/operators/record_count.h"
#include "swbam/cpe_bam_read_steps.h"

#include <climits>

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

namespace {

class CountAction : public SwbamCpeScratchRecord {
public:
    explicit CountAction(swbam::operators::RecordCountPara *para)
        : SwbamCpeScratchRecord(para->scratch_data, para->scratch_capacity) {}
    void Process(bam1_t *, int) {}
};

} // namespace

extern "C" void slave_swbam_record_count(
        swbam::operators::RecordCountPara paras[64]) {
    swbam::operators::RecordCountPara &para = paras[_PEN];
    if (!para.input_block) return;
    if (!para.scratch_data || !para.scratch_capacity ||
        swbam_cpe_decode_bam_block(
            para.input_block, para.un_comp_block, _PEN,
            nullptr, nullptr, nullptr) != 0) {
        para.status = -1;
        return;
    }
    CountAction action(&para);
    const SwbamCpeWalkResult result =
        swbam_cpe_walk_bam_records(para.un_comp_block, INT_MAX, action);
    para.records = result.count;
    para.status = result.read_result == -1 ? 0 : -2;
}
