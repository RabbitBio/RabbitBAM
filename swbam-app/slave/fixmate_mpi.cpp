#include <cstring>
#include <stdint.h>

#include "algorithm_types.h"

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

#include "fixmate_steps.h"
using namespace swbam_fixmate_steps;

extern "C" void slave_mpi_fixmate_plan(MpiFixmatePlanPara paras[64]) {
    const int id = _PEN;
    MpiFixmatePlanPara *para = &paras[id];
    if (para->status != 0 || !para->records || !para->plans ||
        !para->groups) {
        return;
    }

    para->paired_groups = 0;
    para->singleton_groups = 0;
    para->secondary_records = 0;
    para->supplementary_records = 0;
    para->mq_updates = 0;
    para->mc_updates = 0;
    para->ms_updates = 0;

    for (int g = para->group_begin; g < para->group_end; ++g) {
        const int begin = para->groups[g].begin;
        const int end = begin + para->groups[g].count;
        int pre = -1;
        int cur = -1;

        for (int i = begin; i < end; ++i) {
            bam1_t *record = para->records[i];
            if (!record) {
                para->status = -2;
                para->record_index = i;
                return;
            }
            MpiFixmateRecordPlanShared *plan = &para->plans[i];
            memset(plan, 0, sizeof(*plan));
            plan->core = record->core;
            plan->output_data_len = (uint32_t)record->l_data;
            plan->mate_index = -1;
            if (record->core.flag & BAM_FSECONDARY) {
                para->secondary_records++;
                continue;
            }
            if (record->core.flag & BAM_FSUPPLEMENTARY) {
                para->supplementary_records++;
                continue;
            }
            if (pre < 0) {
                pre = i;
                continue;
            }
            cur = i;
            fm_pair_records(para->records, para->plans, pre, cur);
        }

        if (cur < 0 && pre >= 0) {
            MpiFixmateRecordPlanShared *plan = &para->plans[pre];
            plan->core.mtid = -1;
            plan->core.mpos = -1;
            plan->core.isize = 0;
            plan->core.flag &=
                (uint16_t)~(BAM_FPAIRED | BAM_FMREVERSE | BAM_FPROPER_PAIR);
            para->singleton_groups++;
        } else if (cur >= 0) {
            para->paired_groups++;
        }

        for (int i = begin; i < end; ++i) {
            if (fm_plan_output_length(
                    para->records[i], &para->plans[i],
                    para->records, end) != 0) {
                para->status = -3;
                para->record_index = i;
                return;
            }
            para->mq_updates += para->plans[i].update_mq;
            para->mc_updates += para->plans[i].update_mc;
            para->ms_updates += para->plans[i].update_ms;
        }
    }
    para->status = 0;
}

extern "C" void slave_mpi_fixmate_rewrite(
        MpiFixmateRewritePara paras[64]) {
    const int id = _PEN;
    MpiFixmateRewritePara *para = &paras[id];
    if (para->status != 0 || !para->records ||
        !para->output_records || !para->plans ||
        !para->output_offsets || !para->output_data) {
        return;
    }

    for (int i = para->record_begin; i < para->record_end; ++i) {
        const uint64_t offset = para->output_offsets[i];
        const uint64_t next = para->output_offsets[i + 1];
        const MpiFixmateRecordPlanShared *plan = &para->plans[i];
        // The next CPE's slice may start after alignment padding.
        if (next < offset ||
            next - offset < plan->output_data_len ||
            next > para->output_capacity) {
            para->status = -3;
            para->record_index = i;
            return;
        }
        bam1_t *src = para->records[i];
        bam1_t *dst = para->output_records + i;
        memset(dst, 0, sizeof(*dst));
        dst->core = plan->core;
        dst->id = src->id;
        dst->data = para->output_data + offset;
        dst->l_data = (int)plan->output_data_len;
        dst->m_data = plan->output_data_len;
        dst->mempolicy = BAM_USER_OWNS_DATA;
        if (fm_copy_aux_and_append(
                src, dst, plan, para->records) != 0) {
            para->status = -2;
            para->record_index = i;
            return;
        }
    }
    para->status = 0;
}
