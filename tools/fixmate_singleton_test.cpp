#include "algorithm_types.h"
#include <cassert>
#include <cstdio>
#include <cstring>

extern "C" void slave_mpi_fixmate_plan(MpiFixmatePlanPara paras[64]);

static void Record(bam1_t *record, unsigned char *data,
                   uint16_t flags, int position) {
    memset(record, 0, sizeof(*record));
    memset(data, 0, 10);
    memcpy(data, "read", 4);
    data[3] = 0;
    data[4] = 1u << BAM_CIGAR_SHIFT;
    data[8] = 0x10;
    data[9] = 30;
    record->core.flag = flags;
    record->core.tid = 0;
    record->core.pos = position;
    record->core.mtid = 0;
    record->core.mpos = position + 100;
    record->core.isize = 101;
    record->core.qual = 60;
    record->core.l_qname = 4;
    record->core.n_cigar = 1;
    record->core.l_qseq = 1;
    record->data = data;
    record->l_data = record->m_data = 10;
}

int main() {
    bam1_t records[3];
    unsigned char data[3][10];
    bam1_t *pointers[] = {records, records + 1, records + 2};
    MpiFixmateRecordPlanShared plans[3] = {};
    MpiFixmateGroupShared group = {0, 3};
    MpiFixmatePlanPara paras[64] = {};
    paras[0].records = pointers;
    paras[0].plans = plans;
    paras[0].groups = &group;
    paras[0].group_end = 1;

    const uint16_t flags = BAM_FPAIRED | BAM_FREAD1 | BAM_FPROPER_PAIR |
        BAM_FMREVERSE | BAM_FREVERSE | BAM_FMUNMAP | BAM_FDUP;
    Record(records, data[0], flags, 100);
    Record(records + 1, data[1], flags | BAM_FSECONDARY, 200);
    Record(records + 2, data[2], flags | BAM_FSUPPLEMENTARY, 300);
    slave_mpi_fixmate_plan(paras);
    assert(paras[0].status == 0 && paras[0].singleton_groups == 1);
    assert(paras[0].paired_groups == 0);
    assert(plans[0].core.flag == (flags & ~(BAM_FPAIRED | BAM_FMREVERSE | BAM_FPROPER_PAIR)));
    assert(plans[0].core.mtid == -1 && plans[0].core.mpos == -1 && plans[0].core.isize == 0);
    assert(plans[0].mate_index == -1 && plans[0].output_data_len == 10);
    assert(!plans[0].update_mq && !plans[0].update_mc && !plans[0].update_ms);
    assert(records[0].core.flag == flags);
    assert(plans[1].core.flag == records[1].core.flag);
    assert(plans[2].core.flag == records[2].core.flag);

    Record(records, data[0], BAM_FPAIRED | BAM_FREAD1 | BAM_FPROPER_PAIR, 100);
    Record(records + 1, data[1], BAM_FPAIRED | BAM_FREAD2 | BAM_FREVERSE | BAM_FPROPER_PAIR, 200);
    group.count = 2;
    slave_mpi_fixmate_plan(paras);
    assert(paras[0].status == 0 && paras[0].paired_groups == 1);
    assert(paras[0].singleton_groups == 0);
    assert((plans[0].core.flag & BAM_FPAIRED) && (plans[1].core.flag & BAM_FPAIRED));
    assert(plans[0].mate_index == 1 && plans[1].mate_index == 0);
    assert(plans[0].core.mpos == 200 && plans[1].core.mpos == 100);
    assert(plans[0].update_ms && plans[1].update_ms);
    puts("Fixmate singleton/paired planner regression PASS");
}
