#include "BamWrite.h"

BamWrite::BamWrite(int queue_size) {
    int total_bams = 4096; //需要合适的大小
    pool_size = total_bams + 1 ; // 多申请一个，防止边界问题
    bamPool = new bam1_t*[pool_size];
    pool_bg = 0;
    pool_ed = total_bams - 1;
    for (int i = pool_bg; i <= pool_ed; i++) {
        bamPool[i] = bam_init1();
    }

    q_size = queue_size + 5;
    groupQueue = new std::vector<std::vector<bam1_t*>>[q_size];
    q_bg = 1;
    q_ed = 0;

    write_complete = false;
}

BamWrite::~BamWrite() {
    for (int i = 0; i < pool_size-1; i++) bam_destroy1(bamPool[i]);
    delete[] bamPool;
    delete[] groupQueue;
}

bam1_t* BamWrite::getEmpty() {
    while ((pool_ed + 1) % pool_size == pool_bg) {
        usleep(10);
    }
    int num = pool_bg;
    pool_bg = (pool_bg + 1) % pool_size;
    return bamPool[num];
}

void BamWrite::backBam(bam1_t *b) {
    bamPool[(pool_ed + 1) % pool_size] = b;
    pool_ed = (pool_ed + 1) % pool_size;
}

void BamWrite::inputGroup(const std::vector<std::vector<bam1_t*>>& group) {
    groupQueue[(q_ed + 1) % q_size] = group;
    q_ed = (q_ed + 1) % q_size;
}

std::vector<std::vector<bam1_t*>> BamWrite::getGroup() {
    while ((q_ed + 1) % q_size == q_bg) {
        usleep(10);
        if (write_complete && (q_ed + 1) % q_size == q_bg) return {};
    }

    int num = q_bg;
    auto group = groupQueue[q_bg];
    q_bg = (q_bg + 1) % q_size;
    return group;
}

bool BamWrite::isComplete() const {
    return write_complete;
}

void BamWrite::markComplete() {
    write_complete = true;
}
