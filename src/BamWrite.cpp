#include "BamWrite.h"

BamWrite::BamWrite(int queue_size) {
    int total_size = queue_size * MAX_RECORDS_PER_BLOCK;
    pool_size = total_size + 1 ; // 多申请一个，防止边界问题
    bamPool = new bam1_t*[pool_size];
    pool_bg = 0;
    pool_ed = total_size - 1;
    for (int i = pool_bg; i <= pool_ed; i++) {
        // bamPool[i] = bam_init1();

        //分配足够大小的空间
        bam1_t *b = bam_init1();
        // 分配 1KB 对齐内存
        b->data = (uint8_t*)aligned_alloc_custom(64, INIT_DATA_SIZE);
        b->m_data = INIT_DATA_SIZE;
        b->l_data = 0;
        b->mempolicy = BAM_USER_OWNS_DATA;
        bamPool[i] = b;    
    }

    q_size = queue_size + 5;
    groupQueue = new std::vector<std::vector<bam1_t*>>[q_size];
    q_bg = 1;
    q_ed = 0;

    write_complete = false;
}

BamWrite::~BamWrite() {
    for (int i = 0; i < pool_size-1; i++){
        if (bamPool[i]->data) {
            aligned_free_custom(bamPool[i]->data);
        }
    }
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
