#include "BamWrite.h"

BamWrite::BamWrite(int queue_size) {
    int total_size = queue_size * MAX_RECORDS_PER_BLOCK;
    pool_size = total_size + 1;
    bamPool = new bam1_t*[pool_size];
    pool_bg = 0;
    pool_ed = total_size - 1;

    // 批量分配所有 bam1_t 结构体（1 次 calloc 替代 total_size 次 bam_init1）
    bam1_struct_pool_ = (bam1_t*)calloc(total_size, sizeof(bam1_t));
    if (!bam1_struct_pool_) { perror("calloc bam1_t pool"); exit(1); }

    // 批量分配所有 data 缓冲区（1 次对齐大分配替代 total_size 次 aligned_alloc_custom）
    // INIT_DATA_SIZE = 1024 = 16×64，基地址 64 字节对齐后每个子缓冲区自然也是 64 字节对齐
    if (posix_memalign((void**)&data_pool_, 64, (size_t)total_size * INIT_DATA_SIZE) != 0) {
        perror("posix_memalign data pool"); exit(1);
    }

    for (int i = 0; i < total_size; i++) {
        bam1_t *b = &bam1_struct_pool_[i];
        b->data     = data_pool_ + (size_t)i * INIT_DATA_SIZE;
        b->m_data   = INIT_DATA_SIZE;
        b->l_data   = 0;
        b->mempolicy = BAM_USER_OWNS_DATA;
        bamPool[i]  = b;
    }

    q_size = queue_size + 5;
    groupQueue = new std::vector<std::vector<bam1_t*>>[q_size];
    q_bg = 1;
    q_ed = 0;

    write_complete = false;
}

BamWrite::~BamWrite() {
    free(data_pool_);
    free(bam1_struct_pool_);
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

void BamWrite::getEmptyBatch(bam1_t** dst, int n) {
    // 等待池中积累足够槽位
    while (true) {
        int avail = (pool_ed - pool_bg + pool_size) % pool_size;
        if (avail >= n) break;
        usleep(10);
    }
    // 批量 memcpy，处理环形缓冲区绕回
    int next = pool_bg;
    int part1 = pool_size - next;
    if (part1 >= n) {
        memcpy(dst, bamPool + next, (size_t)n * sizeof(bam1_t*));
    } else {
        memcpy(dst,        bamPool + next, (size_t)part1       * sizeof(bam1_t*));
        memcpy(dst + part1, bamPool,       (size_t)(n - part1) * sizeof(bam1_t*));
    }
    pool_bg = (pool_bg + n) % pool_size;
}

void BamWrite::backBam(bam1_t *b) {
    bamPool[(pool_ed + 1) % pool_size] = b;
    pool_ed = (pool_ed + 1) % pool_size;
}

void BamWrite::backBamBatch(bam1_t** const records, int n) {
    if (n <= 0) return;
    int next = (pool_ed + 1) % pool_size;
    int part1 = pool_size - next;
    if (part1 >= n) {
        memcpy(bamPool + next, records, (size_t)n * sizeof(bam1_t*));
    } else {
        memcpy(bamPool + next, records,        (size_t)part1       * sizeof(bam1_t*));
        memcpy(bamPool,        records + part1, (size_t)(n - part1) * sizeof(bam1_t*));
    }
    pool_ed = (pool_ed + n) % pool_size;
}

void BamWrite::inputGroup(std::vector<std::vector<bam1_t*>> group) {
    // move 进队列：避免深拷贝和 ~64 次子 vector 堆分配
    groupQueue[(q_ed + 1) % q_size] = std::move(group);
    q_ed = (q_ed + 1) % q_size;
}

std::vector<std::vector<bam1_t*>> BamWrite::getGroup() {
    while ((q_ed + 1) % q_size == q_bg) {
        usleep(10);
        if (write_complete && (q_ed + 1) % q_size == q_bg) return {};
    }

    // move 出队列：避免深拷贝，queue 槽位变为空状态，下次 inputGroup 复用
    auto group = std::move(groupQueue[q_bg]);
    q_bg = (q_bg + 1) % q_size;
    return group;
}

bool BamWrite::isComplete() const {
    return write_complete;
}

void BamWrite::markComplete() {
    write_complete = true;
}
