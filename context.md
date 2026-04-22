## 项目交接文档（RabbitBAM / 神威平台）

### 0. 会话来源与上下文
- 主要历史会话可参考：[SAM-BAM性能优化交接](3fe13974-290b-4851-960f-b53f85056b95)
- 当前项目路径核心在 `RabbitBAM`，重点文件长期集中于：
  - `src/swbam.cpp`
  - `slave/slave.cpp`
  - `src/BamTools.cpp`
  - `src/BamWrite.cpp`
  - `src/BamWriteComplete.cpp`
  - `src/BamRead.cpp`
  - `src/BamComplete.cpp`
  - `include/BamTools.h`
  - `include/swbam.h`
  - `src/main.cpp`

---

## 1. 核心研究目标

- 在**神威平台（1 MPE + 64 CPE）**上优化 `SAM -> BAM` 与 `BAM -> SAM` 转换性能。
- 优化重点不是单条算法微调，而是：
  - 流程组织（pipeline/fused）
  - 主从核时间掩盖（overlap）
  - 内存分配与对象池
  - 队列与批处理开销
- 保证正确性前提下，逼近从核阶段理论下限。
- 额外功能目标：
  - 支持以子命令方式检查 BAM 是否存在**跨 BGZF 块记录**（用于判定 bam2sam 快速路径是否可用）。

---

## 2. 技术栈与架构

### 2.1 技术栈
- 语言：C++（含 C 风格 HTSlib 接口）
- 依赖：
  - `htslib`（sam/bam/bgzf）
  - `libdeflate`
  - 神威并行运行时（`athread`）
- 构建：CMake + make
- 命令行：CLI11（`main.cpp` 的 `add_subcommand`）

### 2.2 架构模型
- 平台特征：MPE 主核单线程控制 + 64 CPE 从核并行计算。
- 关键约束：
  - 主核与从核内存访问模型受限
  - 同时活跃的从核任务受 `spawn/join` 控制
- 典型设计：
  - 对象池（`BamWrite` / `BamComplete`）
  - 环形队列（读、写、完成队列）
  - 批处理接口（`getEmptyBatch/backBatch`）
  - fused 单线程编排（把 producer/consumer/writer 协同进一个控制流，减少线程切换与锁等待）

### 2.3 当前总体流程（高层）
- `SAM -> BAM`：
  - copy/count -> parse -> pack -> compress -> write/back
- `BAM -> SAM`：
  - read(memcpy模拟磁盘) -> decompress -> parse/read_bam -> format -> write

---

## 3. 已解决的关键技术难点

### 3.1 正确性问题（已修复）
1. **9 条 `bam1_t` 变空（`l_data==0`）**
- 根因：`active_chunks < 64` 时，未清零 chunk 被从核误处理，旧指针被覆盖，污染已入队对象。
- 修复：
  - 未使用 chunk 显式 `text_len=0, count=0`
  - 解析失败导致多分配对象显式回收（`valid_count < pre_alloc_count` 回收尾部）

2. **池泄漏/死锁**
- 根因：批量预分配和回收不平衡，或池规模不足导致 `getEmpty` 阻塞。
- 修复：
  - 引入 batch get/back 接口
  - 调整池容量（并根据飞行中对象数估算）

### 3.2 性能优化（已落地）
1. **批量内存分配**
- `BamWrite`、`BamComplete` 使用大块分配替代逐条分配，减少初始化与碎片开销。

2. **批处理接口**
- `getEmptyBatch` / `backBamBatch` 等，降低函数调用与锁开销。

3. **`bam_len` 预计算**
- 在从核 parse 阶段预计算，减少 pack 阶段额外访存与 cache miss。

4. **fused 流程**
- 引入 `FusedSamToBam`、`FusedBamToSam`，减少多线程流水线中主核调度/同步损耗。

5. **BAM->SAM 的 memcpy 隐藏策略验证**
- 将 read(memcpy) 从 decomp 窗口移到 format 窗口，规避带宽竞争，整体显著提升（该方向验证有效）。

### 3.3 边界功能：跨块记录检测（新增）
- `main.cpp` 新增子命令（`check_cross_block`）
- `BamTools` 增加检测函数
- 现阶段已发现并修复过一次误判逻辑（基于 `block_address` 的边界误判），后续实现应坚持“按块内字节布局扫描记录边界”的判定方法。

---

## 4. 当前瓶颈与结论（重要给后续 Agent）

### 4.1 `SAM -> BAM` 难再提升的主要原因
- 关键路径基本为串行依赖：`copy_count -> parse -> pack/compress`
- `compress` 与 `parse` 都吃从核窗口，难以相互掩盖
- 主核做大规模 memcpy 会明显拖慢（已实验验证）

### 4.2 为什么 `BAM -> SAM` 看起来更容易优化
- 某些阶段（如 format）计算密集、访存竞争相对低，更适合作为 read(memcpy) 掩盖窗口
- 其阶段结构更容易实现有效 overlap

### 4.3 三代测序数据支持的可行性（当前结论）
- 现实现偏向二代数据且假设“记录不跨块”。
- 若支持三代（长读长，跨多块高频）需重构：
  - 跨块拼接状态机
  - 多块连续记录装配
  - 读解码路径接口重写（`read_bam`/block reader）
- 串行拼接是硬开销，并行收益可能被吞噬；如无明确业务刚需，不建议优先投入。

---

## 5. 编码规范与用户个人偏好（高优先级）

### 5.1 你（用户）的明确偏好
- 重视**真实性能数据**，每次改动要能解释“为什么快/慢”。
- 偏好先做**流程级优化**，不是局部计算微优化。
- 强调神威特性：把工作尽量放从核，主核做调度与掩盖。
- 对“省略 read”这类不符合实验建模的方案明确拒绝；`memcpy` 读用于模拟磁盘 I/O，不能删。
- 要求保留旧实现并用编译宏切换新实现（可回退、可对比）。
- 出现“卡住/阻塞”时优先定位池容量、队列背压、批量分配时机。

### 5.2 已形成的代码风格惯例
- 优先 batch API，避免逐条 get/back。
- 在从核阶段尽量顺带计算辅助元数据（如 `bam_len`）。
- 强调“先 correctness，再 perf”，每次性能优化后都要校验输出一致性。
- 关键流程保留明确计时输出（各阶段耗时）。

### 5.3 后续 Agent 的工作原则
- 不要随意改动已验证正确的路径；先做最小侵入改动。
- 每次改动后至少给出：
  - 功能正确性结论
  - 分阶段时间变化
  - 与神威带宽/并行模型一致的解释
- 对可能引入死锁的地方（池容量、预分配时机）优先做容量推导。

---

## 6. 建议下一步（给后续 Agent）

- 在不改变算法语义前提下，优先做：
  1. 阶段粒度计时标准化（固定打印格式，便于横向对比）
  2. 对象池“高水位”统计（验证是否仍有隐性背压）
  3. `SAM->BAM` 中 copy/count 与 parse/compress 的更细粒度重叠实验（仅在不引入主核大 memcpy 的前提下）
- 若推进三代支持，应先做“只检测/只告警”模式，不直接接管主路径。