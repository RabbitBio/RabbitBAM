# SWBAM 论文架构与创新设计草案

## 1. 文档定位

本文用于向导师说明 SWBAM 的研究问题、系统设计、核心创新、实验方案和论文组织。
它不是最终论文，也不把尚未充分验证的能力写成确定结论。当前最适合的论文定位是：

> 面向神威异构众核平台，设计一个以 BGZF 批次为扩展边界的分布式 BAM 处理库，
> 在统一数据通路上同时支持易扩展的通用算子和高性能专用融合算子，并通过统计、
> 转换、全局重排和有状态去重等多类应用验证其通用性与性能。

建议暂用题目：

- 中文：`SWBAM：面向神威异构众核平台的可扩展高性能 BAM 处理库`
- 英文：`SWBAM: An Extensible High-Performance BAM Processing Library for Sunway Many-Core Systems`

## 2. 论文要回答的核心问题

BAM 处理不仅包含 BGZF 解压缩，还包含记录解析、过滤统计、全局重排、跨块状态和
MPI 数据交换。若每个命令分别实现输入、缓冲、CPE 调度和输出，会产生三类问题：

1. **代码不可复用**：不同工具重复实现 BGZF 扫描、MPI 分区、双缓冲和输出布局。
2. **通用性与性能冲突**：逐记录通用回调便于开发，但容易产生中间数据写回、对象
   构造和主核二次扫描；完全专用的融合 kernel 很快，却难以复用。
3. **复杂算法难以扩展**：sort、collate、fixmate、markdup 具有全局通信或跨批状态，
   不能只依靠“并行解压每个 BGZF block”解决。

论文的核心研究问题应写成：

> 如何在神威 MPE/CPE 异构体系上建立一个可复用的 BAM 数据处理框架，同时避免通用
> 抽象进入逐记录热路径，并支持全局重排与有状态算法获得端到端加速？

## 3. 整体故事如何讲

### 3.1 从矛盾出发

论文不要按开发时间叙述“先写命令，后来抽库”。开发顺序不影响论文逻辑。论文应从
最终系统设计反向组织：

```text
现有问题
  通用 BAM I/O 易用但难以充分利用 CPE
  专用 kernel 性能高但代码重复、扩展困难
  复杂 BAM 算法还需要 MPI 重排和跨批状态
          ↓
关键观察
  BGZF block batch 是稳定的粗粒度扩展边界
          ↓
核心设计
  backend-neutral 数据平面
  + batch 级 MPE/CPE 双缓冲运行时
  + generic/fused 双执行路径
          ↓
算法实例
  统计、转换、sort、collate、fixmate、streaming markdup
          ↓
实验结论
  库化开销很小，融合路径更快，复杂应用获得加速并降低部分内存开销
```

### 3.2 一句话主张

全文可以围绕下面一句话展开：

> SWBAM 将抽象成本限制在 BGZF batch 边界，而把逐记录热路径留给 CPE 融合算子，
> 因而能够兼顾库级复用、算法可扩展性和高性能。

这句话必须由三个层次的实验支撑：公共库微基准、设计消融、应用与完整流程。

## 4. 核心创新点

建议把贡献稳定为三项。完整去重流程是系统效果，不单独拔高为第四项核心算法创新。

### 4.1 创新一：以 BGZF batch 为边界的分布式 BAM 数据平面

#### 要解决的问题

传统命令容易把存储读取、BGZF 扫描、MPI 分区、缓冲管理和业务算法写在一起。
这样更换内存、POSIX 或 MPI-IO 后端时需要修改算法，也难以控制大文件内存占用。

#### SWBAM 的实现

数据平面由以下公共组件组成：

```text
BamInputBackend
  ├─ MemoryBamInput
  ├─ PosixBamInput
  └─ MpiIoBamInput
          ↓
BgzfBlockSpan + BgzfBlockBatch
          ↓
MpiBamInputPlan：rank0 扫描、索引广播、rank 分区
          ↓
算法或 CPE batch pipeline
          ↓
RankBodySink
  ├─ MemoryRankBodySink
  ├─ SpoolRankBodySink
  └─ AdaptiveRankBodySink
          ↓
DistributedBamOutput：rank0 gather 或 MPI-IO offset write
```

具体机制包括：

- 输入后端统一输出 BGZF block spans 和对齐 batch，不向算法暴露存储细节。
- rank0 建立全局 block 索引，再广播并按连续 block 区间划分给 MPI ranks。
- POSIX scanner 采用 8 MiB 顺序窗口解析 `BSIZE`，减少大量小 `pread`。
- `RankBodySink` 将 rank 本地压缩 body 的生成与保存位置解耦，可选择内存或临时文件。
- `DistributedBamOutput` 统一计算 body size、prefix、header 和 EOF 布局，支持 rank0
  汇总或各 rank 直接 MPI-IO 写出。
- 算法外排临时文件与最终 rank body spool 分离，分别管理和计时。

#### 为什么不只是工程封装

创新重点不是类和接口本身，而是**选择 BGZF batch 作为系统边界**：

- 比文件级接口更细，可以流水和负载划分；
- 比逐记录接口更粗，不会产生数百万次虚调用或回调；
- 能同时承载 CPE kernel、MPI block 分区和多种存储后端；
- 复杂算法可以保留自己的全局状态，不被迫改写成通用 record callback。

#### 论文中建议这样表述

> We design a backend-neutral distributed BAM data plane that exposes aligned BGZF batches
> rather than records as the extension boundary. This design unifies block indexing, MPI
> partitioning, bounded rank-local output, and distributed file layout without introducing
> per-record abstraction overhead.

#### 支撑证据

- flagstat 库化前后约 `0.052073s -> 0.051994s`，性能持平。
- stats --basic 库化前后约 `0.098962s -> 0.099174s`，回退约 0.21%。
- memory、POSIX、MPI-IO、spool 路径已经形成完整闭环，但高速真实存储上的正式性能
  仍需补测。

### 4.2 创新二：通用与融合共存的 MPE/CPE 批流水执行模型

#### 要解决的问题

如果 CPE 只负责解压，MPE 再逐条解析和计算，会产生额外内存流量；如果每个功能都
编写完整专用流水，又会重复实现双缓冲、launch/join、错误处理和 codec cache。

#### 公共批流水

`RunCpeReadPipeline` 管理两套 compressed batch 和 decoded batch：

```text
时间  ──────────────────────────────────────────────>

CPE   compute batch 0      compute batch 1      compute batch 2
MPE       read batch 1         read batch 2         read batch 3
          <---- overlap ---->  <---- overlap ---->
```

公共 runtime 统一处理：

- 双缓冲分配和轮换；
- CPE kernel launch/join；
- CPE 处理当前批时由 MPE 读取下一批；
- 参数准备、状态校验、结果消费和计时；
- BGZF codec、CRC、parser 与每个 CPE 的 compressor/decompressor cache。

扩展接口 `CpeBatchOperator` 每批调用一次，而不是每条记录调用一次。

#### 通用路径

```text
BGZF batch
  → CPE 通用解压
  → MPE 扫描 record 边界
  → 零拷贝 RawBamRecordView
  → 用户 consumer
```

适用于快速开发记录计数、简单过滤和 SDK 新功能。`RawBamRecordView` 指向 decoded
batch 内存，不构造 `bam1_t`，也不复制 payload。

#### 专用融合路径

```text
BGZF batch
  → CPE 解压 + BAM 解析 + filter/count/extract
  → 小型计数或元数据结果
  → MPE 合并
```

例如 flagstat 的 CPE kernel 一次完成解压、解析和 flags 统计，仅返回计数 slice；
stats 同样在 CPE 上完成 SN、insert size、orientation 和局部有序性统计。这样避免将
完整 decoded records 写回后由 MPE 再扫描。

#### 设计价值

这不是两套互不相关的实现：通用路径与融合路径共享 backend、batch、双缓冲、codec、
parser、生命周期和错误处理，只替换 batch operator。它形成了明确的性能演进路线：

1. 先用 generic path 快速构建正确功能；
2. profiling 后只为热点实现 fused operator；
3. 不修改输入、输出和 MPI runtime。

#### 论文中建议这样表述

> SWBAM provides a dual-path execution model. A generic zero-copy record-view path improves
> programmability, while fused CPE operators combine decompression, parsing, and domain
> computation for performance-critical tools. Both paths share the same batch runtime, and
> dispatch occurs only once per batch.

#### 支撑证据

- 通用 BGZF 解压约 `0.0448s`，聚合解压吞吐约 `17.60 GB/s`。
- 通用 raw record 扫描约 `0.1255s`，约 `21.27 M records/s`。
- 库化后的融合 flagstat/stats 未出现结构性性能回退。
- 尚需正式补做 overlap on/off、generic/fused 同工作量消融实验。

### 4.3 创新三：在公共数据平面上实现全局重排和有状态 BAM 算法

#### 要解决的问题

BGZF block 可以独立解压，但 sort、collate、fixmate 和 markdup 的语义不能局限于
单个 block：

- sort 需要全局坐标顺序和 MPI 数据重分布；
- collate 需要相同 QNAME 连续且 READ1/READ2 顺序正确；
- fixmate 需要跨 block、跨 rank 的 mate group；
- markdup 需要坐标窗口、候选 ownership 和跨 rank 重复判定。

因此论文不能只声称“并行解压”，还要说明算法层如何映射到 MPI + CPE。

#### 应用分类

| 类型 | 实例 | 公共库负责 | 算法实例负责 |
|---|---|---|---|
| 线性流式 | flagstat、stats、filter、转换 | batch、CPE 调度、I/O | 逐记录语义 |
| 全局重排 | sort、collate | block 输入、sink、输出布局 | partition、exchange、merge |
| 有状态处理 | fixmate、markdup | 多遍 batch、spool/MPI-IO | 边界状态、ownership、窗口 |

#### 重点案例 A：collate

建议将 collate 作为“通信与数据布局优化”的代表案例，流程为：

```text
CPE extract QNAME/hash/meta
  → 按 bin/hash 预分组
  → MPI rank ownership 与 ring exchange
  → remote receive 直接写最终 raw arena
  → rank-local loser-tree merge
  → CPE compression + RankBodySink
```

可重点讲：

- 默认 bins 自动对齐 `comm_size`，降低 rank 间负载不均；
- raw arena 避免 `std::vector::resize()` 零初始化；
- remote receive 直接追加到 segment 最终位置，减少接收后二次拷贝；
- bucket/radix 预分组减少全局 QNAME 比较和 cache miss；
- loser tree 及 cursor cache 降低多路归并开销；
- memory 和 external 模式共享输入输出框架，算法外排独立于 body spool。

已有 `WES_0.25G` 初步结果为 `0.708s` 对 SAMtools `4.538s`，约 `6.41x`，但正式
论文必须统一计时和存储环境后重测。

#### 重点案例 B：streaming markdup

建议将 markdup 作为“跨批状态与有界内存”的代表案例，流程为：

```text
首坐标探测
  → CPE 批量提取 pair/single candidate
  → coordinate-owner sparse remote exchange
  → contiguous flat hash + sliding coordinate window
  → duplicate ordinal bitmap
  → CPE rewrite/pack/compress
```

可重点讲：

- coordinate owner 代替全量 hash-owner all-to-all，大多数 candidate 留在本 rank；
- sparse remote pack 只发送真正跨 rank 的 candidate；
- contiguous open-addressing flat hash 避免 `unordered_map` 节点分配和链式访问；
- 滑动坐标窗口淘汰不再可能匹配的 candidate，使 workspace 不随全部记录线性增长；
- duplicate ordinal 使用 bitmap 保存，最终 rewrite 时修改 duplicate flag；
- 双缓冲使下一批 CPE extract 与上一批 MPE merge 重叠；
- 与默认路径输出一致的测试中，workspace 从约 333.12 MB 降到 106.63 MB，降低约
  68%，核心时间还从约 0.815s 降到 0.787s。

这组“更低内存且没有性能代价”的结果比全内存 dedup-pipeline 更适合作为论文亮点。

#### sort 与 fixmate 的定位

- sort 用于证明公共框架可以承载 sample partition、all-to-all 和 k-way merge，并
  支持内排与严格外排。
- fixmate 用于证明框架可以处理跨 block/rank 的小规模边界状态，以及 prefix、middle、
  suffix 三段输出。
- 二者是完整性很强的案例，但不必与 collate、markdup 平均分配篇幅。

## 5. 系统架构图设计

### 图 1：SWBAM 总体架构图

这是全文最重要的图，建议分成四层，并使用不同颜色：

```text
┌──────────────── Applications ────────────────┐
│ flagstat stats filter convert                │
│ sort collate fixmate markdup                 │
└─────────────────────┬────────────────────────┘
                      │ generic / fused operator
┌──────────── Heterogeneous Runtime ───────────┐
│ CPE read pipeline │ CPE write pipeline       │
│ double buffering  │ codec/parser/cache       │
└─────────────────────┬────────────────────────┘
                      │ aligned BGZF batches
┌──────────── Distributed Runtime ─────────────┐
│ block index │ MPI rank plan │ body prefix    │
│ collective status      │ distributed output  │
└─────────────────────┬────────────────────────┘
                      │ backend interface
┌──────────────── Storage Backends ────────────┐
│ memory │ POSIX │ MPI-IO │ memory/spool sink  │
└──────────────────────────────────────────────┘
```

图中应特别画出两条竖向路径：generic path 和 fused path。不要把全部类名放进去，
图注再给出类名映射。

### 图 2：MPE/CPE 双缓冲时间线

横轴表示时间，纵轴分别画 MPE 和 CPE。至少展示三个 batch，并标出：

- MPE read batch `i+1`；
- CPE decode/parse/operator batch `i`；
- join、validate、buffer swap；
- overlap 开启和关闭时的关键路径差异。

这张图对应 overlap 消融实验，不能只画图而不测数据。

### 图 3：generic 与 fused 数据流对比

左侧 generic：

```text
decode → decoded batch → record views → MPE consumer
```

右侧 fused：

```text
decode + parse + operator on CPE → compact result → MPE reduction
```

箭头宽度可表示主存流量，突出 fused path 避免大规模 decoded payload 二次搬运。
若没有硬件计数器，图中不要标具体带宽节省百分比，只描述数据路径。

### 图 4：collate 分布式数据重排

画 3 个或 4 个 MPI ranks 即可：输入 records 经 QNAME hash/bin 映射到 owner rank，
self-owned records 直接保留，remote records 打包交换，最终在 owner 上 merge。
图中标注 raw arena direct receive，解释该优化为什么减少二次拷贝。

### 图 5：streaming markdup 状态窗口

横轴画坐标推进，展示 active window、expired candidates、rank boundary 和 remote
candidate。旁边画 flat hash 与 duplicate bitmap。该图应回答：

- 为什么旧路径需要保存大量 candidate；
- 为什么窗口淘汰不破坏正确性；
- 跨 rank 边界怎样处理；
- 为什么只发送 sparse remote candidates。

### 图 6：实验结果图

建议正文保留四张结果图：

1. 各应用 SWBAM 与 SAMtools 时间柱状图，并在柱顶标加速比；
2. `np=1/2/4/6` 吞吐量或 strong scaling 折线图；
3. generic/fused、overlap、库化前后消融柱状图；
4. 默认/stream markdup 的时间与峰值内存双轴图。

完整去重流程可使用一张 stacked bar，分别显示 collate、fixmate、sort、markdup。

## 6. 论文结构与各章内容

### 1 Introduction

按四段写：

1. BAM 是生物信息流程的重要中间格式，数据规模和处理种类持续增长。
2. 神威提供大量 CPE，但层次化存储、主从核协同和 MPI 通信使通用 BAM 工具难以
   直接获得高性能。
3. 现有 block-parallel 方法难以同时解决代码复用、逐记录抽象开销和复杂全局算法。
4. 引出 SWBAM，并列出三项贡献。

Introduction 中不要放过多命令实现细节，也不要把“实现了很多功能”本身当创新。

### 2 Background and Motivation

建议包括：

- BAM record 与 BGZF block 基础；
- 神威 MPE/CPE、LDM、DMA/主存和 MPI 进程模型；
- 三类 BAM workload：线性、全局重排、有状态；
- motivation experiment：展示解压、MPE record scan、通信和排序等热点。

本章最后提出设计目标：可复用、低抽象开销、可融合、后端可替换、有界内存。

### 3 SWBAM Design

依次说明：

1. 总体分层和 BGZF batch 边界；
2. 输入 backend 与 MPI block plan；
3. MPE/CPE 双缓冲 read/write pipeline；
4. generic/fused 双路径；
5. RankBodySink 与 DistributedBamOutput；
6. 错误传播、资源生命周期和内存策略。

这一章回答“库是什么、为什么这样分层”，不展开某个算法的 hash key 细节。

### 4 Application Mapping and Optimizations

先用表格列出所有应用，再重点展开：

- collate：bin ownership、exchange、raw arena、loser tree；
- streaming markdup：coordinate owner、flat hash、window、bitmap；
- sort：sample/partition/exchange/merge；
- fixmate：name group、rank boundary、rewrite。

flagstat/stats/filter/转换可合并成“线性算子实例”一节，用于说明 API 如何使用。

### 5 Experimental Methodology

必须明确：

- 机器型号、MPE/CPE 数量、MPI ranks、编译器和优化参数；
- SWBAM、SAMtools、htslib、libdeflate 版本；
- 数据集来源、压缩大小、解压大小、records、BGZF blocks、read length；
- 每个实验重复次数，建议至少 5 次，报告中位数及离散程度；
- core time 与 end-to-end wall time 的定义；
- page cache、tmpfs、高速盘、压缩等级和线程数控制；
- 正确性检查方法。

### 6 Evaluation

按问题组织，而不是按命令流水账组织：

- RQ1：公共数据平面的吞吐与扩展性如何？
- RQ2：batch 抽象是否引入额外开销？
- RQ3：generic/fused 与 overlap 分别贡献多少？
- RQ4：复杂应用相对 SAMtools 能获得多少加速？
- RQ5：streaming 和 external 模式能否控制内存及临时空间？
- RQ6：四阶段完整去重流程的总体收益如何？

### 7 Related Work

至少覆盖：

- SAMtools/htslib 的 BGZF、sort、collate、fixmate、markdup；
- RabbitBAM、RabbitBAM-QC、RabbitBAM-SORT；
- 多核/众核上的压缩、基因组格式处理与异构流水；
- MPI-IO 和并行文件处理系统。

重点比较执行模型和抽象边界，不要只比较支持的命令数量。

### 8 Limitations and Conclusion

主动说明当前边界：

- 当前实现假设一条 BAM record 不跨 BGZF block；
- 默认性能基线仍使用 memory backend；
- POSIX/MPI-IO 在真实高速并行存储上的实验尚不充分；
- dedup-pipeline 使用完整内存 BAM 阶段交接，不适合作为无限规模方案；
- SDK 目前主要服务项目内部和示例，还不是完整发布级软件包。

诚实说明限制不会削弱论文，反而能明确结论适用范围。

## 7. 实验设计

### 7.1 Library Microbenchmarks

| 实验 | 自变量 | 指标 | 目的 |
|---|---|---|---|
| BGZF decode | ranks、batch size | GB/s、blocks/s | 测公共读流水上限 |
| Raw BAM scan | ranks、batch size | records/s | 测 generic record view 成本 |
| BGZF encode | level、ranks | GB/s、压缩率 | 测公共写流水上限 |
| Read-write roundtrip | backend、ranks | records/s、wall time | 测完整库数据通路 |
| Input backend | memory/POSIX/MPI-IO | core、I/O、wall | 区分计算与存储瓶颈 |
| Output backend | memory/spool/MPI-IO | wall、峰值内存 | 验证有界输出 |

### 7.2 Library Design Evaluation

必须补齐以下受控消融：

| 消融 | 保持不变 | 只改变 |
|---|---|---|
| 重构前 vs 重构后 | kernel、数据、节点 | 公共库抽象 |
| overlap off vs on | operator、batch、backend | MPE/CPE overlap |
| generic vs fused | 过滤/统计语义、输入输出 | operator 路径 |
| vector/raw arena | collate 算法 | 初始化和数据布局 |
| full remote vs sparse remote | markdup 语义 | candidate pack |
| node hash vs flat hash | key、结果 | hash 结构 |

generic/fused 必须做同工作量对比。例如都实现相同的 MAPQ 过滤并生成同一 BAM，不能
把“只解压”的 generic 与“完整 flagstat”的 fused 时间直接相除。

### 7.3 Application Case Studies

建议至少包含：

- flagstat、stats --basic；
- BAM-to-SAM、SAM-to-BAM、BAM-to-BAM/filter；
- sort、collate、fixmate、markdup；
- 内存允许时的 memory 模式；
- 内存受限时的 sort/collate external 与 streaming markdup。

当前初步结果见 `docs/swbam初步性能结果.md`。正式实验必须统一双方的 I/O 和时间
口径，不能把 SWBAM core time 与 SAMtools 完整 wall time 直接作为最终加速比。

### 7.4 End-to-End Deduplication

主实验建议直接依次调用：

```text
collate → fixmate -m → sort → markdup
```

比较三种配置：

| 配置 | 定位 | 是否放正文 |
|---|---|---|
| SAMtools 四阶段独立命令 | 基线 | 是 |
| SWBAM 四阶段独立命令 | 主要结果 | 是 |
| SWBAM 全内存 dedup-pipeline | 消除中间 I/O 的上界 | 可选/附录 |

主表同时报告：

- 每阶段 wall time；
- 四阶段总 wall time；
- SWBAM 内部 core time breakdown；
- 中间文件总字节数；
- 峰值内存与临时空间；
- 最终输出正确性。

当前已有阶段时间求和约为 SWBAM `2.799s`、SAMtools `10.235s`，初步为 `3.66x`；
但双方计时口径可能不同，只能作为趋势，不能直接成为论文最终数字。

### 7.5 Correctness and Resource Evaluation

正确性至少检查：

- record count 和 mandatory fields；
- flagstat 全部计数、stats SN 和 sort state；
- collate QNAME continuity、READ1/READ2 order、record multiset；
- sort coordinate order 和 record multiset；
- fixmate mate fields/tags；
- markdup duplicate count 与 duplicate key set；
- header 忽略合理的 `@PG` 差异后比较；
- 正常输入、空 body、block 少于 rank 和损坏 BGZF 的错误行为。

资源指标至少包括：

- 每 rank 与全局峰值内存；
- output body memory/spool 占用；
- sort/collate external 临时空间；
- 随输入规模增长的内存曲线；
- streaming markdup 的 active candidate 数和 bitmap 大小。

## 8. dedup-pipeline 的论文定位

### 当前实现本质

当前 pipeline 是阶段级内存集成，而不是算子级深度融合：

```text
完整输入 BAM 广播
  → collate 生成完整压缩 BAM
  → 广播/刷新
  → fixmate 重新解压并生成完整 BAM
  → sort 重新解压并生成完整 BAM
  → markdup 重新解压
```

它减少了中间磁盘文件和进程启动，但没有消除阶段间压缩、解压和完整 BAM 复制。
内存峰值可能包含当前 BAM、下一阶段 BAM、广播刷新副本和 stage workspace。

### 为什么不建议作为主要贡献

- 应用范围受单节点/每 rank 可用内存限制；
- 阶段耦合后错误定位困难；
- 当前核心时间与四个独立阶段求和接近，深度融合收益不明显；
- 会分散论文对公共库、generic/fused 和 streaming markdup 的叙述。

### 推荐处理

保留代码和正确性回归，但在论文中定位为：

> An optional in-memory stage integration that estimates the performance upper bound when
> intermediate storage overhead is removed.

正文完整流程使用四个独立 SWBAM 命令。这样更符合真实使用方式，也能逐阶段验证和
定位性能。若版面不足，pipeline 可放附录或完全不报告，不影响论文核心贡献成立。

只有在后续实现“分布式中间态直接交接，并真正消除阶段间压缩解压”后，才适合把融合
pipeline 升级为核心创新。

## 9. 创新主张与证据对应表

| 主张 | 必须提供的机制证据 | 必须提供的实验 |
|---|---|---|
| 库抽象开销很低 | batch 级接口、不进入 record 热路径 | 重构前后 A/B |
| generic/fused 兼顾易用与性能 | 两条路径共享 runtime 的代码结构 | 同语义 generic/fused |
| MPE/CPE overlap 有效 | 双缓冲时间线和计时点 | overlap on/off |
| 支持复杂分布式算法 | sort/collate/fixmate/markdup 映射 | 应用性能与扩展性 |
| markdup 内存有界 | sliding window、bitmap、spans 多遍读取 | 内存随规模曲线 |
| 后端可替换 | memory/POSIX/MPI-IO、sink/output 解耦 | 高速存储 backend 对比 |
| 输出正确 | 语义 fingerprint 与边界测试 | correctness 表格 |

论文中每个“高效、可扩展、低开销、内存有界”的形容词，都应该能指向该表中的具体
机制和实验。

## 10. 向导师展示的建议材料

建议准备 8 页以内的讨论稿：

1. **问题与目标**：为什么神威上的 BAM 处理不能只靠并行 BGZF 解压。
2. **总体架构图**：四层架构和 BGZF batch 边界。
3. **创新一**：backend、MPI plan、sink、distributed output。
4. **创新二**：generic/fused 双路径与 MPE/CPE overlap。
5. **创新三**：collate 和 streaming markdup 两个算法案例。
6. **已有结果**：应用初步加速比、markdup 内存下降、库化开销。
7. **正确性与限制**：已验证内容、block-aligned 假设、慢盘问题。
8. **待讨论事项**：论文定位、目标会议/期刊、需要补哪些数据集和机器实验。

三分钟口头介绍可以这样组织：

> 我不是只实现了若干 BAM 命令，而是抽象出了一套面向神威的 BAM 批处理库。它以
> BGZF batch 为边界统一内存、POSIX、MPI-IO 输入，MPI block 分区，CPE 双缓冲和
> 分布式输出；对新功能提供零拷贝 record-view 通用路径，对热点功能保留解压、解析、
> 计算融合的 CPE 路径。sort、collate、fixmate 和 streaming markdup 证明这套架构
> 不只适用于简单统计，还能承载全局重排与跨批状态。现有初步结果显示多个应用相对
> SAMtools 获得加速，库化开销低于约 0.3%，streaming markdup 还将 workspace 降低
> 约 68%。下一步需要统一 I/O 口径、补扩展性和消融实验，形成正式论文证据。

## 11. 下一步优先级

### 第一优先级：冻结系统和建立可复现实验

- 固定代码 commit、编译参数、SAMtools/htslib/libdeflate 版本；
- 建立统一 benchmark 脚本和结果 CSV；
- 每项至少运行 5 次，记录中位数、最小值和波动；
- 区分 `core time`、`I/O time` 和 `end-to-end wall time`。

### 第二优先级：补最关键证据

- overlap on/off；
- generic/fused 同语义 A/B；
- `np=1/2/4/6` 扩展性；
- 至少 3 个规模的数据集；
- memory 与真实高速 POSIX/MPI-IO；
- stream markdup 内存随输入规模变化。

### 第三优先级：完善论文表达

- 根据上述设计绘制图 1 至图 5；
- 整理每个命令的算法伪代码和复杂度；
- 明确哪些优化来自公共库，哪些属于算法实例；
- 对照 RabbitBAM 系列工作撰写差异表；
- 与导师确定目标投稿方向后，再决定实验规模和英文写作深度。

### 暂不优先

- 继续大改全内存 dedup-pipeline；
- 为了命令数量增加边缘功能；
- 在没有统一基线前继续追逐单次小幅性能波动；
- 先做动态库、安装包等与核心论文证据关系较弱的工程工作。

## 12. 当前最需要导师确认的问题

1. 论文主定位是否采用“神威 BAM 处理库 + 复杂算法实例”，而不是单个 markdup 工具。
2. 目标投稿更偏高性能计算、并行系统，还是生物信息软件应用。
3. 是否能够获得高速并行文件系统节点完成正式 POSIX/MPI-IO 实验。
4. 是否需要增加更多真实 WES/WGS 数据集及不同 read length、重复率的数据。
5. “BAM record 不跨 BGZF block”的输入假设是否可以作为明确适用范围接受。
6. dedup-pipeline 是否只作为补充结果，正文使用独立四阶段完整流程。

这些问题确定后，再冻结论文标题、实验矩阵和最终贡献措辞。
