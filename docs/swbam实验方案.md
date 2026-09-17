# SWBAM 论文实验方案

## 1. 实验目标

实验围绕以下四个问题组织：

1. SWBAM 公共库的 BAM 读取与读写通路是否高效？
2. MPE/CPE 异构流水、批处理接口和融合算子分别带来多少收益？
3. 基于公共库实现的统计、转换、重排和去重应用相对 SAMtools 表现如何？
4. `collate -> fixmate -m -> sort -> markdup` 完整流程能获得多少总体收益？

论文正文以可复现的独立命令和 file-backed `dedup-workflow` 为主。全内存
`dedup-pipeline` 只作为可选补充，不作为主要性能结论。

## 2. 数据集与实验分工

| 数据集 | 类型 | 目标压缩规模 | 主要用途 |
|---|---|---:|---|
| D1 | WES | 约 9 GB | 全功能、消融和扩展性实验 |
| D2 | WGS | 约 25 GB | 全功能主性能实验和 D1 结论复核 |
| D3 | WGS | 约 85 GB | 大规模流式处理、外排与资源实验 |

正式实验前记录每个数据集的 accession、参考基因组、样本类型、压缩字节数、解压
字节数、记录数、BGZF block 数、平均读长、配对比例及 header MD5。输入数据生成过程
不计入任何命令运行时间。

### 2.1 数据集使用原则

- **D1 和 D2 测试除外排以外的全部项目**：公共库、消融、统计、转换、过滤、
  fixmate、markdup、sort 内排、collate 内排和内排完整流程。
- D1、D2 不运行 sort/collate external，不把自动 spill 的结果混入内排表格；日志必须
  确认实际选择的是 memory mode。
- **D3 测试公共库和其他适用功能，并且 sort/collate 只测试 external mode**。
- D3 不运行 sort/collate memory mode，也不用于主要扩展性实验。
- D3 的完整流程使用 external collate、streaming fixmate、external sort 和 streaming
  markdup，作为大规模、有界内存能力验证。

### 2.2 总体实验矩阵

| 实验项目 | D1 | D2 | D3 |
|---|:---:|:---:|:---:|
| Raw BAM read | 是 | 是 | 是 |
| Raw BAM read-write | 是 | 是 | 是 |
| `bam1_t` read | 是 | 是 | 是 |
| `bam1_t` read-write | 是 | 是 | 是 |
| 主消融实验 | 全部 | 关键项复核 | 否 |
| 统计、转换、过滤 | 是 | 是 | 是 |
| fixmate、streaming markdup | 是 | 是 | 是 |
| sort/collate memory | 是 | 是 | 否 |
| sort/collate external | 否 | 否 | 是 |
| 内排完整去重流程 | 是 | 是 | 否 |
| 外排完整去重流程 | 否 | 否 | 是 |
| 扩展性 | 主要数据集 | 可选复核 | 否 |

## 3. 统一实验口径

### 3.1 软件与参数

- 固定 SWBAM commit、编译器、MPI、HTSlib 和 libdeflate 版本。
- 固定 SAMtools 版本及编译参数，优先与 SWBAM 使用相同 HTSlib/libdeflate 版本。
- 双方统一输出格式、压缩级别 1、过滤条件、排序语义和功能选项。
- 记录节点数、MPI rank 数、每 rank CPE 数、CPU 线程数、内存限制和临时目录。
- SAMtools 的 `-@` 表示附加压缩线程，报告结果时同时写明命令和实际资源配置。

### 3.2 I/O 与时间

SWBAM 当前主要使用 memory backend 模拟高速存储，正式表格应称为
**memory-resident processing** 或“内存驻留处理”，不能直接称为真实磁盘吞吐量。

建议同时保留两类时间：

| 时间 | 用途 |
|---|---|
| core time | 分析算法、CPE kernel、通信、压缩和重排本身 |
| comparable wall time | 与 SAMtools 计算正式加速比 |

公平比较时，SAMtools 输入和输出应放在 `tmpfs` 或同等级高速存储；预加载、首次拷入
内存以及最终归档时间均排除。D3 外排产生的 run/temp I/O 属于算法必要成本，必须计入
wall time，且双方必须使用同一临时存储。

当前 memory backend 可能在每个 rank 保存完整输入。D2/D3 正式测试前必须确认单 rank
内存和节点总内存足够；若不足，应使用按 rank 读取的 POSIX/MPI-IO backend 或实现
rank-local memory preload。不同 backend 的结果必须分表报告，不能静默切换后仍标为
同一 memory 实验。

### 3.3 重复与统计

- 每种配置预热 1 次，不计入结果。
- D1、D2 至少重复 5 次；机器时间紧张时最低为 3 次。
- D3 至少重复 3 次；若成本过高，先跑 1 次功能验证，再完成正式 3 次。
- 主表报告中位数，同时保留最小值、最大值、标准差和原始日志。
- 加速比统一为 `SAMtools median / SWBAM median`。
- 吞吐量同时报告 GB/s 和 M records/s。GB/s 必须注明采用压缩输入字节、解压逻辑
  字节还是输入加输出总字节。

## 4. 公共库微基准

库实验只保留“读”和“读写”两类操作，并分别测试 Raw BAM 与 `bam1_t` 路径。

| 编号 | 测试 | SWBAM 数据路径 | 主要含义 |
|---|---|---|---|
| L1 | Raw BAM read | BGZF decode -> record boundary -> `RawBamRecordView` -> count | 零拷贝可编程读路径 |
| L2 | Raw BAM read-write | L1 -> `RawBamWriter` -> BGZF encode | Raw 身份转换完整闭环 |
| L3 | `bam1_t` read | L1 -> MPE batch materialize -> minimal callback | HTSlib 对象兼容读路径 |
| L4 | `bam1_t` read-write | L3 -> `Bam1Writer` -> BGZF encode | HTSlib 对象兼容读写闭环 |

### 4.1 工作量约束

- read 测试只做最小业务操作，例如累加 records 和 MAPQ，避免业务逻辑掩盖库开销。
- read-write 测试必须保留全部记录和字段，不做过滤；允许 BGZF 压缩字节不同，但解压
  后 record multiset 和 mandatory fields 必须一致。
- Raw 与 `bam1_t` 必须读取同一批输入、使用相同 batch size 和 compression level。
- 库级基线优先编写同语义 HTSlib 小程序：`sam_read1()` 计数，以及
  `sam_read1() + sam_write1()` 回环。SAMtools CLI 结果作为实际工具参考，不代替同语义
  SDK 基线。

### 4.2 报告指标

| 指标 | 说明 |
|---|---|
| wall/core time | 完整可比时间与 SWBAM 内部分解 |
| compressed GB/s | 输入 BAM 压缩字节数除以时间 |
| logical GB/s | 解压 BAM record 字节数除以时间 |
| M records/s | 记录数除以时间 |
| peak memory | 每 rank 最大值和节点聚合值 |
| output bytes | 读写路径压缩后大小 |

### 4.3 `bam1_t` 路径的判定与优化预案

`bam1_t` 当前在 CPE 完成 BGZF 解压和 raw record 边界识别后，由 MPE 批量物化标准
HTSlib 对象。因此它比 Raw 零拷贝路径多一次字段解释和 payload 复制，预期更慢。

正式实验先拆分 `decode/raw scan/materialize/callback/encode/compress`。满足以下任一情况
时，再考虑把部分解析下沉到 CPE：

- `bam1_t` read 明显慢于同语义 HTSlib 基线，例如中位数回退超过 10%～15%；
- `materialize` 占 `bam1_t` read 时间超过约 25%；
- 随 rank 数增加，MPE materialize 成为扩展性瓶颈。

优化顺序建议为：先减少 MPE 对象池和 payload 复制，再让 CPE 批量产生 mandatory-field
descriptor 和规范化 payload offset，最后才评估由 CPE 直接填充数据 arena。不要把
`bam1_t` 的分配、析构及 HTSlib 对象所有权整体搬入 CPE，这会显著增加 LDM 压力和
接口复杂度。无论是否优化，Raw 路径都应保留为最高吞吐路径，`bam1_t` 路径承担生态
兼容和易开发能力。

## 5. 消融实验

正文建议只保留三组主消融，均在 D1、`np=6` 上完成；D2 只复核最关键的一组，D3
不做消融。所有消融必须保证输入、输出语义、压缩级别和资源配置一致。

### A1：MPE/CPE overlap

| 对照 | 唯一变化 | 预期说明 |
|---|---|---|
| overlap off | 串行执行 read、CPE kernel、consume | 双缓冲前基线 |
| overlap on | MPE 读取下一批与 CPE 处理当前批重叠 | 证明异构流水收益 |

建议在 Raw read、Raw read-write 各测一次，并记录 read、kernel、consume、等待和总时间。

### A2：Programmable Batch Path 与 Fusion-Optimized Path

当前代码中的 generic/fused 在论文中建议改称：

- **Programmable Batch Path**：输出 Raw view 或 `bam1_t` batch，由 MPE 后处理；
- **Fusion-Optimized Path**：CPE 内融合 decode、parse 和具体 operator。

选择同一项业务，例如 `MAPQ >= 30` 过滤并写出 BAM，比较两条路径。两边必须输出
相同记录集合，不能用“只解压”对比“完整过滤”。该实验回答可编程接口付出了多少
性能，以及融合算子能收回多少开销。

### A3：Raw view 与 `bam1_t` materialization

在相同 count/MAPQ 统计和相同 read-write 回环上比较 Raw 与 `bam1_t`。报告
materialize/encode 的额外时间和内存。这不是为了证明 `bam1_t` 应该取代 Raw，而是
量化“HTSlib 易用性”的代价，并为是否进行 CPE 解析优化提供依据。

### 5.1 可选应用级消融

版面和时间允许时，再增加以下一项，不必全部执行：

| 可选实验 | 对照 | 价值 |
|---|---|---|
| streaming markdup | 原两遍 stream vs 单遍滑动窗口 | 展示二次解压消除和峰值内存变化 |
| collate receive | 二次 append vs direct arena receive | 展示通信数据布局优化 |
| 库化开销 | 重构前 tag vs 当前库实现 | 证明 batch 抽象没有明显回退 |

已有的历史调优数字可以用于动机说明，但正式消融最好通过同一 commit 的开关重复，避免
把机器波动或其他代码变化当成单项优化收益。

## 6. 应用性能实验

### 6.1 测试功能

| 类别 | 功能 | 独立输入要求 | 主要指标 |
|---|---|---|---|
| 统计 | flagstat、stats `--basic` | 原始 BAM | time、records/s、计数一致性 |
| 转换 | BAM-to-SAM、SAM-to-BAM、BAM-to-BAM | 对应 BAM/SAM | time、GB/s、输出大小 |
| 过滤 | MAPQ 等固定条件 | 原始 BAM | time、kept records |
| 重排 | sort、collate | 原始/坐标 BAM | time、通信量、内存、temp bytes |
| 配对修复 | fixmate `-m` | 固定 name-collated BAM | time、tag/字段一致性 |
| 去重 | markdup `--stream` | 固定 fixmate coordinate BAM | time、内存、duplicate set |

独立命令实验使用预先生成且固定的阶段输入。例如 markdup 不把前面的 collate、fixmate
和 sort 时间算入自身；这些成本只在完整流程实验中累加。

### 6.2 D1、D2 配置

- 测试上表全部功能。
- sort 和 collate 强制 memory mode，并记录估计工作集和实际峰值内存。
- 不创建 external runs，不报告外排数字。
- D1 用于完整功能、消融和扩展性；D2 用于主性能表和大一档规模复核。

### 6.3 D3 配置

- 统计、转换、过滤、fixmate 和 markdup 仍然测试，使用流式或有界内存 backend。
- sort 和 collate 强制 external mode，固定每 rank 算法内存上限和 temp 目录。
- 记录 run 数、merge/consolidation pass、临时读写字节、临时空间峰值和最终峰值内存。
- 不测试 sort/collate memory，也不根据运行时内存情况自动切换模式。

### 6.4 扩展性

为控制实验数量，扩展性集中在 D1。建议选择：

- 库：Raw read、Raw read-write；
- 线性应用：flagstat 或 BAM-to-BAM；
- 全局应用：sort、collate；
- 状态应用：streaming markdup；
- 完整去重流程。

神威侧优先测 `np=1/2/3/6` 的单节点核组扩展性，固定每 rank 使用 64 CPE；如论文要
主张多节点可扩展性，再补 `1/2/4` 节点的代表性项目。SAMtools 测
`1/8/16/24/48` 线程或按机器物理核心数设置的等距档位。扩展性图同时报告 speedup、
parallel efficiency 和 throughput，不能只给最快配置。

应用主性能表统一采用双方各自通过预实验确定的最佳固定配置；配置确定后不得按数据集
临时挑选最快一次。若不做完整扩展性，至少用 D1 做一次资源扫参，说明为什么最终采用
6 ranks/384 CPE 和对应的 SAMtools 线程数。

## 7. 完整去重流程

主流程为：

```text
collate -> fixmate -m -> coordinate sort -> markdup
```

推荐使用 `dedup-workflow` 统一编排四个独立实现，并保留每阶段计时；也可以由脚本逐条
调用。正文不依赖全内存 `dedup-pipeline`。

### 7.1 数据集配置

| 数据集 | collate | fixmate | sort | markdup |
|---|---|---|---|---|
| D1 | memory | streaming | memory | streaming |
| D2 | memory | streaming | memory | streaming |
| D3 | external | streaming | external | streaming |

### 7.2 对比方案

| 方案 | 定位 |
|---|---|
| SAMtools 四阶段流程 | 正式基线 |
| SWBAM 四阶段/dedup-workflow | 论文主要结果 |
| SWBAM dedup-pipeline | 可选的内存阶段交接上界，不进入主加速比 |

### 7.3 报告内容

- 四个阶段各自时间及总时间；
- 总体加速比和各阶段对总体收益的贡献；
- 峰值内存、压缩中间 BAM 总量和 external 临时空间；
- D3 的 run 数、merge pass 与 temp I/O；
- 最终记录数、mandatory-field fingerprint 和 duplicate-key fingerprint；
- 各阶段输出的语义正确性。

如果神威真实磁盘过慢，D1/D2 的主结果可以采用双方等价的 memory-resident 时间，并把
四阶段 core time 相加；但必须明确这是计算路径性能，不是端到端存储性能。D3 external
依赖临时文件，其 I/O 不能从算法时间中删除，应在相同高速临时存储上重新测量。

## 8. 正确性与资源验证

| 功能 | 最低正确性检查 |
|---|---|
| Raw/`bam1_t` read | records、mapped、duplicate、MAPQ、NM 汇总 |
| read-write/转换/过滤 | record multiset、mandatory fields、保留集合 |
| flagstat/stats | 全部计数项、SN 项及 sorted 状态 |
| collate | QNAME continuity、READ1/READ2 order、record multiset |
| sort | coordinate order、record multiset |
| fixmate | mate flags/coordinates/TLEN 及 `MQ/MC/ms` |
| markdup/workflow | duplicate count、duplicate-key set、mandatory fields |

所有实验记录每 rank 峰值内存和全局最大值。D3 额外验证峰值内存不随文件总规模线性
增长，并记录输出 spool 与算法 external temp 是两类独立空间。

## 9. 推荐图表

1. 四项库微基准的 D1/D2/D3 吞吐量分组柱状图。
2. Raw 与 `bam1_t` 的阶段时间堆叠图，突出 materialize/encode 成本。
3. overlap 和 programmable/fusion 两组主消融柱状图。
4. 各应用相对 SAMtools 的加速比图，D1、D2 并列，D3 单独标识 external。
5. D1 的 rank/线程扩展性折线图。
6. D3 external sort/collate 的时间、峰值内存和临时空间图。
7. 完整去重流程的阶段堆叠柱状图。

## 10. 执行顺序

1. 冻结代码、环境、数据集和正确性工具。
2. 用 D1 完成四项库微基准、三项消融和资源扫参。
3. 用 D1、D2 完成除 external 外的全部应用和内排完整流程。
4. 检查 `bam1_t` breakdown，决定是否进行 CPE-assisted materialization 优化。
5. 用 D3 完成流式应用、external sort/collate 和外排完整流程。
6. 汇总原始 CSV、重复运行统计、正确性结果和论文图表。

这一顺序先用较小数据发现实现和口径问题，再投入 D3 的昂贵外排实验；同时确保即使
D3 资源暂时不可用，D1/D2 也能形成一套完整的论文主体结果。
