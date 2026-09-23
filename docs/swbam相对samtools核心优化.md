# SWBAM 相对 SAMtools 的核心优化说明

## 1. 文档目的

本文用于回答两个问题：

1. SWBAM 公共库相对 `HTSlib + SAMtools` 软件栈做了哪些面向神威的优化？
2. `RabbitBAM-MPI` 中各应用为什么可能比对应 SAMtools 命令更快？

这里的“相对 SAMtools”主要表示实现路线和数据通路的差异。SAMtools 本身已经支持
HTSlib BGZF 线程池和多线程压缩；SWBAM 的重点不是简单地把一个串行程序并行化，而是
把 BGZF、BAM 解析、应用算子、MPI 数据重分布和分布式输出一起映射到神威
MPE/CPE 架构。

已有加速比来自早期测试记录，只能作为初步依据。SWBAM 部分数据主要是排除初始文件
加载和最终 dump 的核心时间，已有 SAMtools 数据可能是内存驻留文件上的 wall time，
二者尚未全部统一口径，不能直接作为论文最终结论。

## 2. 一句话概括

> SWBAM 以 BGZF block batch 为并行和扩展边界，在 MPI rank 间分配压缩块，在每个
> rank 内使用 64 个 CPE 执行解压、解析和应用计算，并通过双缓冲、零拷贝 raw record、
> 专用融合算子及分布式输出减少中间数据、内存分配和 MPE 串行工作。

测试配置使用 6 个 MPI rank，即 6 个 MPE 和 384 个 CPE。

## 3. 与 SAMtools 的总体差异

| 方面 | HTSlib/SAMtools 常见路径 | SWBAM 路径 | 主要目的 |
|---|---|---|---|
| 并行层次 | 单进程或 pthread，BGZF 线程池和部分算法多线程 | MPI rank + 每 rank 64 CPE | 利用神威跨核组和核组内并行 |
| 调度粒度 | 以 BAM record 或 HTSlib block/job 为主 | 固定数量 BGZF blocks 组成 batch | 批量启动 CPE，摊薄调度开销 |
| 记录表示 | 通常经 `sam_read1()` 构造/复用 `bam1_t` | 简单路径使用 raw BAM view，复杂路径只抽必要 metadata | 减少对象构造、指针追踪和复制 |
| 应用计算 | 读取记录后在 CPU 上执行命令逻辑 | CPE 可融合 decode、parse、filter/count/extract/rewrite | 避免 decoded 数据回到 MPE 后再次扫描 |
| 全局算法 | 单机内存、线程和临时文件 | MPI sample/owner/bin 分区及 rank 间交换 | 分散排序、聚类和去重工作集 |
| 输出 | HTSlib 顺序写和 BGZF 线程压缩 | CPE 压缩、rank body sink、MPI-IO 分布式布局 | 压缩并行化并避免 rank0 汇总瓶颈 |
| 扩展接口 | HTSlib 是成熟通用格式库 | composable 与 optimized 两条路径共存 | 同时保留易开发性和高性能 |

## 4. 公共库级核心优化

### 4.1 Backend-neutral BGZF 数据面

公共接口把输入介质统一转换成 `BgzfBlockSpan` 和 `BgzfBlockBatch`。应用算法只消费
batch，不直接依赖完整内存文件、POSIX 或 MPI-IO：

```text
Memory / POSIX / MPI-IO
          ↓
BGZF block spans + aligned block batch
          ↓
CPE decode / parse / operator
```

当前实现包括：

- `MemoryBamInput`：完整文件载入内存，用于当前高速 I/O 模拟和核心性能基线。
- `PosixBamInput`：流式读取，不要求完整 BAM 常驻内存。
- `MpiIoBamInput`：MPI-IO 读取 rank 所负责的 blocks。
- `MpiBamInputPlan`：rank0 扫描 BGZF 索引并广播，各 rank 只处理自己的 block 区间。

相对逐条调用 `sam_read1()`，这种接口让 CPE 一次获得一批连续、对齐的压缩块，也让
统计、转换、sort、collate 和 markdup 复用相同的数据入口。

对应代码：`include/swbam/io.h`、`lib/io/`、`include/swbam/mpi_runtime.h`、
`lib/mpi/`。

### 4.2 大窗口 BGZF scanner 与批量读取

POSIX scanner 使用 8 MiB 窗口顺序预读，在内存中解析 BGZF header 的 `BSIZE`，避免
针对约 12100 个 blocks 分别执行一次小 `pread`。处理阶段使用 `preadv` 把连续 span
直接填入最终 block slots；MPI-IO 路径使用派生 datatype 将文件数据散布到同样的
slots，减少中间 scratch buffer。

优化对象是系统调用次数和中间复制，不改变 BGZF 格式语义。Memory backend 本身已经
在内存中，因此主要受益于顺序解析和统一 batch 描述。

### 4.3 MPE/CPE 双缓冲重叠

读流水线维护两套 compressed/decoded batch：

```text
时刻 t：    CPE 处理 batch i
             MPE 同时读取 batch i+1
时刻 t+1：  交换两套 buffer
```

写流水线采用对称设计，让 CPE 压缩当前 payload 时，MPE 消费上一批压缩块。这样减少
“读完再算、算完再写”的串行等待。抽象只在每个 batch 上调用虚接口，不进入逐记录
热循环。

对应代码：`lib/cpe/swbam_cpe_read_pipeline.cpp`、
`lib/cpe/swbam_cpe_write_pipeline.cpp`。

### 4.4 Composable path 与 Optimized path

SWBAM 没有强迫所有应用使用同一种抽象：

- **Composable 路径**：CPE 解压后，`RawBamRecordView` 直接指向 decoded batch 中的原始
  BAM record，不构造 `bam1_t`，适合 SDK 示例、简单扫描和过滤。
- **Optimized 路径**：专用 CPE kernel 在一次启动中完成
  `BGZF decode + BAM parse + application operator`，只回传计数、key、plan 或 bitmap
  等紧凑结果。

Composable 路径强调易扩展；optimized 路径强调性能。二者共享 backend、block batch、codec、
调度和错误处理，因此专用优化不会迫使每个命令重新实现整套 I/O。

这是相对通用 `sam_read1() -> bam1_t -> application loop` 最重要的数据通路区别之一。

### 4.5 Raw record、紧凑 metadata 与 arena

多个热路径不完整构造 `bam1_t`：

- 统计和过滤直接读取 mandatory fields。
- sort 只抽取 `tid/pos/flag/strand/global_order` 和 raw bytes。
- collate 只抽取 QNAME hash、QNAME、flag order、global order 和 raw offset。
- markdup 只交换判重所需的 candidate key、score、ordinal 和 QNAME。

临时 payload 使用 per-block/per-chunk arena 或一次性未初始化 raw arena，避免逐记录
`malloc/free`，也避免 `vector::resize()` 先把大块内存清零、随后又被 CPE 全量覆盖。
collate 的 raw arena 优化曾将 `raw_resize` 从约 2.732 s/rank 降到约 0.003 s/rank。

### 4.6 公共 CPE codec/parser 与对象缓存

公共从核层提供 BGZF inflate/deflate、CRC/ISIZE 校验和 raw BAM fast parser。每个 CPE
缓存 libdeflate compressor/decompressor，避免每批重复创建和销毁。压缩级别 1 还保留
针对 matchfinder 的快速策略。

对应代码：`slave/core/`、`include/swbam/cpe_codec.h`、
`include/swbam/cpe_bam_parser.h`。

### 4.7 RankBodySink 与分布式输出

每个 rank 的压缩 body 可以写入：

- memory：性能基线；
- spool：有界内存临时文件；
- auto：达到内存限制后自动转入 spool。

`DistributedBamOutput` 先计算各 rank body prefix，再选择 gather 到 rank0 或由 MPI-IO
直接写入最终偏移。MPI-IO 路径不要求 rank0 保存完整输出 BAM。

这项优化主要解决大输出的内存可扩展性和集中式汇总瓶颈；默认 memory 模式则保留
现有内存模拟性能。

### 4.8 库化没有替换成熟热路径

公共库只统一 backend、batch runtime 和输出布局。sort 的 bucket exchange、collate
的 QNAME merge、markdup 的 owner 判重等状态仍留在应用层；专用 CPE kernel 也没有
被 composable batch post-processor 替换。已有重构 A/B 结果中：

- flagstat：`0.052073 s -> 0.051994 s`，约快 0.15%；
- stats --basic：`0.098962 s -> 0.099174 s`，约回退 0.21%。

这组结果主要用于说明抽象放在 batch 边界后，库化本身没有明显增加热路径开销。

## 5. 各应用相对 SAMtools 的优化

### 5.1 BAM -> BAM 与过滤

SAMtools `view` 通常通过 HTSlib 解码为 `bam1_t`，在主 CPU 上判断条件，再交给 BGZF
线程池写出。SWBAM 的生产路径主要做了：

1. CPE 融合执行 BGZF 解压、BAM fast parse 和过滤条件判断。
2. passthrough 情况仍走融合 kernel，减少 MPE 二次扫描。
3. 使用 per-BGZF-block arena 保存变长 record data，避免每条记录固定 1 KiB 或独立
   分配。
4. MPE 使用连续 flat pack，CPE 复用 compressor 并行压缩。
5. 输入读取可与当前 CPE 解压重叠，输出可使用 memory/spool/MPI-IO。

Composable SDK 中还有 `RawBamFilterBatchPostProcessor + RawBamWriter`，但论文应用
性能应使用上述 Optimized path，不应把 SDK 示例性能当作生产命令上限。

### 5.2 BAM -> SAM

SAMtools `view` 在 CPU 上逐条解析并格式化 SAM 文本。SWBAM 使用两个 CPE 阶段：

1. CPE 解压并解析 BAM records；
2. CPE 批量执行 SAM 格式化，MPE 只负责收集输出片段和全局布局。

每个 CPE 有独立格式化缓冲；极端记录导致容量不足时只扩展对应 core 并局部重跑，
避免按最坏情况为所有 records 预留空间。初步测试为 `0.234 s vs 1.241 s`，约
`5.30x`，但仍需统一 wall-time 口径复测。

### 5.3 SAM -> BAM

SWBAM 先按换行边界把 SAM body 分给 MPI ranks，再按 chunk/subgroup 处理：

1. CPE 并行解析 SAM 行并生成 BAM raw data；
2. subgroup 根据实际记录量拆分，不受单个固定大池硬上限限制；
3. per-chunk arena 按需增长，超长 SAM 行进入 side buffer；
4. 生成的 raw records 按 BGZF payload 组织，再由 CPE 并行压缩；
5. rank body 最后通过公共分布式输出写出。

相对 SAMtools 通用逐行解析路径，优化重点是多 rank 文本切分、CPE parse 和 CPE
compression。当前尚缺统一的 SAMtools 性能数字。

### 5.4 flagstat

SAMtools `flagstat` 的核心是 `sam_read1()` 后在 CPU 上执行 `flagstat_loop()`。SWBAM
让每个 CPE 对一个 BGZF block 执行：

```text
decode -> raw record scan -> FLAG 分类计数
```

每个 CPE 写自己的计数 slice，MPE 每批只合并 64 份小数组，最后使用 MPI reduction
合并 ranks。完整 records 不回传 MPE。当前结果已与 SAMtools 对齐，核心时间约
`0.052 s`。

### 5.5 stats --basic

SAMtools `stats` 功能很完整，逐记录维护多类统计和分布。SWBAM 当前只对齐 `SN`
Summary Numbers，不能与完整 `samtools stats` 宣称功能等价。

针对 basic 子集，CPE 在 decode/parse 同一遍扫描中计算 mapped/unmapped、MQ0、NM、
bases、quality、insert-size、orientation 和 block 内坐标状态；MPE/MPI 只归并计数、
最大值、直方图和 rank 边界有序性。避免把 266 万条记录传回 MPE。当前核心时间约
`0.099 s`。

### 5.6 coordinate sort

SAMtools sort 在单机内把 `bam1_t` 累积到内存块，多线程排序后形成临时文件，最终
多路归并。SWBAM 使用分布式 sample sort：

1. CPE 解压后直接抽 coordinate key 和 raw bytes，不完整构造 `bam1_t`。
2. 每 rank 本地排序并抽样，rank0 选择 splitters。
3. 已排序 records 用 splitter 单调扫描确定 bucket。
4. CPE 并行 pack remote metadata/raw，MPI `Alltoall`/分块交换目标 bucket。
5. 每个来源 rank 发来的数据段本身有序，因此接收端执行 k-way merge，而不是对全部
   received records 再做一次完整 `std::sort`。
6. 排序后直接把 raw records 填充为 BGZF payload，由 CPE 压缩，避免
   `bam1_t -> raw BAM` 再序列化。
7. 外排模式按每 rank `-m` 形成有序 run，使用有界 ring exchange、匿名临时 extent
   和 loser tree 流式归并。

关键历史收益包括：CPE bucket pack 的总时间约 `13.45 s -> 4.62 s`，final sort 改
为 k-way merge 后约 `9.17 s -> 4.05 s`。WES_0.25G 初步结果为
`0.675 s vs 2.131 s`，约 `3.16x`。

### 5.7 collate

SAMtools 标准 collate 使用 X31-Wang QNAME hash 把 records 分入多个临时 BAM 文件，
再逐 bucket 排序输出；fast mode 使用额外的内存 hash/store。SWBAM 保持相同 hash
语义，但把 bins 分布到 MPI ranks：

1. CPE 融合执行 BGZF decode、QNAME/hash 和 compact metadata 提取。
2. 默认 bins 自动调整为接近 64 且可被 rank 数整除；6 ranks 时为 66，减少 owner
   负载偏斜。
3. 未初始化 raw arena 避免大规模 vector zero-fill。
4. self-owned records 直接 compact 到最终 segment，避免 send buffer 中转。
5. remote exchange 保留 bounded chunked ring，避免单次超大消息和超大临时 buffer。
6. 接收 raw bytes 直接写入预分配的最终 segment arena，去掉
   `recv_raw -> segment.raw` 二次复制；该改动曾将 `recv_append` 从约
   `1.36~1.42 s/rank` 降到约 `0.16 s/rank`。
7. memory merge 使用专用 loser tree、cursor/meta/raw pointer cache，并直接填充压缩
   payload，避免逐记录 `std::function` 调用。
8. external 模式保留 run、extent、consolidation 和文件 cursor，限制算法内存。

输出只要求同 QNAME 连续、READ1 位于 READ2 前以及 record multiset 一致，不要求
不同 QNAME groups 与 SAMtools 具有同一排列。WES_0.25G 初步结果为
`0.708 s vs 4.538 s`，约 `6.41x`。

### 5.8 fixmate -m

SAMtools fixmate 顺序读取 queryname grouped records，在 CPU 上按 template 更新 mate
字段和标签。SWBAM 的优化包括：

1. MPI ranks 分别处理 QNAME groups，显式合并跨 batch、跨 rank 的边界 group。
2. CPE plan kernel 计算字段/tag 修改和输出长度，CPE rewrite kernel 执行实际更新。
3. MPE 只构造 group 描述、prefix/output offsets 和少量边界状态。
4. 输出 records 使用 flat pack、CPE 压缩和公共 RankBodySink。

当前实现对齐 `samtools fixmate -m -z off` 的 mate 坐标、flags、TLEN、`MQ/MC/ms`；
功能选项尚未完整覆盖 SAMtools。WES_0.25G 初步结果为
`0.585 s vs 1.581 s`，约 `2.70x`。

### 5.9 markdup

SAMtools markdup 的主路径是坐标滑动窗口加 pair/single khash，在记录离开窗口时输出；
光学重复和 supplementary 传播可能引入额外状态或额外读写。SWBAM 当前默认高性能
路径采用两遍分布式算法：

1. 第一遍由 CPE 解压并提取 pair/single duplicate candidates，只保存判重所需字段。
2. 默认使用 coordinate owner，将 candidate 交给负责对应坐标区间的 rank；绝大多数
   candidate 留在本 rank。
3. 只对极少数 remote candidates 建立稀疏发送 buffer。已记录的 6-rank 样例中
   `mpi_candidate_bytes` 曾从旧 hash-owner 路径约 360 MB 降到 1793 B。
4. self-owned candidates 原地 compact，避免复制 candidate 和 QNAME。
5. owner 端使用预分配的 flat/open-addressing hash 完成 pair/single 判重，避免 STL
   node hash 的逐节点分配和 cache miss。
6. pair、paired-marker 和 real-single 分流，避免所有 marker 作为完整 candidate 参加
   同一个大排序/判重过程。
7. 判重结果以 bitmap/ordinal 返回；第二遍 CPE 重新解压并直接设置或删除 duplicate
   records，然后走 raw payload 压缩。
8. 普通 `-r` 已融合进第二遍，避免先生成完整 marked BAM 再执行一次 BAM2BAM 过滤；
   `-c -r` 也使用 raw flag/aux 清理和直接跳过 duplicate。
9. `--stream` 提供坐标窗口淘汰的有界候选状态，面向大文件 backend；它和默认两遍
   路径是不同的性能/内存折中。

流式路径测试中，workspace 峰值曾从约 333.12 MB 降到约 106.63 MB，同时核心时间
约 `0.815 s -> 0.787 s`。WES_0.25G 默认 markdup 初步结果为
`0.831 s vs 1.985 s`，约 `2.39x`。

当前 SWBAM markdup 不支持 SAMtools 的全部 optical duplicate、barcode/RG 隔离、
`dt/do` 生成和 supplementary duplicate 传播，论文必须按双方共同功能集比较。

### 5.10 完整去重流程

推荐正式论文主实验使用四个独立阶段：

```text
collate -> fixmate -m -> coordinate sort -> markdup
```

这样双方命令边界一致，也更容易逐阶段验证和定位性能。现有初步核心时间为：

| 实现 | collate/s | fixmate/s | sort/s | markdup/s | 合计/s |
|---|---:|---:|---:|---:|---:|
| SWBAM | 0.708 | 0.585 | 0.675 | 0.831 | 2.799 |
| SAMtools | 4.538 | 1.581 | 2.131 | 1.985 | 10.235 |

按现有数字计算约 `3.66x`，但需要在统一内存驻留输入、压缩级别、线程/核资源和完整
wall time 下重新测量。

`dedup-pipeline` 是实验性的全内存阶段编排：它通过 memory-to-memory BAM 避免真实
中间文件，但阶段间由 rank0 组装完整 BAM 并广播，内存峰值和调试复杂度较高。当前
更适合作为“融合潜力”补充实验，不应替代四个独立命令的主结果。

## 6. 哪些收益来自哪里

可以把性能来源归为四类：

| 收益来源 | 代表功能 | 典型效果 |
|---|---|---|
| CPE 计算并行 | 解压、压缩、格式化、统计、key 提取、rewrite | 减少 MPE 热点 |
| 算子融合 | flagstat、stats、filter、collate extract、markdup candidate | 减少 decoded 中间态和二次扫描 |
| 数据结构与内存布局 | raw view、flat metadata、arena、bitmap、loser tree cache | 减少分配、清零、复制和 cache miss |
| MPI 算法设计 | sort splitter、collate bin owner、markdup coordinate owner | 分散全局工作集并减少通信量 |

因此，论文不能只把结论写成“384 个 CPE 比 SAMtools 快”。更准确的说法是：

> SWBAM 通过异构并行、跨阶段融合、面向 BAM 编码的紧凑数据结构和算法感知的 MPI
> 分区共同获得性能；不同命令的主要收益来源并不相同。

## 7. 当前限制与论文表述边界

1. 当前允许假设单条 BAM record 不跨 BGZF block；必须在论文输入约束中明确说明，
   并在入口检测错误时失败，不能暗示支持任意合法 BAM。
2. 默认 memory backend 是高速 I/O 模拟，不能证明真实并行文件系统上的端到端收益。
3. sort/collate 有外排，流式转换、统计、fixmate、stream markdup 可使用 POSIX/MPI-IO
   和 spool，但这些真实存储路径仍需单独评估。
4. `stats --basic` 只是 SAMtools stats 的功能子集。
5. collate 不支持 SAMtools fast mode；markdup 和 fixmate 也未覆盖全部选项。
6. 输出 BAM 字节级 MD5 不一定一致；应比较 mandatory fields、record multiset、排序
   语义、QNAME continuity 和 duplicate key set。
7. 初步加速比的计时口径还不完全统一，正式论文必须报告双方相同资源下的完整
   wall time，并把 SWBAM core breakdown 作为额外分析。

## 8. 导师常见追问速答

### 这是不是只把 BGZF 解压放到从核？

不是。公共从核层负责 codec/parser；专用 kernel 还执行过滤、flagstat/stats 计数、
sort key 提取、collate QNAME/hash 提取、fixmate plan/rewrite 和 markdup candidate/rewrite。

### SAMtools 也支持多线程，区别是什么？

SAMtools 主要通过 HTSlib BGZF 线程池和部分命令的 pthread 并行利用通用 CPU。
SWBAM 同时使用 MPI 和 CPE，并让应用算子进入 CPE batch；复杂命令还设计了
splitter/bin/coordinate owner 等分布式算法，而不只是并行压缩。

### 为什么一定要做成库？

如果每个命令分别实现扫描、block buffer、CPE launch/join、错误处理和分布式输出，
代码会重复且难以增加 backend。库把这些公共机制放在 batch 边界；新工具可以使用
composable raw view 快速实现，热点明确后再替换为 optimized operator。

### Composable 路径比专用路径慢，库还有意义吗？

有。Composable path 的目标是开发效率和可复用性，Optimized path 承担最终性能。二者共享
同一 I/O/runtime，不要求在易用性与性能之间二选一。正式实验应同时报告 composable、
optimized 和 SAMtools，展示融合的收益及库抽象开销。

### 当前最快的几个关键优化是什么？

- collate：未初始化 raw arena、remote receive 直写最终 arena、专用 loser tree fill。
- sort：CPE raw/key extract、CPE bucket pack、接收有序段 k-way merge。
- markdup：coordinate owner、稀疏 remote candidate、self compact、flat hash 和 bitmap。
- 统计：decode/parse/count 融合，只回传小型统计数组。
- 公共路径：双缓冲 overlap、codec cache 和 batch backend。

### 完整去重实验为什么建议逐阶段执行？

它与 SAMtools 命令边界一致，易于验证每阶段正确性，也能明确加速来自哪一阶段。
全内存 `dedup-pipeline` 可以作为避免中间 I/O 的补充实验，但当前存在完整 BAM 广播
和较高峰值内存，不适合作为唯一主结果。

## 9. 正式实验前需要补齐的数据

- 相同节点、相同核资源、相同压缩级别下的 SAMtools 与 SWBAM wall time。
- `/dev/shm` 或已预热页缓存上的内存驻留对比，以及真实高速存储对比。
- `np=1/2/4/6` 扩展性和 MPE-only/CPE-enabled 对比。
- overlap on/off、composable/optimized、memory/posix/mpiio 的消融实验。
- 峰值内存、MPI 通信字节、sort/collate 临时空间。
- 至少三个数据规模和多种压缩率/记录长度数据集。
- 所有对比的共同功能选项与正确性 fingerprint。

完成这些数据后，本文中的“初步加速”才能升级为论文中的正式性能结论。
