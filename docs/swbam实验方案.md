# SWBAM 论文实验方案

更新：2026-10-12，采用私有分页 RAM 存储与 MPI 区间端到端计时。
论文主实验暂只覆盖 D1/D2；D3 已有结果保留为历史/探索记录，不自动续提。
本方案替代旧 Raw view / MPE `bam1_t` 四项微基准及两路径消融设计。
历史测量保留在 results 中，不作为当前版本的实测结果。
当前实现的构建、正确性、RAM 外排与 D1 单次探测记录见
[RAM 存储与端到端重构验证](results/RAM_存储与端到端重构验证.md)。

**2026-10-11 晚间固定本轮范围：** D1/D2 两项 SDK 加简单 rank 扩展性，
八项独立 app 加普通四步串联。app 固定六节点六活跃核组；每个正式点三次。
SDK 扩展性建议 `1/2/4/6/12`，两数据集使用相同点；完整矩阵须待容量探测后启动。
2026-10-12 已启动 D1 正式半矩阵，66 条性能命令加 8 条输出复读，见
[D1 正式测试](results/RAM_D1_正式性能测试.md)。
D2 六 rank 单次容量/功能探测在 sort 分配失败后停止，状态及配置见
[D1/D2 正式范围与 D2 探测](results/RAM_D1_D2_正式范围与D2探测.md)。
2026-10-12 已按用户要求部署 D2 缺项关卡及有条件正式矩阵，见
[D2 关卡与矩阵](results/RAM_D2_容量关卡与正式矩阵.md)：先验证 sort 的 16G
存储池、markdup、双向转换和 SDK 单 rank 读写；通过且 D1 完成后自动提交。
版本、输入、计数、排序、重复数或内排条件不符均停止，不自动落盘。
03:52，D2 缺项容量关卡已全部通过；10:58 检查时仍等待 D1。
D1 SDK 输出复读 `8506221` 在完成预加载后异常停滞，尚无成功结果，
两份正式矩阵不能视为完成。需先处理此阻塞，不能继续沿用原先结束时间估算。
11:21 已取消旧复读并换节点重试 `8513953`，D1 原地续接且保留已完成成绩；
验证开启既有 stats 阶段诊断并单独设置超时监控，性能命令和冻结程序不变。
11:24 新复读完整通过，D1 已自动续测两 rank SDK count；D2 仍按计划等待
D1 完整矩阵结束，不把本轮 D2 探测样本导入正式三次中位数。
15:00，D1 SDK 两项五点三次已齐，六核组中位数为 9.738284 / 15.483451 s；
整批完成 64/74 命令，但第三轮独立 sort 预加载后疑似停滞，仍需处理。
D2 正式矩阵尚未启动；不能因 SDK 已齐而宣称 app/三轮完整流程也已完成。
15:24，按用户要求取消第三轮停滞 sort、换六个空闲节点重试 `8519141`；
已完成的 64 条记录保持，程序、算法及性能参数不变，独立超时监控防止无限等待。
15:30 新 sort 成功，第三轮 MPI 区间为 57.178582 s、三次中位数 57.240036 s；
D1 已继续第三轮 markdup，D2 将在 D1 全矩阵完成后自动提交。
**16:08:53，D1 正式矩阵全部完成。** 30 条 SDK、24 条独立 app、12 条工作流
阶段及 8 条完整复读通过；三轮四步逐轮合计中位数为 **218.419689 s**。
16:09:17，D2 自动启动正式矩阵；16:55 已完成单 rank SDK count/读写首轮，
MPI 区间分别 48.453229 / 151.375398 s，SDK 输出完整复读通过，正在测试两 rank count。
暂按正式启动后 12-18 小时预留，预计 10 月 13 日 04:10-10:10，调度/存储波动另计。
本轮不增加 app 扩展性、消融或大规模外排性能矩阵；外排仅保留正确性回归。

## 1. 实验要回答的问题

1. 公共库能否高效完成 BAM 读取和读写，并支撑独立 SDK 应用？
2. 批处理与 MPE/CPE 输入预取、输出 flush 重叠是否带来可测收益？
3. 转换、过滤、质检及复杂算法在库的基础上能达到怎样的性能？
4. 普通 `collate -> fixmate -m -> sort -> markdup --stream` 串联效果如何？

架构主线是“四条类型化流水线 + 可组合 operator + app 算法”，不是两套独立
generic/fused 模式。算子的解压、解析和动作可在一次 CPE kernel 中静态组合。
两项 BAM SDK 微基准只覆盖 BAM read/write；SAM 两条流水线通过格式转换实验覆盖。

## 2. 数据集与实验矩阵

| 数据集 | 类型 | 实际压缩字节数 | 记录数 | 用途 |
|---|---|---:|---:|---|
| D1 | WES | 9,891,133,304（9.89 GB） | 177,798,571 | 除外排外的全部实验，主要扩展性 |
| D2 | WGS | 35,255,807,945（35.26 GB） | 365,079,492 | 除外排外的全部实验，规模复核 |
| D3 | WGS | 89,994,058,306（89.99 GB） | 1,210,633,370 | 历史/探索，不纳入本轮主实验 |

数据来源和现有 SAMtools 环境见 [D1/D2/D3 库基准](results/samtools_D1_D2_D3_库基准.md)。
目录名 `D2-WGS-25G`、`D3-WGS-85G` 是历史标签，不代表当前实际 GB 大小。
正式实验固定输入校验和、accession、参考基因组、记录数、BGZF block 数和压缩参数；
数据准备时间不计入命令处理时间。

| 项目 | D1 | D2 | D3 |
|---|:---:|:---:|:---:|
| SDK count / BAM read-write | 是 | 是 | 历史/探索 |
| 转换、过滤、flagstat | 是 | 是 | 历史/探索 |
| fixmate、markdup --stream | 是 | 是 | 历史/探索 |
| sort/collate 内排及内排四步流程 | 是 | 是 | 否 |
| sort/collate 外排及外排四步流程 | 本轮仅小数据正确性 | 本轮不测外排性能 | 历史/探索 |
| 消融 | 本轮不提交 | 本轮不提交 | 否 |
| 两项 SDK 简单扩展性 | 是 | 是 | 否 |

D1/D2 的 sort/collate 日志必须确认 memory mode，不能将自动 spill 混入内排表格。
D3 历史 sort/collate 是 external mode；其他命令不因此自动变成“外排算法”。
2026-10-09 调整：stats 不纳入论文性能比较，不再提交 A4 stats 性能测试；
`stats --basic` 等正确性复读仍保留。collate 及四步流程中的 collate
统一使用 66 个 bins（Samtools `-n 66`，SWBAM `-n 66`），旧 `-n 64`
成绩仅作历史记录，不能改写参数后直接沿用。

## 3. 版本、计时与资源口径

### 3.1 固定配置

- 记录 commit、二进制 SHA-256、编译器、MPI、HTSlib/libdeflate 版本和完整命令。
- 固定 level 1、过滤条件、算法选项；保持 block-aligned record 输入约束。
- 每点至少 3 次，建议升序、降序、升序；正式表格取中位数，保留各次值及范围。
- 同平台扩展性以最小配置为基线，效率为 `T1 / (p * Tp)`。
- 神威节点数、MPI ranks、每 rank 64 CPE 分别列出，不等同于 SAMtools `-@`。

### 3.2 当前主口径与历史结果

| 对象 | 当前主时间 | 排除项 |
|---|---|---|
| SAMtools | A4 bigssd 上进程 wall time | 数据准备；不主动清理 page cache |
| SWBAM app/SDK | 各 rank 的 MPI_Init 前至 MPI_Finalize 后真实耗时，取最大值 | MPI 前 RAM 池预触碰/输入预加载/就绪同步；Finalize 后实际保存 |
| SWBAM 完整处理分项 | 各 rank 的 MPI_Init 返回后至 MPI_Finalize 调用前，独立取最大值 | 上述模拟准备/保存，以及 MPI_Init、MPI_Finalize 两次调用 |

SWBAM 使用 `--storage ram` 显式启用内存模拟，默认真实磁盘。正常 header 解析、
扫描分区、MPI/athread 初始化退出、算法工作区管理、RAM 输入、外排临时读写、
完整输出构建和清理均计入。直接采用 `[swbam-e2e rank=...] seconds=...`，不扣减
真实临时文件时间，也不累加各阶段最大值。SDK 的 `core`、app 的 `4.3` 和已有
核心总时间仅作分项诊断；历史表格不能改名为新端到端成绩。

2026-10-13 增加完整处理分项 `[swbam-processing rank=...]`，
保存为 `processing_max_s`。它包含 athread 初始化/退出、头解析、扫描分区、
通信、算法、完整最终 RAM 输出和工作区清理；不以 `core` 替代它。
各 rank 用同一单调时钟独立取时，不新增测量 collective；结束后打印并由驱动取 max。
MPI 区间和完整处理的最大值可能来自不同 rank，禁止对两个全局最大值套用
rank 0 的初始化时间或统一减去约 7 秒。正式比较同时标明三种口径。
四步按每轮四条命令各自的处理 max 求和，再取三轮中位数；原 MPI 区间流程合计保留。

这是 RAM 存储条件下的 MPI 应用端到端时间，不是 SSD 实测，也不是作业提交到结束
的时间：不含调度排队、launcher 创建进程、预加载和保存。A4 的 bigssd 与神威 RAM
存储条件不同，图表必须标注，不称为相同 I/O 条件下的端到端加速比。

### 3.3 容量与启动

私有 `benchmark-support` 不安装进 SDK。用固定分页池保存本 rank 完整输入分片、
实际外排临时文件和完整输出；`--ram-capacity` 必须显式指定，算法 `-m` 另计。
二者共同占用约 85 GiB 跨核组堆，容量不足报错，禁止落盘或环形覆盖输出。
2026-10-13 起，RAM 与磁盘共用最终输出循环：完整读取 rank body 暂存，
按分布式偏移写入最终文件，不聚集整份文件到 rank 0。RAM 可回收已读完的
源页来降低峰值，但不能以登记页引用替代实际读写。此前页引用版成绩保留为
历史版本，正式比较不得混用两种最终输出路径。

启动脚本提供 `SWBAM_NODE_MEMORY=88000M`，程序按同节点 rank 数、RAM 池预留
和运行时余量自动计算工作预算，日志记录 `[work-memory]`。可选 `--work-memory`
覆盖每 rank 工作预算；sort/collate/markdup 的 `-m` 仍可单独固定，fixmate `-m`
仍是 mate-score 标志。预算是应用级限制而非操作系统 RSS 上限。
正常数据性能测试保持相同算法预算；不为边界缓存增加必填参数或极端数据矩阵。

BAM 依旧按 BGZF 块数分片，SAM 按完整行分片。MPI 前只存原始头、扫描原始信息和
本地数据；MPI 内重新解析、校验、分区，不用准备结果替代正常处理。
六核组配置为六节点、每节点一个 MPI rank、每 rank 64 CPE；其余五核组只贡献内存：

```bash
-N 6 -np 1 -mpecg 6 -cgsp 64 \
-share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
```

启动命令提供匹配的 `--ram-ranks 6` 和唯一绝对 `--ram-session` 路径；父目录需存在，
会话目录必须不存在。准备失败/超时会阻止进入 MPI。默认保存输出发生在 MPI 结束后，
`--no-dump` 仅用于已验证的独立性能命令，仍构建完整输出。

先做小数据和 D1 六核组单次探测，确认容量、正确性及核心时间回退，再决定正式矩阵。
两份完整矩阵已按用户请求授权：D1 已提交，D2 必须先通过缺项容量/正确性关卡，
且等待 D1 完成后再自动提交正式样本。D2 指定关卡现已通过，
但其余 SDK 资源点和新版完整四步的正式成绩尚未取得。
不覆盖历史实验目录或将失败样本记为成功。

## 4. 公共库微基准：两项

| 编号 | 当前示例 | 实际数据路径 | SAMtools 对照 |
|---|---|---|---|
| L1 只读 | `swbam-sdk-optimized-count` | BAM read -> CPE 解压 -> 复用 `bam1_t` 解析 -> count -> MPE 归约 | `view -c` |
| L2 读写 | `swbam-sdk-optimized-bam2bam` | BAM read -> passthrough 解析 -> MPE pack -> BAM write -> CPE 序列化/压缩 -> memory body | `view -O BAM,level=1` |

L1 不是完整 flagstat，也不是只读压缩字节或仅扫描记录长度；L2 不过滤、保留全部记录。
当前 count 不向用户返回 `bam1_t` batch，不能用它证明旧 HTSlib 对象消费接口的性能。
源文件位于 `swbam-sdk/examples/`，默认磁盘模式；写示例必须指定输出路径。RAM 示例命令体如下：

```bash
./swbam-sdk-optimized-count "$INPUT" --storage ram --ram-capacity 2G \
  --ram-ranks 6 --ram-session "$UNIQUE_COUNT_SESSION"
./swbam-sdk-optimized-bam2bam "$INPUT" -o "$OUTPUT" \
  --storage ram --ram-capacity 8G --ram-ranks 6 \
  --ram-session "$UNIQUE_COPY_SESSION"

samtools view --no-PG -@ "$THREADS" -c "$INPUT"
samtools view --no-PG -@ "$THREADS" -O BAM,level=1 -o "$OUTPUT" "$INPUT"
```

报告 `records`、`blocks`、`input_bytes`、`core` 和 MPI 区间主时间；L2 另报完整输出大小。
GB 使用十进制，主图吞吐量以 MPI 区间主时间为分母。核心吞吐量可另列诊断，不能混用。

### 4.1 与 app 的关系

SDK 不依赖 app，复用与 app 相同的流水线和 operator。性能诊断可比较相同
`RunBamTransformPipeline` 的范围，而不能拿 SDK `core` 与 app 全核心时间直接
声称“库化开销为零”。SDK count 与 flagstat 工作量不同，不能当作等价对照。

SDK L2 默认保存完整文件，检查可解析性和记录集合；性能重复可显式 `--no-dump`，
但仍须构建并驻留完整 RAM 输出。正确性另跑 `io-check` 的无过滤
回环（独立解压、CRC/ISIZE、记录边界、总字节数及 MAPQ 校验）和 app 输出检查。
当前 `io-check` 不逐字节对比输入输出记录；完整字段一致性须另做记录比较，不能仅因
`records` 正确就通过。
诊断包含额外校验，不作为库吞吐量成绩。

## 5. 消融：优先验证调度设计

本节仅为后续备选，不进入当前 D1/D2 补测任务；先补齐正式主成绩。
不再测试已删除的 Raw/MPE `bam1_t` 双路径，也不为消融恢复旧实现。

| 优先级 | 实验 | 保持不变 | 唯一改变 |
|---|---|---|---|
| 优先 A1 | BAM read 输入重叠 | count kernel、解析语义、batch 容量 | `overlap_input` off/on |
| 优先 A2 | BAM write 输出重叠 | bam2bam kernel、pack、压缩级别、输出 sink | `overlap_output` off/on |
| 可选 A3 | batch 容量 | 同一算子和相同总输入 | 在现有容量约束内选择批大小 |

A1/A2 选 D1 的同一固定资源点，D2 复核最有效的一项。先固定后处理预取与 source
预取策略，再单独改变所测开关，避免同时改变多个重叠来源。一次核组只能运行一个
CPE kernel；这里是 MPE 预取/flush 覆盖 CPE 工作，不是读写 CPE 同时执行。

**实施状态：尚未开放示例命令行消融开关。** `CpeReadPipelineOptions` 已有
`overlap_input/prefetch_during_post_process`，`CpeWritePipelineOptions` 已有
`overlap_output/overlap_source`；当前高层 count/transform wrapper 未全部透传这些
参数。正式消融前新增小型受控配置入口、验证结果一致，再测量，不能伪造现有 CLI。
A3 也需先确认 operator 的固定容量约束，不是任意改一个参数即可成立。

中间步骤全部沿用同一 kernel；read/kernel/post/flush 的区间可能重叠，主结论看
流水线 wall，不把重叠区间直接相加。可选补充：同输入的单遍/两遍 markdup；只有
当前版本提供可验证选择入口时才执行。跨历史版本比较只能作为重构回归证据。

## 6. 应用与扩展性

| 功能 | 固定输入 | 验证重点 |
|---|---|---|
| flagstat | 原始 BAM | 全部 flagstat 项 |
| BAM2SAM、SAM2BAM | 对应固定 BAM/SAM | header、record multiset、顺序 |
| MAPQ >= 30 过滤 | 原始 BAM | 保留集合、kept/dropped |
| collate | 原始 BAM | QNAME 连续性、记录集合 |
| sort | 同一 fixmate BAM | 坐标顺序、记录集合和 mate 标签保持 |
| fixmate -m | 同一 name-collated BAM | mate 字段和 MQ/MC/ms |
| markdup --stream | 同一 fixmate+sort BAM | duplicate count、标记集合、其他字段 |

stats 仅用于需要的正确性检查，不进入论文性能图表；A4 复读可使用
`quickcheck` 和 `view -c`，SWBAM 保留 `stats --basic` 复读。
sort/collate 明确实际内排或外排。

### 6.1 配置选择

- RAM 小回归：两节点各一个 rank，沿用大共享配置，并显式指定 `--ram-ranks 2`。
- D1 库扩展性沿用已验证的大共享配置，每节点一个 rank，64 个 CPE 参与计算：

```bash
-N "$N" -np 1 -mpecg 6 -cgsp 64 \
  -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
```

SDK 节点数建议 `1/2/4/6/12`，以程序打印的实际 MPI ranks 校验放置。
其余核组提供内存，不计为活跃计算核组。D2 先做六 rank 容量检查，并补单 rank
读写容量探测，再启动与 D1 一致的各扩展性点；D3 已有扩展性仅作历史参考，不续提。
当前版本要求保留完整 RAM 输出，
不能使用旧分块环形覆盖模式扩充容量。
SAMtools 可沿用当前 A4 的 `-@ 0/2/4/8/16/24/32/47/48/64/96/127`
基线；额外线程数不等于实际 CPU 总占用数。

当前八项独立 app 和普通四步流程只测六核组，不增加 app 扩展性矩阵。
主表固定每点三次，不能只挑最快一次；SDK 扩展性另报告加速比、效率和吞吐量。

## 7. 普通四步去重流程

```text
input BAM -> collate -n 66 -> fixmate -m -> coordinate sort -> markdup --stream
```

RAM 主实验使用四条独立命令。两个 dedup 包装命令保留磁盘模式，本轮不扩展其 RAM 模式。
D1/D2 两个重排阶段使用内排；强制外排仅作本轮小数据正确性回归。阶段中间文件是工作流输入输出，
不等同于排序内部 run/temp。

前三步必须在 Finalize 后保存中间文件，下一条命令重新预加载。每阶段保留 MPI 区间
最大值、4.3 诊断、准备/保存耗时和输出大小。每轮将四个 MPI 区间最大值相加，三轮
合计取中位数；同时保留各阶段中位数。不能把历史四步核心时间之和改名为新端到端成绩。

至少验证最终可解析、记录数、坐标有序、duplicate 数和 duplicate 标记集合。
不同 collate 组间顺序可能导致输出排序和标签排列不同；SAM 文本 MD5 不一致时，
进一步核对 record multiset、mandatory fields、AUX 值，不能直接忽略差异。
同参数确定性回归优先用完整 BAM `cmp`，并明确参照的是哪个已验证基线。

## 8. 正式开测前的最小回归

新版已完成一轮小数据回归，见 [RAM 存储与端到端重构验证](results/RAM_存储与端到端重构验证.md)。
当前 RAM 重构单独记录，不能沿用历史回归来宣称新存储路径正确。

使用 `WES_0.25G.bam`，两节点两 ranks，每项一次即可，不据此宣称稳定加速比：

1. 构建当前 app、SDK count 和 SDK bam2bam；记录 commit 和二进制校验和。
2. 两个 SDK 均处理 2,668,351 条记录，读写例保留全部记录。
3. `io-check` 验证无过滤回环及 MAPQ30 输出的可解析性、计数和字节数。
4. 重新依次执行四条命令，不能拿先前运行的阶段输出冒充完整回归。
5. `stats --basic` 检查最终记录数与 sorted 状态；与已验证完整输出比较重复标记。
6. 结果、日志和实际命令记录在 docs/results，失败时先修复再开始大数据实验。

## 9. 图表与后续顺序

1. 两项 SDK 的 D1/D2 RAM MPI 区间吞吐量及 rank 扩展性图。
2. 八项应用性能表，明确平台资源、时间边界和功能范围。
3. 四步流程阶段堆叠图，展示逐轮 MPI 区间合计及各阶段占比。
4. 正确性与容量记录：实际 RAM 临时读写字节数和存储/工作区峰值。
5. overlap 消融与额外外排性能实验仅作后续备选，不进入本轮补测。

顺序：小回归 -> D1 单次代表探测 -> D2 六核组容量探测 -> 冻结版本及容量 ->
D1/D2 两项 SDK 简单扩展性、八项六核组 app 和三轮四步流程
（D1 及通过关卡后的 D2 正式矩阵已授权）。遇到性能瓶颈再按分项优化，不为论文预先恢复
已退役接口或继续大规模重构。
