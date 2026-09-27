# SWBAM 论文实验方案

更新：2026-09-27，按 `a47c668` 的独立 SDK/app 和四流水线架构整理。
本方案替代旧 Raw view / MPE `bam1_t` 四项微基准及两路径消融设计。
历史测量保留在 results 中，不作为当前版本的实测结果。

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
| D3 | WGS | 89,994,058,306（89.99 GB） | 1,210,633,370 | 其他适用功能及外排，不测内排 |

数据来源和现有 SAMtools 环境见 [D1/D2/D3 库基准](results/samtools_D1_D2_D3_库基准.md)。
目录名 `D2-WGS-25G`、`D3-WGS-85G` 是历史标签，不代表当前实际 GB 大小。
正式实验固定输入校验和、accession、参考基因组、记录数、BGZF block 数和压缩参数；
数据准备时间不计入命令处理时间。

| 项目 | D1 | D2 | D3 |
|---|:---:|:---:|:---:|
| SDK count / BAM read-write | 是 | 是 | 容量允许时 |
| 转换、过滤、flagstat、stats --basic | 是 | 是 | 是 |
| fixmate、markdup --stream | 是 | 是 | 是 |
| sort/collate 内排及内排四步流程 | 是 | 是 | 否 |
| sort/collate 外排及外排四步流程 | 否 | 否 | 是 |
| 消融 | 主要实验 | 关键项复核 | 否 |
| 扩展性 | 完整 | 资源允许时复核 | 否 |

D1/D2 的 sort/collate 日志必须确认 memory mode，不能将自动 spill 混入内排表格。
D3 的 sort/collate 只测 external mode；其他命令不因此自动变成“外排算法”。

## 3. 版本、计时与资源口径

### 3.1 固定配置

- 记录 commit、二进制 SHA-256、编译器、MPI、HTSlib/libdeflate 版本和完整命令。
- 固定 level 1、过滤条件、算法选项；保持 block-aligned record 输入约束。
- 每点至少 3 次，建议升序、降序、升序；正式表格取中位数，保留各次值及范围。
- 同平台扩展性以最小配置为基线，效率为 `T1 / (p * Tp)`。
- 神威节点数、MPI ranks、每 rank 64 CPE 分别列出，不等同于 SAMtools `-@`。

### 3.2 两种现有时间分别标注

| 对象 | 当前记录时间 | 排除项 |
|---|---|---|
| SAMtools | bigssd 输入/输出的进程 wall time | 数据准备；不主动清理 page cache |
| SWBAM SDK | `core`，各 rank 流水线墙钟时间的最大值 | load、scan/plan、初始 body reserve、最终统计归约；示例无最终 gather/dump |
| SWBAM app | 各命令现有核心时间及 4.3 | 完整输入 load、最终真实 dump；具体额外排除项以日志和源码为准 |

SDK count 的 `core` 包含内存 batch 读取、CPE 解压/解析/计数及 MPE 批次归并。
SDK bam2bam 的 `core` 包含内存 batch 读取、解析、pack、压缩和本 rank 内存 body 写入。
它不含完整输出 header/EOF、rank 汇集及全文件输出布局，不能等同于 app 的 4.1~4.6。

现阶段按已有部署测量，不要求为了本轮实验切换 SAMtools 到 tmpfs。两种时间边界在
图表中分别写明；如给出 `SAMtools wall / SWBAM core`，标为该口径下的耗时比，
不称为同 I/O 条件下的端到端加速比。也不要混用早期 fat 和当前 A4 的结果。

### 3.3 大数据容量预检查

memory backend 每 rank 保存完整输入，不能用“文件大小 / ranks”估计输入内存。
还需容纳本 rank 输出、算法工作区；完整 memory output 在 rank0 汇集时还需整份输出。
现有大共享配置约 85 GiB 可用堆不意味着约 90 GB 的 D3 能连同输出和工作区装入。

D3 先核算容量，再选择 POSIX/MPI-IO 输入和 spool/MPI-IO 输出，或另行实现并验证
rank-local preload。不同 I/O 策略单列，不能静默切换。外排 run/temp 与 rank body
spool 是独立开销，记录各自字节数、真实时间、模拟时间；不能把必要临时 I/O 当作
初始 load/最终 dump 一并删除。当前计划不宣称全部内存状态与输入规模无关。

## 4. 公共库微基准：两项

| 编号 | 当前示例 | 实际数据路径 | SAMtools 对照 |
|---|---|---|---|
| L1 只读 | `swbam-sdk-optimized-count` | BAM read -> CPE 解压 -> 复用 `bam1_t` 解析 -> count -> MPE 归约 | `view -c` |
| L2 读写 | `swbam-sdk-optimized-bam2bam` | BAM read -> passthrough 解析 -> MPE pack -> BAM write -> CPE 序列化/压缩 -> memory body | `view -O BAM,level=1` |

L1 不是完整 flagstat，也不是只读压缩字节或仅扫描记录长度；L2 不过滤、保留全部记录。
当前 count 不向用户返回 `bam1_t` batch，不能用它证明旧 HTSlib 对象消费接口的性能。
源文件位于 `swbam-sdk/examples/`，两者均只接收输入文件参数，默认 memory 模式：

```bash
./swbam-sdk-optimized-count "$INPUT"
./swbam-sdk-optimized-bam2bam "$INPUT"

samtools view --no-PG -@ "$THREADS" -c "$INPUT"
samtools view --no-PG -@ "$THREADS" -O BAM,level=1 -o "$OUTPUT" "$INPUT"
```

报告 `records`、`blocks`、`input_bytes`、`core`；L2 另报 `body_bytes`。
GB 使用十进制：输入 GB/s = `input_bytes / core / 1e9`；M records/s =
`records / core / 1e6`。L2 可补充 `(input_bytes + body_bytes) / core / 1e9`，
明确后者只统计输出 body，不称为磁盘带宽。

### 4.1 与 app 的关系

SDK 不依赖 app，复用与 app 相同的流水线和 operator。性能诊断可比较相同
`RunBamTransformPipeline` 的范围，而不能拿 SDK `core` 与 app 全核心时间直接
声称“库化开销为零”。SDK count 与 flagstat 工作量不同，不能当作等价对照。

SDK L2 当前不 dump 完整文件，只检查保留数。正确性另跑 `io-check` 的无过滤
回环（独立解压、CRC/ISIZE、记录边界、总字节数及 MAPQ 校验）和 app 输出检查。
当前 `io-check` 不逐字节对比输入输出记录；完整字段一致性须另做记录比较，不能仅因
`records` 正确就通过。
诊断包含额外校验，不作为库吞吐量成绩。

## 5. 消融：优先验证调度设计

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
| flagstat、stats --basic | 原始 BAM | 全部 flagstat 项；对应 SN 项及排序状态 |
| BAM2SAM、SAM2BAM、BAM2BAM | 对应固定 BAM/SAM | header、record multiset、顺序 |
| MAPQ >= 30 过滤 | 原始 BAM | 保留集合、kept/dropped |
| sort、collate | 同一原始 BAM | 坐标顺序 / QNAME 连续性、记录集合 |
| fixmate -m | 同一 name-collated BAM | mate 字段和 MQ/MC/ms |
| markdup --stream | 同一 fixmate+sort BAM | duplicate count、标记集合、其他字段 |

SAMtools 完整 `stats` 工作量大于 SWBAM `stats --basic`；可核对对应 SN，但性能
表必须注明功能范围，不能声称是全功能等价加速。sort/collate 明确实际内排或外排。

### 6.1 配置选择

- 小回归：单节点 `-N 1 -np 6 -cgsp 64 -share_size 14000 -cache_size 32`。
- D1 库扩展性沿用已验证的大共享配置，每节点一个 rank，64 个 CPE 参与计算：

```bash
-N "$N" -np 1 -mpecg 6 -cgsp 64 \
  -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
```

节点数 `1/2/4/8/16/24/32/48`，以程序打印的实际 MPI ranks 校验放置。
其余核组提供内存，不计为活跃计算核组。D2 先做容量检查与代表点，再决定是否复核
全曲线；D3 不做扩展性。SAMtools 可沿用当前 A4 的 `-@ 0/2/4/8/16/24/32/47/48/64/96/127`
基线；额外线程数不等于实际 CPU 总占用数。

应用扩展性优先选择 sort、collate、markdup 和四步流程；线性应用用代表点确认即可，
避免每个命令都铺满资源矩阵。根据 D1 预实验预先确定主表配置，不能只挑最快一次；
最佳吞吐、资源消耗和并行效率分别报告。

## 7. 普通四步去重流程

```text
input BAM -> collate -> fixmate -m -> coordinate sort -> markdup --stream
```

主实验使用四条独立命令或 `dedup-workflow`，不依赖全内存 `dedup-pipeline`。
D1/D2 两个重排阶段使用内排；D3 两个阶段使用外排。阶段中间文件是工作流输入输出，
不等同于排序内部 run/temp。

每阶段保留 4.3、命令核心时间、load/dump 和实际输出大小。四个核心时间之和称为
“四步核心时间合计”，不是端到端 wall；不把 SAMtools 全流程 wall 和这个合计混作
相同边界。完整 wall 可另外记录，但不是本轮 memory-core 主指标。

至少验证最终可解析、记录数、坐标有序、duplicate 数和 duplicate 标记集合。
不同 collate 组间顺序可能导致输出排序和标签排列不同；SAM 文本 MD5 不一致时，
进一步核对 record multiset、mandatory fields、AUX 值，不能直接忽略差异。
同参数确定性回归优先用完整 BAM `cmp`，并明确参照的是哪个已验证基线。

## 8. 正式开测前的最小回归

当前版本已完成一轮，见 [SDK 与四步去重回归结果](results/最终版_SDK与四步去重回归.md)。

使用 `WES_0.25G.bam`，单节点 6 ranks，每项一次即可，不据此宣称稳定加速比：

1. 构建当前 app、SDK count 和 SDK bam2bam；记录 commit 和二进制校验和。
2. 两个 SDK 均处理 2,668,351 条记录，读写例保留全部记录。
3. `io-check` 验证无过滤回环及 MAPQ30 输出的可解析性、计数和字节数。
4. 重新依次执行四条命令，不能拿先前运行的阶段输出冒充完整回归。
5. `stats --basic` 检查最终记录数与 sorted 状态；与已验证完整输出比较重复标记。
6. 结果、日志和实际命令记录在 docs/results，失败时先修复再开始大数据实验。

## 9. 图表与后续顺序

1. 两项 SDK 的 D1/D2 吞吐量及 rank 扩展性图，D3 可运行配置单列。
2. 输入/输出 overlap 消融图及对应流水线时间线；无开关实验前不填写收益。
3. 应用性能表，明确平台资源、时间边界和功能范围。
4. 四步流程阶段堆叠图，展示核心时间合计及各阶段占比。
5. D3 外排的 run/pass、临时空间、实际/模拟 I/O 和峰值内存表。

顺序：小回归 -> D1 两项 SDK 与代表配置 -> 受控消融 -> D1/D2 全功能与四步流程
-> 容量可行时 D3 其他功能及外排。遇到性能瓶颈再按分项优化，不为论文预先恢复
已退役接口或继续大规模重构。
