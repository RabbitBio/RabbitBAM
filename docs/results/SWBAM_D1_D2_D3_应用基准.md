# SWBAM D1/D2/D3 应用基准

## D1 六核组测试（2026-10-01）

03:18 启动 D1 应用矩阵，首项 flagstat 作业为 **8367514**。本轮只测 D1，
计划 10 个应用各三轮，首轮额外执行四项输出统计检查；不修改应用或库代码。

| 项目 | 配置 |
|---|---|
| 版本 | `9a16bb3fbccb2a6593327eae73b8dc03cbc2b723` |
| 二进制 SHA-256 | `6732983a5014e0473728ad0f85f60c97ab9a2902c082be990385699a26932f81` |
| 固定节点 | `33449-33454`，六节点、每节点一个 MPI rank |
| 活跃计算资源 | 6 MPE + 384 CPE；每节点其余核组提供内存 |
| 输入 | `D1-WES-9G.coordinate.bam`，9,891,133,304 bytes、177,798,571 records |
| 输入 / rank body | `memory / memory` |
| 最终验证输出 | `mpiio`，在核心计时之后单独写盘 |
| BAM 压缩级别 | 1 |
| collate / sort | `-m 32G`，日志必须确认 `mode=memory` |
| markdup | `--stream -l 300 -m 8G` |

资源参数：

```bash
-N 6 -np 1 -node 33449-33454 -mpecg 6 -cgsp 64 \
  -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
```

每轮顺序为 flagstat、stats --basic、BAM2BAM、MAPQ >= 30 过滤、collate、
fixmate -m、sort、markdup --stream、BAM2SAM、SAM2BAM。四步去重顺序执行，
上一阶段的 BAM 作为下一阶段输入。collate 使用 `-n 66`；各轮覆盖同名输出，
保留最后一轮文件。若 SAM 大于 72 GiB，则单独跳过全内存 SAM2BAM 并记录原因，
不将其他输入策略混入本轮。

计时：有输出命令取 `Complete the total (4.1~4.6) cost`，另记 `4.3`；
flagstat、stats 取 `444Complete the total body cost`。物理 load、模拟输出
内存分配和最终 MPI-IO 验证写盘按应用现有口径排除，分别保留相关日志。
最终 MPI-IO 不代替核心阶段的模拟内存写，其用途是避免 rank0 汇集完整输出。
对照 A4 时使用 [Samtools 应用基准](samtools_D1_D2_D3_应用基准.md) 中的
bigssd wall time，表中应明确这是不同 I/O 时间边界的耗时比。
Samtools stats 与 SWBAM stats --basic 的统计范围也不同。

首轮对 BAM2BAM、过滤、SAM2BAM（容量允许时）、markdup 的输出执行
`stats --basic`，检查可解析、记录数及坐标排序；过滤输出记录数与 kept 统计核对。
这些检查不代替逐记录字段或 duplicate 标记集合的独立对照。

原始记录与完整提交命令：

```text
本地：build_sunway_sduhpc/benchmarks/apps_D1_n6_20261001_031800/
远端：online/guoshi/wzs/wzs_data/swbam_bench/D1_apps_current/apps_D1_n6_20261001_031800/
服务：swbam-apps-d1-n6-20261001-031800-resume.service
驱动：tools/bench_apps_d1.sh（运行时使用本地目录中的 driver.sh 冻结副本）
```

`results.tsv` 区分 bench 与 validation，保存核心/4.3/load/dump 时间、记录数、
输出大小和提交 wall；`commands.sh` 保存完整命令，`logs/` 保存逐项日志。
失败立即停止后续提交。完成时 `state` 为 `COMPLETE`，若有容量跳过则为
`COMPLETE_WITH_SKIPS`。正式性能表待三轮完成后汇总。

### 启动状态

首项 flagstat 作业 8367514 正常结束，记录数为 177,798,571：

| core / s | 4.3 / s | load / s | 进程总时间 / s | 提交到结束 / s |
|---:|---:|---:|---:|---:|
| 3.598642 | 3.354887 | 111.400876 | 122.087559 | 141 |

首项结束后，驱动因后台 PATH 不含 rg 以及汇总行的 CRLF 解析失败而暂停；
已修正脚本，恢复这条有效结果，未重跑或修改应用。03:22 继续提交 stats --basic
（8367591）。按首项 load 与已有同规模测试估计，其余三轮测试约需 2～4 小时，
主要不确定性是共享存储持续读写速度和 SAM 转换的容量；这不是核心性能预测。

03:23 stats --basic 已正常完成并被驱动自动收录：core=3.171133 s、
4.3=2.927941 s、load=70.184053 s，记录数正确且 `is sorted: 1`。
驱动随后自动进入 BAM2BAM，确认修正后的日志解析和连续提交正常。

### 03:30 Collate 失败后暂停

截至 2026-10-01 09:48 检查，驱动已退出，`state=FAILED rc=1 current_job=8367733`，
本批没有仍在运行的作业。正式矩阵仅完成首轮 **4/30** 项；四项额外输出验证尚未开始。

| 首轮应用 | 作业 | core / s | 4.3 / s | load / s | 最终写盘 / s |
|---|---:|---:|---:|---:|---:|
| flagstat | 8367514 | 3.598642 | 3.354887 | 111.400876 | 0 |
| stats --basic | 8367591 | 3.171133 | 2.927941 | 70.184053 | 0 |
| BAM2BAM | 8367624 | 9.267312 | 8.519630 | 71.405070 | 44.332816 |
| MAPQ >= 30 | 8367682 | 9.525545 | 8.807564 | 70.219963 | 41.325690 |

上述应用输入记录数均为 177,798,571，BAM2BAM 全部保留；过滤保留
170,359,591 条。输出大小分别为 11,467,682,880 和 10,899,175,329 bytes。
这些是首轮计数及运行自检结果，尚无三轮中位数或独立输出验证结论。

collate 作业 8367733 于 03:30 失败：

```text
MPI BAM collate selected mode=memory memory_estimate_max=10220195878 memory_limit=34359738368
ERROR: collate extract failed at global block 400359 status=-2 record=0.
```

该作业输入 load=70.102850 s，随后进入提取阶段并显式返回失败，未产生有效性能结果。
从核 `slave_mpi_collate_extract` 的 `-2` 同时表示解压器分配/解压/CRC 失败，以及
QNAME 长度或 NUL 终止符校验失败；`record=0` 不能单独判定是哪条分支。
因此当前不能直接断言服务器故障、输入损坏或已确定代码根因，也不是外排被误选。
下一步应隔离该块，核对压缩内容与从核错误分支，再决定修复和续跑；此前 2～4 小时
完工估计因失败暂停已失效。

### 10:12 根因诊断：Collate raw arena 并发写竞争

已在原出错节点 33451、使用原报错批次的数据复现代码缺陷：内排提取阶段
让不同 CPE 直接写入按 ISIZE 紧密拼接的 raw arena，块间没有对齐隔离。
神威编译器对部分未对齐写生成 `ldl_u` / 掩码合并 / `stl_u` 的机器字读改写。
相邻块逻辑上不重叠，但边界会共享一个机器字；并发写可能把邻块尾部覆盖成旧值。

对应代码为 `swbam-app/algorithms/collate/swbam_collate_mpi.cpp` 中
`CollateExtractOperator::Prepare()` 的 `raw_->data + raw_pos` 和
`raw_pos += isize`。这与 BAM record 跨 BGZF 块无关。

诊断证据：

- 原始块 400359 位于文件偏移 4,858,789,689，压缩长度 16,804 bytes，
  ISIZE=65,195，含 220 条记录。它及同批 64 块均通过独立 gzip CRC、
  QNAME 和完整 record 边界检查，没有发现输入块损坏。
- 原二进制在原节点处理这 64 块的小样本成功，说明错误不是必现的格式错误。
- 仅细分错误日志的全量 D1 诊断作业 8372379 也成功：records/received 均为
  177,798,571，4.3=58.297711 s，4.1~4.6=59.363980 s，load=101.056937 s，
  dump=71.770070 s。但该次没有独立验证输出，不作为正式成绩或正确性证明。
- 小样本压力测试只调整 CPE 启动时序及输出布局，不改变输入或提取算法。
  主核将从核产物逐字节与独立解压结果比较，能检出 CRC 检查之后才发生的覆盖。

| 诊断作业 | 紧凑并发布局 | 对齐并发布局 | 串行对照 |
|---|---|---|---|
| 8372439：遍历地址对齐偏移，扰动相邻 CPE 时序 | 4096 批中检出 2213 个坏块 | 4096 批、0 错误（独立 64 KiB 槽） | 8 批、0 错误 |
| 8372472：保持原批次起始偏移 `%64=16`，针对块 400359 | 512 批中检出 332 个坏块 | 512 批、0 错误（逐块 `align_up(ISIZE,64)`） | 8 批、0 错误 |

压力测试刻意放大竞争窗口，上述数字不是生产失败概率。第二组明确复现
块 400359 尾部从 offset 65193 开始被覆盖为预填的 `0xa5`，CRC 从预期
`faa92014` 变为 `81e1b265`。诊断版返回 `status=-23 record=0`（CRC 分支）；
原版同一分支返回 `status=-2 record=0`，与凌晨故障现象一致。
部分坏块的从核 status 仍为 0，说明邻核可能在本核完成 CRC 后才破坏尾部，
不能以偶然重跑成功代替修复。这两项 probe 检出坏块后主动返回非零，
调度器显示的失败是诊断预期结果，不是新的节点故障证据。

建议修复：raw arena 容量按逐块 `align_up(ISIZE,64)` 求和，各块起始地址按
64 字节对齐；保留实际 ISIZE / raw_len，元数据使用含 padding 的 raw_offset，
后续 pack 仅复制 record 有效字节。这样保留直接解压快路径，不新增整批复制。
D1 每 rank 最多增加约 8.2 MiB padding；性能影响需修复后实测。
同时应细分解压、CRC、QNAME 错误码，避免以后只看到通用 `-2`。

本次仅诊断，**未修改正式 SDK/App 源码或原基准二进制，未恢复矩阵**。
正式继续前先修复该竞争，再做输出正确性检查。

复现脚本、样本、诊断二进制和日志：

```text
本地：build_sunway_sduhpc/benchmarks/collate_D1_diagnosis_20261001/
远端：online/guoshi/wzs/wzs_data/swbam_bench/D1_apps_current/collate_diagnosis_20261001/
inspection.json / snapshot_probe.log / full_D1_diag.log
arena_probe.log / arena_stress.log / exact_probe.log
probe_alignment_sweep.cpp（4096 批）/ probe.cpp（512 批、原偏移）
```

### 修复与回归（2026-10-01）

已修改正式代码：

- `CollateEstimateRawBytes()` 的 memory/backend 两个入口均按
  `align_up(ISIZE,64)` 累计容量，`Prepare()` 按相同跨度分配块位置，包含最后一块的
  padding；实际解压容量、raw_used 和 record raw_len 仍保持原始长度。
- raw arena 继续一次分配、CPE 直接写入，不增加中间复制、清零或同步。
  MPI/self pack 按元数据取有效 record，不传输 padding，不改变输出数据和分组规则。
  standalone 及 memory-to-memory helper 复用同一修复；外排原有独立 decoded 槽不变。
- 将 BGZF 输入、解压器分配、inflate、CRC、QNAME 长度、QNAME NUL 错误分别标为
  -20/-21/-22/-23/-24/-25，并在主核报错中打印 reason 和实际/预期 raw 长度。
- 在流水线契约文档补充 CPE 并行输出区域必须隔离未对齐机器字写的约束。

Sunway `RabbitBAM-MPI` 构建与 `git diff --check` 通过。修复二进制 SHA-256：
`3376c05a9297fa4fbd4cfa200e4cac1e782db2a6d831ad387794209c654e8b6b`。

| 验证 | 作业 | 结果 |
|---|---:|---|
| 正式 Prepare + CPE + pack 压力回归 | 8372635 | 512 批、每批 14,090 records，逐字节及 padding/arena guard 检查通过 |
| 原故障 64 块，memory 输入 | 8372644 | 精确 raw record multiset 一致，continuity/read_order errors 均为 0 |
| 同一样本，POSIX 输入 | 8372650 | 同上；同时覆盖 backend 容量估算 |
| 全 D1，六节点六核组内排 | 8372661 | 正常结束；records=received=177,798,571，groups=88,498,907 |
| 全 D1 输出重新读取 flagstat | 8372847 | 正常结束；16 项统计与原始 D1 完全一致，177,798,571 条记录均可解析 |

压力回归直接调用当前应用源码的容量估算、`CollateMemoryExtractOperator::Prepare()`、
正式 CPE kernel、Validate/PostProcess/Pack；仅测试 wrapper 在邻核启动前加延迟，
未把测试延迟加入产品代码。还检查两批起始 offset、步长边界、metadata offset、
padding 不进入 pack。64 块有效长度为 4,169,347 bytes，两个估算入口均得到
4,171,072 bytes（增加 1,725 bytes padding）。
两个小样本输出均含 14,090 条记录、7,130 个 QNAME group。

全 D1 仍使用本页原来的资源/压缩/内存模式，输出为独立修复验证文件，不覆盖原矩阵。

| 单次对照 | 修复前成功的诊断运行 | 修复后 |
|---|---:|---:|
| 4.3 / s | 58.297711 | 58.116073 |
| 4.1~4.6 核心 / s | 59.363980 | 59.183538 |
| 输入 load / s（排除） | 101.056937 | 71.568922 |
| 最终 MPI-IO dump / s（排除） | 71.770070 | 75.572216 |
| 算法 tracked_peak_max / bytes | 20,441,062,578 | 20,445,103,522 |

单次未观察到性能回退，不作为三轮中位数结论。算法追踪峰值仅增加 4,040,944 bytes
（约 3.85 MiB）；此值不等于进程 RSS，也不包含完整预加载输入等所有内存。
最终输出 `fixed.name.bam` 为 18,755,398,492 bytes。重新读取输出的 flagstat 检查通过，
16 项统计与原始 D1 的 `bench_flagstat_r1.log` 完全一致；4.3 为 4.267056 s，
核心 body 为 4.511705 s，输入 load 194.053542 s 单独排除。
本轮精确 record multiset、QNAME 连续性及 read order 验证覆盖上述小样本；
完整 D1 验证的是全量可解析性和 flagstat 统计一致性，未做全量逐记录或 QNAME 分组比较。

本地回归资料位于 `build_sunway_sduhpc/benchmarks/collate_fix_20261001/`，
SWLS 对应目录为 `online/guoshi/wzs/wzs_data/swbam_bench/D1_apps_current/collate_fix_20261001/`。
`build_regression.py` / `regression.cpp` / `delay_slave.cpp` 复现压力测试；
`run_smoke.sh` / `verify_sample.py` 复现小样本的精确比较。
`full_D1_fixed.log` 保存完整 D1 运行日志，`verify_D1_flagstat.log` 保存全量输出重读结果。
原矩阵曾因 collate 错误暂停；下节记录使用修复版的续跑结果。

### 修复版续跑（2026-10-01 10:37）

已恢复原 D1 矩阵，从首轮 collate 开始，作业 **8373003**，节点仍为
`33449-33454`。远端二进制 SHA-256 已核对为上述修复版 `3376c05a...e8b6b`。
首轮已完成的 flagstat、stats、BAM2BAM、过滤四项保留；这些功能未受本次
collate 修改影响。剩余 26 项性能测试和 4 项输出验证自动顺序执行，失败即停。
不把前面的修复诊断作业当作正式三轮成绩。

原二进制、失败日志及续跑前结果在本地运行目录的
`archive_before_collate_fix_20261001_1036/` 保留，远端旧二进制另存为
`RabbitBAM-MPI.pre-collate-fix-20261001`。上述四项既有结果使用旧版，后续
所有作业使用修复版。首次替换遇到冻结文件只读，短暂启动的旧版
作业 8372992 已主动取消，未计入结果；修正文件权限并核对二进制后才正式续跑。

后台服务为 `swbam-apps-d1-n6-20261001-collate-fixed.service`，
继续记录原目录的 `results.tsv`、`commands.sh`、`submit.log` 和 `logs/`。
按当前 load/dump 与完整 D1 collate 的实测耗时，预计剩余约 **2～4 小时**，
即北京时间约 12:40～14:40；排队、共享存储波动或新的失败会延长此估计。
load 与验证写盘仍不计入核心成绩，但会影响实际等待时间。

### D1 正式结果（2026-10-01 13:03 完成）

后台状态 `COMPLETE`，10 个应用均完成三轮，首轮另做 4 项输出验证；
`skipped.tsv` 无跳过项。下表单位为秒，使用各命令现有核心计时口径，
即有输出命令的 `4.1~4.6`、flagstat/stats 的 `444 total body`。

| 应用 | 第 1 次 | 第 2 次 | 第 3 次 | 中位数 |
|---|---:|---:|---:|---:|
| flagstat | 3.598642 | 3.614396 | 3.614975 | **3.614396** |
| stats --basic | 3.171133 | 3.167709 | 3.166262 | **3.167709** |
| BAM2BAM | 9.267312 | 9.247038 | 9.249369 | **9.249369** |
| MAPQ >= 30 过滤 | 9.525545 | 9.510121 | 9.508204 | **9.510121** |
| collate，内排 | 59.269599 | 59.256867 | 59.243102 | **59.256867** |
| fixmate -m | 36.282483 | 36.287498 | 36.288087 | **36.287498** |
| sort，内排 | 49.113062 | 49.097309 | 49.095607 | **49.097309** |
| markdup --stream | 49.242175 | 49.233443 | 49.242145 | **49.242145** |
| BAM2SAM | 13.916938 | 13.915934 | 13.916317 | **13.916317** |
| SAM2BAM | 11.957723 | 11.959146 | 11.956962 | **11.957723** |

四步去重各阶段中位数之和为 **193.883819 s**，仅表示四次命令的核心时间相加；
三轮四步各自的完整时间仍可从 `results.tsv` 逐轮复核。
过滤保留 170,359,591 / 177,798,571 条记录。首轮 BAM2BAM、过滤、
markdup、SAM2BAM 的输出均能被 `stats --basic` 完整读取，记录数与预期一致，
`is sorted: 1`；markdup 输出统计有 50,084,605 条重复记录。
这组检查没有逐记录比较，也未与 Samtools 核对重复标记集合。

与 A4 Samtools 基准对照时，SWBAM 是排除物理 load/dump 的内存模拟核心时间，
Samtools 是 bigssd 上含磁盘读写的 wall time；因此两者之比仅代表**当前实验口径**，
不能解释为同一端到端条件下的系统加速比。`stats --basic` 与 Samtools `stats`
统计范围不同，不做直接倍数对比。

## D2 六核组容量探测与应用矩阵（2026-10-01）

D2 输入为 `D2-WGS-25G.coordinate.bam`，实测 35,255,807,945 bytes、
365,079,492 条记录。继续固定节点 `33449-33454`，每节点 1 rank、64 CPE，
`-mpecg 6 -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc`。
内存输入、rank body 内存、level 1 压缩及 MPI-IO 最终输出均与 D1 相同，
物理 load/dump 单独记录。

D2 六核组 collate 的 `6 * rank_input_bytes + 256MiB` 策略阈值约为 36.64 GB。
第一次用 `-m 40G`（42.95 GB）已选择内排，但实际 tracked workset 达
47.61 GB（六个 rank 的最大值，作业 8375501），触发内存上限检查，未产生有效成绩。
因此保留日志 `logs/bench_collate_r1_m40_failed.log`，改用 `-m 48G`
重新探测（作业 8375600）。第二次输入 load=338.457074 s，仍选中内排，
但核心阶段未打印完成摘要就以 exit code 1 结束，日志未报告具体错误分支；
目前无法判定是内存分配、MPI 交换还是其他原因，也没有有效输出或成绩。
第二次日志保存在 `logs/bench_collate_r1.log`，四步容量探测因此未通过。
sort 原计划保持 `-m 32G`，但未被提交；四步流程先单列为待诊断。
首轮容量探测顺序为 collate、fixmate -m、sort、markdup --stream，
必须正常退出、记录数匹配，且 collate/sort 为 `mode=memory`。
本次两次探测均未得到有效 collate 输出，因此四步链路没有继续提交。
后续独立应用矩阵跳过这四步，每个可独立执行的应用计划三次。

BAM2SAM 放在其他应用三轮之后，避免它的内存需求影响四步流程复测。
SAM2BAM 仅在实际生成的 SAM 不超过 72 GiB 全内存安全阈值时执行，
否则将实际字节数写入 `skipped.tsv`。阈值是当前全文件 memory backend 的
预防性限制，不表示 SAM2BAM 算法本身不支持大文件。

| 项目 | 位置或状态 |
|---|---|
| 本地结果目录 | `build_sunway_sduhpc/benchmarks/apps_D2_n6_20261001_140017/` |
| SWLS 结果目录 | `online/guoshi/wzs/wzs_data/swbam_bench/D2_apps_current/apps_D2_n6_20261001_140017/` |
| 二进制 SHA-256 | `3376c05a9297fa4fbd4cfa200e4cac1e782db2a6d831ad387794209c654e8b6b` |
| 首次容量探测 | 8375501，`-m 40G` 工作集超限，已停止 |
| 容量重探 | 8375600，`-m 48G` 核心阶段返回失败，待定位 |
| 独立应用矩阵服务 | `swbam-d2-apps-simple-20261001-140017.service` |
| 独立应用首项 | flagstat 作业 8375701，14:26 完成 |

驱动为 [D2 测试脚本](../../tools/bench_apps_d2.sh)，续跑闸门为
[探测后续跑脚本](../../tools/bench_apps_d2_continue.sh)。
`pilot48_submit.log`、`matrix48_submit.log`、`simple_submit.log`、
`results.tsv`、`skipped.tsv`、
`commands.sh` 与 `logs/` 保存命令、原始结果和错误。矩阵时间表待完成后填写。

独立矩阵在 16:07 完成 flagstat、stats --basic、BAM2BAM、MAPQ >= 30
过滤各三轮；两项输出首轮使用 `stats --basic` 重新读取，均记录数匹配、
`is sorted: 1`。过滤保留 328,059,484 / 365,079,492 条记录。
下表为各命令现有核心时间口径，单位秒：

| 应用 | 第 1 次 | 第 2 次 | 第 3 次 | 中位数 |
|---|---:|---:|---:|---:|
| flagstat | 9.220707 | 9.213445 | 9.214041 | **9.214041** |
| stats --basic | 8.545855 | 8.541871 | 8.545384 | **8.545384** |
| BAM2BAM | 25.055052 | 25.044349 | 25.056666 | **25.055052** |
| MAPQ >= 30 过滤 | 25.038824 | 25.035280 | 25.042066 | **25.038824** |

BAM2BAM 和过滤三次输出大小分别稳定为 39,508,070,956 bytes 和
34,885,964,113 bytes。首轮重读只验证可解析性、数量及排序状态，
没有做逐记录字段比对。

随后首个 BAM2SAM 作业 8377360 在 16:11 失败：输入 load=248.798189 s，
`4.3 OptimizedBamToSamMPI`=20.255727 s，但 4.5b 报
`[rank 5] ERROR: failed to allocate simulated output memory. size=26897303854`。
本路径同时保留每 rank 的完整 D2 输入和 SAM rank body，在模拟写阶段还要
分配约 26.90 GB 的第二份输出，超出当次可用内存。它没有完成 4.1~4.6，
因此不能将 4.3 当作正式 BAM2SAM 成绩；SAM2BAM 也未被提交。
当前 `state=FAILED`，后台服务已停止，SWLS 无遗留测试作业。

### D2 四步内排容量重探

六核组 collate 进一步设为 `-m 64G`（作业 8378534，固定节点
`33449-33454`）：输入 load=423.417244 s，内排策略已生效，但仍在
4.3 阶段以 exit code 1 退出，未产生输出。`-m` 只是算法工作集上限，
不会增加实际可分配内存。
原始日志位于 `build_sunway_sduhpc/benchmarks/apps_D2_pipeline_n6_m64_20261001_172518/`。

诊断重跑作业 8378979 确认六个 rank 均完成 extract、sort 和全部 exchange，
tracked workset 最大为 47,607,489,490 bytes，且通过 `-m64G` 检查；
随后全部在 merge 的 rank body `Append` 报 `errno=12 (Cannot allocate memory)`。
失败时每 rank 已写约 5.79–6.06 GB 压缩 body。collate 初始化内存 sink
时以本 rank 压缩输入大小作为初始容量，而公共 `EnsureMemoryCapacity()`
超出容量后按 2 倍 `realloc`；同时每 rank 仍保留约 35.26 GB 完整输入和
上述 collate 工作集。因此是输出缓冲扩容的实际内存峰值问题，不是
`-m` 阈值问题。诊断日志见
`build_sunway_sduhpc/benchmarks/apps_D2_collate_diag_n6b_20261001_181004/`。

随后提交 12 核组容量探测（作业 8378648，每节点 1 rank、64 CPE，
`-health 0` 自动选择健康节点，collate 72 桶，collate/sort/markdup
`-m 64G`）。每个 rank 仍加载完整输入，但分到的 collate 数据和工作集约减半。
四步命令依次为 collate、fixmate `-m`、sort、markdup `--stream`；
只有前一步成功且记录数匹配才提交下一步。此配置的结果须单独标注为
12 核组，不能与六核组成绩合并。
原始日志目录为
`build_sunway_sduhpc/benchmarks/apps_D2_pipeline_n12_m64_20261001_173546/`。
首次指定的 12 节点范围被调度器拒绝，改用 `-health 0` 后正常运行。
输入 load=319.191308 s，collate 选择内排；但 4.3 阶段持续超过 12 分钟
仍无完成摘要，17:55 主动终止作业。
诊断重跑作业 8379073 显示，rank 0–9 已从 exchange step 1 进入 step 2，
但 rank 10、11 仍停在 step 1 的第 16 个 chunk，后续不再推进，已终止。
当前 `CollateExchangeMemory()` 以各 rank 局部的 `send_done/recv_done`
决定何时退出该 step，随后立即进入下一 step；环形 `MPI_Sendrecv`
却需要相互配对的 rank 在同一轮继续参与。不同 rank 提前退出会造成
协议死锁。原始日志见
`build_sunway_sduhpc/benchmarks/apps_D2_collate_diag_n12_20261001_181835/`。
本次同样没有有效 collate 成绩或 `name.bam`；依赖它的 fixmate、sort、
markdup 均未提交。修复需分别处理内存输出峰值与 ring exchange
收发完成协议，并重新验证性能与输出正确性。

### 2026-10-01 Collate 修复与四步续测

本轮不改变排序、分桶、压缩或模拟 I/O 口径，修复两处生命周期问题：

- **交换协议**：每条有向连接只发一次终止 header，空连接也发送一次。
  收、发方向分别完成后使用 `MPI_PROC_NULL`；payload 为零的方向不再
  发送仍需对端匹配的零长度消息。内排、外排共用 header/bytes helper。
  这允许各 rank 的 chunk 数不同，不需要逐 chunk 增加全局 barrier/allreduce。
- **内存峰值**：交换完成后，所有 segment（包括 self）已独立拥有数据。
  在 merge/压缩输出前释放提取阶段的 records capacity 和 raw arena，
  D2 预计每 rank 减少约 23 GB 的无用驻留，给 rank body 扩容留出空间。
  原 `tracked_peak` 仍在释放前统计，不把峰值人为减小；释放耗时计入
  `extract.free` 和核心时间。完整输入、memory body 和最后模拟 memcpy 保持不变。

验证：`RabbitBAM-MPI` 构建、`git diff --check` 通过。定向 MPI 回归分别在
6、12 rank 上测试空数据、不对称多 chunk、稀疏连接和 metadata 容量边界，
内排/外排均通过；接收数据在释放源缓冲后逐字节检查。
回归源码为 `tools/collate_exchange_regression.cpp`，可用
`python3 tools/build_collate_exchange_regression.py` 在现有构建基础上生成。
作业 8379642、8379656 的日志在
`build_sunway_sduhpc/benchmarks/collate_exchange_fix_20261001/`。
随后作业 8379722 增加“部分 rank 没有输入 run、其他 rank 仍传多个 chunk”的
12-rank 回归，也全部通过；日志为同目录的 `ring_empty12.log`。

真实 BAM smoke（作业 8379657，6 rank、66 bins）通过：14,090 条 records、
7,130 个 QNAME groups；与输入的完整 record multiset 一致，
`continuity_errors=0`、`read_order_errors=0`。这是小数据正确性回归，
不是 D2 正式性能成绩。

D2 续测使用固定节点 `33449-33454`，每节点 1 rank/64 CPE，大共享配置不变，
collate/sort/markdup 的算法上限设为 `-m 64G`，诊断关闭。计划四步各三次，
前一步失败即停止，第一轮最终输出用 `stats --basic` 重读检查。

| 项目 | 位置或状态 |
|---|---|
| 本地目录 | `build_sunway_sduhpc/benchmarks/apps_D2_fixed_n6_20261001_191100/` |
| SWLS 目录 | `online/guoshi/wzs/wzs_data/swbam_bench/D2_apps_current/apps_D2_fixed_n6_20261001_191100/` |
| 后台服务 | `swbam-d2-fixed-n6-20261001-191100.service` |
| 二进制 SHA-256 | `bb5590df816ea637d7a03074b3b2c140047a83189846b478cbb59f0386487100` |
| 首项作业 | 8379671，collate 第 1 次，19:19 成功退出 |

19:15 首次 collate 已通过此前失败的 merge/body 扩容点，完成核心阶段：

| 指标 | 第 1 次（完整作业已成功退出） |
|---|---:|
| 输入 load，排除在核心外 | 370.709865 s |
| 4.3 算法 | 134.939913 s |
| 4.6 模拟输出 memcpy | 2.235203 s |
| 4.1～4.6 核心 | **137.733860 s** |
| input / received records | 365,079,492 / 365,079,492 |
| groups | 181,559,971 |
| BGZF blocks | 1,927,266 |
| 交换阶段 tracked peak（非进程 RSS） | 47,607,489,490 bytes |

### 22:57 状态核查：前两步成功，sort 内存不足停止

collate 首轮在 19:19 成功退出，fixmate 首轮在 19:33 成功退出。
二者记录数均为 365,079,492。下表为首轮有效结果，尚未形成三次中位数：

| 命令 | 作业 | 4.1～4.6 核心/s | 4.3/s | input load/s | 真实 dump/s | 输出 bytes |
|---|---:|---:|---:|---:|---:|---:|
| collate | 8379671 | **137.733860** | 134.939913 | 370.709865 | 218.555723 | 52,603,854,132 |
| fixmate -m | 8379833 | **85.017577** | 82.056637 | 545.414502 | 206.111531 | 53,772,662,489 |

随后 sort 作业 8380023 读取上述 fixmate 输出，load=560.107761 s，
扫描到 2,017,608 个数据块，选择 `mode=memory`，策略估计的最大工作集
为 27,155,366,938 bytes。进入 4.3 后六个 rank 均抛出未捕获的
`std::bad_alloc`，被 signal 6 终止，19:43 后台矩阵自动停止。
没有有效 sort 核心时间或排序输出，markdup、最终 stats 重读和后两轮均未提交。
22:57 核查后台服务为 failed，SWLS 无遗留测试作业。

代码核查显示 sort 内排存在完整本地 raw、完整 packed send raw 和完整
received raw 同时驻留的阶段；raw arena 还按倍增扩容，metadata 及辅助数组
另占内存。每 rank 的完整输入也已从原始 D2 的 35.26 GB 增至 53.77 GB。
因此 `3 * rank_input_bytes + 256MiB` 的模式估计不能代表此时的进程峰值，
增大 `-m` 无法增加实际可分配内存。日志尚未标出抛异常的具体分配点，
不能仅凭此日志断言失败发生在某个特定 vector 的 reserve/assign。

下一步需先处理 sort 的驻留峰值：检查 raw 预分配与扩容、发送缓冲/源数据的
释放顺序，并评估每 rank 仅加载所属压缩输入的模拟方式，再恢复六核组内排。
改变 rank 数或切换外排的结果须另列，不能作为当前六核组内排的补测成绩。
此前 2～3 小时的完成估计已因 sort 容量失败失效。

记录数相符不代替 D2 逐记录正确性验证；本轮已做的小数据精确比对与
尚未完成的大数据最终重读检查需分别报告。

### Rank-local 输入适配（2026-10-01 晚）

新增 `--io-backend rank-memory`，只预加载各 rank 所属的连续 BGZF 块；
原 `memory` 全文件复制模式保留，四步算法、CPE kernels 和 MemReader 快速路径不改。
此模式也可用于 `flagstat`、`stats --basic` 和文件串联 `dedup-workflow`。
`tools/bench_apps_d2.sh` 增加 `SWBAM_APPS_IO_BACKEND=rank-memory`，用于
`SWBAM_APPS_PIPELINE_ONLY=1` 或容量 pilot；旧结果目录/二进制不覆盖。

时间口径：4.1 与 444 body wall 扣除 `simulation_prepare_cost()`，即 root
扫描的真实 pread 时间及本 rank 分片 allocation/load/wait。扫描解析、索引构建、
MPI 分区/广播、正常内存批读、算法和模拟输出复制仍计入；最终真实 dump 仍排除。
日志逐 rank 列出 offset/bytes/blocks 和 excluded_prepare_s，后置 222 行汇总准备时间。

六 rank 下的输入驻留估计（按块数切分，实际不保证字节等分）：

| 输入 | 原每 rank 全文件 bytes | 新每 rank 平均 bytes |
|---|---:|---:|
| D2 initial | 35,255,807,945 | 5,875,967,991 |
| D2 fixmate 输出 | 53,772,662,489 | 8,962,110,415 |
| D3 initial | 89,994,058,306 | 14,999,009,718 |

这只减少模拟磁盘输入副本，不减少排序 raw/交换状态和输出驻留。D2 内排须重新
做容量探测；D3 仍按原方案考虑外排，不能据此保证全部内存峰值小于大共享上限。

构建 `RabbitBAM-MPI`、`git diff --check` 和 benchmark shell 语法检查通过。
小数据回归尝试遇到环境阻碍：上传 online 文件系统报 `Disk quota exceeded`；
固定节点 33451 提交被拒绝（`The specified node to be allocated is NOT ok`）。
改用 base 临时保存二进制、健康节点自动分配的探测作业 8383154 也未进入程序：
计算节点不能访问 base 路径（`not found`）。这些不是新适配器的运行测试结果。
截至本次记录，尚无 rank-memory 性能或正确性成绩，未启动新的 D2/D3 矩阵。

待 online 空间清理后：先小数据对照 `memory` / `rank-memory` 的完整四步输出，
再从已保留的 D2 `apps_D2_fixed_n6_20261001_191100/out/fixmate.bam` 续测 sort。
历史 collate/fixmate 与新版本结果需标明各自输入模式和二进制版本。

### 文件转换接入与服务器核查（2026-10-02）

`run_all --io-backend rank-memory` 已接入 BAM→BAM/filter、BAM→SAM、SAM→BAM。
BAM 输入复用公共 rank-local BGZF 计划和 MemReader 快速路径；输出 header 单独
生成，不能把不含完整 header 的 rank 分片当成完整 BAM 文件前缀。
SAM 输入仅支持未压缩 `.sam`，按完整行切分，每 rank 预加载自己的连续范围。
两类输入均保持原来的解析、过滤、序列化、压缩 kernel 和批次调度。

计时继续以 4.1～4.6 为主：排除真实输入预加载、BAM 扫描/SAM 边界查找的
物理读盘时间，以及最终真实写盘；正常扫描解析、索引、MPI 分区/广播、
内存批读、算法和模拟写出复制仍计入。SAM 的分片 allocation/load/wait 与
边界读盘明细逐 rank 打印，后置 222 行汇总模拟环境准备耗时。

BAM→SAM 的 rank-memory 模式改为最多预分配 64 MiB 输出，再按实际输出增长，
不再预留压缩输入大小的五倍。完整 rank body、扩容峰值和模拟输出副本仍占
内存，因此不能据此保证 D3 转换都能运行，仍须容量探测。最终输出建议用
`--io-output-backend mpiio`，避免 root gather 完整输出；真实 MPI-IO 写盘不计入
上述核心时间。`tools/bench_apps_d2.sh` 的 rank-memory 模式现可覆盖完整矩阵，
但脚本的 SAM 平均分片大小检查仅作输入准入，不代表工作集峰值保证。

示例命令（六节点/六计算核组的大共享资源参数沿用前文）：

```bash
RabbitBAM-MPI run_all -i input.bam -o output.bam \
  --io-backend rank-memory --rank-body-backend memory \
  --io-output-backend mpiio --compress-level 1
```

将输出改为 `.sam` 即 BAM→SAM；输入改为未压缩 `.sam`、输出 `.bam` 即 SAM→BAM；
BAM→BAM 可继续使用原过滤参数。旧 `memory`、POSIX、MPI-IO 输入模式保留。

本地验证：`RabbitBAM-MPI`、`swbam-sdk-optimized-count` 和
`swbam-sdk-optimized-bam2bam` 构建通过，`git diff --check` 与
`bash -n tools/bench_apps_d2.sh` 通过。新 app 二进制 SHA-256：
`737b31b919a2d700bc1b59503bd8ca37956649777e023e6299bb50c31dc460ff`。

服务器健康探测使用**已部署的旧二进制**，不能作为新 rank-memory 实现的回归：

| 项目 | 核查结果 |
|---|---|
| 作业 | 8383507，自动健康节点，单节点 6 ranks/384 CPE |
| 时间 | 2026-10-02 00:23:52 提交，00:24:12 成功结束 |
| 输入 | 之前 collate 小数据 `batch.bam`，64 个数据块 |
| 命令 | `flagstat --io-backend memory` |
| 记录数 | 14,090；mapped 14,060；duplicates 0 |
| 结论 | 当前已有程序的读取、MPI 和 CPE 执行可用，不是正式性能成绩 |

此前探测作业 8383464 使用 login 的逻辑 base 路径，被计算节点报 `not found`；
改成 `readlink -f` 得到的 online1 物理路径后，上述作业成功。
今后提交二进制/输入路径需使用计算节点可见的物理路径。

新二进制上传到
`online/guoshi/wzs/wzs_data/swbam_bench/conversion_rank_memory_20261002/`
仍报 `Disk quota exceeded (122)`，rsync 退出码 11。文件系统整体约有 1.1 P
剩余空间，故不是整个文件系统写满；具体用户/组/项目配额层级尚未查明。
当前可以运行已部署的只读程序，但新版本上传和大输出测试仍受配额阻塞。
未删除原始数据或旧 benchmark BAM，也未提交新的 D2/D3 矩阵。
本次结束时 `bjobs` 无本用户作业；新转换路径的运行正确性与性能验证待空间
问题解决后完成，不将构建成功或旧版健康探测记作新版本测试通过。

### 2026-10-02 小数据回归与 D2 续测

04:00 核查 SWLS 可连接、可提交作业，新二进制可以上传，原 D2 的
`apps_D2_fixed_n6_20261001_191100/out/name.bam` 和 `fixmate.bam` 均保留。
采用 `-health 0` 自动选择健康节点，避免此前被拒绝的固定节点；不会取消
同账户其他项目的作业。小文件/程序写入成功不代表大输出配额足够，后续仍
以实际作业是否成功写完为准，未清理旧 benchmark 或原始数据。

先在六节点六 rank 上用 14,090 条记录的 `batch.bam` 对比完整 `memory`
和 `rank-memory`。两者均在 markdup 解码时失败，作业分别为 8386463、8386500。
本地 gzip CRC 校验及 record 布局检查表明：collate 和 fixmate 输出可解析，
排序输出从第 437 条记录起出现长度字段为零的坏记录。新旧输入模式均复现，
不是新 rank-local 预加载特有的错误；sort 的 records/received 计数正确也
不能证明输出内容正确。

定位到内排 `slave_mpi_sort_bucket_pack` 的并行变长 raw bytes 打包边界：
相邻 CPE 的 record 区间可能共享一个非对齐机器字，字节写存在 read/modify/write
竞争。本轮只改变边界复制：CPE 继续生成全部 metadata 和打包区间内的主体
records；每个区间的首尾 record 由 MPE 在 join 后复制。最多每 rank 128 条
边界 records，源/目标布局、MPI counts、排序语义和其余并行打包路径不变。
边界复制计入 `bucket_pack` 和核心时间，不改外排的 range-pack kernel。

修复后构建 `RabbitBAM-MPI` 和 `git diff --check` 通过，二进制 SHA-256：
`00fa056a33e22a39a9aa11bb0569e83beed4f0539f94c65f498e36042e9231aa`。

| 小数据检查 | 作业或结果 |
|---|---|
| 完整 memory sort / markdup | 8386534 / 8386538，均成功 |
| rank-memory sort / markdup | 8386542 / 8386544，均成功 |
| 排序输出精确检查 | 14,090 条，coordinate key 非递减，与 fixmate 输入完整 record multiset 一致 |
| 两模式输出对照 | sort、markdup 的完整有序 raw record 列表分别逐条一致；markdup 有 4,270 条 duplicate records |
| rank-memory markdup 重读 | 8386559，14,090 条，`is sorted: 1` |
| rank-memory BAM→SAM→BAM | 8386568 / 8386581，均为 14,090 条 |
| SAM 往返输出重读 | 8386589，14,090 条，`is sorted: 1` |

小数据资料保存在
`build_sunway_sduhpc/benchmarks/apps_D2_rank_smoke_{memory,rank}_n6_20261002_040600/`
及 `apps_D2_rank_smoke_fixed_{memory,rank}_n6_20261002_041300/`。
后两个目录包含修复后日志/结果/输出；对应 SWLS 目录位于 `D2_apps_current/`。
以上不是 D2 正式性能结果，也未据单次小数据时间推断性能回退比例。

04:17 启动独立后台续测，先复用已保留的 D2 fixmate BAM，顺序运行 sort、
`markdup --stream` 各三次，第一轮额外重读 markdup 输出；随后补做 BAM2SAM、
SAM2BAM 各三次并重读最后的往返输出。先前已完成的 flagstat、stats、
BAM2BAM、过滤三轮及 collate/fixmate 首轮不重跑。
本轮不补 collate/fixmate 的第 2、3 次，因此目前仍不能给出四步各三次的
完整中位数表；旧完整 memory 成绩与新 rank-memory 成绩单列。

| 项目 | 本轮配置 |
|---|---|
| 本地目录 | `build_sunway_sduhpc/benchmarks/apps_D2_rank_memory_n6_20261002_041800/` |
| SWLS 目录 | `online/guoshi/wzs/wzs_data/swbam_bench/D2_apps_current/apps_D2_rank_memory_n6_20261002_041800/` |
| 后台服务 | `swbam-d2-rank-memory-n6-20261002-041800.service` |
| sort 输入 | 原 `apps_D2_fixed_n6_20261001_191100/out/fixmate.bam`，53,772,662,489 bytes |
| 预期记录数 | 365,079,492 |
| 资源 | `-N 6 -np 1 -health 0 -mpecg 6 -cgsp 64 -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc` |
| 输入/输出 | `rank-memory` / memory rank body / MPI-IO 最终验证输出，BAM level 1 |
| 算法上限 | sort、markdup 各 `-m 64G`；sort 日志必须为 `mode=memory` |
| 首项作业 | 8386664，04:17 提交并 RUN；截至 04:18 正在扫描/准备输入，容量结果待确认 |

驱动增加尾部续测、测试输入/记录数和重复次数配置，保留完整矩阵默认行为。
同时修正 222 准备时间的解析：按 `cost` 后的数字读取，不再把日志最后的
`(rank-local preparation)` 当成零。实际提交命令记录在运行目录的
`commands.sh`，每次日志独立保存，容量、写盘或结果检查失败即停。

启动配置：

```bash
SWBAM_APPS_TAIL_ONLY=1
SWBAM_APPS_SORT_INPUT=/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data/swbam_bench/D2_apps_current/apps_D2_fixed_n6_20261001_191100/out/fixmate.bam
SWBAM_APPS_NODE_RANGE=auto
SWBAM_APPS_IO_BACKEND=rank-memory
SWBAM_APPS_COLLATE_MEMORY=64G
SWBAM_APPS_SORT_MEMORY=64G
SWBAM_APPS_MARKDUP_MEMORY=64G
```

核心计时保持 4.1～4.6，预加载/扫描物理读盘和最终真实 dump 排除；
正常内存读写与处理仍计入。开启本地用户 linger 以维持后台驱动，SSH/IDE
断开不依赖当前会话继续；服务器停机仍会中断。按历史 load/dump 与任务数量，
初步估计剩余约 2～4 小时，排队、共享存储波动或容量/配额失败会改变估计。
本次启动记录不代表首项大数据已通过，后续结果需从实际日志确认。

### 2026-10-02 09:57 核查：分片加载成功，sort 仍分配失败

作业 8386664 于 04:35:30 因六个 rank 的未捕获 `std::bad_alloc`、signal 6
失败，后台服务已自动停止，`results.tsv` 只有表头。sort 没有有效核心成绩
或输出文件，后续 markdup、文件转换和验证均未提交；此前 2～4 小时的完成
估计失效。本次核查没有提交新的作业，也没有取消同账户其他项目作业。

每 rank 预加载成功，输入 bytes 为 8,961,925,271～8,962,310,494，约 8.35 GiB，
证明该次运行不再每 rank 保存完整的 53.77 GB fixmate 输入。准备时间中，
root 真实扫描读盘为 780.337549 s，分片 allocation/load/wait 最大
211.197199 s，合计排除约 991.534748 s。4.1 正常扫描/分区等成本为
0.696826 s，4.2 为 0.000046 s，随后选择 `mode=memory`、`-m 64G` 并失败。
这些分项仅为失败诊断，不计为正式 sort 性能结果。

本轮日志没有磁盘配额错误或坏节点证据，直接失败原因是内存分配失败。
分片输入已减少驻留，但 sort 内排的 raw arena、发送/接收 raw 和 metadata
仍会同时占用空间；具体失败分配点及当时实际峰值尚未记录，不能从
27,155,366,938 bytes 的模式估计推断进程峰值，更不能仅增大 `-m` 就保证成功。
下一步应定位分配点和缓冲生命周期，评估降低 sort 驻留峰值后再续测六核组内排。

### 2026-10-02 10:19：sort 交换阶段提前释放源缓冲，续测 D2

此次仅调整 app 内排 sort 的缓冲生命周期，不修改 SDK、CPE kernel、
排序键、分桶规则、压缩级别或 MPI payload 交换方式。输入仍使用
`rank-memory`，不为节省内存而改回实际磁盘读取。

旧实现的接收数组在 CPE pack 完成前分配，本地 raw/metadata 一直到整个
交换结束才释放，因此本地源、发送副本、接收副本及分桶索引同时驻留。
现在保留 CPE pack 与 MPI counts 交换的重叠，但顺序改为：

```text
CPE pack + MPI counts exchange
→ join → MPE 补齐区间首尾 records
→ 释放 local raw、local metadata 和三个逐记录索引
→ 分配 receive arrays → self copy / MPI Sendrecv
→ 释放 send buffers → k-way merge → 原 BAM write 流水线
```

没有增加完整 record 拷贝或增加 payload 通信轮数；原边界记录复制修复保留。
源缓冲释放计入 `cleanup`，接收分配仍计入 `bucket_pack`，整体均包含在 4.3。
代价是接收数组分配/初始化不再与 CPE pack 重叠，另外增加三次阶段级
MPI readiness 归约，避免某个 rank 分配失败后其他 rank 继续等待 payload。
分桶索引、send/receive 分配失败会报告具体阶段及申请数量。

新增默认关闭的 `SWBAM_SORT_MEMORY_DIAG=1`，每 rank 仅打印 setup、extracted、
allocate-send、allocate-receive 四个阶段的缓冲容量/请求 bytes。它不是进程 RSS，
不包含输入驻留、rank body 和运行时内存；不能用这些值宣称已测得完整进程峰值。
小数据日志中六个 rank 的 allocate-receive 均确认 local raw、metadata 和索引
已归零。完整 D2 的实际容量是否通过仍待首轮作业结果，不以模式估算代替实测。

构建 `RabbitBAM-MPI`、`bash -n tools/bench_apps_d2.sh`、`git diff --check` 通过。
新二进制 SHA-256：
`59a95b06a43ebd8a9eeaa4b9c4c5b6424255df682640325474d910ce2ad8fa65`。

#### 小数据正确性

输入沿用此前的 14,090-record fixmate BAM，六节点各一个活跃核组，
rank-memory 输入、memory rank body、MPI-IO 最终输出、BAM level 1。

- sort 8390800、markdup 8390807、`stats --basic` 重读 8390812 均完成。
- sort 输出 coordinate key 非递减，完整 record multiset 与输入一致。
- 新旧 sort、markdup 的完整有序 raw record 列表分别逐条一致。
- markdup 仍为 4,270 条 duplicate records；复读为 14,090 条、`is sorted: 1`。
- 日志、命令、输出保存在 `build_sunway_sduhpc/benchmarks/apps_D2_sort_peak_smoke_n6_20261002_101300/`。

#### 0.25 GB 性能对照

使用 `WES_0.25G.bam`（256,845,825 bytes、2,668,351 records），修改前二进制为
本节之前的 `00fa056a...`，修改后为 `59a95b06...`。均为六节点六活跃核组、
`-health 0` 自动分配节点、`sort -m 8G`、rank-memory 输入、memory rank body、
MPI-IO 输出、压缩 level 1，诊断关闭。每版三次，执行顺序前/后、后/前、前/后。
物理 preload 和最终 MPI-IO 输出排除在下表时间之外；三次独立分配节点，不固定节点。

| 版本 | 作业号（r1/r2/r3） | 4.3 三次 / s | 4.3 中位数 / s | 4.1～4.6 三次 / s | 4.1～4.6 中位数 / s |
|---|---|---|---|---|---|
| 修改前 | 8390822 / 8390838 / 8390844 | 0.678720 / 0.676017 / 0.675215 | 0.676017 | 0.697443 / 0.694022 / 0.693776 | 0.694022 |
| 修改后 | 8390826 / 8390832 / 8390846 | 0.689157 / 0.688496 / 0.689585 | 0.689157 | 0.707745 / 0.707012 / 0.708129 | 0.707745 |

4.3 中位数增加 **1.94%**，4.1～4.6 中位数增加 **1.98%**；本轮小数据对照
低于约 5% 的目标，但不能推断完整 D2/D3 的退化比例。最后一轮新旧 BAM
文件 `cmp` 返回 0，逐字节一致。完整命令、六次日志和驱动位于
`build_sunway_sduhpc/benchmarks/sort_peak_wes_n6_20261002_101500/`；
对应 SWLS 目录均在 `wzs_data/swbam_bench/D2_apps_current/` 下。

#### D2 后台续测

10:19 提交首项 sort 作业 **8390860**，节点 **24430–24435**，截至提交后核查为
RUN，正在扫描/准备分片输入，尚无核心时间或容量通过结果。随后由驱动串行
运行 sort、`markdup --stream` 各三次及 BAM2SAM、SAM2BAM 各三次，并复读验证。
复用原 D2 fixmate 输入，不重复已完成的 collate/fixmate 首轮和简单应用三轮。
首轮 sort 开启上述四阶段诊断，后两轮关闭；任一阶段/结果失败即停止矩阵。

| 项目 | 配置 |
|---|---|
| 本地目录 | `build_sunway_sduhpc/benchmarks/apps_D2_sort_peak_n6_20261002_101900/` |
| SWLS 目录 | `online/guoshi/wzs/wzs_data/swbam_bench/D2_apps_current/apps_D2_sort_peak_n6_20261002_101900/` |
| 后台服务 | `swbam-d2-sort-peak-n6-20261002-101900.service` |
| 输入 | 原 `apps_D2_fixed_n6_20261001_191100/out/fixmate.bam`，53,772,662,489 bytes |
| 预期记录数 | 365,079,492 |
| 资源 | `-N 6 -np 1 -health 0 -mpecg 6 -cgsp 64 -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc` |
| 后端/预算 | `rank-memory` / memory rank body / MPI-IO；sort、markdup 各 `-m 64G` |
| 计时 | 仍使用既有核心口径，物理准备和最终真实输出排除 |

实际展开的全部 `bsub` 命令写入本地 `commands.sh`，逐项日志在 `logs/`，
成功结果追加 `results.tsv`，运行/失败状态在 `state`。驱动以用户 systemd 服务
后台运行，linger 已开启，不依赖本次 SSH/IDE 会话。依据前次 D2 准备时间，
首项预计约 20～30 分钟、整批约 4～6 小时；这只是无失败时的估计，
存储波动、排队、配额或新的容量失败会使它失效。

### 2026-10-02 11:57：修复 sort 分块 MPI 消息匹配，重新续测

上一轮作业 8390860 的日志停在 10:34:51 的 `allocate-receive`，直到
11:45 用户同意停止时仍无输出或有效成绩。已先停止
`swbam-d2-sort-peak-n6-20261002-101900.service`，再仅取消该作业，
确认它不再 RUN；本地 `state` 标记为 CANCELLED，旧日志保留。
此前完成时间估计失效，该轮不得计入 sort 的性能结果。

代码检查发现 `MpiSortSendRecvBytes` 在一个方向传输完毕后，仍向真实邻居
发送或接收 count=0 消息。MPI 零长度消息仍需要匹配；各 rank 的收发长度
不同，循环次数也不同，已退出循环的邻居可能不会再发出匹配消息。
本轮改为已完成方向使用 `MPI_PROC_NULL`，并检查 `MPI_Sendrecv` 返回值。
这个内部函数同时用于内排和外排；非空 payload 的 64 MiB 分块、标签、
数据复制、交换顺序和 CPE kernel 不变，之前的提前释放源缓冲优化保留。

新增 `tools/sort_exchange_regression.cpp`，直接包含生产实现并调用内部
byte-exchange 函数，使用与 app 相同的构建/链接配置，不加入 SDK 或应用热路径。
原回归构建脚本增加可选算法参数，默认 collate 行为不变：

```bash
python3 tools/build_collate_exchange_regression.py build_sunway_sduhpc sort
```

| 检查 | 结果 |
|---|---|
| 修复前回归，8392313 | 全空 case=0 通过，稀疏收发 case=1 不完成；随后仅终止该回归作业 |
| 修复后回归，8392325 | 四个 case 全通过：全空、稀疏/空方向、不等块数且超过两个 64 MiB 块、重复使用 metadata/raw tags |
| 数据核对 | 每个接收 byte 的来源、目标与序号 pattern 均验证通过 |
| 构建/静态检查 | `RabbitBAM-MPI` 构建、回归构建脚本语法、`git diff --check` 通过 |

回归资源为 `-N 1 -np 6 -health 0 -cgsp 64 -share_size 14000 -cache_size 32`。
修复前后回归的独立二进制和日志保存在
`build_sunway_sduhpc/benchmarks/sort_byte_exchange_n6_20261002_115000/`，
对应 SWLS 目录在 `D2_apps_current/`。这确定了函数中存在可复现的消息匹配
错误；尚未取得上一轮 D2 停滞时的运行栈，所以不将其具体停滞点宣称为已直接定位。

应用小回归继续使用六节点六活跃核组、rank-memory 输入、memory rank body、
MPI-IO 输出和 BAM level 1：sort 8392333、markdup 8392335、重读 8392344 均完成。
14,090 条记录 coordinate key 非递减，与原 fixmate 完整 multiset 一致；
新旧 sort/markdup 有序 raw records 逐条一致，duplicate records 仍为 4,270，
复读 `is sorted: 1`。目录为
`build_sunway_sduhpc/benchmarks/apps_D2_sort_mpi_smoke_n6_20261002_115400/`。

同一目录额外保留 WES_0.25G 的单次 sort（8392342），诊断关闭，核心 4.3 为
**0.687269 s**，4.1～4.6 为 **0.705865 s**，2,668,351 条记录，`mode=memory`。
与上节提前释放优化后的三次 4.3 中位数 0.689157 s 相近，本次没有观察到
额外明显回退；仅为一次探测，不是新的三次中位数结论。新输出与上节最后一轮
完整 BAM `cmp` 返回 0。命令在 `wes_command.sh`，日志在 `logs/wes_sort.log`。

新增默认关闭的两个阶段标记：`receive-ready` 在所有 rank 接收数组分配完成后
打印，`exchange-done` 在 payload 交换完成后打印。这样可以区分内存分配阶段
和 MPI 交换阶段，而不改变核心计时或增加逐记录诊断。

新二进制 SHA-256：
`85637c33868736a560889aa9d276922711db743df8e026b22f6ac18af8de7cc9`。
11:57 独立后台重提 D2，首项 sort 为 **8392376**，已 RUN，初始节点
**24430–24431、21641–21642、21647、24298**，截至启动核查正在准备输入，
尚未确认完整 D2 容量/交换通过或产生核心成绩。

| 项目 | 新轮次 |
|---|---|
| 本地目录 | `build_sunway_sduhpc/benchmarks/apps_D2_sort_mpi_fix_n6_20261002_115700/` |
| SWLS 目录 | `online/guoshi/wzs/wzs_data/swbam_bench/D2_apps_current/apps_D2_sort_mpi_fix_n6_20261002_115700/` |
| 后台服务 | `swbam-d2-sort-mpi-fix-n6-20261002-115700.service` |
| 资源 | 六节点各一个活跃核组，`-mpecg 6 -cgsp 64 -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc`；`-health 0` |
| 输入/计时 | 复用原 D2 fixmate BAM，rank-memory；真实准备和最终 MPI-IO 输出仍排除 |
| 后续安排 | sort、markdup 各三次；首轮验证后继续；再测 BAM2SAM、SAM2BAM 各三次及重读 |
| 诊断/失败策略 | 首轮 sort 六个阶段标记，后两轮关闭；容量、输出或结果验证失败即停 |

原简单应用三轮及 collate/fixmate 首轮不重测。新轮次命令、逐项日志、
成功结果和状态仍分别记录在 `commands.sh`、`logs/`、`results.tsv`、`state`。
首项按历史准备/输出速度粗估 20～30 分钟，整批无失败时约 4～6 小时；
应先以本轮首项的实际完成确认继续运行，不能把旧轮次长期 RUN 当作正常进度。

### 2026-10-02 13:02 核查：sort 首轮完成，markdup 长时间无进度日志

sort 作业 **8392376** 于 **12:16:40** 完成，提交至完成共 1,146 s。
六个 rank 均打印 `receive-ready` 和 `exchange-done`，随后成功压缩并输出，
该次完整 D2 内排实际越过了此前的分配/交换停滞阶段，未见 `bad_alloc` 或 MPI 错误。

| sort 首轮指标 | 实测 |
|---|---|
| 资源/后端 | 六节点六活跃核组，rank-memory 输入、memory rank body、MPI-IO 最终输出，level 1、`-m 64G`、`mode=memory` |
| 输入/收到记录数 | 均为 365,079,492 |
| 4.3 核心算子时间 | 110.944666 s |
| 4.1～4.6 核心总时间 | **113.292516 s** |
| 输入扫描/分片加载等真实准备 | 845.352271 s，排除在核心时间之外 |
| 最终真实 MPI-IO 输出 | 161.424594 s，排除在核心时间之外 |
| 输出文件大小 | 40,691,022,479 bytes（约 40.69 GB） |

结果已追加本轮 `results.tsv`，输出保留于 `out/fixmate.sorted.bam`。
这是一次有效完成的性能记录，不是三次中位数；完整输出的独立
`stats --basic` 重读和后续去重输出验证仍待完成，不将算子汇总等同于全部正确性验证。

驱动于 **12:16:40** 继续提交 markdup 作业 **8392543**。截至 **13:02**，
驱动和该作业仍为 RUN，但日志最后更新时间为 **12:24:51**，停在
`boundary halo=fallback blocks=64`。六个 rank 的输入准备均成功，真实准备
最大 483.826498 s，4.1 为 0.542801 s，4.2 为 0.000022 s；没有 4.3 完成信息、
markdup 输出文件或正式性能成绩，验证、其余轮次和文件转换均尚未开始。

约 38 分钟没有新进度，不能据 RUN 状态声称正常推进。当前日志未报告内存
分配、配额或解码错误，但也不足以区分计算缓慢与阻塞。静态检查确认
markdup 的 `MdSendrecvBytes` 已使用 `MPI_PROC_NULL`，不能直接将此轮停滞
归因于刚修复的 sort 零长度消息问题；需要进一步阶段诊断或运行栈定位。
原整批时间估计暂不可靠。本次核查只记录状态，没有取消、新提交或修改测试程序。

### 2026-10-02 13:21 诊断：无映射文件尾触发不必要的全局 fallback

markdup 作业 **8392543** 仍为 RUN，六个任务均未退出。`cnload -j 8392543 -l`
显示六个节点均为 busy，两个网络接口和文件系统均为 OK；未发现节点故障报告。
输入加载已成功，日志约 56 分钟没有新的阶段信息，不能作为有效性能成绩。

为核对 halo 判定，只读取 `out/fixmate.sorted.bam` 中六个 rank 各自的首尾
32 个 BGZF 块，使用 Python 3 的 `gzip` 和 `struct` 验证块内 record 边界并提取
`tid/pos`。没有扫描整份文件、提交新作业或修改 BAM；抽查耗时 **1.222 s**。
坐标使用 BAM 内部的零基 `tid/pos`：

| rank | 首记录 | 首 32 块末记录 | 尾 32 块首记录 | 尾 32 块未映射/总记录 |
| ---: | --- | --- | --- | ---: |
| 0 | 0:9995 | 0:47495 | 2:755182 | 0/5958 |
| 1 | 2:800238 | 2:844607 | 4:96673760 | 0/5972 |
| 2 | 4:96718700 | 4:96763385 | 7:85406080 | 0/5961 |
| 3 | 7:85454109 | 7:85505131 | 11:57719502 | 0/5942 |
| 4 | 11:57782367 | 11:57838319 | 16:64553676 | 0/5957 |
| 5 | 16:64598328 | 16:64645918 | -1:-1 | 6511/6511 |

全部实际存在的左右相邻边界都满足 300 bp 的覆盖要求。只有 rank 5 的尾部
坐标检查失败：`MdBuildOnePassBoundaryPlan` 无条件要求
`probe.last_first_tid >= 0`，即使该 rank 根本没有右邻居。随后 `MPI_Allreduce`
将这一局部失败变为全局 `coverage_ok=0`，因此所有 rank 都使用 fallback。
这是已由输入抽查和代码共同确认的、不必要的优化退出条件，不是输入损坏。

静态检查还发现 fallback 的放大风险：边界候选由 `MdWindowBoundaryKey`
保留到最终 finalize，按顺序写出会在首个 unresolved block 处停下，后续
decoded batches 可能持续积压；`MdOnePassFindRecord` 又从 pending 队列头
线性搜索 ordinal，积压时逐记录处理成本会随队列增长。D1 首轮为
`halo=enabled`，其核心时间为 49.242175 s，不能用该结果估计 D2 fallback。

以上确认了 fallback 的触发原因及代码层面的性能风险，**尚未取得运行栈或
批次进度，不能断言此次实际停在某个 MPI 调用，也不能确认风险就是全部停滞原因**。
尝试 `bdebug 8392543` 及针对节点 24430 的取样，管理器接受了请求，但未取得
可解析的栈输出。下一步应按是否存在相邻 rank 来检查相应 halo，并增加默认
关闭的低频阶段/批次/pending 诊断，先验证 fallback 行为，再恢复矩阵。
本次没有修改算法、取消作业或提交新作业。

### 2026-10-02 14:11 修复末尾 halo 判定并恢复 D2

本次停止旧驱动并仅取消停滞的 markdup 作业 **8392543**（13:51:46 退出），
保留 sort 首轮 **8392376** 的有效成绩和完整输出。旧二进制、markdup 日志及
本地驱动/命令记录归档在本轮目录的
`archive_markdup_halo_20261002_135700/`，SWLS 也保留旧二进制和失败日志。

`MdBuildOnePassBoundaryPlan` 现在仅在存在左邻居时检查首部 halo，仅在存在
右邻居时检查尾部 halo。最后一个 rank 的无映射文件尾不再导致不必要的
全局 fallback；真实相邻边界不足、空 rank 等情况仍保留原有回退。
没有修改 duplicate 判定、CPE kernel、MPI payload 或逐记录处理路径。

增加默认关闭的 `SWBAM_MARKDUP_DIAG=1`：首批、每 512 批及补齐回调打印
候选交换、window 处理、决策交换和压缩阶段，以及 pending/captured/retired
块数。诊断用于确认运行进度，不把 `RUN` 状态当成成功；恢复矩阵仅首轮
markdup 开启，后两轮关闭。

| 修复验证 | 结果 |
|---|---|
| 构建/静态检查 | `RabbitBAM-MPI` 和边界回归程序构建成功；`git diff --check`、驱动 `bash -n` 通过 |
| 生产 halo 探测回归，8393887 | 正常 mapped halo、最后 rank 的 unmapped tail 均 enabled；真实左/右覆盖不足均 fallback，四项通过 |
| 六节点小数据 markdup，8393889 | 14,090 条记录，halo enabled，六个 rank 均完成；首轮诊断开启 |
| 独立完整复读，8393895 | 14,090 条记录，4,270 条 duplicate，`is sorted: 1` |
| 逐条对比 | 排除 header 后，新旧有序 raw BAM records 完全相同，包括 duplicate flags |

回归资源与正式实验相同：六节点、每节点一个 rank，
`-mpecg 6 -cgsp 64 -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc`。
初次使用单节点六 ranks 的测试工具作业 **8393854** 在进入 fixture 写入前
发生段错误，没有有效测试结论；改用上述正式资源配置并补充测试阶段日志后
通过。未取得首次测试工具失败的栈，不将其归因为服务器或生产 halo 修复。
回归二进制、完整命令、输出和日志在
`build_sunway_sduhpc/benchmarks/markdup_halo_fix_n6_20261002_135700/`，
SWLS 对应目录位于 `D2_apps_current/`。

小数据单次 4.3 为 0.053786 s、4.1～4.6 为 0.061539 s，诊断开启且数据很小，
仅用于正确性验证，不用于推断 D2 性能或宣称性能没有回退。

新二进制 SHA-256：
`842f013a608a92ae3b7a98c8a39fe77cd4fb1417f7003cf5a1100ca7f018a67b`。
后台服务为 `swbam-d2-markdup-halo-fix-n6-20261002-141100.service`，复用
`apps_D2_sort_mpi_fix_n6_20261002_115700/` 的驱动和 `results.tsv`，新增启动
日志为 `halo_fix_submit.log`。配置为 resume/tail-only、rank-memory、sort 和
markdup 均 `-m 64G`；简单应用与 collate/fixmate 不重复测试。

14:11:29 首项恢复为 D2 markdup **8393989**，已启动并准备各 rank 输入。
后续依次独立复读首轮去重输出、完成其余两轮 sort/markdup，再测 BAM2SAM、
SAM2BAM 各三次和往返复读。物理输入准备与最终 MPI-IO 输出仍排除在正式
4.1～4.6 时间之外，正常内存批读/计算/模拟写复制仍计入。此刻尚无新 D2
markdup 核心成绩；整批耗时需结合首项实际运行和共享存储速度更新。

截至 **14:26**，新作业 8393989 和后台驱动仍运行，日志尚未打印分片加载
完成或进入 halo/core 的信息。六个节点的网络和文件系统健康检查均为 OK，
远端二进制 SHA-256 与新构建一致；没有新的容量、解码或配额错误报告。
输入准备已比上一轮慢，但当前证据不能直接确定是共享存储延迟还是其他等待，
也不能宣称完整 D2 已恢复正常。仅小数据和定向边界回归已经通过。

旧驱动已停止，新驱动在独立后台服务中继续，当前会话结束不影响它；任何
提交、容量或汇总验证失败都会停止后续矩阵。若后续核心处理和存储正常，
按剩余 sort/markdup、六项文件转换和两项复读粗估还需 **4～8 小时**，
输入扫描/加载及大型 SAM 输出波动可能延长该范围。续跑的
`logs/bench_markdup_r1.log` 本地日志仅包含新作业；SWLS 同名 `-o` 日志由调度器
追加，包含旧作业尾部，分析时以新作业或本地独立提交日志为准。

### 2026-10-02 14:43 核查：D2 markdup 首轮和完整复读通过

修复后的 D2 markdup **8393989** 于 **14:34:00** 完成，提交至完成共
1,351 s。日志明确打印 `boundary halo=enabled blocks=64`，六个 rank 均持续
越过批次阶段，直到 5,120 批及最终汇总；本次完整 D2 没有再次出现此前的
长期 fallback 停滞，也没有容量、解码或 MPI 错误报告。

| markdup 首轮指标 | 实测 |
|---|---|
| 资源/后端 | 六节点六活跃核组；rank-memory 输入、memory rank body、MPI-IO 最终输出；level 1、`--stream -l 300 -m 64G` |
| 输入/输出记录数 | 均为 365,079,492，removed=0 |
| duplicate 标记记录数 | 8,894,488 |
| 4.3 核心算子时间 | 164.979447 s |
| 4.1～4.6 核心总时间 | **167.253638 s** |
| 真实输入准备 | 1,018.876642 s，排除 |
| 最终真实 MPI-IO 输出 | 142.177528 s，排除 |
| 输出大小 | 40,718,592,510 bytes |
| 诊断状态 | 首轮低频批次诊断开启；后两轮关闭，最终性能仍需三次结果 |

随后独立 `stats --basic` 作业 **8394244** 于 **14:42:23** 完成：完整解析
365,079,492 条记录，`is sorted: 1`，`reads duplicated: 8894488`，与 markdup
汇总一致。该验证确认可解析性、记录数、coordinate 排序及 duplicate bit
计数，不等同于已与 Samtools 完整逐记录核对 duplicate 决策。

有效成绩与验证已经追加 `results.tsv`。本轮 sort 首次 4.1～4.6 为
113.292516 s，markdup 首次为 167.253638 s；当前都不能替代三次中位数。
后台驱动已自动提交第二轮 sort **8394298**（14:42:23 启动），截至本次
核查仍为 RUN，正在准备输入，尚无第二轮核心成绩。剩余为第二/三轮
sort、markdup，BAM2SAM/SAM2BAM 各三次及最终往返复读，粗估仍需
**4～7 小时**，主要受真实扫描/加载和大型 SAM 写出速度影响。后台服务
继续运行，本次核查没有额外提交、取消或修改算法。

### 2026-10-02 21:45 核查：三轮 sort/markdup 完成，补跑文件转换尾项

sort、markdup 三轮全部成功，各轮均为 365,079,492 条记录；markdup 三轮
均为 `halo=enabled`、marked=8,894,488、removed=0。核心总时间如下，仍采用
4.1～4.6，不计真实输入准备与最终 MPI-IO 输出：

| 功能 | 第一次/s | 第二次/s | 第三次/s | 三次中位数/s |
|---|---:|---:|---:|---:|
| coordinate sort | 113.292516 | 113.269820 | 113.293716 | **113.292516** |
| markdup --stream | 167.253638 | 167.166642 | 167.182997 | **167.182997** |
| BAM→SAM | 35.243717 | 35.231319 | 待补测 | 尚未形成 |
| SAM→BAM | 29.548895 | 29.545193 | 待补测 | 尚未形成 |

sort 和 markdup 的极差/中位数分别约 0.0211% 和 0.0520%，核心成绩稳定。
markdup 首轮低频诊断开启，后两轮关闭。先前首轮去重输出的独立完整复读
通过；SAM 往返的独立复读尚未执行，不能把转换算子的汇总当作最终验证。
两轮 BAM→SAM 输出均为 162,143,967,115 bytes，两轮 SAM→BAM 输出均为
39,508,302,575 bytes；算子汇总的记录数均正确。

后台于 **17:10:43** 停止，触发点是第三轮 BAM→SAM 的**提交失败**：

```text
job submit failed, ret = -19, reason: No enough compute nodes
```

该次没有分配 job ID，没有运行转换程序，不是新的算法、容量或解析错误。
成功结果保留在 `results.tsv`。失败日志及此前启动日志已另存于
`archive_submit_no_nodes_20261002_214500/`。

21:42 核查时 `q_share` 为 Open:Active，调度器报告有空闲节点；这仅表示
当前资源状态，不能反推 17:10 的拒绝原因。使用相同参数和二进制恢复原
resume 驱动，所有成功项目均跳过。第三轮 BAM→SAM **8399353** 于
**21:44:27** 提交成功，已进入 rank-memory 输入准备；随后只补第三轮
SAM→BAM 和 `stats --basic` 往返复读，没有修改代码或扩大实验范围。

原后台服务 `swbam-d2-markdup-halo-fix-n6-20261002-141100.service` 已恢复运行。
按前两轮转换（BAM→SAM 提交至完成约 20～24 分钟、SAM→BAM 约 12 分钟）
与复读准备耗时估计，剩余约 **45～90 分钟**；最终输出与完整验证成功前
仍不将整个 D2 矩阵标为完成。

### 2026-10-03 01:38 核查：三轮性能结果齐全，补交最终复读

第三轮 BAM→SAM **8399353** 已于 10 月 2 日 22:08:41 完成。第三轮
SAM→BAM **8399690** 的远端日志也有完整算法汇总、MPI-IO 写出和
`TOT TIME 666.774413`，无 ERROR 记录，输出文件修改时间为 22:19:48。
各轮转换均为 365,079,492 条记录，BAM→SAM 输出均为 162,143,967,115 bytes，
SAM→BAM 输出均为 39,508,302,575 bytes。当前三次核心总时间汇总为：

| 功能 | 第一次/s | 第二次/s | 第三次/s | 三次中位数/s |
|---|---:|---:|---:|---:|
| coordinate sort | 113.292516 | 113.269820 | 113.293716 | **113.292516** |
| markdup --stream | 167.253638 | 167.166642 | 167.182997 | **167.182997** |
| BAM→SAM | 35.243717 | 35.231319 | 35.268094 | **35.243717** |
| SAM→BAM | 29.548895 | 29.545193 | 29.537109 | **29.545193** |

计时仍为六节点六活跃核组、rank-memory 的 4.1～4.6，物理输入准备与最终
输出排除。第三轮 SAM→BAM 的 4.3=27.829587 s，真实输入准备=407.112743 s，
真实最终输出=222.982986 s。

本地后台服务在 22:20:06 以 status=1 退出，本地 SAM→BAM 日志截断在输入
切分阶段，`state` 文件为空，第三轮结果及最终复读没有写入。远端完整日志
证明转换已执行到最终输出和进程计时，不需要重跑这一轮。原始本地截断日志
保留，完整远端日志另存为 `logs/bench_sam2bam_r3.remote.log`，核对单一汇总、
记录数、ranks、完成标记和输出大小后补回 `results.tsv`；缺失的 SSH 提交至
返回耗时记为 `NA`，不使用程序自身 TOT TIME 冒充该指标。

当前本地 Btrfs 文件系统显示使用率 100%，仍有约 15 GiB 空闲。日志记录
不完整是已确认的事实，但没有完整错误消息，**不能仅凭当前磁盘状态确定
22:20 的驱动退出原因**；本次没有删除数据、修改算法或构造未测结果。

01:37:07 恢复原 resume 服务，三轮性能项目全部跳过，仅提交最终
`stats --basic -i out/sam2bam.bam --io-backend rank-memory` 复读 **8401608**。
作业已启动，正在输入准备，验证要求仍为完整 365,079,492 条记录且
`is sorted: 1`。首轮 markdup 独立复读已通过；最终 SAM 往返复读尚待结果，
预计约 **10～20 分钟**，主要取决于真实扫描/加载，不计入新的性能中位数。

### 2026-10-03 02:58 核查：D2 补测矩阵完成

最终 SAM 往返复读 **8401608** 于 **01:48:31** 完成，作业日志包含
`Job 8401608 has been finished.`，驱动随即打印 COMPLETE，`state` 也为
`COMPLETE`；后台服务正常退出，不再运行。

`stats --basic` 完整解析 **365,079,492** 条记录，`is sorted: 1`，无错误报告。
其中主记录 363,119,942 条、supplementary 1,959,550 条；不能把排除补充比对的
`raw total sequences` 与所有记录总数混淆。该输出来自初始 coordinate BAM 的
SAM 往返，不是 markdup 输出，`reads duplicated: 0`；先前独立 markdup 复读的
8,894,488 条 duplicate 标记仍单独记录。

复读实际提交至完成 684 s，其中物理扫描/分片加载 653.002692 s；排除准备后
body 时间 9.440242 s，4.3 为 8.902860 s。这是正确性验证，不加入三次性能统计。
它确认可解析性、完整记录数和 coordinate 排序，尚未进行整份 BAM 逐记录
字段等价检查，不能将该复读称为完整内容逐字节验证。

本轮 `results.tsv` 经结构化核对，共 **12 项性能结果 + 2 项独立复读**，
sort、markdup、BAM→SAM、SAM→BAM 均各三次，没有重复行或跳过项。
三次核心中位数保持为 **113.292516、167.182997、35.243717、29.545193 s**。
第三轮 SAM→BAM 的提交至返回耗时仍为 NA，恢复的核心成绩来源保持注明。

因此这批 D2 sort/markdup 与文件转换补测已结束；不扩大为 collate/fixmate
已经各测三轮或 D3 已完成的结论。本次仅更新结果记录，没有新增或取消作业、
修改算法、重跑已完成测试。

### 2026-10-03 03:25 启动 D2 collate/fixmate 的 rank-memory 三次补测

本轮仅补齐 collate、fixmate，各三次，不重跑已经完成的简单应用、sort、
markdup 或文件转换，也不启动 D3。旧整文件 `memory` 的单次成绩仍保留在
原目录，不与本轮 `rank-memory` 数据混算中位数。

| 项目 | 配置 |
|---|---|
| 本轮目录 | `apps_D2_collate_fixmate_rank_n6_20261003_032500` |
| 本地根目录 | `build_sunway_sduhpc/benchmarks/` 下同名目录 |
| SWLS 目录 | `online/guoshi/wzs/wzs_data/swbam_bench/D2_apps_current/` 下同名目录 |
| 输入 | `D2-WGS-25G/initial/D2-WGS-25G.coordinate.bam`，35,255,807,945 bytes |
| 预期记录数 | 每次输入、输出均为 365,079,492 |
| 计算资源 | 六节点、每节点一个 MPI rank 和一个活跃核组，共 6 MPE + 384 CPE |
| 大共享资源 | `-mpecg 6 -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc`；其余核组仅提供内存 |
| 输入/输出 | `--io-backend rank-memory --rank-body-backend memory --io-output-backend mpiio` |
| collate | `-n 66 -m 64G -T tmp/collate`，必须选择 memory 模式 |
| fixmate | `-m`，使用同一轮的 collate 输出；这里的 `-m` 是生成 mate score 的开关 |
| 压缩 | BAM level 1 |
| 诊断 | 关闭 `SWBAM_COLLATE_DIAG` |
| 正式计时 | `Complete the total (4.1~4.6) cost`，同时记录 4.3 |
| 排除项 | 真实输入扫描/分片加载、最终 MPI-IO 写出；正常内存批读、算法和模拟输出复制仍计入 |

冻结使用当前已验证的 `RabbitBAM-MPI`，没有修改或重构算法。二进制
SHA-256 为 `842f013a608a92ae3b7a98c8a39fe77cd4fb1417f7003cf5a1100ca7f018a67b`，
与上一批 D2 sort/markdup/转换实验一致。只给 `tools/bench_apps_d2.sh` 增加
`SWBAM_APPS_COLLATE_FIXMATE_ONLY=1` 模式，并修正这两份 QNAME 分组输出的
复读检查，不再错误地要求 coordinate 排序；其他输出的排序检查保持不变。
脚本语法、`git diff --check`、隔离的模拟提交测试均通过，模拟测试还确认
只执行六项性能测试和两项复读、错误记录数会停止、原 coordinate 检查仍有效，
resume 不会重复提交完成项。模拟结果仅验证脚本，未当作性能数据保存。

完整提交命令逐项写入本轮 `commands.sh`，冻结的脚本、配置与 SHA 在
`driver.sh`、`launch.sh` 和 `metadata.txt`。每次提交采用以下资源前缀，
`EXE` 为本轮目录中二进制的 SWLS 绝对路径，`LABEL` 为独立日志名：

```bash
bsub -q q_share -I -b -K -J "$LABEL" \
  -N 6 -np 1 -health 0 -mpecg 6 -cgsp 64 \
  -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc \
  -timelimit 01:30:00 -o "logs/$LABEL.log" \
  -x SWBAM_COLLATE_DIAG=0 "$EXE" [command and arguments]

"$EXE" collate \
  -i ../../../D2-WGS-25G/initial/D2-WGS-25G.coordinate.bam \
  -o out/name.bam -n 66 -m 64G -T tmp/collate \
  --io-backend rank-memory --rank-body-backend memory \
  --io-output-backend mpiio --compress-level 1

"$EXE" fixmate -m -i out/name.bam -o out/fixmate.bam \
  --io-backend rank-memory --rank-body-backend memory \
  --io-output-backend mpiio --compress-level 1

"$EXE" stats --basic -i out/name.bam --io-backend rank-memory
"$EXE" stats --basic -i out/fixmate.bam --io-backend rank-memory
```

顺序为 collate 1 → fixmate 1 → 两份输出独立复读 → collate 2 → fixmate 2 →
collate 3 → fixmate 3。复读确认完整可解析性与总记录数，不要求 coordinate
排序，也不据此声称全部 QNAME 分组和 mate 字段已逐条验证。三轮复用两个
输出文件，只保留最后一轮文件，避免额外复制；不把 BAM 下载到本地。

后台服务为 `swbam-d2-collate-fixmate-rank-n6-20261003-032500.service`，
独立于当前会话运行。首项 collate **8402980** 于 **03:25:03** 启动，节点为
`43664-43665,43654,43666,43671-43672`，远端二进制 SHA 核对一致，程序已打印
`rank-memory` 和 `6 MPE + 384 CPEs`，当前仍在输入准备，没有本轮有效成绩。
提交或汇总检查失败时脚本停止后续矩阵，不取消其他实验的作业。

按之前完整 D2 collate/fixmate 的实际运行耗时估算，本轮含复读约需
**2～3 小时**，即暂估 10 月 3 日 **05:25～06:25** 完成；排队和共享存储
波动可能延长该范围。这是估计，不是完成承诺；`results.tsv` 齐全且
`state=COMPLETE` 前不将这两项的三次补测标为完成。

### 2026-10-03 12:08 核查：D2 collate/fixmate 三轮补测完成

本轮于 **05:24:41** 正常结束，`state=COMPLETE`，后台服务为 inactive/dead，
`ExecMainStatus=0`。核对八份独立日志，均有正确的单次算法汇总、`TOT TIME`
及调度器完成标记，无 ERROR 或非零退出报告；SWLS 当前账户没有运行中的
作业。本次仅核查与记录结果，没有新提交或取消作业，也没有修改算法。

| 功能 | 第一次核心/s | 第二次核心/s | 第三次核心/s | 核心中位数/s | 4.3 中位数/s |
|---|---:|---:|---:|---:|---:|
| collate，内排 | 136.460361 | 136.563449 | 136.543614 | **136.543614** | 133.677977 |
| fixmate -m | 84.926698 | 85.016266 | 85.030173 | **85.016266** | 81.983670 |

资源和计时口径沿用本轮配置：六节点六活跃核组、rank-memory 输入、memory
rank body、MPI-IO 最终输出、BAM level 1，主成绩为 4.1～4.6，不含物理输入
准备和最终写盘。collate 三轮均为 `mode=memory`，`-m 64G`。极差占中位数
分别为 **0.0755%**、**0.1217%**，没有再次出现容量失败或交换停滞。

| 功能/轮次 | 作业号 | 输入准备/s（排除） | 最终 MPI-IO/s（排除） | 提交至完成/s |
|---|---:|---:|---:|---:|
| collate 1 | 8402980 | 472.288577 | 213.470599 | 850 |
| fixmate 1 | 8403084 | 635.533891 | 222.617644 | 971 |
| collate 2 | 8403404 | 471.341523 | 300.498908 | 933 |
| fixmate 2 | 8403524 | 668.508628 | 305.674719 | 1081 |
| collate 3 | 8403647 | 434.297829 | 297.428471 | 894 |
| fixmate 3 | 8403766 | 671.649396 | 191.319788 | 970 |

从首次提交 03:25:02 至 COMPLETE 共 **1 小时 59 分 39 秒**，包含两次
正确性复读及作业间的驱动开销。物理 I/O 有波动，核心成绩保持稳定；这些
准备/写盘时间未混入正式核心中位数。

所有性能轮次的记录数均为 **365,079,492**，collate 的 received 数也与
输入总数一致。三轮 collate 输出均为 **52,603,854,132 bytes**，fixmate
均为 **53,772,662,489 bytes**，远端最终文件大小与日志一致。

首轮输出独立复读 **8403215 / 8403300** 均完整解析 365,079,492 条记录，
`is sorted: 0`，符合 QNAME 分组输出不保证 coordinate 顺序的预期。
这确认完整可解析性与记录数；不将其扩大为逐条 QNAME 分组、mate 字段或
全部内容等价验证。第三轮输出保留在本轮 SWLS 的 `out/name.bam`、
`out/fixmate.bam`；本地未下载大型输出。

按既定混合计时口径，A4 Samtools collate/fixmate 的端到端中位数
429.26 / 178.85 s，除以本轮 SWBAM 内存模拟核心时间，数值比为
**3.14 / 2.10**。这不是双方都包含真实磁盘 I/O 的端到端加速比。

#### D2 独立应用矩阵汇总

三批 `results.tsv` 合并核对，共 **10 个功能，各三次，30 项性能记录**，
以及 **6 项输出复读**，没有重复的功能/轮次键。旧 collate/fixmate 单次
整文件 memory 成绩不加入本表；其他功能保持其已完成的原始成绩。

| 功能 | SWBAM 主成绩中位数/s | 本轮输入模式 |
|---|---:|---|
| flagstat | 9.214041 | memory |
| stats --basic | 8.545384 | memory |
| BAM→BAM | 25.055052 | memory |
| MAPQ ≥ 30 过滤 | 25.038824 | memory |
| collate，内排 | 136.543614 | rank-memory |
| fixmate -m | 85.016266 | rank-memory |
| coordinate sort，内排 | 113.292516 | rank-memory |
| markdup --stream | 167.182997 | rank-memory |
| BAM→SAM | 35.243717 | rank-memory |
| SAM→BAM | 29.545193 | rank-memory |

原始来源分别为 `apps_D2_n6_20261001_140017`（前四项），
`apps_D2_collate_fixmate_rank_n6_20261003_032500`（collate/fixmate），以及
`apps_D2_sort_mpi_fix_n6_20261002_115700`（后四项）。所有项目使用六节点
六活跃核组、内存模拟 I/O；只读统计主时间为 body，写出项目为 4.1～4.6。
原有首轮诊断开关和第三轮 SAM→BAM 提交 wall 缺失的说明仍有效，不将本表
称为全部使用同一输入驻留方式、同一诊断状态的重测。

至此 D2 的独立应用三次矩阵已完成；四步去重是分阶段独立计时，不据此
声称已有新的整体端到端流程三次成绩。D3 应用外排实验尚未启动。
