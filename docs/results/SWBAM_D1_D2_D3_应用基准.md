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

### 2026-10-03 14:17 启动 D3 可行项目与外排回归

使用已保存的 **0309267** app/SDK 源码及与 D2 相同的冻结二进制，SHA-256
`842f013a608a92ae3b7a98c8a39fe77cd4fb1417f7003cf5a1100ca7f018a67b`。
本次仅扩展测试驱动，没有修改 SDK、app 算法、CPE kernel 或计时实现。

D3 输入为 `D3-WGS-85G/initial/D3-WGS-85G.coordinate.bam`，远端实测
**89,994,058,306 bytes**，预期 **1,210,633,370 records**。六节点、每节点
一个 MPI rank 和一个活跃核组（6 MPE + 384 CPE），其他核组仅提供内存。
统一使用 rank-memory 分片预加载、memory rank body、MPI-IO 最终输出、
BAM level 1，诊断关闭；保留真实输出以进行阶段串联和独立复读。

```bash
-N 6 -np 1 -health 0 -mpecg 6 -cgsp 64 \
  -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
```

本轮根目录为
`build_sunway_sduhpc/benchmarks/apps_D3_n6_20261003_141600/`，后台服务为
`swbam-d3-apps-n6-20261003-141600.service`。三个阶段使用独立子目录，本地
和 SWLS 名称均为根目录名加 `_external_smoke`、`_simple`、`_external`，
SWLS 统一放在 `wzs_data/swbam_bench/D3_apps_current/` 下。日志、实际完整
提交命令、13 列原始结果分别保存在各阶段的 `logs/`、`commands.sh`、
`results.tsv`；冻结驱动和全部配置位于根目录。

| 阶段 | 内容 | 正式轮次 | 每项时限 |
|---|---|---:|---|
| external_smoke | WES_0.25G 的 collate → fixmate → sort → markdup，四份输出完整复读 | 各一次，仅回归，不加入 D3 成绩 | 15 分钟 |
| simple | D3 flagstat、stats --basic、BAM→BAM、MAPQ≥30；首轮回环/过滤输出复读 | 各三次，12 项性能 + 2 项复读 | 2 小时 |
| external | D3 collate → fixmate → sort → markdup；首轮四份输出复读 | 各三次，12 项性能 + 4 项复读 | 6 小时 |

小回归使用 **256,845,825 bytes、2,668,351 records**，collate/sort 均
`-m 128M`，确保实际走外排并产生临时读写。正式 D3 collate/sort 均
`-m 8G`，六 ranks 的合计算法预算为 **48 GiB**；collate 使用 `-n 66`，
fixmate 的 `-m` 为 mate-score 开关，markdup 使用 `--stream -l 300 -m 8G`。
不测试 D3 内排，不复测已经完成的 SDK 矩阵。

`tools/bench_apps_d2.sh` 保留原 D2 默认配置，只新增 D3 数据集、external
模式检查、转换跳过、时间上限和外排原始指标支持；`tools/bench_apps_d3.sh`
负责按阶段编排。语法及 `git diff --check` 通过。隔离模拟提交测试确认：
原 D2 矩阵仍为 30 性能项 + 4 复读；D3 回归/基础/外排分别为 4+4、12+2、
12+4，未执行 SAM 转换。错误模式或没有真实 temp write 会拒绝结果，回归
失败时仅阻止 D3 外排，基础矩阵仍运行。模拟测试不作为性能数据。

真实执行顺序为小外排回归 → D3 基础矩阵 → D3 外排矩阵。仅当小回归全部
通过时启用 D3 外排；D3 外排第一轮四个功能和输出检查全部成功后才开始
第二、三轮。任一正式阶段失败停止后续提交，不自动改 ranks、放大内存或
退回内排。每次提交使用独立日志，输出文件在该阶段内按轮复用，保留最后
一轮文件；不下载大型 BAM 到本地，也不删除之前 D1/D2 的输出。

只读主指标为 body，写出主指标为 4.1～4.6，同时记录 4.3、物理输入准备与
最终 MPI-IO 写盘。输入准备和真实最终输出排除；正常内存读取、算法和模拟
写复制仍计入。外排沿用现有模拟高速临时 I/O 的报告模型，**不是全部阶段
都使用单一实测墙钟边界**：collate 排除实际临时读写，sort 使用源码中
已有的分项模拟时间模型。完整实际 wall、临时读写细分仍在原始日志保存，
`external_metrics.tsv` 另列 actual wall、runs、segments、临时读写字节数和
受控算法内存峰值。驱动要求 `mode=external`，且 runs、临时读和写字节数
均大于零，不能仅以指定 `-m 8G` 就宣称实测了外排。

暂不提交 BAM↔SAM：A4 的 D3 SAM 输出为 **532,628,139,096 bytes**，
六 ranks 的平均份额约 **82.7 GiB**，尚未加输入、body 扩容、工作区或模拟
输出副本就接近大共享空间。因此仅切换 rank-memory 不能使这两项在六核组
配置下安全运行；本轮不通过改 app、改成 spool 或增至更多 ranks 混入口径。

14:17 首项小回归 collate **8408851** 已启动，使用节点 `21992-21997`，
远端 SHA 与冻结二进制一致。它在 14:17:32 完成：mode=external，
2,668,351 条记录、12 runs、72 segments，临时写/读分别为
1,401,725,890 / 1,402,164,407 bytes，受控峰值最大 128 MiB。随后 fixmate
**8408854** 在 14:18:02 完成，记录数一致；剩余回归及正式 D3 尚待完成。
这些是回归状态，不是 D3 的正式性能结果。

基础四项及复读按 D2 准备/输出耗时与 D3 规模粗估约 **6～10 小时**。
完整选定矩阵暂按 **1～2 天**预留，需依据首轮 D3 外排实际临时 I/O 更新；
共享存储变慢可能延长或触发时限。外排仍须实际写/读临时文件，排除其报告
时间不代表可以省略真实执行时间。文件系统显示约 1.1 PiB 空闲，但 `quota`
没有返回账户额度，不能据此保证账户可写全部临时/输出数据；遇到配额错误
驱动停止并保留失败记录。会话关闭不影响独立后台服务。

### 2026-10-03 14:27 调整 D3 本批范围

小数据 collate、fixmate 以及两份输出的 `stats --basic` 完整复读均通过，
记录数均为 2,668,351。复读作业分别为 **8408859**、**8408868**。
外排 sort 探测 **8408872** 在 14:19 启动，六 rank 的输入加载完成，确认
`mode=external -m 128M`，但数分钟内没有核心完成摘要；约 14:24 只取消了
这一项探测。静态检查发现其环形交换结束消息有匹配风险，但没有运行时
调用栈，不能据此断言确切停点，也不能直接归因于服务器异常。

因此上述八项正式矩阵计划暂不全部执行，**本批实际安排六项各三次，共
18 项性能运行和 4 项 D3 首轮输出复读**：

| 顺序 | 实际安排 | 状态 |
|---|---|---|
| 1 | flagstat、stats --basic、BAM→BAM、MAPQ≥30，各三次及两份输出复读 | 基础矩阵已启动 |
| 2 | external collate `-m 8G`、fixmate `-m`，各三次及两份输出复读 | 等基础矩阵成功后自动续跑 |
| 暂缓 | external sort、markdup --stream | 小回归未通过，不提交 D3 作业 |
| 不提交 | BAM↔SAM | 六 ranks 内存容量不足，理由见上文 |

D3 第一项 `flagstat` **8408916** 于 **14:24:27** 提交，在节点
`21992-21997` 运行；14:27 日志显示 rank-memory 输入及 header 初始化完成，
尚处于输入分片准备/扫描阶段，没有正式核心成绩。collate/fixmate 独立续跑
服务为 `swbam-d3-collate-fixmate-after-simple-20261003-141600.service`，
它还会检查小回归四项成功记录；原完整外排阶段仍保持阻止状态。
续跑子目录名为 `apps_D3_n6_20261003_141600_collate_fixmate`，冻结驱动、
续跑脚本和阻止原因均保留在本轮根目录中。

本批时间预算调整为：基础矩阵约 **6～10 小时**，含外排 collate/fixmate
约 **12～24 小时**。这是按既有结果和规模推算的预留时间，不是已测结果；
首轮 D3 大文件加载、临时 I/O 或排队变慢时仍可能延长。两个后台服务独立
于 SSH 会话，关闭本地终端不会停止实验。只更新测试脚本与记录，没有改动
app/SDK 算法，也没有将取消的 sort 探测计入性能结果。

### 2026-10-04 00:30 D3 基础矩阵完成，collate 输出模拟失败

基础矩阵于 **2026-10-03 19:36:15** 全部完成：四项各三次，共 12 项性能
运行及 2 项首轮输出复读。仍为六节点、六活跃核组，rank-memory 输入、
memory rank body、level 1；没有将物理加载或最终 MPI-IO 写盘计入主成绩。

| 项目 | 第一次/s | 第二次/s | 第三次/s | 主成绩中位数/s | 4.3 中位数/s |
|---|---:|---:|---:|---:|---:|
| flagstat | 28.346675 | 28.605849 | 28.572994 | **28.572994** | 26.552981 |
| stats --basic | 26.276818 | 26.500211 | 26.528699 | **26.500211** | 24.531829 |
| BAM→BAM | 75.457844 | 75.695684 | 75.749789 | **75.695684** | 69.341413 |
| MAPQ≥30 | 75.336684 | 75.653971 | 75.558507 | **75.558507** | 69.654980 |

只读主成绩为 `444Complete the total body cost`；写出主成绩为
`Complete the total (4.1~4.6) cost`，不是单独 4.3，也不是提交 wall。
各项三次极差/中位数均小于 1%，核心成绩稳定。单次物理输入准备约
749～1,258 秒，BAM→BAM 最终写盘约 412～587 秒、过滤约 346～513 秒，
均单列并排除。全部基础项含复读的实际历时约 5 小时 12 分。

性能作业号：flagstat 为 8408916/8410391/8410443；stats 为
8409104/8410405/8410448；BAM→BAM 为 8409293/8410408/8410529；过滤为
8409595/8410431/8410630。三个回环输出大小均为 101,753,253,749 bytes，
保留记录数均为 1,210,633,370；三个过滤输出大小均为 91,171,148,189 bytes，
保留记录数均为 1,108,854,424。

首轮复读 8410152、8410361 完整解析成功，分别确认上述两种记录数，且
`SN is sorted = 1`。这里的总记录数取 `OptimizedStatsMPI total_records`，
不是排除 secondary/supplementary 后的 `SN raw total sequences`。
这验证可解析性、记录数和坐标排序，不据此声称逐字节与输入完全相同。

后续 D3 外排 collate **8410634** 于 10 月 3 日 19:37 提交，21:06:39
报退出码 1。六个 rank 的输入分片加载均成功，确实进入 external 模式，
`-m 8G`，实际 run arena 为 8,543,787,008 bytes；核心 4.3 已返回
**895.530010 秒**（实际算法 wall 为 4,111.443868 秒），rank body 为
6 个 memory、0 个 spool。随后 4.4、4.5a 完成，rank 0 打印 4.5b 后，
没有 4.6、总主时间、统计归约或最终输出完成行。

日志与 `CollateGatherOutput` 中某个非 0 rank 分配完整模拟输出副本失败
后 collective 返回的分支一致，疑似输出副本造成内存峰值超限；没有记录
失败 rank、请求大小或 errno，故暂不把具体内存上限作为已证实结论。
当前日志没有配额或临时文件 I/O 错误，也不支持直接归因于坏节点。
**这一失败运行的 4.3 仅作诊断，不纳入正式性能成绩。**

collate/fixmate 续跑服务已失败停止，fixmate 尚未提交；D3 sort/markdup
仍因小数据外排 sort 回归未通过而暂缓，SAM 转换仍未安排。00:30 查询
账户 `bjobs` 无作业，本次没有重新提交，也没有修改算法代码。建议先处理
输出模拟副本的内存峰值，再续测 collate/fixmate，并单独修复外排 sort。
本轮远端基础及失败 collate 的原始日志已同步到对应本地 benchmark 子目录。

### 2026-10-04 输出模拟副本容量修复与 D3 续测

本次只修改 app 的 collate/fixmate 输出模拟阶段，两处共用
`swbam-app/include/simulated_write_buffer.h`。SDK、CPE kernel、记录解析、
排序、交换、压缩及实际 rank body 未修改；不重新提交 D3 已完成的基础四项。

正常情况保留完整副本分配及一次 `ReadAt` 复制，D1/D2 不强制切换为窗口。
**仅当全量副本分配失败时**，复用 64 MiB 输出窗口，按源 offset 连续读完
本 rank 的全部 body，包括短尾块；不跳过复制，不改为丢弃 sink，也不在
计算期间回退读取磁盘。真实最终输出仍从原 rank body 写出，窗口只用于
模拟 SSD 写入，不作为最终 BAM。增加失败 rank、body/window 字节数和 errno
诊断；窗口也分配失败或 body 读取失败时仍报错并由驱动停止。

分配及回退尝试计入 4.5b，仍不计入主成绩；所有窗口复制计入 4.6，输入
物理准备及最终 MPI-IO 写盘仍排除。窗口重复使用会改变目的端缓存/页面
行为，不能把 D3 的窗口 4.6 当作与全量副本完全相同的内存访问模式；因此
日志明确标识 fallback，继续保留 4.3/4.6 分项及完整主时间。核心算法路径
没有改变，不以窗口模式更短的模拟复制时间声称算法优化。

`RabbitBAM-MPI` 构建及 `git diff --check` 通过。辅助类的独立回归测试为
`tools/simulated_write_buffer_regression.cpp`，使用有限分配器覆盖全量成功、
全量失败后窗口成功、两个分配均失败、全部字节/尾块、前后 envelope 边界、
空 body、读取失败、大小溢出及 reset；O2 与 ASan/UBSan 两种构建均通过。
它仅测试 app 辅助类，不新增 x86 应用构建或修改 SDK 示例。

```bash
g++ -std=c++11 -O2 -Wall -Wextra -Werror \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  tools/simulated_write_buffer_regression.cpp \
  -o /tmp/swbam-simulated-write-buffer-regression
/tmp/swbam-simulated-write-buffer-regression
cmake --build build_sunway_sduhpc --target RabbitBAM-MPI -j 4
```

新冻结二进制 SHA-256 为
`49999a9cad48bd0eeba9291fb77ae13f5d652cacbea11d531944ae972b47b475`；
源码基于已提交 0309267 加上述 app-only 修复，本次没有创建新提交。
运行根目录为 `build_sunway_sduhpc/benchmarks/`
`apps_D3_simwrite_fix_n6_20261004_005100/`，保留冻结 driver、二进制、配置及
smoke/formal 两个独立子目录，远端仍位于 `D3_apps_current/`。

小回归继续使用 WES_0.25G、六节点六活跃核组、rank-memory、memory body、
MPI-IO 输出、level 1，外排 collate `-m 128M` 接 fixmate `-m`。
00:52:58 提交的 collate **8411637** 已成功，主时间/4.3/4.6 分别为
1.185345/1.164424/0.012353 秒；fixmate **8411638** 主时间/4.3 为
0.550946/0.530439 秒。两者均为 2,668,351 records，未触发窗口回退。
与上一轮小回归的 1.470861/0.551128 秒相比，没有观察到主成绩回退，但
这是不同运行时段的各一次回归，不是受控三次性能对比，不能据此声称加速。
两份输出完整复读和旧版逐字节比较通过后，才批准 D3 正式续跑。

00:54:58 小回归全部完成，输出复读作业为 **8411640/8411642**，记录数
均为 2,668,351。随后远端 `cmp` 对 name.bam、fixmate.bam 均返回 0，与
上一版冻结小回归输出逐字节一致，输出大小分别为 281,864,351 和
290,796,959 bytes；批准记录保存在运行根目录的 `smoke_approved`。

D3 正式 collate 首轮 **8411648** 于 **00:57:58** 提交，00:58 开始在
节点 `23820-23825` 运行；远端二进制 SHA 与冻结版本相同，日志确认
rank-memory 和 6 MPE + 384 CPE，目前仍在输入准备阶段，尚未验证 D3
窗口回退和最终输出能够完整完成。后台服务为
`swbam-d3-simwrite-formal-20261004-005100.service`。

正式范围仍为 external collate `-m 8G`、fixmate `-m` 各三次，六节点
六活跃核组，首轮两份完整输出复读通过后才提交第二、三轮。全部命令和
结果位于 `_formal/commands.sh`、`results.tsv`、`external_metrics.tsv`，
任一失败自动停止，不改变 ranks 或退回磁盘核心 I/O。sort/markdup 未提交。
参考前次 collate 实际 wall，首项约 **1.5～2 小时**，这批预计 **8～12
小时**，需按首轮 fixmate 与存储状态更新；后台服务不依赖当前会话。

### 2026-10-04 02:39 D3 collate 首轮成功，fixmate 正在运行

正式首轮 collate **8411648** 于 **02:33:05** 被驱动确认成功，实际提交
wall 为 **5,706 秒**（约 95 分钟），完成 external 模式 `-m 8G`：

| 指标 | 首轮结果 |
|---|---:|
| 主成绩 4.1～4.6 | **841.383074 s** |
| 4.3 外排模拟算法时间 | 832.717749 s |
| 4.6 全部 body 的模拟写复制 | 6.872714 s |
| 物理输入准备，排除 | 1,164.564523 s |
| 实际算法 wall，含真实临时 I/O | 3,903.427085 s |
| 最终 MPI-IO 输出，排除 | 604.296278 s |
| total/received records | 1,210,633,370 / 1,210,633,370 |
| 输出大小 | 160,004,190,693 bytes |
| runs / segments | 54 / 324 |
| 合计 temp write / read bytes | 859,720,107,303 / 859,986,407,590 |
| tracked 算法内存峰值最大值 | 8,589,934,592 bytes |

rank 1～5 的全量模拟副本分配失败后，明确记录了 **67,108,864-byte**
窗口回退；每个 rank 仍复制约 26.66 GB 的全部 body。rank 0 保留全量
快路径。所有 rank body 均为 memory，没有切换为 spool；4.6、统计归约、
最终 MPI-IO 和 `TOT TIME` 均成功完成，证明本轮已越过原先失败的输出
模拟阶段。4.6 约占主成绩 **0.82%**，但窗口的缓存行为说明仍按上文保留。
这些是一次成功成绩，尚不替代三次中位数，亦不使用上次失败运行作加速基线。

fixmate 首轮 **8411876** 于 **02:33:05** 接续提交，02:39 查询仍在
节点 `23820-23825` 运行；日志确认 rank-memory、6 MPE + 384 CPE，尚处
输入准备阶段，没有错误或核心完成行。首轮 fixmate 和两份输出的完整复读
尚未完成，因此暂不宣称 D3 全部正确性验证已通过。

后台服务仍 active，驱动将在首轮及复读成功后续跑第二、三轮，本次查询
没有重复提交或修改代码。按首轮 collate wall 与原计划，剩余暂预留
**6～10 小时**，需据首轮 fixmate 和共享存储耗时更新。

### 2026-10-04 09:27 fixmate 内存 body 写入失败，续测停止

本地后台服务现为 failed，退出码 1。已保存日志显示 fixmate 首轮
**8411876** 在 **03:08:09** 报失败，实际总时间 2,096.282737 秒：
六 rank 的分片加载均成功，每 rank 输入约 26.667 GB，物理准备
1,837.577251 秒排除；4.1/4.2 完成后，六 rank 均报
`ERROR: BAM record writer failed to append BGZF block.`，随后 collective
清理退出。没有 4.3 完成、4.5b 窗口分配、4.6 或最终输出完成行。

因此此次失败发生在 **fixmate 核心中的压缩块追加到内存 rank body**，
不是已经修复的额外模拟副本分配阶段。错误由
`swbam-sdk/lib/cpe/swbam_cpe_record_write_adapter.cpp` 的批次输出检查打印，
可能来自块长度检查或 sink Append；日志未记录具体块、body 已用/容量或
errno，尚不能据此唯一确定根因。

静态检查发现 `MemoryRankBodySink` 超出 capacity 时按两倍 `realloc`。
fixmate 当前初始输出容量约等于输入分片 26.667 GB，输出增加 mate 标签后
可能触发向 53.334 GB 扩容；如果分配器需要先分配新区域再搬移旧数据，
输入 + 旧输出 + 新输出的瞬时需求约 106.668 GB，超过共享堆预算。
这是优先核查的内存峰值风险，不把该估算当作已测峰值或确认的 errno。
后续优先考虑仅在 app 初始化时预留输出余量，避免在核心阶段大规模搬移，
并增加一次失败诊断；不因此把 memory body 改成磁盘 spool。

目前 D3 collate 仍只有首轮成功成绩 **841.383074 秒**，两份 D3 输出完整
复读及第二、三轮均未提交。fixmate 没有有效性能成绩，上一节剩余耗时
估计已不适用。09:27 的两次 SWLS SSH 查询均连接超时，不能确认当前远端
队列或服务器整体状态，也不据此把凌晨的程序错误归因于网络/服务器。
本次只更新状态记录，没有修改算法或重新提交作业。

### 2026-10-08 12:03 连接恢复，修正 fixmate 预留后续测

SWLS SSH 和调度已恢复，提交前本账号队列为空。上次 collate 首轮
`out/name.bam` 仍为 160,004,190,693 bytes，因此保留作业 **8411648** 的
有效成绩，不重跑首轮，也不重跑已完成的 D3 四项基础测试。

仅修改 standalone fixmate 的 app 层初始化：当输入驻留内存且显式选择
memory rank body 时，初始输出容量由本 rank 压缩输入大小改为其两倍。
该值等于原先第一次扩容后的容量；预留仍计入 **4.2**，没有把正常分配
移出主成绩。对当前 D3 分片，输入加预留输出约 80 GB，避免第一次
`realloc` 搬移时输入、旧输出、新输出同时存在的约 106.7 GB 风险；这是
容量估算，不是实测峰值或上次 ENOMEM 的确证。若输出超过预留，仍沿用
原扩容行为。追加失败时在批次层记录 rank、body 已用/容量和 errno，成功
路径没有增加逐记录诊断。SDK、CPE kernel、算法、CLI、auto/spool 行为及
物理 preload/dump 排除口径均未改变。

构建 `RabbitBAM-MPI`、`git diff --check`、host 模拟写窗口回归通过。
新二进制 SHA-256：
`ecd6df4fa34645ac8bb1d1db06fdb27c5e2b5618b310a109589959006883a0e9`。
本次不创建 Git 提交，已有修改保留。冻结源码 patch、二进制、driver、
配置、小回归及失败前的完整 formal 日志/结果备份保存在
`build_sunway_sduhpc/benchmarks/apps_D3_fixmate_reserve_20261008_115900/`。

WES_0.25G、六节点六 ranks 的小回归于 **12:02:44** 完成：

| 项目 | 作业 | records | 主时间 / 4.3 (s) |
|---|---|---:|---:|
| external collate，`-m 128M -n 66` | 8463448 | 2,668,351 | 1.292027 / 1.271885 |
| fixmate，`-m` | 8463449 | 2,668,351 | 0.542349 / 0.522687 |
| collate 输出完整复读 | 8463454 | 2,668,351 | 0.086475 / 0.069774 |
| fixmate 输出完整复读 | 8463464 | 2,668,351 | 0.080225 / 0.072245 |

两份输出与 10 月 4 日冻结小回归输出远端 `cmp` 均返回 0，大小仍为
281,864,351 / 290,796,959 bytes。fixmate 日志确认初始化容量为两倍输入，
没有追加错误；其主时间与此前 0.550946 秒接近。但两者是不同时间各一次
运行，不据此宣称加速或已验证 D3 大容量分配。

**12:03:07** 以 `RESUME=1` 启动正式续跑，驱动明确跳过 collate 首轮，
首先提交 fixmate **8463470**；查询状态 RUN，六节点
`21768-21771,21972,21982`，日志确认 rank-memory、6 MPE + 384 CPE。
当前尚在物理输入准备阶段，没有 4.2/4.3 完成，不能宣称大数据修复已通过。
资源仍为 `-N 6 -np 1 -health 0 -mpecg 6 -cgsp 64 -share_size 100
-cross_size 88000 -cache_size 32 -xmalloc`，memory body、MPI-IO 最终输出、
level 1，collate 仍为 external `-m 8G -n 66`，诊断开关关闭。

正式结果与完整命令继续写到原 `_formal/results.tsv`、`commands.sh` 和
`external_metrics.tsv`。旧失败 fixmate 日志已改名保留为
`bench_fixmate_r1.failed_job8411876.log`，旧二进制在本地备份及远端
`RabbitBAM-MPI.before-reserve-20261008` 保留。collate 首轮使用旧 SHA，后续
使用新 SHA；本次未修改 collate 代码，版本区别在此明确记录。

后台服务 `swbam-d3-resume-20261008-115900.service` 在本地 user systemd
运行且启用 linger，不依赖当前会话。fixmate 首轮及两份 D3 输出完整复读
均通过后才进入第二、三轮，任一失败自动停止。sort/markdup 未提交。
按之前实际 wall 暂估首项 **45～75 分钟**、剩余整批 **6～10 小时**，
共享存储、排队和新容量分配结果仍有不确定性，核心计时排除物理 I/O，
但等待任务完成的实际时间仍包含这些 I/O。

### 2026-10-08 15:05 fixmate 首轮成功，输出复读异常等待

D3 fixmate **8463470** 于 **13:02:47** 被驱动确认成功，提交 wall 为
3,579 秒，日志有最终 `TOT TIME` 和正常结束：

| 指标 | 首轮结果 |
|---|---:|
| 主成绩 4.1～4.6 | **265.171089 s** |
| 4.3 核心算法 | 256.187366 s |
| 4.2 输出缓存初始化 | 0.000446 s |
| 4.6 模拟输出复制 | 7.193328 s |
| 物理输入准备，排除 | 2,620.142041 s |
| 最终 MPI-IO 输出，排除 | 664.670761 s |
| records | 1,210,633,370 |
| 输出大小 | 163,882,611,815 bytes |

六 rank 的初始 body 容量均成功预留到约 53.33 GB，实际 body 约
27.31 GB，没有再次出现压缩块追加失败；六 rank 的额外模拟副本均回退
到 64 MiB 窗口，仍复制全部 body 字节。说明该版本已完成 D3 fixmate
全流程，但这只是一次有效运行，不能替代三次中位数或输出完整复读。

**13:02:49** 开始的 collate 输出复读 **8463955** 截至 **15:05** 仍为
RUN，节点 `21768-21771,21621,21972`。本地及远端日志均只有
`111Complete the initialization cost`，没有 input open、header、分片准备
或核心完成行。结合当前代码位置，只能定位到输入打开及其 MPI 同步附近，
不能把这两小时描述为正常分片加载，更不能直接归因于算法、磁盘或坏节点。
登录节点对同一文件 `stat` 正常，首 64 KiB 的 `dd` 在 0.0317791 秒完成；
这只证明登录节点的小范围访问可用，不证明计算节点或全量存储吞吐正常。

后台服务仍 active，第二、三轮尚未提交，fixmate 输出复读也尚未开始。
已有 collate/fixmate 首轮成绩均保留，本次只查询和记录状态，没有取消或
重复提交作业。原整批耗时估计暂不能保证，需要先确认复读等待原因。

### 2026-10-08 15:19 异常验证停止，固定节点排队续测

进一步检查后，复读 **8463955** 仍只有初始化完成行。按源码位置，它
停在 `MpiBamInput::Open` 或附近同步，尚未进入 BGZF 全局扫描、分片加载
与 CPE 计算，因此不能解释为 D3 正常大文件读取耗时。`bjobinfo` 显示
六任务仍为 RUN，`cnload` 的计算网络和 GFS 状态正常；没有取得可用的
计算节点栈，尚不能唯一归因于某个 syscall、MPI 或某个坏节点，也没有
证据指向内存不足。与成功 fixmate 的节点组合相比，异常验证将 `21982`
替换为 `21621`，这仅是相关线索，不将 `21621` 标记为已确认坏节点。

**15:09:10** 仅终止异常验证 **8463955**，驱动随之停止。完整日志归档为
`validation_verify_collate_r1.stopped_job8463955.log`，collate **8411648**
和 fixmate **8463470** 的成功成绩及输出全部保留。没有取消其他用户任务，
也没有将未验证结果视为正确性检查通过。

新增的代码仅是 app 的 `stats` 进度诊断：`SWBAM_STATS_DIAG=1` 时逐 rank
打印 `input_open_begin/done`、`input_prepare_begin/done`、`core_begin/done`，
默认关闭。不改 SDK、CPE kernel、算法或计时排除项。本轮仅在**验证命令**
开启诊断，collate/fixmate 性能命令保持关闭。新二进制 SHA-256 为
`8a60cf8e52f395553131042cd1b50138a56b57f669a055ad1ef956bc042746aa`，
构建、`git diff --check` 和 Bash 语法检查通过。

当前原成功节点 `21768-21771,21972,21982` 被其他任务占用。另外尝试在
实时监控显示 idle 的 `44030-44031` 提交小 stats 检查，默认 health 和
health 0 均被调度器以 **ret=-34 / specified node is NOT ok** 拒绝，
没有生成作业 ID。其 `bqueues` 配置状态也与实时空闲显示不一致。这说明
实时 idle 不等于调度可分配，不把拒绝归因于 SWBAM 执行错误；这两次
提交未执行程序，不计入任何测试成绩。

为继续实验，driver 新增可选批作业模式 `SWBAM_APPS_BATCH=1`，默认仍
保留原交互式提交；明确节点列表支持逗号组合。批作业每 60 秒查询自身
作业状态，结束后读取远端完整日志并沿用原 records/core/`TOT TIME`
检查，不跳过正确性门槛。每项排队加运行上限为 6 小时 5 分钟，超限或
任一失败会停止并终止仅当前作业。mock 回归已覆盖 PEND→RUN→结束、异常
调度响应和缺失完成摘要，均通过；这不是实际硬件回归通过的替代。

**15:19:06** 开始新的 WES_0.25G 小回归，首项 external collate
**8465063** 已提交并处于 **PEND**，固定上述六节点。资源仍为六节点、
每节点一个 MPI rank、每 rank 一个活跃核组/64 CPE、`mpecg 6`、
`cross_size 88000`、memory body、MPI-IO 最终输出、level 1。小回归先跑
collate `-m 128M -n 66`、fixmate `-m` 及两份完整输出复读，全部通过
后才自动 `RESUME=1` 原 D3 formal：跳过两个成功首轮，重试两份 D3
输出验证，再补齐 collate/fixmate 第二、三轮。D3 collate 仍为 `-m 8G`。
sort/markdup、转换及已完成基础四项不在此次范围内。

冻结 driver、新二进制、源码 patch、原 formal 完整备份、拒绝提交日志和
mock 测试位于 `build_sunway_sduhpc/benchmarks/`
`apps_D3_stats_retry_20261008_151200/`。后台服务为
`swbam-d3-fixed-retry-20261008-151200.service`，不依赖当前会话；旧服务已
停止，新服务已确认 active。现阶段等待资源，不能声称已排除原等待问题。
资源等待时长未知；获得节点并通过小回归后，D3 剩余预计 **6～10 小时**，
需要根据后续两份复读及共享存储情况调整。

### 2026-10-08 15:30 用户要求不等固定旧节点，换组恢复运行

`bjobs` 默认只显示本账号的作业，不是全机节点占用清单。15:27 的
`cnload` 明确显示固定六节点中 `21982` 仍被 **8464522** 占用，其余五个
已经 idle，因此 **8465063** 仍需等待该特定节点。按用户要求终止这条
排队任务；它没有开始执行，不产生性能成绩。旧后台服务已随之停止。

进一步核对本机 `bsub -h`：`-health 0` 仅最佳节点，`1` 为最佳加正常
节点，`2` 为全部节点。实时 idle 不能证明节点满足提交时的健康等级等
条件；上节 ret=-34 的观察不应单独用于认定服务器故障。此次用
**health 1** 扩大选择范围，不用 health 2，也不对任何节点进行复位。

新组 `43872-43877` 的小数据 stats 检查 **8465160** 于 15:28 成功结束：
六节点、六 ranks、384 活跃 CPE，仍用 `mpecg 6`、`cross_size 88000`、
`share_size 100`、`cache_size 32`、`xmalloc` 和 rank-memory。六 rank 的
input open / prepare / core done 均齐全，records 为 **2,668,351**，
4.3 为 0.066557 秒，body 为 0.075221 秒，`TOT TIME` 10.008667 秒。
该项只是小数据正确性和节点可用性检查，不混入正式性能结果。
本次仅 driver 加入 `SWBAM_APPS_HEALTH_LEVEL`（默认 0），二进制仍为
前述 **8a60cf8e…**，没有改算法或内存模拟口径。Bash 语法、diff check
及 health 1 配置下的 driver mock 回归均通过。

**15:30:18** 直接接续原 D3 formal，跳过已经成功的 collate 和 fixmate
首轮；新的 collate 输出复读 **8465172** 已在 **43872-43877** 实际 RUN，
不是排队。日志显示六 rank 均 `input_open_done`，header 成本仅
**0.023529 秒**，已进入 `input_prepare_begin`。说明当前已越过前一次
停滞点，但仍未完成 D3 全量加载、计算和复读，不能宣称全部验证通过。
换节点后恢复是环境相关线索，不将它当作已确定原故障根因。

此次通过新组的小 stats 检查批准继续，无需重跑旧组尚未启动的四项
小回归。后续顺序仍为 D3 collate/fixmate 两份输出复读，再补齐两个
命令第二、三轮；任一失败停止。保持六节点六活跃核组、memory body、
MPI-IO 最终输出、level 1、collate `-m 8G -n 66`。验证诊断开启，性能
命令诊断关闭。明确记录首轮旧节点/health 0 与续跑新节点/health 1 的
差异；各轮实际节点见独立日志和 `commands.sh`，不能声称节点完全一致。

新服务为 `swbam-d3-new-nodes-20261008-153000.service`，已确认 active。
本次改用交互式提交，后续若固定新组被占用或运行失败则停止，不再无限
等待旧组。冻结 driver、配置、二进制、probe 完整日志及源码 patch 位于
`build_sunway_sduhpc/benchmarks/apps_D3_new_nodes_20261008_153000/`。
正式结果仍写原 `_formal/results.tsv`，预计剩余 **6～10 小时**，依共享
存储扫描/加载和最终写出耗时调整；这些实际等待仍排除在模拟核心成绩外。

### 2026-10-09 02:15 核查：D3 collate/fixmate 三轮完成

续跑于 **2026-10-08 22:28:18** 正常结束，`state=COMPLETE`。后台服务
`swbam-d3-new-nodes-20261008-153000.service` 已 inactive/dead，
`ExecMainStatus=0`；本次查询当前账户 `bjobs` 无作业。原 formal
`results.tsv` 共 **6 项性能结果和 2 项输出复读**，功能/轮次键没有重复。
八份成功日志均有算法汇总、`TOT TIME` 和调度器完成标记，没有程序
ERROR 或非零退出报告；旧失败、取消日志保留但不纳入成绩。

| 功能 | 第一次核心/s | 第二次核心/s | 第三次核心/s | 核心中位数/s | 4.3 中位数/s | 极差/中位数 |
|---|---:|---:|---:|---:|---:|---:|
| collate，外排 | 841.383074 | 829.791935 | 820.021205 | **829.791935** | 821.119045 | 2.5744% |
| fixmate -m | 265.171089 | 265.027532 | 265.062352 | **265.062352** | 256.124491 | 0.0542% |

继续采用六节点六活跃核组、rank-memory 输入、memory rank body、BAM
level 1 和最终 MPI-IO 输出。主成绩取 `Complete the total (4.1~4.6)
cost`，排除物理输入准备、4.5b 模拟输出缓冲分配及最终真实写盘。
collate 为 `-m 8G -n 66`，三轮均实际进入 external 模式；其 4.3 继续
沿用每 rank 实际算法时间减真实临时读写时间的模拟模型，临时内存复制
仍计入。不是磁盘端到端 wall，也不是所有项目使用同一原始墙钟边界。

| 功能/轮次 | 作业号 | 输入准备/s（排除） | 最终 MPI-IO/s（排除） | 提交至完成/s |
|---|---:|---:|---:|---:|
| collate 1 | 8411648 | 1164.564523 | 604.296278 | 5706 |
| fixmate 1 | 8463470 | 2620.142041 | 664.670761 | 3579 |
| collate 2 | 8466084 | 1204.479121 | 867.868804 | 6424 |
| fixmate 2 | 8467294 | 2134.234816 | 847.296486 | 3267 |
| collate 3 | 8467743 | 1319.029928 | 918.345152 | 6317 |
| fixmate 3 | 8468720 | 2073.366222 | 896.435544 | 3260 |

collate 三轮实际算法 wall 分别为 3903.427085、4314.593508、
4048.791867 秒。每轮共 54 runs、324 segments，临时写入
859,720,107,303 bytes、读取 859,986,407,590 bytes，受控算法峰值
最大 8,589,934,592 bytes；该峰值不包括输入、rank body 及全部运行时
内存。真实临时文件仍需执行，因此排除其报告时间不等于省略实际等待。

所有性能轮次 records 均为 **1,210,633,370**，collate 的 received
也一致；fixmate 均为 601,203,069 paired groups、0 singleton。
三轮输出大小分别固定为 collate **160,004,190,693 bytes**、fixmate
**163,882,611,815 bytes**，远端最后一轮文件大小与记录一致。六个 rank
的实际 body 均驻留内存，0 spool；collate 的五个非 root rank 和 fixmate
的六个 rank 使用 64 MiB 模拟输出窗口，全部 body 字节复制仍计入 4.6。

首轮输出完整复读 **8465172 / 8465743** 分别于 16:15:38、17:07:07
被驱动确认完成，两者 records 均为 1,210,633,370、`is sorted: 0`，
符合 QNAME 分组输出不保证 coordinate 顺序的预期。确认的是完整可解析
和记录数，不替代逐条 QNAME 连续性、mate 标签或全部内容的独立验证；
第三轮仅核对算法汇总和文件大小，没有新增全量复读。

各轮节点和二进制差异按前述记录保留，不能将此表称为同一节点组、完全
相同二进制的受控重测。A4 Samtools collate/fixmate 端到端中位数
1830.58 / 571.70 秒除以上述模拟核心中位数，数值比约 **2.21 / 2.16**，
只代表既定不同 I/O 边界的耗时比。

原始结果、完整命令、八份成功日志和外排指标继续保存在
`build_sunway_sduhpc/benchmarks/apps_D3_simwrite_fix_n6_20261004_005100/`
下的 `apps_D3_simwrite_fix_n6_20261004_005100_formal/`，远端对应
`wzs_data/swbam_bench/D3_apps_current/` 同名目录。最终 name/fixmate BAM
保留，可供后续 sort 接续，不重新生成已有阶段。

至此 D3 SDK 扩展性、基础四项以及 collate/fixmate 各三轮已有结果。
外排 sort、markdup 和 BAM↔SAM 仍未完成；没有完整四步 D3 核心合计。
下一步优先定位小数据外排 sort 的停滞，验证交换终止消息、空方向及最终
记录数和坐标顺序，再探测 D3 sort 的容量并接 markdup。暂不补跑独立
三轮全流程，也不为 BAM↔SAM 静默改成磁盘 spool 或增加核组数。
本次仅核查并更新结果文档，没有提交、取消作业或修改算法。

### 2026-10-09 02:49 外排 sort 修复与 D3 后两阶段续测

本次确认并修复两个程序缺陷，不将已复现的停滞直接归因于服务器：

- `MpiSortStrictExchangeRuns`：每条方向必须发送一次终止 header，包括
  空方向；已完成的发送/接收改用 `MPI_PROC_NULL`，接收完成标志保持
  单调。真实和模拟两套 header 均使用同一规则，避免不等长分块产生
  无匹配消息。原有 metadata/raw 字节交换及计时模型不变。
- `MpiSortStrictRunArena::can_append`：先检查新的 metadata 起点不早于
  已使用 raw 末端，再做无符号减法，避免回绕后错误允许写入、覆盖
  已存记录。空间不足时沿用现有 spill，不新增批次复制或改变排序键。

仅修改主核外排 sort 的容量判定、消息结束协议和可关闭的阶段诊断；
未修改 SDK、库流水线、CPE kernel、逐记录解析、压缩级别和过滤语义。
增加的容量比较不涉及排序/压缩热路径；本轮尚无 D3 修复前后完整成绩，
不能声称已量化证明性能无退化。

`RabbitBAM-MPI` 构建和 `git diff --check` 通过。
`tools/sort_exchange_regression.cpp` 增加实际外排 run 交换用例和 arena
边界检查，使用生产交换函数和 CPE 打包，覆盖空方向、非对称多块、
不同 rank 的 run 数及 payload/ordinal 完整性：

| 验证 | 作业号 | 结果 |
|---|---:|---|
| 旧终止协议、六节点 | 8470898 | 字节交换通过，在非对称 external case 1 停滞，人工取消 |
| 旧容量检查负例 | 8471109 | 六 rank 均触发 `arena metadata overlap` 断言，随后取消清理 |
| 修复版完整定向回归 | 8471129 | arena 检查、四项字节交换、四项外排交换共九项通过，退出 0 |
| 14,090 条小 BAM 外排 sort / markdup | 8470958 / 8470960 | 两份输出完整有序 raw record 列表与旧正确基线一致，4,270 条重复标记一致 |

14,090 条对照还检查 sort 前后 record multiset 一致及坐标顺序；输出
复读为 8470959 / 8470977。该小对照在 arena 边界补丁之前执行，后续
0.25G 复测则包含两个修复。中间单节点回归启动失败以及取消的 WES
探测均保留日志，不作为成绩；定向六节点负例不是正式性能结果。

关闭诊断后，用最终冻结二进制再次测试 WES_0.25G 的已有 fixmate BAM
（290,796,959 bytes、2,668,351 records），于 **02:48:53** 全部完成：

| 项目 | 作业号 | 核心/s | 4.3/s | 验证 |
|---|---:|---:|---:|---|
| sort，external，`-m 128M` | 8471181 | 0.659449 | 0.640691 | 18 runs、108 segments，真实临时文件确实执行 |
| 完整复读 sort | 8471185 | 0.076676 | 0.068474 | 2,668,351 records，`is sorted: 1`，duplicate 0 |
| markdup，`--stream -l 300 -m 8G` | 8471218 | 0.784874 | 0.766244 | marked 15,695，未删除记录 |
| 完整复读 markdup | 8471223 | 0.076334 | 0.068183 | 2,668,351 records，`is sorted: 1`，duplicate 15,695 |

WES sort 的输入准备 2.566985 秒、最终 MPI-IO 0.940475 秒，markdup
分别 2.224940 / 1.007359 秒；这些不计入核心成绩。此表仅为小数据回归，
不纳入 D3 汇总，也不能根据单次小数据数值推断 D3 正式中位数。

#### D3 续测配置

复用已生成的第三轮 fixmate 输出，不重新生成 collate/fixmate：

```text
/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data/swbam_bench/D3_apps_current/apps_D3_simwrite_fix_n6_20261004_005100_formal/out/fixmate.bam
```

输入大小 **163,882,611,815 bytes**，预期 **1,210,633,370 records**。
固定节点 `43872-43877`、`-health 1`、六节点六 ranks，每 rank 一个
活跃核组，共 6 MPE + 384 CPE；其余核组仅提供共享内存。
每项三轮，BAM level 1。完整实际提交命令保存在 formal `commands.sh`，
资源和应用参数为：

```bash
bsub -q q_share -I -b -K -N 6 -np 1 -node 43872-43877 -health 1 \
  -mpecg 6 -cgsp 64 -share_size 100 -cross_size 88000 \
  -cache_size 32 -xmalloc -timelimit 06:00:00 \
  RabbitBAM-MPI sort -i "$D3_FIXMATE" -o out/fixmate.sorted.bam \
  -m 8G -T tmp/sort --io-backend rank-memory \
  --rank-body-backend memory --io-output-backend mpiio --compress-level 1

# markdup 使用相同资源参数，各轮独立提交；不是上式命令的 shell 管道。
RabbitBAM-MPI markdup --stream -l 300 -m 8G \
  -i out/fixmate.sorted.bam -o out/markdup.bam \
  --io-backend rank-memory --rank-body-backend memory \
  --io-output-backend mpiio --compress-level 1
```

顺序为 sort1 → sort 输出完整复读 → markdup1 → markdup 输出完整复读
→ sort2 → markdup2 → sort3 → markdup3 → 最终 markdup 输出完整复读，
共 **6 项性能测试 + 3 项复读**。复读使用 `stats --basic`，要求完整
可解析、记录数一致和 `is sorted: 1`；重复标记数留在各独立日志。
只对冻结 benchmark driver 增加最后一轮复读，不改项目通用 driver。
此轮仅测两个后续阶段，不新增独立三轮完整四步流程。

主成绩继续取 `Complete the total (4.1~4.6) cost`。输入采用 rank-memory，
每 rank 只预加载自己的数据，排除物理扫描/加载和等待；正常 scan/split、
计算、通信和模拟内存读写按现有实现计入。排除 4.5b 模拟输出分配和
最终真实 MPI-IO。外排 sort 保留原有模拟临时 I/O、模拟 MPI 和读写重叠
计时模型，真实临时文件仍需执行；不能将含真实临时 I/O 的 `444` 墙钟
当作本轮主成绩，也不能把模拟成绩写成实测 SSD 端到端时间。
算法 tracker 的 `-m 8G` 峰值不包含输入、rank body 和全部运行时 RSS。

冻结二进制 SHA256：
`2dbc8f71194c0b000a82e18d5de6e9967e3f0b76daf9f84eec8a95396ee7abac`。
Git 基线 `0309267c6ef7c1b4809ab8a04ad37eb6d7c53c6d`，本次修复尚未提交；
源码 patch、二进制、驱动、批准记录、所有小回归日志均归档于：

```text
build_sunway_sduhpc/benchmarks/apps_D3_sort_tail_n6_20261009_022500/
```

正式结果位于其 `apps_D3_sort_tail_n6_20261009_022500_formal/`；远端
对应 `online/guoshi/wzs/wzs_data/swbam_bench/D3_apps_current/` 同名
formal 目录。BAM 输出留在 SWLS，不下载大文件到本地。

后台服务 `swbam-d3-sort-markdup-n6-20261009-022500.service` 于
**02:49:29** 启动，核查为 active/running。首轮 sort **8471267**
已调度到固定六节点并进入输入准备，尚无 D3 sort 核心结果；它同时
承担当前配置下的首轮容量探测。任何提交、内存、模式、记录数或排序
验证失败都会停止后续项，保留日志，不将失败结果计入中位数。
预计全轮预留 **8～14 小时**，即 10 月 9 日约 **11:00～17:00**，
随真实临时 I/O 和共享存储状态变化；物理等待被排除在成绩外并不意味着
不用等待。后台服务已脱离当前会话，本次确认启动后不持续轮询。

### 2026-10-09 11:41 核查：D3 sort 两轮、markdup 一轮完成，调度中断

后台服务于 **08:13:09** 以 `ExecMainStatus=1` 停止，formal
`state=FAILED`。不是程序仍在后台运行，也不是完整三轮已经完成。
已保存三项性能结果、两项首轮输出复读；成功功能/轮次键不重复。

| 项目 | 作业号 | 核心/s | 4.3/s | 输入准备/s（排除） | 最终写盘/s（排除） | 提交至结束/s |
|---|---:|---:|---:|---:|---:|---:|
| sort 1，external | 8471267 | **412.420477** | 405.949475 | 2485.551822 | 403.202436 | 7225 |
| sort 输出完整复读 | 8472142 | 29.355777 | 27.555239 | 1442.368085 | 0 | 1495 |
| markdup 1，stream | 8472267 | **623.526667** | 617.152828 | 965.390935 | 394.346491 | 2011 |
| markdup 输出完整复读 | 8472351 | 29.373421 | 27.576973 | 1323.452782 | 0 | 1401 |
| sort 2，external | 8472437 | **413.267653** | 406.817305 | 2414.656076 | 533.063014 | 7282 |

两轮 sort 核心极差相对两次均值约 **0.2052%**，目前仅有两轮，不能
当作正式三轮中位数。两轮均为 96 runs、576 segments，local/received
records 均为 **1,210,633,370**。实际算法 wall 分别 4302.216198 和
4299.467294 秒，临时写入 933,825,003,934 bytes、读取
934,106,945,464 bytes。算法受控峰值最大 8,589,934,592 bytes，
run arena 最大 5,138,861,248 bytes；不代表完整进程或节点 RSS。
两轮最终输出大小一致，为 **105,760,974,710 bytes**，核查远端现存
第二轮文件大小相符，但没有新增对第二轮文件的完整复读。

markdup 首轮完整处理 **1,210,633,370** records，marked 为
**42,383,102**，未删除记录，输出 **105,926,943,620 bytes**。
算法 tracker 峰值最大 1,008,910,694 bytes，pending 峰值最大
492,554,114 bytes。两项首轮复读均解析全部 1,210,633,370 records，
`is sorted: 1`；sort 输出 duplicate 为 0，markdup 为 42,383,102，
与算法 marked 汇总完全一致。该核对不替代与独立参考实现逐条 duplicate
决策的完整比对。所有成功性能日志均确认 6 memory body ranks、0 spool，
有 `TOT TIME` 和调度完成标记，未发现程序 ERROR；没有改变模拟计时口径。

**08:13:08** 第二轮 markdup 在提交时被拒绝：

```text
job submit failed, ret = -34, reason: The specified node to be allocated is NOT ok
```

没有生成作业 ID，也没有开始运行 markdup；此行不构成第二轮性能结果。
驱动按失败即停策略退出，没有继续提交 sort3、markdup3 或最终复读。
因此本次不能把服务失败解释为上述已完成算法的运行失败。

11:41 后核查 `cnload -c 43872-43877 -l`，六节点均 busy，正在执行
其他作业 **8474363**，两条高速网络均 OK、GFS 为 `1,2ok`，没有显示
down/softerr/hardft。本次实时占用不能倒推出 08:13 拒绝时的具体原因，
但说明现在不能立即在原固定组启动交互式作业；也不能仅凭 ret=-34 认定
节点硬件损坏。当前账户 `bjobs` 没有本轮 D3 作业；未干预其他任务。

尚需 **markdup2 → sort3 → markdup3 → 最终 markdup3 完整复读**。
沿用原二进制、参数和输出可 `RESUME=1` 跳过全部成功键，无需重跑
collate/fixmate 或成功的 sort 两轮。按现有实际提交耗时估算，剩余约
3.53 小时纯作业等待，建议预留 **4～6 小时，另加资源等待**；此前
11:00～17:00 完工窗口已因调度中断失效，不再表示当前预计结束时间。
本次仅核查、记录结果，没有修改程序、取消其他作业或重新提交实验。

### 2026-10-09 12:27 更换空闲节点，断点续测 D3 后两阶段

按用户要求继续提交。原组 `43872-43877` 仍被其他任务占用，选择同一
机柜内六个实时 idle、两条高速网络 OK、GFS `1,2ok` 的节点：
**43844-43845,43880-43881,44006,44017**，仍用 `-health 1`。
没有取消或干预其他任务，也没有使用 health 2、复位节点或更改程序。

新组先用原冻结二进制提交小 stats 检查 **8475329**，成功退出 0：
WES_0.25G fixmate 输入完整解析 **2,668,351** records，六 rank 均
完成 input open、rank-memory prepare 和 core。4.3 为 0.072194 秒，
body 为 0.080136 秒，`TOT TIME` 11.031122 秒。该项仅验证新组可用性，
不计入 D3 性能成绩。

正式二进制、driver SHA256 与上节冻结记录一致；sort `-m 8G`、markdup
`--stream -l 300 -m 8G`、BAM level 1、六节点六 ranks/活跃核组、
cross-size 88000、rank-memory 输入和 memory body 均不变。物理准备、
最终 MPI-IO 与模拟输出分配的排除口径不变，性能诊断关闭。唯一配置
变化为固定节点列表及 `SWBAM_APPS_BATCH=1`：按批作业提交，在两项
之间节点暂时不可用时允许等资源；每项资源等待加运行仍受原六小时
上限和 driver 超时保护，程序/验证失败则停止后续项。

启用 `RESUME=1`，确认五条已成功结果与续跑前快照逐字段一致，没有
重复功能/轮次键。首轮 sort、markdup、两份验证及 sort2 全部跳过。
续跑顺序仍为 **markdup2 → sort3 → markdup3 → 最终完整复读**，继续
写原 formal `results.tsv`、`external_metrics.tsv` 和独立日志。
原失败提交及其他日志、结果、二进制完整本地快照位于原 root 下的
`before_resume_20261009_122600/`，不将失败日志覆盖后误认为未曾中断。

后台服务 `swbam-d3-tail-resume-n6-20261009-122600.service` 于
**12:27:03** 启动，核查 active/running。第二轮 markdup **8475352**
已于 12:27 调度到新组六节点，状态 **RUN**，不是等待原忙碌节点。
输入复用已有第二轮 `fixmate.sorted.bam`，大小核对为
105,760,974,710 bytes；尚无 markdup2 核心成绩。

续跑配置、节点检查、probe 完整日志和批准记录归档于同一 root 的
`resume_20261009_122600.sh`、`resume_20261009_122600.txt`、
`logs/new_nodes_probe_20261009_122600.log`；服务输出为
`resume_20261009_122600.service.log`。配置及批准记录也已上传远端归档。
新旧节点组的差异明确保留，不将最终三轮表声称为完全同一节点组的
受控重测。

依据已完成的实际等待时间，剩余四项约 3.53 小时，预留 **4～6 小时**，
预计 10 月 9 日约 **16:30～18:30** 完成；不含额外排队或异常存储
延迟。后台运行不依赖当前会话，本次确认启动后不持续监测。

### 2026-10-09 15:34 核查：markdup 第二轮完成，sort 第三轮运行中

续跑服务仍为 active/running，原 formal `results.tsv` 已有四项性能
结果和两项首轮完整复读结果。第二轮 markdup **8475352** 于
**13:19:41** 完成，没有程序 ERROR；新增成绩如下：

| 项目 | 核心/s | 4.3/s | 输入准备/s（排除） | 最终写盘/s（排除） | 提交至结束/s |
|---|---:|---:|---:|---:|---:|
| markdup 2，stream | **623.514428** | 617.104938 | 1881.386726 | 584.033538 | 3155 |

该轮 records/kept 均为 **1,210,633,370**，marked 为 **42,383,102**，
输出 **105,926,943,620 bytes**，均与首轮一致；6 memory body ranks、
0 spool。两轮核心时间相差 0.012239 秒，约 0.002%，尚不作为三轮
中位数，也不替代第三轮输出的最终完整复读。

第三轮 sort **8475750** 于 13:19:42 提交，排队后 **13:49:00** 开始，
核查状态 RUN，仍使用上节新组六节点。六 rank 均已成功完成输入准备
（`ok=1`），排除的最大准备时间为 **2936.181341 秒**；扫描
6,606,356 个 BGZF blocks，4.1 为 1.864731 秒，4.2 为 0.000064 秒，
选中 external，`-m 8G`。当前日志尚未输出 4.3 或最终核心成绩，未见
ERROR；关闭诊断时算法阶段没有持续进度打印，不能据此认定停滞。

剩余 **当前 sort3 → markdup3 → 最终 markdup3 完整复读**。参考前两轮
sort 实际算法约 72 分钟，以及本轮准备耗时，预计从本次核查起还需
**约 1.5～3 小时，约 17:00～18:30 完成**，额外排队或存储变慢会顺延。
本次仅核查和更新文档，没有新提交或取消作业；后台驱动继续运行。

### 2026-10-09 18:08 三轮性能完成，更换节点运行最终复读

第三轮 sort **8475750** 于 **16:03:17** 被驱动确认完成，第三轮
markdup **8476555** 于 **17:21:04** 完成。原 formal `results.tsv`
已有六条性能记录与两条首轮完整复读记录，未发现重复功能/轮次键：

| 功能 | 第一次核心/s | 第二次核心/s | 第三次核心/s | 核心中位数/s | 极差/中位数 |
|---|---:|---:|---:|---:|---:|
| sort，外排 | 412.420477 | 413.267653 | 413.132094 | **413.132094** | 0.2051% |
| markdup，stream | 623.526667 | 623.514428 | 623.691306 | **623.526667** | 0.0284% |

第三轮 markdup 的 records/kept 均为 **1,210,633,370**，marked 为
**42,383,102**，removed 为 0；输出 **105,926,943,620 bytes**，
与前两轮一致。以上仍采用原内存模拟计时口径；外排 sort 为模型时间，
不能将三轮波动小解释为时间模型误差同样小。

最后的 `stats --basic` 完整复读 **8478709** 从 17:21 起处于 PEND。
核查 `cnload` 后，原固定组中的 `43844-43845,43880-43881` 已被
**8476617** 占用，调度器报告 `RMS_NODE_UNAVAIL_ERROR`。
按用户要求停止旧续跑服务并仅取消此未开始的验证作业；没有重跑性能
命令，也没有取消其他用户或当前账号的无关任务。

换节点过程中，实时空闲的 `22336-22341` 在提交前被 **8478888**
抢先分配，小检查 **8478889** 因此仍排队，已取消且未执行程序。
随后全队列只显示三个 idle 节点 `43238,44030-44031`，三节点小检查
又被调度器以 `ret=-19 / No enough compute nodes` 拒绝，没有作业 ID。
这些探测均不作为正确性通过或性能成绩；idle 显示不能代替实际分配。

最终保持原六节点、六 ranks/活跃核组、`health 1`、跨核组大共享内存
配置与 `rank-memory` 输入，去掉固定节点限制，由调度器选择空闲健康
节点。没有使用 health 2、复位节点或改动算法/二进制/计时实现。
新验证 **8478910** 于 **18:04:57** 提交、**18:04:58** 实际开始，
状态 **RUN**，分配到 **22338-22343**。截至 18:08 六 rank 均已
`input_open_done` 并进入 `input_prepare_begin`，未见程序 ERROR。
当前仅确认已进入正常输入准备阶段，尚未完成最终全量解析验证。

新后台服务为 `swbam-d3-final-validation-auto-n6-20261009-175600.service`，
核查 active/running；`RESUME=1` 明确跳过全部八条已有成功结果，只执行
最后的 `validation_verify_markdup_r3`。六小时排队/运行上限、输出记录数
与排序状态检查、失败即停策略保持不变。后台运行不依赖当前会话。

原八条结果、命令、日志、state 与二进制已在原 root 的
`before_validation_node_switch_20261009_175600/` 快照保留；新续跑脚本、
批准记录和服务日志分别为 `resume_validation_20261009_175600.sh`、
`validation_resume_approved_20261009_175600`、
`resume_validation_20261009_175600.service.log`。切换前快照与新配置已上传
原远端 formal 下的 `audit_validation_node_switch_20261009_175600/`。
新最终验证仍使用原 `logs/validation_verify_markdup_r3.log` 并在通过后
归入原 `results.tsv`。A4 补测未受此次操作影响。

参考首轮 markdup 完整复读约 23 分钟，当前任务预计从启动起约
**25-40 分钟**，主要等待共享存储的物理扫描/加载；相应准备耗时继续
排除在模拟核心时间外。存储变慢会延长，最终验证通过前不将整批标记
为 COMPLETE。BAM/SAM 转换不在此次续跑范围内。

### 2026-10-09 18:54 核查：D3 sort/markdup 整批完成

最终完整复读 **8478910** 于 **18:34:15** 被驱动确认成功，提交至结束
为 **1,758 秒**（29 分 18 秒）。原 formal `state=COMPLETE`，续跑服务
`swbam-d3-final-validation-auto-n6-20261009-175600.service` 已正常退出，
`Result=success`、`ExecMainStatus=0`，不是后台仍在计算或等待资源。

| 最终 markdup 输出复读检查 | 结果 |
|---|---:|
| 完整解析 records | **1,210,633,370** |
| coordinate 排序状态 | **is sorted: 1** |
| duplicate records | **42,383,102** |
| 复读核心 body/s | 29.353531 |
| 4.3/s | 27.569111 |
| 物理输入准备/s（排除） | 1675.242235 |
| 进程总时间/s | 1711.682093 |

完整日志有唯一算法汇总、`TOT TIME`，未发现程序 ERROR 或非零退出。
记录数与 duplicate 数均与三轮 markdup 汇总一致；确认的是完整可解析、
记录数、坐标排序与重复标记数量，不将其表述为全部 record 内容、mate
标签或重复判定语义与 Samtools 逐条一致的独立验证。

原 formal `results.tsv` 现有 **9 条**唯一记录：sort/markdup 各三轮，
首轮 sort/markdup 输出复读，以及最终第三轮 markdup 输出复读。切换前
八条结果逐字段保持不变，没有重复提交性能轮次，最终验证不作为新的
性能样本。此前各算法三轮中位数保持不变。

至此，本批 D3 外排 sort/stream markdup 性能与计划内输出复读均结束；
先前基础四项与 collate/fixmate 也各有三轮结果。**BAM/SAM 双向转换
仍未测试**，也没有新增同一轮次完整四步串联的统一计时成绩，不能将
此次 COMPLETE 标记解释为 D3 实验方案全部项目均已完成。
本次仅核查、校验现有结果和更新文档，没有新增或取消 SWLS/A4 作业。

### 2026-10-09 D3 BAM/SAM 转换补测启动

本批补测 BAM→SAM 和 SAM→BAM，各三次，保持六节点、六活跃核组
（6 MPE + 384 CPE）。不重跑此前完成的八项应用，也不修改 SDK 或 CPE
kernel。应用的普通转换路径保留；只有设置
`SWBAM_CONVERSION_MEMORY=<chunk-size>:keep|discard` 才启用分块适配。

运行根目录：
`build_sunway_sduhpc/benchmarks/apps_D3_conversions_n6_20261009_193400/`。
后台服务：`swbam-d3-conversions-n6-20261009-193400.service`。
冻结二进制 SHA-256：
`71acde4f9eb128b50fb801b022873c4bac88d882cf1491faa0363592b96e9941`。
根目录保存冻结驱动、编排脚本、转换源码快照及 `submit.log`；各子批次
独立保存 `commands.sh`、`results.tsv`、日志和逐 rank 的
`conversion_chunks.tsv`。远端子批次在
`online/guoshi/wzs/wzs_data/swbam_bench/D3_apps_current/` 下。

#### 有界内存和真实执行

- BAM→SAM 仍只预加载本 rank 的 BAM 范围；按完整 BGZF block 将该范围
  切成不超过 1 GiB 的输入块。每块完整调用现有 BAM read→SAM write
  流水线，保留本块 SAM，不保留整个 rank 的 SAM。
- SAM→BAM 每个 rank 只读取自己完整行边界内的 SAM 范围；每次物理加载
  不超过 1 GiB，非末块截至最后一个换行符，尾部字节由下一次重新加载。
  每块完整调用现有 SAM read→BAM write 流水线，BAM 压缩级别仍为 1。
- 实际输出 body 使用内存 writer，跨块复用其容量；每块 CPE 全部结束后，
  将本块全部输出字节复制到 64 MiB 循环窗口，累计真实 memcpy 时间和字节数。
  不用抽样或规模外推代替完整 D3 转换。
- 首轮 `keep` 在本块计算和模拟写完成后保存验证 body，并最终 MPI-IO
  写出 SAM/BAM；实际存储 I/O 不与本块 CPE 计算重叠。后两轮 `discard`
  仍执行完整转换和全部输出 memcpy，但只保留逻辑大小，不保存/写盘。
  不落盘轮次必须与首轮实际文件大小一致。

分块会增加流水线初始化次数，并可能改变 SAM→BAM 的压缩块打包边界；
这些初始化、打包和压缩开销均保留在成绩中，不能称为完全没有性能影响。
普通 D1/D2 历史结果不替换为本批分块成绩。

#### 计时和验证门槛

主成绩继续取 `Complete the total (4.1~4.6) cost`：全局分区/广播、本地
reader/writer 初始化、每 rank 各块转换时间之和的最大值、输出布局和
每 rank 全部输出 memcpy 时间之和的最大值。不是逐块取最大值再相加。
输入物理扫描/加载、验证 body 保存/刷新、循环模拟窗口分配以及最终
MPI-IO 写盘排除。输出复制是在 CPE 已结束后测量的，避免简单减去与
计算重叠的真实 I/O 时段而同时删掉计算时间。正常 CPE 操作、MPE
批次处理、内存 writer 增长和块边界选择仍计入。

`4.3`、逐 rank 的 core/sim-write/preload/persistence、最大输入/输出块、
输入/输出/复制字节数另行保存。较宽的 body/TOT 诊断仍会受实际存储及
rank 等待影响，不替代上述主成绩；不能据此声称与真实 SSD 端到端等价。

正式矩阵前依次执行：WES_0.25G 普通转换、8 MiB 分块转换、往返 BAM
再导出 SAM；三份完整 SAM SHA-256 必须相同。再运行小数据两轮以验证
第二轮不落盘，并用 D1 普通/1 GiB 分块各一次检查性能；若任一转换
回退超过 20%，停止正式提交，此阈值只是防止严重回退，不表示已达成
5% 的性能等价结论。小数据与 D1 均不加入 D3 正式中位数。

D3 首轮 SAM→BAM 后完整运行 `stats --basic`，要求 1,210,633,370 条
记录及 `is sorted: 1`，通过后才运行第二、三轮。每次正式转换还要求
完整记录数、六 rank 汇总、每 rank copied-bytes=output-bytes 和正常退出。
失败则保留日志并停止后续提交；小数据内容一致及 D3 解析/计数/排序
不等同于已对全量 D3 做逐 record 内容比对。

#### 正式命令

统一启动参数仍为：

```bash
-N 6 -np 1 -health 1 -mpecg 6 -cgsp 64 \
  -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
```

首轮转换用 `-x SWBAM_CONVERSION_MEMORY=1G:keep`，后两轮用
`-x SWBAM_CONVERSION_MEMORY=1G:discard`；命令体分别为：

```bash
RabbitBAM-MPI run_all -i "$D3_BAM" -o out/initial.sam \
  --io-backend rank-memory --rank-body-backend memory \
  --rank-body-temp-dir tmp --io-output-backend mpiio --compress-level 1

RabbitBAM-MPI run_all -i out/initial.sam -o out/sam2bam.bam \
  --io-backend rank-memory --rank-body-backend memory \
  --rank-body-temp-dir tmp --io-output-backend mpiio --compress-level 1
```

`$D3_BAM` 为 D3 初始 coordinate BAM，不使用 markdup 输出代替。
首轮验证存储预计需要约 1.1 TB 的额外峰值（SAM 验证 body 与最终 SAM
同时存在），物理读写可能主导等待时间。当前文件系统空闲不代表已确认
账户配额足够；配额失败时停止，不自动删除既有数据。首轮输出是后续
SAM→BAM 的真实输入，不能为了省磁盘而在首轮就跳过 SAM 保存。

本批已通过 Sunway 构建、native 适配层回归、脚本语法检查和
`git diff --check`。19:33 已开始真实小数据普通 BAM→SAM；正式 D3
尚待前置检查，不将启动状态当作六项性能结果均已完成。

小数据初次复读 **8479311** 正常退出、完整记录数正确，但驱动按
coordinate 输入要求 `is sorted: 1`，而 WES_0.25G 本身为未排序输入，
实际返回 0。因此初次自动提交在验证门槛处停止，不是转换或服务器
报错。已将驱动的预期排序状态参数化：小数据要求 0，D1/D3 仍要求 1。
原驱动/日志保存在 `archive_initial_driver/`，二进制不变；保留两项已有
普通转换结果，19:36 从复读继续，**8479318** 通过。分块两项
**8479322 / 8479324** 也已完成，记录数均为 2,668,351，逐 rank 全部
输出字节均经过实际模拟写复制；其完整复读 **8479325** 通过。8 MiB
只是小数据边界回归配置，不能将其核心时间并入 1 GiB 正式成绩。

#### 20:12 前置检查与续跑

小数据普通输出、8 MiB 分块输出和往返 BAM 再导出的完整 SAM
SHA-256 均为
`df653f8f35d2d251325c10cf75d0723941ed28cda77cf756b5f51b8308220752`。
对应复读均通过，记录数为 2,668,351；第二轮不落盘测试
**8479363 / 8479365** 也通过，实际执行了全部转换及输出复制，逻辑大小
与首轮文件相同，dump 时间为 0。

D1 普通对照 **8479366 / 8479406** 两项完成，其往返输出复读
**8479419** 确认 177,798,571 条记录且 `is sorted: 1`。

| D1 单次对照 | 普通核心/s | 1 GiB 分块核心/s |
|---|---:|---:|
| BAM→SAM | 15.868810 | 13.362721 |
| SAM→BAM | 11.893100 | 11.997753 |

分块 BAM→SAM **8479426** 于 20:06:50 正常完成，SAM 文件为
68,705,660,624 字节，与普通对照大小相同。驱动因交互日志末尾的 CR
误判 `keep_output=1`，在登记结果时停止；已修正解析、从完整原日志
恢复该结果和六 rank 诊断，没有重新执行该性能轮次，也未更换二进制。
恢复过程保存在 `recovery_d1_bam2sam_8479426.txt`；修正版源码/驱动归档为
`conversion_sources_driver_v3.tar.gz`。修正后的冻结驱动 SHA-256 为
`ff4ef6f2299f35a11a7de09f189ab81eea8d80c282d6c49e33e84f5d0dafe415`。

20:12 续跑作业 **8479533** 为 D1 分块 SAM→BAM，不是正式 D3 样本；
20:16:52 完成，记录数正确，单次核心时间比普通对照增加约 0.9%。
分块 SAM→BAM 的 `results.tsv/load_s` 只保存原 222 阶段，本次为
4.066837 秒；另有逐块物理加载最大值 116.777998 秒，单独保存在
`conversion_chunks.tsv` 和 excluded-preparation 日志中，两者均排除在
核心时间外，不能把 `load_s` 一列解释为全部物理加载。
通过其完整输出复读及性能门槛后，后台自动提交正式 D3 矩阵。D1 前置
输出暂保留以支持断点续跑，额外约 160 GB，本批整体存储预留约 1.3 TB。
当前预留约 **6-12 小时**完成后续全部补测，主要受真实输入准备、首轮验证保存
和输出复读限制；物理 I/O 变慢或排队会延长，不计入主性能成绩。

#### 20:19 正式 D3 矩阵已启动

D1 分块往返输出完整复读 **8479554** 于 20:19:32 通过：
177,798,571 条记录、`is sorted: 1`。性能门槛文件
`d1_timing_gate.txt` 显示 BAM→SAM 比值 0.8421、SAM→BAM 比值 1.0088，
均未触发回退停止条件；这只是单次配对探测，不是新的 D1 三轮中位数。

正式首轮 D3 BAM→SAM **8479563** 于 **20:19:35** 提交并运行，当前
节点为 `21634-21635,21623-21624,21631,21636`。日志已确认
`rank-memory` 输入、BAM→SAM 模式及 6 MPE + 384 CPE，暂无程序报错；
仍在输入准备阶段，不能声称这条 D3 命令已完成。

正式目录：
`apps_D3_conversions_n6_20261009_193400_formal/`。后台总状态为
`RUNNING formal`，服务仍活动；后续自动执行首轮 SAM→BAM、完整复读，
验证通过再进行两项第二、三轮 `discard`。任一失败会停止后续提交，
完整六项结果及正式三轮中位数尚待完成。可结束本次聊天，后台继续运行，
不需要持续监测；当前按约 **6-12 小时**预留，约至次日凌晨/上午。

#### 2026-10-10 00:24 D3 转换三轮全部完成

最后一项 SAM→BAM **8481158** 已正常退出，驱动于
**2026-10-10 00:24:05 CST** 标记 `COMPLETE`。顶层与 formal 子批次
状态均为 `COMPLETE`，后台服务已退出，`Result=success`、
`ExecMainStatus=0`。正式矩阵从 10 月 9 日 20:19:35 开始，
实际等待约 **4 小时 4 分 30 秒**，包括真实输入准备、首轮保存及完整复读。

主成绩仍为内存模拟 SSD 的 `Complete the total (4.1~4.6) cost`，
保持六节点、六活跃核组、1 GiB 输入块；首轮 keep、后两轮 discard。

| 项目 | 三次主时间/s | 中位数/s | 三次 4.3/s | 逻辑输出字节 |
|---|---|---:|---|---:|
| BAM→SAM | 88.818279 / 88.825900 / 88.813831 | **88.818279** | 68.420337 / 68.414813 / 68.428195 | 532,628,139,096 |
| SAM→BAM | 92.870567 / 92.858387 / 92.858332 | **92.858387** | 89.151889 / 89.156350 / 89.155758 | 101,753,708,292 |

| 项目/轮次 | 作业号 | 222 输入准备/s（排除） | 最终写盘/s（排除） | 提交至结束/s |
|---|---:|---:|---:|---:|
| BAM→SAM 1 | 8479563 | 1272.044829 | 2870.706983 | 5484 |
| SAM→BAM 1 | 8480629 | 36.909973 | 526.204263 | 2288 |
| BAM→SAM 2 | 8480829 | 1379.156367 | 0.000000 | 1497 |
| SAM→BAM 2 | 8480958 | 0.008734 | 0.000000 | 1405 |
| BAM→SAM 3 | 8481086 | 1078.381501 | 0.000000 | 1193 |
| SAM→BAM 3 | 8481158 | 0.009197 | 0.000000 | 1451 |

SAM→BAM 的 222 列不是全部物理加载，其逐块 preload 每 rank 累计
最大值分别为 **1373.443062 / 1291.395022 / 1334.892111 s**，也排除
在主成绩外。BAM→SAM 已一次预加载 rank 输入，逐块 preload 为 0；
首轮验证 body 的 persistence 另外记录，不计入主成绩。

正确性核查：六项均汇总 **1,210,633,370** 条 records/kept、六 ranks
且正常退出；36 条逐 rank 诊断均满足 copied-bytes=output-bytes，
keep/discard 标记与轮次一致。同一转换的三轮逻辑输出大小完全一致，
且与首轮实际 SAM/BAM 文件大小相同。首轮往返 BAM 的完整复读
**8480741** 通过：全部记录数正确、`is sorted: 1`；验证不并入性能中位数。
此结论仍不代表完成了全量 D3 的逐 record 内容一致性比对。

完整七条结果（六条性能、一条验证）、36 条逐 rank 诊断、实际命令及
各项日志已归档到本地 formal 目录的 `results.tsv`、
`conversion_chunks.tsv`、`commands.sh` 和 `logs/`。此前 D3 应用成绩
不重测、不覆盖；本次核查没有修改算法、取消或重新提交作业。

### 2026-10-10 D2/D3 普通四步串联重测启动

按用户要求新增两组真正逐轮执行的普通串联：
`initial BAM -> collate -n 66 -> fixmate -m -> sort -> markdup --stream`。
每个数据集三轮，每轮都从其初始 coordinate BAM 重新生成中间文件，
不使用之前独立应用的阶段输出，不使用全内存 `dedup-pipeline`。

运行根目录：
`build_sunway_sduhpc/benchmarks/apps_D2_D3_workflows_n6_20261010_003937/`。
当前后台服务：`swbam-d2-d3-workflows-n6-20261010-003937-resume.service`。
顺序为小数据 memory 串联、小数据 external 串联、D2 三轮、D3 三轮；
任一阶段失败停止后续提交，避免两个数据集同时争抢共享存储。
正式 D2/D3 不因对方失败自动跳过验证或改变算法模式。

| 参数 | D2 | D3 |
|---|---|---|
| 初始记录数 | 365,079,492 | 1,210,633,370 |
| 节点/活跃核组 | 6/6，6 MPE + 384 CPE | 同左 |
| collate | `-n 66 -m 64G`，memory | `-n 66 -m 8G`，external |
| sort | `-m 64G`，memory | `-m 8G`，external |
| markdup | `--stream -l 300 -m 64G` | `--stream -l 300 -m 8G` |
| 输入/body/最终输出 | rank-memory / memory / mpiio | 同左 |
| BAM level | 1 | 1 |
| 预期 marked（已验证基线） | 8,894,488 | 42,383,102 |

统一资源前缀为 `-N 6 -np 1 -health 1 -mpecg 6 -cgsp 64
-share_size 100 -cross_size 88000 -cache_size 32 -xmalloc`，诊断关闭。
之前转换使用的两个节点已被其他任务占用，因此本轮由调度器选择满足
health 1 的空闲节点，各命令实际节点以其作业记录为准，不假称固定同组。
二进制复用刚完成转换回归的冻结版本，SHA256 为
`71acde4f9eb128b50fb801b022873c4bac88d882cf1491faa0363592b96e9941`。
不修改 app 算法、SDK、CPE kernel 或现有阶段计时实现。

主指标为每轮四个阶段 `Complete the total (4.1~4.6) cost` 相加，
再对三轮合计取中位数，**不是四个独立阶段中位数相加**。
继续排除物理输入扫描/加载、模拟输出窗口分配及最终真实 MPI-IO；
实际内存读取、模拟写 memcpy、算法/通信仍按已有实现计入。
D3 外排仍使用既定的临时 I/O/通信模拟模型，真实临时文件必须执行，
且验证 `mode=external`、runs 和临时读写字节为正。不把真实磁盘
读写改称内存实测，也不把核心合计称为真实 SSD 端到端时间。

完整原始阶段结果保存在各子目录的 `results.tsv`，每轮合计和三轮
中位数另保存为 `workflow_results.tsv`、`workflow_summary.json`。
额外 wall 诊断是四个阶段提交至结束时间之和，不含验证，不冒充包含
所有驱动间隔的整个流程 wall。输出仍真实写盘用于下一阶段输入，
最后一轮四份 BAM 留在本轮独立目录，不删除先前实验或原始数据。

首轮和末轮完整复读最终 markdup 输出，要求全部记录数正确、坐标有序、
duplicate flag 计数相同；每轮 marked/removed 还须与已验证基线一致。
D3 首轮另复读 collate、fixmate、sort 中间输出，正确性复读不计成绩。
这些检查不等于已验证全量逐记录 duplicate 标记集合。
正式矩阵共 24 项性能和 7 项复读；小数据 4+1/4+4 仅作回归。

脚本语法、合计统计自检及 `git diff --check` 通过；自检专门区分
median(sum) 与 sum(median)，并拒绝缺项、重复项或错误 marked 数。
本次仅补测试编排、末轮复读与汇总脚本，冻结在 `driver.sh`、
`orchestrator.sh`、`summary.py`，配置/哈希为 `plan.json`、`metadata.sha256`。
初始服务于 00:40 启动。两种小数据串联回归均已通过：

| 小数据回归 | collate/sort 上限 | markdup 上限 | 四阶段核心合计/s | 最终完整复读 |
|---|---|---|---:|---|
| memory | `1G` | `1G` | 2.774881 | **8481371**：2,668,351 records，sorted=1，duplicate flags=15,695 |
| external | `128M` | `8G` | 3.296204 | **8481398**：同上；首轮中间输出复读也通过 |

external 小回归最初误把用于强制 collate/sort 外排的 `128M` 也传给
markdup，作业 **8481382** 在基础窗口内存检查处正常报错退出。
这是测试编排参数错误，不是服务器故障；当时尚未启动正式 D2/D3。
失败日志及初始配置保存在 `archive_smoke_markdup_8481382/`。
仅将小回归 markdup 恢复为 `8G` 后断点续跑，作业 **8481395** 通过；
已成功的小回归阶段没有重测，失败样本不计成绩，正式 D2/D3 参数不变。
上述小数据合计只是单轮回归门槛，不是正式 D2/D3 性能结果。

正式 D2 首轮 collate **8481404** 于 **00:50:43 CST** 提交，核查为
`RUN`，实际节点 `26980-26983,21623-21624`。命令确认 `-n 66 -m 64G`、
rank-memory 输入、memory body、mpiio 输出及 6 MPE + 384 CPE。
当前尚在输入准备/分区阶段，不能据此宣称正式首轮已经完成。
后台续跑服务为 active/running，D2 三轮与首末轮验证完成后自动启动 D3；
D3 此时已安排在编排中，但尚未提交到调度队列。

按已有阶段真实等待粗估 D2 约 4-6 小时、D3 约 16-24 小时，
整批约 20-30 小时；排队、共享存储或失败会使估计顺延或失效。

#### 2026-10-10 09:55 控制机重启后恢复串联编排

本次核查发现本地控制机已于 **09:14:15 CST** 重启：journal 的前一
boot 最后记录为 01:05，原 transient user service 不再存在。
因此本地状态仍留在 D2 第一轮 fixmate，后续 sort 没有自动提交；
不能把失去 service 的默认 `Result=success` 当作实验完成。
A4 驱动在 A4 本机持续运行，没有受到这次本地重启影响。

远端 fixmate **8481582** 已完成，日志包含完整六 rank 统计、唯一全局
汇总、MPI_Finalize 和 TOT TIME，没有 ERROR；总记录数 365,079,492。
输出 BAM 为 **53,772,662,489 字节**，BGZF EOF 正确，远端二进制
SHA256 与本批冻结版本一致。调度器已清理历史记录，不能再从 bjobs
取得该作业退出状态。最终流程输出仍按原计划做完整复读，不把日志
计数与 EOF 检查代替全量正确性验证。

| D2 第一轮已完成阶段 | 作业 | 主核心/s | 4.3/s | 输入准备/s | 最终 MPI-IO/s |
|---|---:|---:|---:|---:|---:|
| collate | 8481404 | 136.690927 | 133.909455 | 488.719045 | 205.674971 |
| fixmate | 8481582 | 83.223580 | 80.298961 | 679.375747 | 190.511290 |

从远端完整日志恢复 fixmate 样本，并保存原本地截断日志和旧元数据到
`archive_reboot_20261010/`，恢复依据为 `recovery_fixmate_8481582.json`。
该条 `submit_wall_s` 因控制机中断无法确知，明确留空；程序自身
TOT TIME=965.494009 s 不冒充“提交至结束”时间。
汇总脚本因此只对缺失的整轮 wall 诊断输出 null，四阶段 core 合计
与其中位数仍按原样计算。新增缺失 wall 自检通过，未改计时实现、
算法、kernel、数据流、资源参数或已完成性能样本。

新的后台服务
`swbam-d2-d3-workflows-n6-20261010-003937-resume2.service`
于本地 **09:55:25 CST** 启动，PID **34741**，核查 active/running。
已有小回归和正式 collate/fixmate 均直接跳过，没有重复运行性能命令。
首轮 sort **8486308** 于 **09:55:34 CST** 提交，核查为 RUN，节点
`33240-33245`，六节点六活跃核组、`-m 64G`、rank-memory/memory/mpiio
及 level 1 保持不变。当前仍在输入准备阶段，D2 正式三轮串联尚未
完成，D3 本轮完整串联还未开始，将在 D2 完成验证后自动接续。
先前独立应用与 D3 两项转换的已完成结果不受此次编排恢复影响。

#### 2026-10-10 12:32 第二次本地重启，迁移编排至 SWLS

本地控制机于 **12:19:49 CST** 再次重启，前一 boot 最后 journal
记录为 11:09:12，原 resume2 transient service 丢失。D2 第二轮 sort
**8486922** 实际已完成：完整日志含六 rank 汇总、memory 模式、
365,079,492 local/received records、TOT TIME 和作业 finished 标记。
输出为 40,691,022,479 字节；核心 **112.645974 s**、4.3 **110.310764 s**，
物理准备 647.606735 s、MPI-IO 158.002873 s。没有程序错误。

从完整日志恢复该条性能样本，不重复运行 sort；丢失的提交 wall 留空，
不以 TOT TIME=925.317782 s 代替。原元数据/CSV 快照保存在
`archive_reboot2_20261010/`，恢复依据为 `recovery_sort_8486922.json`。
D2 首轮四阶段核心合计为 **498.858882 s**，完整最终输出复读
**8486775** 已通过。三轮中位数仍需全部运行完成后才能计算。

为避免继续依赖本地控制机，将同一批次的编排和结果登记迁到 SWLS
登录节点后台；计算作业仍在原调度队列，不在登录节点做 BAM 计算。
远端控制目录为：

`/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data/swbam_bench/workflow_control/apps_D2_D3_workflows_n6_20261010_003937/`。

后台 PID/PGID **29858**，于 **12:32:08 CST** 通过 nohup/setsid 启动，
PPID=1；以后本地断线或重启不再丢失这轮编排。SWLS 登录节点自身
退出、重启或作业失败仍可能中断，不能承诺无条件完成。
远端 `state`、`submit.remote.log` 和子目录 CSV 为最新依据，本地文件
只是快照；位置另记于本地 `controller_location.json`，不要同时重新
启动旧本地 systemd 服务。

`driver_remote.sh` 仅将 SSH 分发/上传改为在 SWLS 本机执行/复制，
`orchestrator_remote.sh` 调用该驱动。共享文件系统不支持 flock，
首次迁移启动 PID 29688 在锁检查处退出，未提交新计算作业；已改为
原子 mkdir 互斥锁，启动/互斥自检、完整冻结文件 SHA 校验及汇总
自检通过，初次失败日志保留在 `submit.remote.lock_probe.log`。
算法、二进制、六节点六活跃核组、memory/external 模式、`-m`、
rank-memory/memory/mpiio、压缩级别及核心计时口径全部不变。
控制侧提交 wall 的记录主机变化需注明，不把它当作新的核心性能口径。

核查 **12:32:30 CST**：D2 第二轮 markdup **8487264** 为 RUN，节点
`23654-23655,43590-43593`，已进入 rank-memory 初始化；已有阶段全数
跳过，没有覆盖旧阶段原始日志或重复性能命令。后续为第三轮四步、
末轮复读和 D3 三轮串联。D2 预计还需 **2-3 小时**，随后 D3 约
**16-24 小时**，整批约 **18-27 小时**，即正常情况下 10 月 11 日
上午至下午完成；长排队、真实存储波动或失败会顺延或中断。

#### 2026-10-10 15:49 D2 第三轮 sort 异常，更换节点续测

D2 前两轮四阶段均已完成；第二轮四阶段核心合计为 **498.797216 s**。
第三轮 collate **8487335** 为 **136.490567 s**，fixmate **8487371** 为
**83.303610 s**，记录数均为 365,079,492。尚未取得第三轮 sort、markdup
和末轮完整复读，因此还不能公布本批 D2 三轮流水线中位数。

第三轮 sort **8487444** 从 **13:21:03 CST** 运行至 **15:46:57 CST**，
约 146 分钟，明显超过前两轮的 15-19 分钟。取消前远端应用日志为零字节，
sort 输出仍是第二轮文件，没有本轮 preload/core/finished 信息；不能将
它当作有效性能样本，也不能把这段时间直接归因于真实磁盘加载。
六个任务显示 RUN，节点网络和 GFS 检查为 OK；`bdebug` 请求被接受，
但没有取得可解析的运行栈。计算节点 SSH 不允许当前账号登录，因此没有
直接确认 MPI_Init、文件读取或具体算法中的阻塞位置，也未确认硬件损坏。

小数据节点检查 **8489079** 在 `32890-32895` 正常完成，使用同一冻结
二进制和六核组大共享资源：rank-memory 的六 rank 加载均成功、记录数
**2,668,351**，程序 TOT TIME **11.244217 s**。此检查只验证可用性，不计入
正式性能表；`stats --basic` 仍仅作正确性检查，不作为论文性能对照项目。
该节点组随后被其他任务 **8489088** 占用，未取消或干预其他任务。
两次指定节点提交返回 `ret=-34`，没有获得作业 ID，也没有性能样本。
尝试复查原节点时，两节点已被 **8489137** 占用，未能进行原节点小数据对照。

仅取消异常作业 **8487444**，其旧编排随后正常失败退出并释放互斥锁。
取消前后日志、结果快照、任务/节点状态和检查日志归档于同一 benchmark
root 下的 `recovery_sort_20261010_154800/`；已有成功样本全部保留。
实际 sort 日志含旧作业取消标记和新作业输出，旧标记不是新作业错误；
原始异常日志另外保存，后续仍要求唯一全局汇总、六 rank 和完整结束标记。

于 **15:49:46 CST** 断点提交第三轮 sort **8489135**，使用当时空闲且
网络/GFS 正常的 **43628-43633**，已 RUN；**15:51 CST** 核查已打印
初始化、rank-memory header 和 `6 MPE + 384 CPEs` 信息。相比原作业恢复
了程序级进展，是环境相关线索，不是已确定原故障根因；仍须等待完整运行
和正确性复读，不能据 RUN 或初始日志宣布通过。

**16:06 CST** 再核查，六 rank 的 preload 均为 `ok=1`，总数据块
2,017,608；准备时间最大 **845.416474 s**，其中 root 物理扫描读取
**751.762588 s**、分片加载及等待最大 **93.653886 s**。已进入
`mode=memory` 内排，未报告错误。该次准备时间接近首轮的 864.844517 s，
不作为核心 sort 成绩；此时仍须等待核心结果、输出和完整结束标记。

**16:09:59 CST**，重试 sort **8489135** 完整结束并登记为有效第三轮：

| 指标 | 第三轮 sort |
|---|---:|
| 4.1-4.6 核心时间 | 112.719856 s |
| 4.3 | 110.376840 s |
| 排除的物理输入准备 | 845.416474 s |
| 排除的最终 MPI-IO | 233.064603 s |
| 提交至完成 wall | 1,213 s |
| local/received records | 365,079,492 / 365,079,492 |
| 输出大小 | 40,691,022,479 bytes |

日志含六 rank 统计、唯一全局汇总、MPI_Init/Finalize 和 TOT TIME；
第三轮文件的修改时间为 16:09:43，大小与前两轮一致。三次 sort 核心为
112.664029 / 112.645974 / 112.719856 s，中位数 **112.664029 s**；
第三次相对该中位数仅慢约 **0.05%**。这证明同一冻结版本换节点后能
正常完成这一输入，不能据此唯一证明原节点硬件故障或排除偶发程序问题。
第三轮输出的完整复读仍在后续最终验证中，不以汇总和大小相同代替它。

第三轮 markdup **8489315** 已自动接续，核查 **16:10 CST** 为 RUN，
仍使用 `43628-43633`。待其结束和末轮复读后，编排自动启动 D3 三轮；
此时 D2 剩余预计 **25-40 分钟**，D3 仍预计另需 **16-24 小时**。

远端后台 PID/PGID **12763**、PPID=1，继续第三轮 sort、markdup、最终
复读后自动进入 D3 三轮普通串联。当前日志改为
`submit.recovery_20261010_154941.log`，位置和 PID 已更新到
`controller_location.json`；不要同时重启旧本地编排。
后续固定上述六节点，health 1；各轮节点变化以 `commands.sh` 和调度日志
为准，不声称整批使用完全相同节点。没有使用 health 2 或复位节点。

二进制 SHA-256 仍为 `71acde4f9eb128b50fb801b022873c4bac88d882cf1491faa0363592b96e9941`。
本次没有修改算法、库、SDK、kernel 或计时实现。D2 仍为每 rank `-m 64G`
内排，D3 仍为 `-m 8G` 外排；rank-memory/memory/MPI-IO、level 1、
六节点六活跃核组和物理 I/O 排除口径不变。预计恢复正常后 D2 剩余
**约 40-70 分钟**，随后 D3 **约 16-24 小时**；这是包含真实准备、输出
和验证的等待估计，不是报告的核心时间，排队或再次异常会使估计失效。

#### 2026-10-10 21:50 核查：D2 完整串联结束，D3 首轮外排 sort 运行中

D2 第三轮 markdup **8489315** 于 **16:25:02 CST** 完成，核心为
**166.227385 s**。末轮完整复读 **8489456** 于 **16:33:52 CST** 完成，
记录数为 **365,079,492**，排序状态为 1，重复 flag 数为 **8,894,488**，
与首轮及预期基线一致。三轮 marked 均为 8,894,488，removed 均为 0；
本批 D2 三轮普通串联及首/末轮完整输出验证已全部通过。

| D2 轮次 | collate/s | fixmate/s | sort/s | markdup/s | 四阶段核心合计/s |
|---|---:|---:|---:|---:|---:|
| 1 | 136.690927 | 83.223580 | 112.664029 | 166.280346 | 498.858882 |
| 2 | 136.435232 | 83.185024 | 112.645974 | 166.530986 | 498.797216 |
| 3 | 136.490567 | 83.303610 | 112.719856 | 166.227385 | 498.741418 |

逐轮核心合计的中位数为 **498.797216 s**，约 **8.31 分钟**；不是独立
应用中位数之和，也不是实测 SSD 端到端 wall。前两轮因本地重启缺失部分
提交 wall，其整轮 wall 诊断仍留 null；第三轮四阶段提交 wall 合计为
4,086 s，不包括末端复读，不将它当作主性能指标。

编排已自动进入 D3，后台 PID 12763 仍正常存在，`state=RUNNING D3`。
目前第一轮完成如下，使用同一冻结版本、六节点六活跃核组、`-n 66`、
collate/sort `-m 8G` 外排和原内存模拟口径：

| D3 项目 | 作业号 | 核心/s | 排除的输入准备/s | 排除的最终写盘/s | 提交至结束/s |
|---|---:|---:|---:|---:|---:|
| collate 1，external | 8489489 | 866.830909 | 1233.804073 | 611.437932 | 5894 |
| fixmate 1 | 8490185 | 265.719968 | 1912.100457 | 632.552083 | 2829 |
| collate 输出完整复读 | 8491178 | 33.335827 | 2152.706200 | 0 | 2217 |
| fixmate 输出完整复读 | 8491448 | 33.016609 | 2791.862375 | 0 | 2846 |

两项性能汇总及两项完整复读的记录数均为 **1,210,633,370**。
collate 输出为 160,004,190,693 bytes，fixmate 输出为
163,882,611,815 bytes；复读只用于正确性，不作为 stats 性能对照。

首轮 sort **8491845** 于 **20:23:42 CST** 在 `43628-43633` 启动。
六 rank 已成功完成分片预加载，准备最大 **2003.293404 s**，其中物理
扫描读取 **1727.149069 s**、分片加载及等待最大 **276.144335 s**；
20:57 后进入 external 算法阶段。截至 21:50 未报告错误，也未产生
最终核心成绩。其总运行时长约 1.5 小时，仍在此前完整 D3 外排 sort
约两小时的范围内；关闭阶段诊断后，没有连续进度日志，不能仅以 RUN
状态或节点健康状态证明算法持续前进，也不能据此判定停滞。

本次只核查和同步结果，没有取消、新提交或修改实验程序。D3 尚需首轮
sort/markdup 及相应验证、第二/三轮四阶段和最终复读；按已完成阶段和
历史外排耗时粗估仍需 **12-18 小时**，即正常情况下 10 月 11 日上午至
下午。实际临时文件 I/O、真实输入/输出等待被排除在核心成绩外，但仍
影响完成时间；排队、存储波动或运行失败会改变该估计。

#### 2026-10-11 晚间核查：D2/D3 普通串联全部完成

22:21 左右查询 SWLS，远端顶层和 D3 子目录 `state` 均为 `COMPLETE`，
控制 PID 12763 已退出。恢复日志于 SWLS 时钟 **10:35:50 CST** 打印
`D2/D3 ordinary workflows complete`。D2 已完成结果不变；D3 也已完成
三轮普通四步串联以及首轮各阶段、末轮最终输出完整复读。

| D3 轮次 | collate/s | fixmate/s | sort/s | markdup/s | 四阶段核心合计/s |
|---|---:|---:|---:|---:|---:|
| 1 | 866.830909 | 265.719968 | 417.696993 | 623.157072 | 2173.404942 |
| 2 | 859.625863 | 265.726986 | 417.488205 | 623.231092 | 2166.072146 |
| 3 | 759.552180 | 265.816599 | 418.047763 | 623.290652 | 2066.707194 |

逐轮合计再取中位数为 **2166.072146 s**，约 **36.10 分钟**。
D2 相同合计口径中位数仍为 **498.797216 s**。上述均属于冻结旧版
`4.1~4.6` 内存模拟核心时间，不是新版 `--storage ram` 的实测 MPI 区间。
不能通过加固定 7 s/命令或更改列名将历史成绩转换成新版端到端成绩。

D3 每轮 records/kept 均为 1,210,633,370，marked 均为 42,383,102，
removed 均为 0。首轮最终复读 8493610 和末轮 8498591 均完成，记录数、
`is sorted: 1`、duplicate flag 数均与算法汇总一致。末轮日志含正常
MPI_Finalize、TOT TIME 及 `Job 8498591 has been finished.`。
各轮最终文件为 105,926,943,620 bytes。此验证不等于全量字段一致性证明；
新版 RAM 回归中发现并修复的外排边界保护问题未进入这批冻结旧二进制。

第二轮 collate/fixmate/sort/markdup 作业分别为
8493754 / 8494301 / 8494596 / 8495943；第三轮分别为
8496165 / 8496751 / 8497733 / 8498414。全部已有完整性能登记，
没有待运行的本批 benchmark。D3 两项转换的三轮结果已于 10 月 10 日
00:24 完成，见前文；此次核查也确认末轮 SAM->BAM 日志正常结束。

顶层原恢复日志、D2/D3 `results.tsv`、workflow 汇总、命令和阶段日志
保留在远端原目录。本地新快照为原控制归档目录下
`snapshot_20261011_222300/{D2,D3}/`，不覆盖此前的运行状态快照。
本次没有取消、新提交或修改任何算法/二进制；其他可见作业未干预。

后续建议冻结新版 RAM 二进制，先做 D2 六 rank 容量与输出探测，再重测
D1/D2 两项 SDK、六核组独立应用和普通四步串联。RAM 容量与算法 `-m`
工作区另计但共同占用节点内存，不能直接将旧 `-m 64G` 配上任意 RAM 池。
D3 旧成绩继续作为历史/探索记录，不自动扩展本轮论文主实验范围。
