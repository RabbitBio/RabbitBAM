# RAM D1 正式性能测试

日期：2026-10-12。本轮补新版 RAM 条件下 MPI 端到端成绩，历史核心时间不改名。
驱动：`tools/ram_d1_matrix.sh`；冻结应用和两个 SDK 二进制与 D2 探测相同。
不修改算法、四条流水线或 CPE kernel。

## 启动状态

SWLS 时钟 01:07:31 启动，原后台 PID 22798。
第一项 `sdk_count_n1_r1` 作业 **8506089** 已成功，节点 **22456**：
177,798,571 条记录，核心 14.765757 s，MPI 区间 max **21.903565 s**，
MPI 外准备 195.962940 s，12G 池的文件页 peak=9,927,917,568 bytes。
这是单核组首轮样本，不是六核组成绩或正式中位数。

其后驱动遇到 SWLS Bash 4.2 在 `set -u` 下展开空数组的兼容问题。
已改为兼容的空数组展开，并在服务器 Bash 4.2 实测验证；
01:12:43 以后台 PID **24264** 续接，复用已验证 count，不重复提交。
`sdk_bam2bam_n1_r1` 作业 **8506152** 已成功，节点 **22456**，
MPI 区间 max **56.433480 s**、核心 49.286531 s，默认保存完整输出。
这里的单节点代表一个 rank、一个 64 CPE 活跃计算核组，其他核组仅提供内存。
10:58 检查：仅上述两条正式样本完成，随后 SDK 输出完整复读
`validation_sdk_n6_r1` 作业 **8506221** 仍显示 `RUN`，已经超过九小时。
六个 rank 的预加载/ready 标记均于 01:22 前后完成；日志随后未更新，
没有复读结果或成功退出总结。原因尚未定位，不能仅凭这些日志归因于服务器。
驱动 PID **24264** 因此未继续提交。D2 容量关卡已全部通过，但按串行规则等待
D1，见 [D2 关卡与矩阵](RAM_D2_容量关卡与正式矩阵.md)。

按用户要求，现已取消旧复读 `8506221`，调度状态确认为 `EXIT`。
旧日志、ready 标记清单、任务/节点信息、成绩、作业映射和原驱动均保存在
`retries/validation_sdk_n6_r1.stopped_job8506221/`，没有将失败复读记为成功。
11:21 提交重试 **8513953**，固定替代节点 **22447-22452**；
**11:24:44，驱动确认完整复读通过**：177,798,571 records、sorted=1，
六 rank 均 `status=0`，MPI 区间 max=10.552032 s、核心=3.267513 s、
MPI 外准备 max=162.184295 s。相同冻结程序和同一输出完整复读成功，
不支持“这份输出必然不可解析”的判断；原停滞的具体原因仍未定位。
此前成功的 D2 节点当时已被其他任务占用，因此改选当前 idle、计算网络和 GFS
状态正常的上述六节点；没有沿用旧复读的 `24528-24531,22456,24532`。

D1 驱动以 PID **15364** 续接，仍使用原冻结脚本和程序，复用原先两条成功成绩。
本次只在验证命令通过 `-x SWBAM_STATS_DIAG=1` 开启已有阶段诊断，不改变性能命令。
新日志/会话为 `validation_sdk_n6_r1.retry1`；原日志名链接到新日志，活动作业映射
指向新 ID，旧映射完整留档。复读成绩仍归类为 validation，不纳入正式中位数。

独立 watchdog 只监控本次新复读：开始 RUN 后最多 30 分钟，提交至结束最多一小时，
超时仅取消 `8513953`。六项 mock 检查在本地及 SWLS Bash 4.2/Python 3.6 通过，
覆盖 DONE/EXIT/KILLING 不重复取消、运行/排队超时只取消指定 ID、无效 ID 拒绝。
watchdog 已在本次作业退出时结束，没有超时取消新作业。
维护窗口内短暂挂起 D2 等待驱动，现已恢复；D2 仍等 D1 全矩阵完成再正式开测。
11:24:44 自动提交下一项两 rank SDK count **8513986**；11:26 检查为 `RUN`，
节点 `22448-22449`，D1 状态已推进到 `RUNNING sdk_count_n2_r1`。
11:29 检查，该 count 已通过，MPI 区间 max=14.634874 s、核心=7.393490 s、
准备 max=156.180891 s；下一项两 rank SDK BAM 回环 **8514059** 已 `RUN`，
同为 `22448-22449`。原先两条成绩保留，当前已有三条正式性能样本和一条成功验证。

远端唯一目录：

```text
/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data/swbam_bench/ram_D1_matrix_20261012_010800
```

驱动在 SWLS 登录节点后台串行提交计算作业，断开本地会话不停止测试。
任一命令错误、记录数不符、内排项目 spill 或正确性复读失败即停止后续；
不自动调整资源、落盘或重复提交已完成的样本。
`BENCH_RESUME=1` 仅用于核对已有状态后的显式续接，会校验二进制 SHA-256。

## 矩阵与容量

共 **66 条性能命令 + 8 条验证命令**：

- 两个 SDK，各 1/2/4/6/12 ranks，每点三次，共 30 条。轮次资源顺序升序、降序、升序。
- flagstat、过滤、collate、fixmate、sort、markdup、BAM->SAM、SAM->BAM，
  固定六节点六活跃核组，每项三次，共 24 条。
- 三轮普通 `collate -> fixmate -> sort -> markdup --stream`，每轮四条独立 MPI 命令，
  共 12 条。每轮重新从原始 BAM 开始，并保存全部阶段输出。
- SDK 首次完整输出、SAM->BAM 首次完整输出各复读一次；
  每轮 workflow 的 sort 和 markdup 输出分别完整复读，共 8 条。
  stats 仅为验证，不计入性能图或四步合计。

| SDK ranks | count 池/rank | BAM 回环池/rank |
|---:|---:|---:|
| 1 | 12G | 24G |
| 2 | 6G | 16G |
| 4 | 4G | 8G |
| 6 | 2G | 8G |
| 12 | 2G | 4G |

| 六 rank app | RAM 池/rank | 固定参数 |
|---|---:|---|
| flagstat | 2G | 完整 flagstat |
| filter | 8G | MAPQ>=30，level 1 |
| collate | 8G | `-n 66 -m 24G`，必须 memory |
| fixmate | 12G | `-m` 添加 ms 标签 |
| sort | 8G | `-m 32G`，必须 memory |
| markdup | 8G | `--stream -l 300 -m 8G` |
| 两项转换 | 24G | 完整 SAM/BAM，level 1 |

RAM 池全部预触碰；日志的文件页 peak 不是整体内存占用。算法工作区另计。
所有 rank 仍使用 `-mpecg 6 -cgsp 64 -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc`。
自动选择健康节点，不固定历史占用节点；每条命令核对实际 MPI rank 数。

## 输入与验证

原始 BAM：`D1-WES-9G/initial/D1-WES-9G.coordinate.bam`，
9,891,133,304 bytes，177,798,571 条记录，814,215 个 BGZF 数据块。

独立 app 复用已验证的固定中间输入目录：

```text
swbam_bench/D1_apps_current/apps_D1_n6_20261001_031800/out/
```

其中原始 SAM 为 68,705,660,624 bytes；name、fixmate 和 fixmate.sorted BAM
分别为 18,755,398,492 / 19,289,393,734 / 12,041,427,629 bytes。
这些固定文件只用于独立命令。三轮串联不能借用它们，必须用本轮前一阶段输出。
输入绝对路径、大小和 mtime 保存在 `inputs.txt`；实际命令保存在 `commands.sh`。

每个成功作业必须有完整且不重复的 rank ID、全部 `status=0`、RAM 存储总结及
全局记录数 177,798,571。过滤应保留 170,359,591；markdup 应 marked=50,084,605。
完整复读要求可解析、记录数正确、坐标 sorted=1，markdup duplicate bit 计数与汇总一致。
这些检查加上已通过的小数据独立字段回归，不等同于 D1 全量跨工具字段/标记集合比较。

除 SDK 首次回环和 SAM->BAM 首次输出外，独立写性能命令使用 `--no-dump`。
它仍保存完整 RAM 输出，绝非环形覆盖或只统计字节数。每轮 workflow 四个阶段
都默认在 Finalize 后保存，下一条命令重新预加载，中间准备和实际保存不计成绩。

## 计时与等待

`results.tsv` 主列 `e2e_max_s` 是各 rank `MPI_Init` 前至 `MPI_Finalize` 后真实耗时最大值。
`core_s`、`prepare_max_s`、`pool_peak_max_bytes` 和 `output_bytes` 仅用于诊断。
`workflows.tsv` 每轮累加四个 MPI 区间，三轮合计中位数保存为 `workflow_summary.txt`。
验证命令不会加入合计，历史独立阶段中位数不会导入本轮。

按已有 D1 每条物理准备约 2-3 分钟，以及 SAM 扫描、串联输出保存和验证开销，
最初预留 **4-6 小时**，排队和共享存储波动另计。
原 SDK 输出复读异常停滞，令先前 05:10-07:10 的结束估算失效。
该停滞已通过取消旧任务、换节点重试恢复；从 11:25 起，若持续正常，
剩余 D1 暂预留约 4-6 小时，即约 15:30-17:30，调度和共享存储波动另计。
不是保证；旧卡住作业不补入成绩，新复读也只用于 validation。

冻结 SHA-256：

- app：`5fcc7b5e6758b8c503d75019c2a14e03db49db96bdf8c031687200f1bd9a8399`。
- count：`925af1c5e20d79cb57bddf56ab0a5c570056b27feaab467291d54420d1b80b0a`。
- writer：`8c6ee35ca1cb7b8ce34881384de5bfbb24fcf89303b730bd1759b0968096bd8f`。
- 当前 driver：`71e515d1e8d13cef0b5f106be33f3cc536053feb3ae49462c7915f27322a05ec`。
- 原 driver：`71dfb4b8962f66a05aa688ae02a2410eab6d34ffb680ca460da976c0df23c0f9`，
  留存于 `bin/ram_d1_matrix.submitted_v1.sh`。程序二进制及首项计时不变。

驱动的七个 mock 调度测试覆盖矩阵完整性、轮次顺序、rank 去重、错误退出、
容量错误、禁止混入 external sort、无重复提交续接和逐轮合计；不是真实性能测试。
`bash -n` 与 `git diff --check` 通过。A4 与 D2 状态见
[正式范围与 D2 探测](RAM_D1_D2_正式范围与D2探测.md)。
