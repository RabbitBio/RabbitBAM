# Markdup 固定探测与部分预处理

日期：2026-11-08。承接 `markdup_POS_halo预处理验证.md`，本轮只修改 markdup 的边界预处理，不替换原冻结程序或正式性能表。

## 实现

每个 rank 固定探测首 32 + 尾 32 个 BGZF 块（不足时取并集），不再为覆盖 `b±L` 扩读。解压、解析和候选提取后，仍按源记录 POS 选择目标 halo 内已观察到的候选，再复用现有 owner 路由和 pair/single hash 判重。

| 模式 | 预处理 | 主循环 |
| --- | --- | --- |
| `halo=enabled`，覆盖完整 | 与上一版相同：仅闭合 key 的赢家最终 KEEP，已知败者 DUP，其余 Pending；闭合 marker 可以过滤 | 保留正常快速路径及现有安全检查 |
| `halo=partial`，覆盖不足 | **只确认已知 DUP**；不产生最终 KEEP；全部未决赢家和 paired-marker 保留 | 保守窗口、边界进度、候选通信与判定回传均保留 |
| `halo=fallback`，探测元数据不可用 | 不做提前判定 | 原保守处理，例如空 rank 或无法保留探测元数据 |

`enabled` 作为 plan 内部布尔值表示预处理已执行，不等于覆盖完整；新增 `partial` 显式区分。日志和汇总工具也区分三种模式。

部分预处理只过滤 source ID 已确认 DUP 的普通候选。未决赢家不会被当作 KEEP 跳过，而是在主循环按原输入顺序重新提取、进入窗口。partial 分支不向窗口预置另一套赢家状态，避免重复计数或 ordinal 映射变化；CPE key、score、QNAME tie-break、窗口和压缩热路径不变。

逻辑 `[b-L,b+L]` 不是实际数据覆盖证明：partial 分支即使某个 key 的位置支撑区间位于逻辑 halo 内，也不能据此提前 KEEP。此版选择全局保守的 DUP-only partial 分支，不新增逐 key 的实际覆盖证明。

## 对 D1 / D2 的影响

上一版六 rank 实验的 D1、D2 都在初始 64 块内完整覆盖，因此新方案在相同输入分区和参数下仍应走原快速分支，不进入 partial 分支。预处理规模不扩大，CPE、窗口和主循环候选筛选规则不变，仅在预处理阶段多一次可用性状态检查及模式字段。

这是代码路径分析，不是新版本全量 D1/D2 性能保证。本轮不自动重跑 D1/D2；覆盖完整的 20 万条样本用于简单兼容性检查。

## 验证

| 检查 | 结果 |
| --- | --- |
| Sunway APP / halo regression 构建 | 通过 |
| 主机元数据测试 | partial 不提前 KEEP、marker 存在性及原 key 角色检查通过 |
| 通信汇总工具测试 | 4 项通过，校验 partial、缺失 rank、错误字节合计及普通候选/marker 分离 |
| 两 rank CPE/MPI 六项回归，作业 `9077239` | 全部通过；覆盖完整保持旧路径，覆盖不足固定 64 块并启用 DUP-only 预处理 |
| partial 单端 / 双端实际输出 | 两 rank 测试中逐条 raw BAM 字节及顺序与不做预处理的保守窗口参考一致；每个 BGZF block 解码并验证 CRC |
| 六 rank 20 万条样本，作业 `9077240` | 200,000 条、22,259 条 DUP；全部记录字节、顺序、CRC 与已验证基准一致 |
| D3 全量，作业 `9077269` | 六 rank 成功；165,373,770 条、33,171,753 条 DUP，均与历史计数一致；完整 RAM 输出大小一致 |

20 万条样本仍为完整覆盖，预处理最终 KEEP=124、DUP=21、普通 Pending=123、未闭合 marker=128，与上一版一致。预处理候选 payload=3,093 B，主循环普通候选、marker 和候选 payload 全为零；4.3 为 0.160753 s，仅作小样本 smoke，不视为正式加速比。

## D3 实验

同输入、同六节点六 rank、每 rank 64 CPE；其他核组只供内存。继续使用现有分布式 RAM storage：输出完整写入 RAM，仅跳过最终磁盘 dump。preload/dump 不计 processing 或 4.3。

```text
markdup --stream -l 300 -m 8G --compress-level 1
--storage ram --ram-capacity 24G --ram-ranks 6 --no-dump
-N 6 -np 1 -mpecg 6 -cgsp 64
-share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
```

输入为原 D3 小鼠 `fixmate + sort` 输出，16,546,778,000 B、165,373,770 条记录；预期 DUP=33,171,753。每项单次诊断，不与旧正式中位数混为同轮对照。

| 指标 | 上一版：扩读后 fallback | 本版：固定探测，允许 partial |
| --- | ---: | ---: |
| 作业 | 9076526 | 9077269 |
| 最后一次探测块数 / rank | 512 | 64 |
| 预处理最终 KEEP | 0 | 0 |
| 预处理 DUP | 0 | 13,267 |
| 主循环普通候选发送数 | 548,126 | 541,474 |
| 主循环 paired-marker 发送数 | 547,945 | 547,945 |
| 主循环候选 payload/B | 89,254,397 | 88,650,061 |
| 预处理候选 payload/B | 0 | 1,208,091 |
| 候选 payload 总量/B（两阶段） | 89,254,397 | 89,858,152 |
| 主循环判定 payload/B（辅助指标） | 8,770,016 | 8,663,584 |
| 预处理判定 payload/B（辅助指标） | 0 | 53,216 |
| 4.3/s | 102.156094 | 103.156585 |
| processing max/s | 103.745166 | 104.747855 |

用户关心的候选记录通信：普通候选数为 `main_candidates-main_markers`，paired-marker 单列；候选 payload 含 candidate structs、marker 及必要 QNAME，不是整条 BAM record。预处理赢家/marker 在主循环可能再次发送，因此同时报告两阶段总 payload，不能仅以主循环下降宣称总通信下降。MPI count/status/progress collective 不纳入候选 payload。

### D3 结果分析

六 rank 都为 `partial=1 conservative=1 probe_blocks=64 final_keep=0`；预处理保留普通 Pending=1,793、未闭合 marker=15,044，主循环继续处理。DUP=13,267 是所有 source rank 的合计，不全是远程候选，所以主循环远程普通候选减少量为 6,652（1.21%），不是 13,267。

主循环候选 payload 减少 604,336 B（0.68%），但预处理多发 1,208,091 B，最终两阶段候选总量增加 603,755 B（0.68%）。marker 必须保留，部分赢家也需要重进主循环，因此本轮没有获得候选总通信减少，更没有实现 D3 主循环零候选通信。

本版 4.3 单次比上一版慢 0.98%，processing 慢 0.97%。没有同版本重复样本或交错旧版对照，不能认定为确定的回归，也不能直接称为噪声；目前只能结论为性能接近、未观察到提速。候选 MPI 传输及交换同步计时 `candidate_mpi` 的 rank 最大值约 0.117 s，而 `group` 最大值约 46.764 s；减少这部分候选 payload 并不等于减少主要窗口/分组开销。

两个作业的完整 RAM 输出大小都是 16,614,394,500 B；RAM 池最大峰值同为 6,211,764,224 B/rank，各 rank 的临时数据读写总量也一致。不能把本轮性能变化解释成省略输出或减少暂存。全量只完成生成阶段计数和大小核对，未对本版 D3 全量逐记录复读比较。

## 文件与复现

- 新构建：`build_bounded_io/app_markdup_partial_halo/`。
- APP SHA-256：`2f9a0f206b3c3fb3a55e46d77ecd476d96577f2a2ba33a6a553775970fd2824e`。
- 本地原始日志和源码快照：`build_bounded_io/verification/markdup_partial_halo_20261108/`。
- SWLS 目录：`/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data/swbam_bench/markdup_partial_halo_20261108/`。
- 准确集群命令由 `tools/submit_markdup_pos_halo.sh` 写入各轮 `*.command.sh`。
- 汇总：`tools/markdup_halo_metrics.py LOG --records 165373770 --duplicates 33171753`。

旧版 `markdup_pos_halo_20261108/bin/RabbitBAM-MPI-final` 及其全部结果保留；未修改论文。D3 若仅核对生成阶段总数，不能宣称本版全量逐记录一致。
