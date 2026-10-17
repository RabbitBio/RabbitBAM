# markdup 回退索引与重复更新优化

## 2026-10-17 13:52 最终更新

本页索引版验证链于 12:43:23 COMPLETE，7/7 样本完成。D3 292.167561 s 的完整输出复读通过；同期 D1 old/new 为 48.693476 / 49.364532 s（新版 +1.38%），D2 new/old 为 170.370591 / 166.283629 s（新版 +2.46%）。D1/D2 均 halo enabled、linear lookup index_peak=0，记录数与生成阶段 duplicate 总数正确；本轮是单次对照，不替换三次正式成绩，也不宣称零回退。

独立增量版验证链也于 12:54:42 完成，D3 单次 241.256594 s，完整输出复读通过；详见 [增量边界与安全释放最终探测](markdup_增量边界与安全释放验证.md)。两版日志、冻结指纹和输出目录分开保存，没有混用成绩。13:50 bjobs 无本轮活动作业。下文旧状态保留为历史。

## 2026-10-17 12:13 恢复更新

原复读 `8617619` 在六 rank 预加载后没有后续处理日志，07:00 达到六小时驱动上限被取消；未记作通过。新会话重试 `8624788` 已于 12:09:17 完整验收：165,373,770 条、sorted=1、duplicate bit=33,171,753，processing=3.209754 s。原 292.167561 s 生成样本复用，不重跑；D1/D2 新旧对照继续执行，目前不能提前宣布无回退。

恢复脚本为 `tools/resume_markdup_lookup.py`，原冻结驱动和程序保持不变；原失败状态与日志保留，重试有独立作业号及会话。另已独立实现、构建并启动后续增量边界验证，详见 [增量边界与安全释放](markdup_增量边界与安全释放验证.md)，不混入本页索引版成绩。

更新时间：2026-10-17 01:02 CST。状态：实现、主机检查及六核组小数据逐记录比较通过；D3 单次 processing 已降至 292.167561 秒，完整输出已保存，复读 `8617619` 进行中，随后执行 D1/D2 新旧对照。仅做针对性复验，不启动正式矩阵。

## D3 首次实测

作业 `8617469` 于 01:00:40 经驱动完整验收，六 rank 均 status=0，records=165,373,770、marked=33,171,753，完整输出 16,614,394,500 bytes 已在 MPI_Finalize 后保存。仍触发原实际数据的 halo fallback，没有跳过回退或减少输入。

| 指标 | 前一版边界优化 | 本轮索引/更新优化 |
| --- | ---: | ---: |
| processing max/s | 388.807467 | 292.167561 |
| MPI 区间 max/s | 395.921945 | 299.301961 |
| rank 2 group/s | 185.994514 | 155.769095 |
| rank 3 group/s | 134.659185 | 69.624752 |
| 算法跟踪峰值/bytes | 1017512304 | 1017326352 |
| RAM 临时文件累计读、写各 bytes | 19003959792 | 19003959792 |

processing 减少 96.639906 秒，即 24.86%，约 4.87 分钟。为独立单次探测，不是三轮正式中位数，也不是仅根据分项计时推算收益。输出大小和 RAM 临时 I/O 未减少。六 rank 均确认 `lookup=indexed`，索引峰值分别为 40/18/239/322/54/78 个节点，最大约 15 KiB。

完整复读 `8617619` 于 01:00:41 提交，记录数、排序、重复标记数仍待复读验收。D1/D2 新旧对照尚未完成，不能宣布其新版性能无回退。D3 当前仍慢于 a4 markdup 的 132.50 秒三次 wall 中位数；剩余开销仍主要在密集边界的反复判重处理，是否进一步增量维护应在本链结束后判断。

## 范围

上一轮边界存在标记合并和按需更新，将 D3 单次 processing 从 3543.573882 秒降到 388.807467 秒。其完整输出复读 `8616499` 已通过：165,373,770 条、坐标有序、duplicate bit=33,171,753。不是全量逐记录跨工具比较。

本轮继续解决两项冗余工作，不修改判重规则、从核、codec、批次大小、压缩/暂存调度或 SDK 公共接口：

1. 回退路径每次记录决策原来线性扫描 pending 批次链表，再线性扫描批内至多 64 个块。新增 app 私有非拥有式区间索引 `markdup-pending-index.h`，批次查找 O(log B)，最近批次命中 O(1)，批内按起始 ordinal 二分定位。保留对已封闭块的迟到决策检查。
2. 已保留的回退 carry 候选在首次进入窗口时已经 pending，后续继续保留时不重复调用 pending 更新。新候选照常标记，最终 keep/duplicate 决策不省略。

## 生命周期与预算

- 只有 halo 预处理未成功时启用索引，正常路径继续原线性查找，不创建索引节点；旧候选免重复 pending 也仅在有回退边界时启用。
- 每个非空 pending 批次一个节点，payload 转存不影响索引；批次完全压缩暂存后，先移除索引和最近命中缓存，再释放批次。空批次不占节点，不能误删相同起始 ordinal 的非空批次。
- 节点存储估算为 `sizeof(value_type) + 4*sizeof(void*)`，64 位构建为 48 bytes/批。纳入原决策 metadata 的 `-m/4` 限制和总体工作内存跟踪，不新增用户参数，不侵占 RAM 存储池计数。
- 正常路径的额外分支和编译布局仍可能造成微小性能变化，因此不声称未经实测的绝对零回退。D1/D2 用同输入、新旧程序单次交错对照；超过约 5% 标记需复核，不重跑完整矩阵，也不替换旧成绩。

## 已完成检查

- `tools/markdup_pending_index_test.cpp`：区间命中/缺口、重叠拒绝、空批次、乱序回收、最近命中失效、地址复用和预算计数；普通及 ASan/UBSan 通过。
- `tools/markdup_window_test.cpp`：多批判重对照、标记合并、无新边界批次、optimistic 本地/远端与 EOF、正常路径无边界表分配、旧回退 winner 只标记一次 pending；普通及 ASan/UBSan 通过。
- 独立 Sunway RAM 支持版和无私有支持的磁盘版 app 构建、`git diff --check` 通过。验证驱动四项配置检查和底层九项日志验收检查通过。
- 六节点小数据 `8617454`：200,000 条记录、duplicate bit=22,259，与上一轮输出逐记录字节/顺序一致，坐标排序与 gzip CRC 检查通过。该小数据走 halo enabled，回退索引由 D3 全量实际数据另验。

本地主机检查命令：

```bash
g++ -std=c++11 -O2 tools/markdup_pending_index_test.cpp -o build_bounded_io/verification/markdup_pending_index_test
mpicxx -std=c++11 -O2 -I swbam-app/include -I swbam-sdk/include -I ext/htslib-1.20 tools/markdup_window_test.cpp -o build_bounded_io/verification/markdup_window_test_lookup
python3 -m unittest discover -s tools -p test_markdup_lookup_recheck.py -v
python3 -m unittest discover -s tools -p test_d3_probe.py -v
```

## 冻结配置与实测链

旧控制程序为 D1/D2 原正式 markdup，SHA-256：
`3dc30d0a9ac3ba2a384fec617b66f32d7801a53d389288db0f51f110c7726a51`。

新版独立目录 `build_bounded_io/app_markdup_lookup/`，SHA-256：
`b6d8469064cb7b0f9b0b357012ef73f577b4c38847e018104b9e0eb8924b16e4`。

六节点、每节点一个活跃核组，64 CPE；`--stream -l 300 -m 8G --compress-level 1`，`SWBAM_NODE_MEMORY=88000M`。池容量为小数据 2G/rank、D1 8G/rank、D2/D3 24G/rank，复读 8G/rank。主成绩为 processing 各 rank 最大值，MPI 区间另记；预加载和 MPI_Finalize 后保存不计入这两个区间。

顺序为小数据逐记录比较、D3 新版完整输出、D3 完整复读、D1 old/new、D2 new/old。D1/D2 四个单次控制样本使用 `--no-dump`，仍生成完整 RAM 输出并检查生成阶段记录数/重复标记数，不将其表述为新的全量逐记录对照。D3 重新保存、预加载、完整解析最终 BAM。

日志必须确认新版 D3 为 `lookup=indexed`；D1/D2 为 halo enabled 且六 rank 均 `lookup=linear index_peak=0`。验证驱动遇到错误即停；单个提交有 6 小时总上限，调度时限 3 小时，不无限等待。

```text
本地：/home/wzs/work/data/D3-mouse-SRR7768905/markdup_lookup_20261017_005000
远端：/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data/swbam_bench/markdup_lookup_20261017_005000
驱动：tools/markdup_lookup_recheck.py
tmux：markdup_lookup_20261017
进度：state.json、driver.log、jobs.jsonl
样本：results.tsv、logs/
汇总：summary.json（所有关卡完成后才写）
```

D3 全量作业 `8617469` 于 00:49:55 提交，01:00:40 验收通过。按历史预加载、保存和处理耗时，本链从 00:49 启动起暂估约 45--75 分钟，约 01:35--02:05，排队与共享存储波动另计。没有自动串接 D3 正式三轮矩阵。
