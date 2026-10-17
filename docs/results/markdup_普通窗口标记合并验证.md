# markdup 普通窗口存在标记合并

更新时间：2026-10-17 15:22 CST。状态：D3 同指纹三次 markdup 确认、完整输出复读及 D1/D2 单次新旧控制全部完成。主链 15:11:00 COMPLETE，当前 bjobs 无用户作业。D3 processing 中位数 **103.548938 s**，MPI 区间中位数 **110.647584 s**；D1/D2 新版单次分别慢约 2.15% / 2.84%，未超过预设 5% 复核线，但不能宣称正常路径完全无回退。原正式成绩、源码与程序指纹均保留。

## 三次确认结果

| 样本 | 作业 | processing max/s | MPI 区间 max/s |
| --- | --- | ---: | ---: |
| 首轮，完整保存并复读 | 8625242 | 103.554015 | 110.647584 |
| 第二次，完整 RAM 输出 | 8625269 | 103.505648 | 110.615678 |
| 第三次，完整 RAM 输出 | 8625280 | 103.548938 | 110.667069 |
| 中位数 | | **103.548938** | **110.647584** |

三次确认链于 14:17:22 COMPLETE，processing 范围为 103.505648--103.554015 s，全部健康样本均保留。参考 a4 samtools markdup 三次 wall 中位数 132.50 s，processing 对应约 1.28x、MPI 区间约 1.20x；SWBAM 使用 RAM 模拟存储而 a4 使用真实 NVMe，不能把这一比较称为同存储介质的程序 wall 比较。这是同版本 markdup 的三次确认，不等于已完成全部 D3 应用和四步串联正式矩阵。

生成阶段均为 165,373,770 条、marked=33,171,753、六 rank status=0，完整输出大小均为 16,614,394,500 bytes。首轮保存输出完整复读通过数量、排序及 duplicate bit 检查，小数据另有 200,000 条逐记录字节/顺序比较；不据此宣称全量跨工具逐字段一致。

## 正常路径对照完成

| 数据 | 旧版作业 | 旧版 processing/s | 新版作业 | 新版 processing/s | 新/旧变化 |
| --- | --- | ---: | --- | ---: | ---: |
| D1 | 8625268 | 48.746952 | 8625279 | 49.796452 | +2.15% |
| D2 | 8625384 | 166.206185 | 8625363 | 170.925445 | +2.84% |

D2 new 于 15:01:17 验收，D2 old 于 15:11:00 验收，主链随即 COMPLETE。两版均为 365,079,492 records、marked=8,894,488、完整 RAM 输出 40,718,592,510 bytes；六 rank 均成功。两版 RAM 池峰值均为 14,168,358,912 bytes/rank，累计 RAM 临时数据读、写各 13,766,224,386 bytes。生成日志对照不等于 D2 全量输出逐记录比较。

D1/D2 均 halo=enabled；新版六 rank 均 lookup=linear、index_peak=0，没有增量边界或普通窗口合并计数。D2 rank 1 仍有 33 批原始未决暂存、读写各 127,409,080 bytes，这是正常路径的有界延后缓存，不是边界 fallback，也不能将其与全部 RAM 临时 I/O 混同。

本次为单次配对控制，不能据此把约 2--3% 差异认定为随机噪声或断言无性能影响；目前没有超过预设阈值的性能异常。没有重跑或覆盖 D1/D2 正式矩阵。建议保留冻结版本与这组兼容性对照，下一步再启动 D3 正式应用/串联矩阵；本轮未自动提交该矩阵。以下为恢复过程与当时的等待估计，均已结束。

14:11:46 原驱动将 D1 new 作业 8625279 的 STARTING 当作失败退出，应用实际继续正常结束。已完成 D1 old 为 48.746952 s，D1 new 原日志为 49.796452 s，单次差约 +2.15%，低于预设 5% 复核线，但不能据此宣称零回退或修改历史正式成绩。两项记录数与 duplicate 总数一致，均使用正常 halo 路径。

14:46:56 使用独立恢复入口 `tools/resume_markdup_markers.py` 续跑：保留原脚本、源码、二进制指纹与失败状态；只对 bjobs 的 STARTING 最多等待五分钟，原始调度状态另存 `scheduler-startup.jsonl`，不改写成成功状态。EXIT 仍由原验收拒绝，DONE/消失作业仍须原完整日志检查；已有作业复用原 job ID，不重复提交。新增六项驱动测试及既有检查共 36 项通过，`git diff --check` 通过。

恢复 tmux 为 `markdup_markers_recovery_20261017`，日志为主链根目录同级 `markdup_markers_20261017_135600.startup-recovery.log`。完成日志重新验收后只剩 D2 new/old 单次控制，预计从恢复起约 20--30 分钟，即约 15:07--15:17，排队及共享存储波动另计；此前 14:35--15:00 的估计已被驱动中断作废。未启动正式矩阵。

14:47:00 已验收并复用 D1 new 原作业 8625279，D2 new `8625363` 随即提交并 RUN，六节点 `26628-26633`；完成后自动接 D2 old，同程序参数与容量不变。没有重复提交 D1 或 D3。

## D3 首个全量样本

| 指标 | 上一版增量边界 | 本版 marker 合并 |
| --- | ---: | ---: |
| processing max/s | 241.256594 | 103.554015 |
| MPI 区间 max/s | 248.389116 | 110.647584 |
| rank 2 group/s | 137.905347 | 47.060118 |
| rank 5 group/s | 78.665821 | 36.972471 |
| 工作集跟踪峰值/bytes | 826929189 | 547254309 |
| RAM 存储池峰值/bytes | 6211764224 | 6211764224 |
| RAM 临时数据累计读、写各 bytes | 16457456316 | 16457456316 |

完整保存 16,614,394,500 bytes，生成阶段记录数 165,373,770、marked=33,171,753、六 rank status=0。该单次 processing 比上一版缩短 57.08%，低于 a4 samtools 三次 wall 中位数 132.50 s，约为其 1.28x 速度；包含 MPI 初始化退出的区间也低于该参考值（约 1.20x）。仍需注明 SWBAM RAM 模拟存储和 a4 真实 NVMe wall 的口径差异；没有把程序全部 wall、物理预加载和保存纳入 SWBAM 的这两个区间。

六 rank 均实际触发 fallback，并打印 marker 合并计数。最重的 rank 2 单批保留 marker 峰值为 707，rank 5 为 703。累计 visits/retained 为批间 metadata 次数，不是输出记录数。原始未决数据转存仍为 594 批、读写各 2,389,566,787 bytes，完整输出和临时 I/O 都未减少；主要收益来自减少窗口重复计算及其同步等待，不能归因于少写输出或删除数据。

原完整复读、D1/D2 正式结果及本版源码/程序冻结均独立保留。14:04:05 本版生成样本验收通过，14:04:06 提交完整复读 `8625261`。下面启动前状态保留为过程记录。

14:08:00 复读验收通过：165,373,770 records、sorted=1、duplicate bit=33,171,753；processing=3.202919 s。该完整复读检查数量、排序及重复总数，不表述为全量跨工具逐记录比较。D1 old 控制作业 `8625268` 于 14:08:01 提交并 RUN，随后自动 D1 new、D2 new/old，不提前宣布本版正常路径无性能回退。

## 定位依据

上一版增量边界的 D3 单次 processing 为 241.256594 s，rank 3 group 已降至 19.649916 s，但 rank 2 仍 137.905347 s、rank 5 仍 78.665821 s。检查发现普通窗口仍把所有未过期 paired-marker 放入下一批，旧去重条件只覆盖边界更新。回退的边界候选已被增量表剥离后，内部密集坐标中的重复存在标记没有合并，会在后续每批反复扫描、写入 carry。

paired-marker 是“有配对记录存在”的 key，不是输出记录，也不是待决 winner。相同完整 single key 的存在标记仅需一份；不同参考、位置、方向和其他 key 字段仍严格区分。不能合并真正的 BAM 记录或任意 winner 决策。

主机重现使用 128 批共 131,072 个存在标记、8 个不同 key，并混入待决 singles。改动前 carry 保存 131,072 个 marker，改动后为 8 个，参考逐记录最终 keep/duplicate 决策一致。主机 O2 的该段诊断为 0.148048 / 0.001932 s，工作集跟踪 29,502,656 / 232,896 bytes；这是结构重现，不是论文样本，也不能据此外推整条 D3 加速比。

## 实现边界

`swbam_markdup_stream_mpi.cpp` 仅当 optimistic one-pass 且存在 fallback 边界配置时，使用已有 flat table 对下一批保留的 marker 按完整 key 合并。复用既有表和工作集统计，不新增用户参数，不扩大窗口、不改变 key/score/QNAME tie-break、不改变 MPI、CPE、codec 或压缩暂存调度。

普通 `halo=enabled` 不启用此合并，非 optimistic 接口保留原有行为。标记失效仍依据原滑动窗口坐标；被标记的 singleton、pair winner、远端 keep/duplicate 和 EOF 决策均完整保留。不存在物理输入/临时 I/O 扣时、删改数据或以不完整输出替代写出。

回退结束时打印 `md-window-reuse`：marker_visits 与 marker_retained 是各批累计 metadata 访问/保留次数，marker_peak 是单批下一轮保留 marker 的峰值；不是 BAM 记录计数，也不是文件或 I/O 字节。

## 检查与冻结

- `tools/markdup_window_test.cpp` 增补普通密集坐标重现及 40 组确定性随机多批判重，与同语义全量参考逐记录决策比较；覆盖 marker 先于/晚于 singleton、配对同分、QNAME、远端路由及窗口过期。现有边界和正常路径回归全部保留。普通与 ASan/UBSan 通过。
- 四项新驱动配置检查、既有八项驱动检查和九项日志验收检查通过，`git diff --check` 通过。
- 独立 RAM 构建 `build_bounded_io/app_markdup_markers/`，SHA-256 `c7fcfe9e5e06aa88550dac2463d3b73b2308ea19fc2c63f58808ab85d5acb297`；磁盘独立构建 `build_bounded_io/app_markdup_markers_posix/` 通过，SHA-256 `3f47bea3c51b170f7078f1a593e86c53b83b1fbd3341ba198504d08a84f0fafc`。markdup 从核四个执行代码段与上一版逐字节一致。
- 原 241.256594 s 的增量程序和 D1/D2 正式程序均未覆盖，源码快照保存在新链 `source-snapshot/`，冻结源码与脚本分别记录指纹。

## 六核组验证链

```text
驱动：tools/markdup_marker_recheck.py
tmux：markdup_markers_20261017
本地：/home/wzs/work/data/D3-mouse-SRR7768905/markdup_markers_20261017_135600
远端：/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data/swbam_bench/markdup_markers_20261017_135600
```

与上一版使用同一完整排序输入、六节点每节点一个活跃核组、`--stream -l 300 -m 8G --compress-level 1`。D3/D2 RAM 池 24G/rank，D1 8G/rank，小数据 2G，复读 8G。依次执行小数据逐记录比较、D3 全量实际输出及保存、完整复读、D1 old/new 和 D2 new/old 单次控制。只有完整输出成功且复读记录数/排序/duplicate 总数正确，才接受性能成绩。

D3 新程序必须有六 rank marker 合并计数且确实减少保留项；D1/D2 必须正常 halo、linear lookup、没有 fallback 合并计数。D1/D2 单次新旧慢于 5% 标记需复核，不替代三轮正式成绩。首个 D3 结果正常执行预算约 10--20 分钟，整条链约 45--90 分钟，排队与共享存储波动另计；不自动启动正式矩阵，也不因目标时间而删除失败或慢样本。

13:56:23 小数据作业 `8625241` 验收，processing=0.151996 s；13:56:29 完成 200,000 条、22,259 个 duplicate bit、逐记录字节/顺序、排序及 gzip CRC 比较。D3 全量 `8625242` 于 13:56:29 提交并 RUN，六节点 `26628-26633`。

## 同版本两次确认

为避免只用最快单次宣布达标，使用 `tools/confirm_markdup_markers.py` 在首轮完整复读通过后，补同输入、同六核组、同参数、同指纹的两次 markdup，仍生成完整 RAM 输出；仅 `--no-dump` 省去 MPI_Finalize 后的物理保存，不省略计时内的完整 RAM 写出。首轮已保存并完整复读，后两次检查生成阶段数量、duplicate 总数、rank 成功和回退路径。保留全部健康样本，最终取三次中位数并报告范围，不仅取最小值。不启动其他应用或四步正式矩阵。

```text
tmux：markdup_marker_confirm_20261017
本地：/home/wzs/work/data/D3-mouse-SRR7768905/markdup_marker_confirm_20261017_140500
远端：/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data/swbam_bench/markdup_marker_confirm_20261017_140500
```

驱动校验首轮 binary SHA、完整生成日志和完整复读结果，最长等待两小时；首链 FAILED 不继续提交。六项汇总/版本隔离/样本数检查通过。两次确认通常约需 12--25 分钟，D1/D2 配对链通常还需 25--45 分钟，排队及共享存储波动另计。完成后分别读取 `confirmation_summary.json` 与主链 `summary.json`，不把单次控制或新性能样本覆盖历史正式成绩。

14:08:09 第二次 D3 `8625269` 已 RUN，节点 `26634-26639`，与 D1 控制的六节点分开。两次确认暂估 14:20--14:35 完成，D1/D2 控制暂估 14:35--15:00，存储波动和排队另计。代码、二进制均冻结，不再为追求更小的单次值改动本轮程序。
