# Markdup POS halo 预处理验证

日期：2026-11-08。独立实验，不替换历史冻结二进制、正式性能表或论文文件。

本文记录前一版“扩读至覆盖完整，否则 fallback”的实现和冻结结果。后续固定 64 块、允许部分预处理的实现与 D3 实验见 [固定探测与部分预处理验证](markdup_固定探测与部分预处理验证.md)；本文末尾的设计讨论保留为当时方案记录。

## 本轮实现

目标是先避免不完整边界组的赢家被提前判为 KEEP，再测量主循环是否仍需传输候选。

1. 对每个 rank 分界坐标 `b`，选择源记录 **POS 位于 `[b-L,b+L]`** 的候选；重叠 halo 合并。传输候选元数据和必要 QNAME，不传整条 BAM。
2. 初始读取首尾各 32 个 BGZF 块，覆盖不足时按 32/64/128/256 块扩读。必须越过坐标端点，或已读完该 rank，才能排除端点上的同 POS 记录遗漏。达到上限或探测预算不足时保留原保守窗口。
3. 复用原有 pair/single hash、质量分数、QC、QNAME 和稳定次序比较。paired-marker 单独处理，不能按普通记录的最高分规则删除。
4. 已知竞争失败的记录可预判 DUP。只有完整位置支撑区间被选中的 POS halo 覆盖，组内最佳记录才可预判 KEEP；其余最佳候选保持 Pending，在主循环重新提取并进入原窗口。
5. 只有最终 KEEP/DUP 候选和已闭合的 marker 才从主循环过滤。未闭合 marker 仍参与判重。不能证明剩余依赖局限于 rank 内时，继续启用保守边界进度机制。
6. 保留候选路由、判定回传，并新增阶段独立通信计数；没有为追求“零通信”删除安全路径。

### 如何避免提前 KEEP

按现有 key 构造中本端使用 start/end 的角色，在位移受 `L` 约束时：

| 本端去重坐标 | 同 key 记录的可能 POS 区间 |
| --- | --- |
| unclipped start `K` | `[K,K+L]` |
| unclipped end `K` | `[K-L,K]` |

双端按完整 pair key 的 orientation/leftmost 区分本端角色，不把远端配偶坐标误作本端 POS。

例如 `b=10000,L=300`：`u: POS=10250,150M,K=10250,score=3500` 在 halo 内，但它的可能竞争区间延伸到 10550，因此不能提前 KEEP。后续 `v: POS=10320,70S80M,K=10250,score=4100` 在 halo 外，仍可在本地窗口替换 u：u 变 DUP，v 继续 Pending，安全进度到达后才最终 KEEP。

**适用条件和保护：** 这里依赖 POS 到键所用 unclipped 坐标的位移界，不能仅凭 query 长度断言所有 CIGAR 都满足。超长 D/N 或 clipping 等可能超出 `L`。若主循环出现已预判闭合 key 的新候选，程序明确失败，不静默沿用旧 KEEP；可增大 `-l` 或设置 `RABBITBAM_MARKDUP_TWO_PASS=1` 重跑。本轮不宣称支持任意 CIGAR 下的零候选通信。

## 验证结果

| 检查 | 结果 |
| --- | --- |
| Sunway APP 与 MPI/CPE halo probe 构建 | 通过 |
| 主机元数据测试 | POS 选择、fringe Pending、后续高分替换、双端 key 角色、partial marker 通过 |
| 原有滑动窗口回归 | `markdup window decisions PASS` |
| 2 rank 真实 CPE/MPI probe | 六项通过：正常 halo、unmapped 尾部、head/tail 覆盖不足回退、fringe 不提前 KEEP、端点同 POS 平台扩读；后者同时验证第二批双端 QNAME 同分比较 |
| 6 rank 小样本完整输出 | 200,000 条记录，22,259 条 DUP；全部 BAM record 字节及顺序与已验证基准一致，CRC 通过 |
| D1 全量 6 rank 首轮及最终版 | 两次均为 177,798,571 条记录、50,084,605 条 DUP，与历史计数一致；六 rank 成功 |
| D2 全量 6 rank | 365,079,492 条记录、8,894,488 条 DUP，与历史计数一致；六 rank 成功；halo enabled |
| D3 全量 6 rank | 165,373,770 条记录、33,171,753 条 DUP，与历史计数一致；六 rank 成功；halo fallback |

小样本首次作业 `9075331` 在算法开始前被自动工作内存预算检查拒绝；显式设置 `SWBAM_NODE_MEMORY=88000M` 后重试 `9075357` 成功。失败日志保留，不计为算法测试通过或性能样本。

### 通信口径

- `[md-halo rank=...]`：预处理传输字节、最终 KEEP/DUP、普通 Pending 和未闭合 marker 数。
- `[md-comm rank=...]`：仅主循环及 EOF 收尾的远程候选数、其中 marker 数、候选/QNAME 字节及判定字节。所有 rank 求和是总发送 payload，不重复累计接收量。
- `main_candidates` 包含 marker；`main_markers` 是其子集。
- 用户本轮关注的是候选记录通信：普通候选数取 `main_candidates-main_markers`，paired-marker 单列。候选 payload 包含两类元数据及必要 QNAME，不是完整 BAM record，也不把判定回传算成候选通信。
- **不包含** MPI count/status/progress collective 和协议开销；payload 为零不等于没有 MPI 调用。

| 数据 | halo | 预处理候选 payload/B | 预处理判定 payload/B | 主循环远程候选 | 主循环候选 payload/B | 主循环判定 payload/B |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 20 万条小样本，6 rank | enabled；非保守窗口 | 3,093 | 8 | 0 | 0 | 0 |
| D1 全量首轮，6 rank | enabled；非保守窗口 | 71,221 | 1,192 | 0 | 0 | 0 |
| D1 全量最终版，6 rank | enabled；非保守窗口 | 71,221 | 1,192 | 0 | 0 | 0 |
| D2 全量，6 rank | enabled；非保守窗口 | 6,821 | 24 | 0 | 0 | 0 |
| D3 全量，6 rank | fallback；保守窗口 | 0 | 0 | 1,096,071（含 marker） | 89,254,397 | 8,770,016 |

小样本各 rank 合计仍有 123 个普通 Pending、128 个未闭合 marker，但其后续处理均为本地，没有候选/判定 payload 通信。这说明“不提前 KEEP”与“主循环零候选传输”在该样本上可以同时成立，不是任意输入的保证。

D1 首轮作业 `9075374` 也是初始首尾各 32 块即覆盖成功。预处理合计最终 KEEP=1,488、DUP=1,545、普通 Pending=1,106、未闭合 marker=1,772。主循环与 EOF 收尾的候选和判定 payload 均为零；这不表示 Pending 被直接当作最终 KEEP。

| D1 诊断版本 | 作业 | processing max/s | 4.3/s |
| --- | --- | ---: | ---: |
| 首轮 POS halo | 9075374 | 52.000255 | 50.965994 |
| 最终版：保护检查先按本 rank 已闭合 key 的坐标范围预筛选 | 9075475 | 51.020365 | 49.994925 |

两版判定与通信计数一致。以上不是同一版本的多轮正式中位数，也没有同轮旧版对照，不能据此宣称算法整体加速。最终版完整 RAM 输出大小为 12,192,952,100 字节。

最终 APP SHA-256：`24cef7db4ea1ca6c111367fe4057f3a35b488fc51b24409d7af1fc9aebc2836d`。首轮 APP SHA-256：`2dd90be832e721ec392463a0b3f4f17764af1e798604b203a6e7fb5de111bac7`。两版均单独保留。最终版小样本 `9075458` 再次通过 200,000 条完整 record 字节、顺序和 CRC 比较。

端点扩读测试最初把 head 索引末端误当作首尾探测并集末端，作业 `9075459` 的测试断言失败；修正测试为检查首尾并集覆盖后，`9075549` 的六项测试通过，未因此修改生产判重逻辑。日志均保留。

进一步将该端点样本改为带 `MC/ms` 的双端候选，并在第二批放置同分竞争者，作业 `9075646` 六项再次通过，确认跨 batch 的 QNAME offset 重定位及 tie-break；生产 APP 指纹未变。

## 构建与复现

生产修改位于 `swbam-app/algorithms/markdup/`；SDK、CPE key 构造和 hash/window 判重规则未改。

```bash
cmake -S swbam-app -B build_bounded_io/app_markdup_pos_halo \
  -DCMAKE_TOOLCHAIN_FILE="$PWD/swbam-sdk/cmake/SunwayToolchain.cmake" \
  -DSWBAM_DIR="$PWD/build_bounded_io/install/lib/cmake/SWBAM" \
  -DCMAKE_BUILD_TYPE=Release -DSWBAM_BUILD_MARKDUP_HALO_TESTS=ON
cmake --build build_bounded_io/app_markdup_pos_halo -j 8
```

`tools/markdup_halo_plan_test.cpp` 是主机元数据测试；`tools/markdup_boundary_regression.cpp` 直接调用生产 probe 和真实 CPE。集群提交参数见 `tools/submit_markdup_pos_halo.sh`。

```bash
mpicxx -std=c++11 -O2 -ffunction-sections -fdata-sections \
  -I swbam-app/include -I swbam-sdk/include -I ext/htslib-1.20 \
  tools/markdup_halo_plan_test.cpp \
  swbam-app/algorithms/markdup/swbam_markdup_stream_mpi.cpp \
  -Wl,--gc-sections -o /tmp/swbam-markdup-halo-test
/tmp/swbam-markdup-halo-test
```

远端独立目录：

```text
/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data/
  swbam_bench/markdup_pos_halo_20261108/
```

小样本使用 `-N 1 -np 6`；D1 使用 `-N 6 -np 1 -mpecg 6 -cross_size 88000`，实际计算仍为 6 rank、每 rank 64 CPE。

本地原始日志、完整提交命令和源码快照：`build_bounded_io/verification/markdup_pos_halo_20261108/`。远端 `bin/RabbitBAM-MPI-final` 是最终 APP；每轮 `*.command.sh` 记录准确参数与使用的二进制。

D1 输入为已生成的 `fixmate.sorted.bam`，12,041,427,629 字节、177,798,571 条记录，不是最初压缩大小约 9.89 GB 的 BAM。D1 使用现有分布式 RAM storage，保留完整内存输出，仅跳过最终磁盘 dump；性能只取 processing/core，不计 preload 和 dump。此次 D1 只核对记录与 DUP 计数，不能替代全量 record multiset 验证。

## D2 / D3 后续实验

两项使用同一最终 APP（`24cef7db...`），不修改算法，各运行一次。资源与 D1 相同，RAM storage capacity 改为 24G/rank，仍为 `--stream -l 300 -m 8G --compress-level 1 --no-dump`。脚本为 `tools/submit_markdup_pos_halo.sh`，汇总入口为 `tools/markdup_halo_metrics.py`；后者校验六 rank 完成状态、记录和 DUP 计数，以及阶段 payload 之和与累计统计一致。

| 数据 | 作业 | 输入大小/B | 预期记录数 | 预期 DUP 数 | 状态 |
| --- | --- | ---: | ---: | ---: | --- |
| D2 | 9076523 | 40,691,022,479 | 365,079,492 | 8,894,488 | 六 rank 成功；计数通过 |
| D3 | 9076526 | 16,546,778,000 | 165,373,770 | 33,171,753 | 六 rank 成功；计数通过 |

主指标为主循环实际发送的普通候选数、marker 数和候选 payload；判定字节仅作辅助项。覆盖不足若走 fallback，将如实记录，不把它列为“预处理成功后零候选通信”。正在讨论的“固定探测、允许部分预处理”方案尚未应用到这两个作业。

| 最终版数据 | 4.3/s | processing max/s | 主循环普通候选 | 主循环 paired-marker | 主循环候选 payload/B |
| --- | ---: | ---: | ---: | ---: | ---: |
| D1 | 49.994925 | 51.020365 | 0 | 0 | 0 |
| D2 | 173.348805 | 177.031260 | 0 | 0 | 0 |
| D3 | 102.156094 | 103.745166 | 548,126 | 547,945 | 89,254,397 |

D2 初始首尾各 32 块即覆盖成功，最终 KEEP=187、DUP=9、普通 Pending=224、未闭合 marker=228；预处理候选 payload=6,821 B，主循环普通候选、marker 及候选 payload 全部为零。完整 RAM 输出为 40,718,592,510 B。与 D1 一样，这里的零只指候选数据传输，不表示没有 MPI 计数或同步调用。

D3 扩读到首尾合计 512 块仍未满足覆盖条件，因此没有执行 owner 预判，全部走原保守窗口。主循环实际发送普通候选 **548,126** 个、paired-marker **547,945** 个；合计候选 payload **89,254,397 B**（约 89.25 MB）。六 rank 发送量如下：

| 源 rank | 普通候选 | paired-marker | 候选 payload/B |
| --- | ---: | ---: | ---: |
| 0 | 1 | 1 | 163 |
| 1 | 2 | 2 | 326 |
| 2 | 545,293 | 545,337 | 88,809,185 |
| 3 | 2,830 | 2,605 | 444,723 |
| 4 | 0 | 0 | 0 |
| 5 | 0 | 0 | 0 |

D3 processing max=103.745166 s，4.3=102.156094 s，完整 RAM 输出为 16,614,394,500 B。这里只核对生成阶段记录和 DUP 计数，不等于本版全量逐记录比较；单次诊断也不构成新性能中位数。汇总 JSON 与完整日志位于本地验证目录。

上一轮失败断言作业 `9075459` 在打印 MPI_Abort 后仍占用节点，本轮已通过 `bkill` 清理；不影响有效成绩和 D2/D3 作业。

## 下一版设计讨论（尚未实现）

建议改为固定首 32 + 尾 32 块探测，允许覆盖不足时的部分预处理，不再因覆盖检查失败就丢弃整个预判机会。最小风险实现保留完整覆盖分支：

| 实际探测覆盖 | 预处理动作 | 主循环 |
| --- | --- | --- |
| 完整覆盖目标 halo | 保留现有逐 key 的 KEEP / DUP / Pending 判定 | 使用现有安全快速路径；不承诺任意输入零候选传输 |
| 覆盖不完整 | 仅最终确认已知 DUP；赢家 Pending，未闭合 marker 保留 | 保留保守窗口与候选路由，继续接收未探测候选 |

覆盖检查由“是否执行预处理”的门槛，改为决定“哪些判定能提前结束”的条件。若完全不检查覆盖，也可以先采用 DUP-only 预判，但不能因此关闭主循环边界守护。以后若需要让部分覆盖中的个别 key 也提前 KEEP，必须另建实际探测覆盖证明；不能把逻辑 `[b-L,b+L]` 当成已经读到的完整数据范围。

解码范围和预判范围应区分：当前已经解压、解析并提取首尾全部 64 块的候选，之后才按源 POS 筛选参与边界预判的集合。建议仍以 `b±L` 内已探测候选为预判重点。全部 64 块参与预判并非必然错误，按现有 owner 路由也不是所有候选都上网络，但会扩大预处理 hash、QNAME 和状态管理范围；D1 主循环候选传输已为零，扩大内区预判未必有收益。

例：`b=10000,L=300`，右侧只完整探测到 POS 10120，而某 start-key `K=9990` 的可能竞争区间到 10290。已知败者可以 DUP，最佳候选仍 Pending；主循环读到 POS 10200 的同 key 新候选时，仍可能需要跨 rank 发送。此时继续通信是正确处理未探测数据，而非预处理逻辑失效。`-l` 的窗口和边界安全作用仍保留。

`git diff --check -- swbam-app tools docs/results/markdup_POS_halo预处理验证.md` 用于本轮文件检查；工作树中用户自己的论文修改未调整，其已有行尾空白不属于本轮变更。
