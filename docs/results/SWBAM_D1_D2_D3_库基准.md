# SWBAM D1：当前 SDK 只读与读写基准

## 1. 本轮配置

首次启动：2026-09-27 19:44（北京时间）。当前状态：20:16 启动关闭诊断的新正式批次，
两次 6 核组读写复测通过，20:19:59 自动恢复 1～12 核组扩展性测试。
本轮测试当前四流水线 SDK，不使用旧 Raw/MPE bam1_t 示例，也不把 flagstat 当作 count。

| 项目 | 配置 |
|---|---|
| 代码版本 | `a47c66850ac0c0141c92c8c85038994fad41d5cb`；本轮不修改算法 |
| 输入 | D1-WES-9G.coordinate.bam |
| 压缩输入 | 9,891,133,304 bytes |
| 记录 / BGZF blocks | 177,798,571 / 814,215 |
| SDK 只读 | `swbam-sdk-optimized-count`，CPE 解压、解析 bam1_t 后计数 |
| SDK 读写 | `swbam-sdk-optimized-bam2bam`，无过滤，BAM level 1 |
| 活跃核组数 | 1 / 2 / 4 / 6 / 8 / 12 |
| 布局 | 每节点 1 MPI rank，1 MPE + 64 CPE 参与计算 |
| 大共享内存 | `-mpecg 6 -share_size 100 -cross_size 88000 -xmalloc` |
| 计划重复 | 新正式批次每项每档 3 次，共 36 个作业，串行提交 |
| 输出 | 仅 rank-local memory body；不做最终 gather/dump，不生成验证 BAM |

其余核组用于提供内存，不计入活跃计算核组；但每个 rank 实际占用一个节点的资源，
不能把六节点布局写成普通单节点六核组。SDK memory backend 每 rank 读取整份 D1，
加载总流量随 rank 数增长，因此排队和 load 可能远长于核心处理时间。

新批次先测 `rw_n6_r1`、`rw_n6_r2`，两次成功、输出 body 大小相同且最大/最小 core
不超过 1.20 时继续。两次计入正式重复，完整矩阵仍为 36 个作业。失败则停止，不自动提交扩展性。
随后第一轮资源顺序 `6,1,2,4,8,12`，第二轮 `12,8,6,4,2,1`，第三轮 `1,2,4,6,8,12`，
跳过已经完成的两次 6 核组读写。第二轮先测读写再测只读，其余轮相反。
没有额外预热作业，不清理文件系统缓存。
每个作业执行时限 30 分钟；失败或记录校验不通过时驱动停止，不自动重试。

## 2. 时间与正确性边界

主指标为 SDK 打印的 `core`，取各 rank 流水线墙钟时间的最大值：

- 不计入磁盘 load、全局 block scan/plan、最初的 body reserve 和最终 MPI 统计归约。
- 只读包括 batch 内存读取、CPE 解压/解析/count 和 MPE 批后归并。
- 读写包括 batch 内存读取、解析、pack、CPE 序列化/压缩及内存 body 写入。
- 读写不计完整 BAM header/EOF、输出汇集；不是 app 的 4.1~4.6，也不是磁盘 wall。

每次检查进程成功退出、实际 ranks、records、blocks、input_bytes、core>0；读写
还需 body_bytes>0，且示例内部要求 kept=total。沿用小数据回归已通过的实现，
本轮为避免 dump 成本不做 D1 全字段对照，不能表述为“已证明 D1 每个输出字节正确”。

`submit_wall_s` 另记整个提交等待时间，包含传输、排队、load、核心与 MPI 生命周期，
只用于估计何时跑完，不用于计算库吞吐量。

## 3. A4 对照与待填结果

沿用 [D1/D2/D3 库基准](samtools_D1_D2_D3_库基准.md) 中 A4 的 SAMtools 1.20：

| 测试 | A4 最佳配置 | A4 wall 中位数 / s | 本轮 SWBAM 6 核组 core 中位数 / s |
|---|---|---:|---:|
| 只读 view -c | `-@8` | 11.36 | 待完成 |
| BAM level 1 读写 | `-@48` | 31.29 | 待完成 |

按用户当前口径对照 A4 bigssd wall 与 SWBAM memory core，耗时比为前者除以后者；
不是同 I/O 边界的跨平台端到端加速比。正式结果报告 3 次数据、中位数、相对 1 rank
的加速比/效率、M records/s、输入 GB/s；读写补充输入+body GB/s，GB 为十进制。

### 首轮探测（不是三次中位数）

| 项目 | 作业号 | ranks | core / s | 提交到结束 / s | A4 wall / 本次 core |
|---|---:|---:|---:|---:|---:|
| 只读，第 1 次 | 8319214 | 6 | 2.468453 | 142 | 4.60x |
| 读写，第 1 次 | 8319246（取消） | 6 | | | |
| 读写，诊断重跑（不计正式重复） | 8319443 | 6 | 8.580724 | | 3.65x |

首个只读作业实际 records=177,798,571，blocks=814,215，input_bytes=9,891,133,304，
均与预期一致。诊断读写也符合这些值；两个单次 core 均低于对应 A4 wall 基线，
但不能替代后续三次中位数，且诊断读写开启了日志。

## 4. 运行与取回结果

脚本：[tools/bench_sdk_d1.sh](../../tools/bench_sdk_d1.sh)。实际命令模板：

```bash
bsub.py swls -q q_share -I -b -K -J "sdk_d1_${MODE}_n${N}_r${R}" \
  -N "$N" -np 1 -mpecg 6 -cgsp 64 \
  -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc \
  -timelimit 00:30:00 -o "logs/${MODE}_n${N}_r${R}.log" \
  -x SWBAM_DIAGNOSTICS=0 \
  "${RUN_DIR}/bin/${EXE}@${REMOTE_DIR}" \
  ../../../D1-WES-9G/initial/D1-WES-9G.coordinate.bam
```

`MODE=read/rw` 对应两个 SDK 可执行文件。脚本保存二进制快照，不受后续重新构建影响。

```text
原批次本地 RUN_DIR（已停止，历史日志保留）：
/home/wzs/work/code/RabbitBAM/build_sunway_sduhpc/benchmarks/sdk_D1_20260927_194500_a47c668

原批次 SWLS REMOTE_DIR（登录目录下相对路径）：
online/guoshi/wzs/wzs_data/swbam_bench/D1_sdk_current/sdk_D1_20260927_194500_a47c668
```

RUN_DIR 下 `results.tsv` 每完成一项追加一行；`state` 显示当前任务或 COMPLETE/FAILED；
`logs/*.submit.log` 保存提交和完整输出，远端 `logs/*.log` 保存调度日志。

原驱动由 sduhpc 的 user systemd 服务托管，当前已停止：

```bash
systemctl --user status swbam-sdk-d1-20260927-194500 --no-pager
```

不要在同一 RUN_DIR 另开驱动或覆盖快照。结束本次 Codex 会话不影响服务；当前用户
`Linger=no`，应保留 sduhpc 登录会话，不保证退出所有系统登录后服务仍存在。重启
sduhpc 会中断后续提交。要终止时还需检查 SWLS 是否仍有已提交作业，停止本地服务不保证
调度器内的作业同时取消。

| 原批次快照 | SHA-256 |
|---|---|
| count | 788c5e60661adb831eeee5844325b76920c659e2e0570aca4038b80272c1a1df |
| bam2bam | c130b9a57f90ac83b3e4c62f13763fcae1fa87216a3f52ba593b3a5ab55c721d |

完整结果尚未产生，不用小数据成绩或旧 app 成绩替代本轮 SDK 的 D1 成绩。

首次目录 `sdk_D1_20260927_194000_a47c668` 的作业 8319188 使用登录节点绝对输入
路径，在读入阶段退出，未产生 SDK 结果。该日志保留；当前批次改用 job 目录下的
共享文件相对路径。失败试投不计入三次有效重复，也不计入性能统计。

## 5. 读写诊断（2026-09-27）

原读写作业 8319246 运行超过 15 分钟仍未输出 SDK 汇总；原程序只在最终打印结果，
无法据此判定停在 load 还是计算。已停止整轮自动提交并取消该作业，不将等待时间当作 core。

新增可选开关 `SWBAM_DIAGNOSTICS=1`，不修改算法、CPE kernel 或数据流：

- 每 rank：从核初始化、输入 load/header、MPI 汇合、scan/plan、body reserve、流水线与退出。
- BAM read：初始化、首批和每 256 批的 kernel/post-processing 进度、最终 flush/清理。
- 日志写 stderr 并主动刷新；未设置开关或设为 `0` 时不打印。
- 诊断日志包含在运行耗时内，不能把本轮时间并入正式性能中位数。

两个 SDK 目标已重新构建成功，`git diff --check` 通过。本次只重跑一个 6 核组读写作业，
作业号 **8319443**，其余资源参数和输入不变，增加 `-x SWBAM_DIAGNOSTICS=1`。

```text
本地诊断目录：
build_sunway_sduhpc/benchmarks/sdk_D1_diag_20260927_201000
日志：submit.log
远端目录：
online/guoshi/wzs/wzs_data/swbam_bench/D1_sdk_current/sdk_D1_diag_20260927_201000
远端日志：logs/rw_diag.log
服务：swbam-sdk-d1-diag-20260927-201000
```

完整诊断命令（在本地诊断目录执行）：

```bash
bsub.py swls -q q_share -I -b -K -J sdk_d1_rw_diag \
  -N 6 -np 1 -mpecg 6 -cgsp 64 \
  -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc \
  -timelimit 00:30:00 -o logs/rw_diag.log -x SWBAM_DIAGNOSTICS=1 \
  ./bin/swbam-sdk-optimized-bam2bam@online/guoshi/wzs/wzs_data/swbam_bench/D1_sdk_current/sdk_D1_diag_20260927_201000 \
  ../../../D1-WES-9G/initial/D1-WES-9G.coordinate.bam
```

诊断二进制 SHA-256：`e617161e4186798901dc910f25556c9dfa26f9103f28fb95e2e615191395babc`。
本次节点为 `33448-33453`；作业正常结束，本地提交服务退出状态为 0。

| 阶段 / 指标 | 本次结果 |
|---|---:|
| 各 rank load | 62.524 / 65.945 / 67.953 / 68.776 / 68.951 / 91.906 s |
| scan/plan | 约 0.236～0.241 s/rank |
| 初始 body reserve 调用 | 均小于 0.001 s（不代表全部页面已物理写入） |
| 各 rank pipeline | 约 8.535～8.581 s |
| 最大 core | 8.580724 s |
| records / blocks | 177,798,571 / 814,215 |
| 总 body bytes | 11,467,657,614 |
| 输入吞吐量 | 1.153 GB/s |
| 输入 + body 吞吐量 | 2.489 GB/s |
| 记录吞吐量 | 20.721 M records/s |

6 个 rank 的 `pipeline.done` 均为 `ok=1`，MPI 归约、内存释放和退出正常。未写磁盘 BAM，
因此这里只确认运行与计数自检通过，不宣称完成 D1 输出内容的逐字段验证。

**结论：** 本次没有复现长时间停滞，主要前置等待来自 load 及等待最慢 rank；
核心读写本身已在约 8.6 s 完成。上一次作业没有阶段日志，不能追溯判定它具体卡在何处，
也不能仅凭本次成功排除间歇性程序问题。后续先关闭诊断重复 6 核组读写，稳定后再恢复
1～12 核组扩展性。诊断结束时未自动提交后续作业；新正式批次见下一节。

## 6. 关闭诊断后恢复正式测试

2026-09-27 20:16 启动新批次，显式传入 `-x SWBAM_DIAGNOSTICS=0`。
使用包含可选诊断代码的新二进制快照；核心算法与诊断重跑相同。此前只读探测和诊断
读写都保留为探测结果，不混入本批次的三次中位数。

```text
本地 RUN_DIR：
/home/wzs/work/code/RabbitBAM/build_sunway_sduhpc/benchmarks/sdk_D1_20260927_201700_diagoff
SWLS REMOTE_DIR：
online/guoshi/wzs/wzs_data/swbam_bench/D1_sdk_current/sdk_D1_20260927_201700_diagoff
服务：swbam-sdk-d1-20260927-201700
```

| 当前快照 | SHA-256 |
|---|---|
| count | 2848f114e30322d6fc462df9e7490b58692f5aa259b34b5e5a1596a654dd4ce2 |
| bam2bam | e617161e4186798901dc910f25556c9dfa26f9103f28fb95e2e615191395babc |

| 6 核组读写正式复测 | 作业号 | core / s | 提交到结束 / s | body bytes |
|---|---:|---:|---:|---:|
| 第 1 次 | 8319501 | 8.579890 | 105 | 11,467,657,614 |
| 第 2 次 | 8319520 | 8.575233 | 101 | 11,467,657,614 |

两次记录数、block 数和输入字节数符合 D1 预期，无诊断日志，作业正常退出。
两次 core 差约 0.054%，输出 body 大小一致，稳定性检查通过。按本轮比较口径，
两次读写与 A4 最佳 wall 31.29 s 的耗时比均约 3.65x；正式中位数仍待第 3 次完成。
稳定性判定与后续提交由脚本自动完成；`scaling_ready` 文件出现表示两次检查通过，
`results.tsv` 持续追加有效结果，`state` 显示当前作业或 `COMPLETE/FAILED`。

20:19:59 已开始 `read_n6_r1`，剩余 34 项串行运行。依据两次 6 核组作业整体等待
105/101 s，以及更多 rank 同时加载完整 D1 的额外流量，预计从 20:20 起还需约
**60～100 分钟**，即约 21:20～22:00 完成。排队或文件系统再次停滞会延长该估计。
本地服务持续提交，失败或超时立即停止，不自动重试；本次会话结束后继续运行。
