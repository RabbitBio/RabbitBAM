# 流式 Markdup 公共 BAM 读写流水线迁移验证

日期：2026-09-26。WES 小数据单次 smoke，不作为稳定加速比或完整边界覆盖。

2026-09-27 补充：完整四步串联触发了旧 halo 候选被批回调清空的漏标问题，已修复。
本文保留迁移时历史结果；最新根因、性能及完整 BAM 回归见
[去重串联与边界回归](去重串联_边界回归验证.md)。

## 改动范围

默认单遍 `--stream` 接入：公共 BAM read -> app 候选/窗口/决策 -> 公共 BAM write。
CLI、kernel、判重语义、边界 halo、输出顺序及 block-aligned 约束不变。

- `MdExtractCandidatesPass` 改用 `RunCpeReadPipeline` 和 app `MdCandidateReadOperator`。
  `DuringKernel()` 在当前槽解压时归并上一槽候选，保留两槽覆盖；窗口后处理仍在
  decode join 后执行。候选提取 CPE 保留为 app 专用阶段。
- 单遍写端使用 `CpeWritePipelineSession::Start/Complete`，压缩期间执行 MPE
  窗口计算和 MPI 交换，完成后 Flush，再退休 pending blocks；异常清理先 join。
- decoded arena 继续由 app 持有并转移到 pending 环，不新增整批 raw 拷贝。
  超过 CPE 安全压缩尺寸的合法块仍按完整 record 边界重新打包。
- 公共 `Submit/SubmitWithPrefetch` 复用 Start/Complete，同步调用的生命周期不变。

候选读端同时被两遍和默认非流式算法复用。`-r/-c`、强制两遍 fallback 的第二遍
rewrite/write 以及首坐标预探测尚未迁移；`swbam_markdup_stream_mpi.cpp` 本轮未改动。
完整 dedup-pipeline 没有在本轮重新测试。

## 环境与性能

SWLS 单节点、6 MPI ranks、每 rank 64 CPE；memory 输入/输出，压缩 level 1，
`-m 8G --stream -l 300`。输入 `../performance_test/fixmate.sorted.bam`：
241,865,073 bytes、12,675 个 BGZF 数据块、2,668,351 条记录。

| 测试 | 作业号 | 4.3 / s | 4.1~4.6 / s |
| --- | ---: | ---: | ---: |
| 迁移前单遍 stream | 8303098 | 0.724634 | 0.740173 |
| 迁移后单遍 stream | 8303222 | 0.722752 | 0.740504 |
| 迁移后 stream -r，两遍 fallback | 8303231 | 0.769285 | 0.785518 |
| BAM2BAM 复读回环，共同写接口回归 | 8303240 | 0.183315 | 0.197702 |

单遍 4.3 约 -0.26%，核心总时间约 +0.04%，视为基本持平，不宣称稳定加速。
`-r` 未另跑迁移前性能基线，表中不是它的迁移前后对比。
整文件输入 load、最终 gather/dump 不计核心时间，4.1~4.6 还排除输出 malloc 4.5b。

### 内存与分项

| 单遍指标 | 迁移前 | 迁移后 |
| --- | ---: | ---: |
| tracked_peak_max / bytes | 126,461,103 | 139,046,575 |
| tracked_peak_max / MiB | 120.60 | 132.61 |
| stream_pending_peak_max / bytes | 12,619,908 | 12,619,908 |
| rewrite_decomp / s | 0 | 0 |
| extract_merge max / s | 0.160 | 0.158 |
| group max / s | 0.191 | 0.192 |
| write max / s | 0.026 | 0.026 |

增加约 12 MiB/rank 固定 workspace：公共 read 额外的两套 decoded buffer 约 8 MiB，
公共 write 比旧单输出槽多约 4 MiB。app 自有 decoded arena 的零拷贝移交保留，
pending 峰值未增加；新增容量已纳入原有 `-m` 检查。
`tracked_peak` 不是进程 RSS，不含 backend 的完整输入与最终输出驻留。

读端 `candidate_decomp` 改为 runtime kernel wall 减重叠 candidate merge，包含
少量 launch/调度开销，不能把约 0.002 -> 0.003 s 全部解释为解压变慢。
`compress` 仍是 Complete 暴露等待及校验时间；公共 write `timing.kernel` 则包含
Start 到 Complete 之间的窗口计算/通信，不能重复相加为核心时间。

## 正确性

- 迁移前后单遍完整 BAM `cmp` 返回 0，包括 header、record/AUX 和压缩字节。
- 输入/输出 records=2,668,351，marked=15,695；pair/single duplicates 分别为
  15,550 / 145；输出 BGZF blocks=12,675，与基线一致。
- BAM2BAM 回环保留全部 2,668,351 条记录。解压后逐条比较完整 BAM record，
  包括 AUX、FLAG 和顺序，均与单遍输出一致；不是只比较聚合计数。
- `--stream -r` 删除 15,695 条，保留 2,652,656 条。逐条比较确认输出恰好等于
  单遍基线按 `FLAG & 1024 == 0` 筛选后的有序记录序列，包括全部 AUX。
- Sunway 构建和 `git diff --check` 通过，四个计算任务正常退出。

未测试所有 backend/spool 组合、空 rank、异常输入、压缩级别、halo fallback 或
全部 `-c/-r` 组合。公共同步 writer 的错误回调分支进行了代码检查，没有故障注入。

## 复现命令

在 `build_sunway_sduhpc` 运行，基线先用迁移前二进制，并将 `after` 换成 `before`：

```bash
bsub.py swls -q q_share -I -b -J markdup_pipe_after \
  -N 1 -np 6 -cgsp 64 -share_size 14000 -cache_size 32 \
  -o markdup_pipe_after_20260926.log \
  -x LD_LIBRARY_PATH=../ext/htslib-1.20:../ext/libdeflate-1.20/build \
  ./RabbitBAM-MPI@online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc \
  markdup -i ../performance_test/fixmate.sorted.bam \
  -o ../performance_test/markdup_pipe_after_20260926.bam \
  --compress-level 1 -m 8G --stream -l 300
```

两遍回归增加 `-r`，输出和日志使用 `markdup_pipe_remove_20260926`。
BAM2BAM 回环沿用同样资源参数，命令部分替换为：

```bash
run_all -i ../performance_test/markdup_pipe_after_20260926.bam \
  -o ../performance_test/markdup_pipe_roundtrip_20260926.bam --compress-level 1
```

SWLS `performance_test` 目录：

```bash
cmp markdup_pipe_before_20260926.bam markdup_pipe_after_20260926.bam
```

回环和删除验证使用只读、常量内存的 gzip/BAM 顺序解析：分别跳过 BAM header，
读取每条长度前缀及完整 record 字节；回环逐条相等，删除版本仅与基线非 DUP
记录比对，并检查两边同时 EOF。该校验不忽略 AUX，也不依赖指纹的统计一致性。

日志在 SWLS `online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc/`：
`markdup_pipe_{before,after,remove,roundtrip}_20260926.log`。
