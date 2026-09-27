# Sort 公共 BAM 读写流水线迁移验证

日期：2026-09-26。以下为小数据单次 smoke，不作为稳定加速比或完整边界测试。

## 改动范围

- 内排、严格外排均使用 `RunCpeReadPipeline`，app 的提取 operator 绑定原 CPE
  raw/key kernel，并将结果保存到持久 raw/metadata 或有界 run arena。
- 最终输出使用 `RunCpeWritePipeline`；归并/打包 source 和统计适配留在 sort app，
  原 CPE 压缩 kernel 通过公共 `CpeWriteKernelOperator` 调度。
- 写流水线新增可选 `overlap_source` / `SubmitWithPrefetch()`：当前批压缩期间，
  MPE 写出上一批并准备下一批。默认关闭 source overlap，已有调用继续原行为。
- 排序键、稳定次序、采样、MPI exchange、外排文件和归并算法未变。
- 公共读端四套 batch buffer 纳入外排 `-m` 跟踪，较旧读 workspace 多约
  8 MiB/rank，必要时减少 run arena。该预算不是进程总 RSS，也不包括完整输入驻留。

## 配置与结果

SWLS 单节点、6 MPI ranks、每 rank 64 CPE；memory 输入/输出，level 1。
输入 `../performance_test/WES_0.25G.bam`：12,100 个 BGZF 数据块，2,668,351 条记录。

| 测试 | 作业号 | 4.3 / s | 4.1~4.6 / s | 说明 |
| --- | ---: | ---: | ---: | --- |
| 迁移前内排，`-m 8G` | 8302402 | 0.664794 | 0.680522 | 实测核心 wall |
| 迁移后内排，`-m 8G` | 8302469 | 0.663378 | 0.680881 | 实测核心 wall |
| 迁移后外排，`-m 128M` | 8302480 | 0.669081 | 0.695606 | 原模拟高速临时 I/O 模型 |

内排 4.3 单次差约 -0.21%，本轮未观察到明显回退。上述核心时间不包含完整文件
load、最终 gather/dump；4.1~4.6 还排除模拟输出内存分配 4.5b。

外排实际 4.3 wall 为 **32.652825 s**，不能把模拟的 0.669081 s 当成真实磁盘性能。
本轮没有补跑旧版外排性能基线。外排每 rank 生成 3 个 run，总计 18 个；
`tracked_peak_max=134217728`，没有超过 `128 MiB/rank` 算法预算。

## 正确性

- 迁移前后内排 BAM：`cmp` 返回 0，完整文件逐字节一致。
- 迁移后内排与外排：压缩块划分不同，解压后的完整 BAM 字节流 `cmp` 返回 0。
- `stats --basic` 复读两份输出（8302512、8302518）：均为 2,668,351 条记录，
  `is sorted=1`，primary=2,656,392，supplementary=11,959，NM mismatches=878,007。
- 公共 `io-check`（8302520）通过 memory/POSIX/MPI-IO 输入、128-block 合成压缩、
  完整 BAM 回环及 MAPQ >= 30 过滤回环；过滤保留 2,481,669 条。
- Sunway `RabbitBAM-MPI` 构建和 `git diff --check` 通过。

未覆盖所有压缩级别、空 rank、跨块畸形输入或外排多级 consolidation；本轮保持
block-aligned 输入约束。dedup-pipeline 的 sort helper 复用同一核心，但未另跑完整去重回归。

## 复现命令

在 `build_sunway_sduhpc` 执行。迁移前后分别运行下列命令，使用 `before` / `after`
文件名区分输出；下例为迁移后：

```bash
bsub.py swls -q q_share -I -b -J sort_pipeline_after \
  -N 1 -np 6 -cgsp 64 -share_size 14000 -cache_size 32 \
  -o sort_pipeline_after_20260926.log \
  -x LD_LIBRARY_PATH=../ext/htslib-1.20:../ext/libdeflate-1.20/build \
  ./RabbitBAM-MPI@online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc \
  sort -i ../performance_test/WES_0.25G.bam \
  -o ../performance_test/sort_pipeline_after_20260926.bam \
  --compress-level 1 -m 8G
```

外排使用相同资源参数，改为 `-m 128M -T ../performance_test/tmp/sort_pipeline_20260926`，
输出/日志分别为 `sort_pipeline_external_20260926.bam` / `.log`。
验证任务沿用同样资源参数，命令分别是：

```bash
stats -i ../performance_test/sort_pipeline_after_20260926.bam --basic
stats -i ../performance_test/sort_pipeline_external_20260926.bam --basic
io-check -i ../performance_test/WES_0.25G.bam
```

在 SWLS 的 `performance_test` 目录进行字节比较：

```bash
cmp sort_pipeline_before_20260926.bam sort_pipeline_after_20260926.bam
cmp <(gzip -cd sort_pipeline_after_20260926.bam) \
    <(gzip -cd sort_pipeline_external_20260926.bam)
```

第二条需 Bash。日志保存在 SWLS
`online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc/`：
`sort_pipeline_{before,after,external,stats,external_stats,io_check}_20260926.log`。
