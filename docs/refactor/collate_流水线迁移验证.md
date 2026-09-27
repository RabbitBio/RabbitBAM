# Collate 公共 BAM 读写流水线迁移验证

日期：2026-09-26。WES 小数据单次 smoke，目标是确认接口迁移后的正确性和明显性能回退，
不作为稳定加速比或完整边界测试。

## 改动范围

```text
公共 BAM read -> collate app 算法 -> 公共 BAM write
```

- 内排、外排均通过 `CollateReadSource` / `CollateExtractOperator` 接入
  `RunCpeReadPipeline`。复用 `BamReadBatchOperator` 和 `BindBamReadBatch`，
  保留原 CPE 解压、QNAME/hash/meta 提取 kernel。
- 内排继续直接解压到持久 raw arena；外排在批后处理里追加到 run arena，满后 spill。
  未引入先解压到临时 buffer、再复制整个内排 raw 的路径。
- `CollatePayloadSource` 将 app loser tree 归并结果逐批交给 `RunCpeWritePipeline`。
  开启 `overlap_source`，保留压缩当前批时 flush 上一批、填充下一批的覆盖关系。
- hash/bin、稳定记录序号、owner、self/remote exchange、分桶排序、cursor cache、
  外排 consolidation 和最终输出布局不变。内排不增加逐记录虚调用。
- 公共读端持有两套 compressed、两套 decoded batch，读 workspace 较旧版增加约
  8 MiB/rank；外排固定预算取读、交换、压缩阶段的最大值，并据此分配 run arena。
  `tracked_peak` 不是进程总 RSS，也不包含完整输入驻留和最终输出内存。

实现位于 `apps/mpi/algorithms/collate/swbam_collate_mpi.cpp`，没有改动 CPE kernel、
CLI 或公共流水线实现；写端复用上一轮 sort 迁移新增的 source overlap 接口。

## 配置与性能

SWLS 单节点、6 MPI ranks、每 rank 64 CPE；memory 输入/输出、压缩 level 1、`-n 66`。
输入 `../performance_test/WES_0.25G.bam`：12,100 个 BGZF 数据块、2,668,351 条记录。

| 测试 | 作业号 | 4.3 / s | 4.1~4.6 / s | 计时说明 |
| --- | ---: | ---: | ---: | --- |
| 迁移前内排，`-m 8G` | 8302605 | 0.690722 | 0.708037 | 核心 wall |
| 迁移后内排，`-m 8G` | 8302676 | 0.688869 | 0.705829 | 核心 wall |
| 迁移后外排，`-m 128M` | 8302682 | 1.416973 | 1.440407 | 模拟高速临时 I/O 时间 |

内排 4.3 单次变化约 **-0.27%**，本次未观察到明显性能回退。核心时间不包含完整
输入 load、最终 gather/dump；4.1~4.6 还排除模拟输出内存分配 4.5b。

外排实际 wall 为 **27.651204 s**，不能把模拟值 1.416973 s 当成真实磁盘性能。
每 rank 2 个 run，共 12 个；`run_arena_bytes=88070144`，
`tracked_peak_max=134217728`，consolidation_passes=0。未补跑旧版外排性能基线。

### 分项解读

| 内排 max / s | 迁移前 | 迁移后 |
| --- | ---: | ---: |
| extract CPE | 0.048 | 0.049 |
| extract merge | 0.051 | 0.050 |
| local sort | 0.124 | 0.124 |
| exchange | 0.215 | 0.214 |
| merge wall | 0.199 | 0.200 |
| compress pipeline | 0.194 | 0.194 |

这些 max 可能来自不同 rank，不能相加还原核心时间。`compress_pipeline` 已包含
覆盖期间的 MPE fill/flush；`compress_fill` 和 `write` 是嵌套观察值。
`extract_alloc/free` 现在只统计 operator metadata 池，公共 batch 分配释放仍计入
核心总时间。本轮关闭 read overlap，保持 read 与 CPE 提取分项不重叠。

外排沿用 `simulated = actual - temp_write_actual - temp_read_actual`；模拟 memcpy
已包含在 actual 内，不再额外加一遍。旧 `accounted/unaccounted` 使用了未扣真实
临时 I/O 的 exchange/merge 分项，本次外排也出现负 unaccounted；这不是负执行时间，
也不能用于解释 CPU 时间构成。本轮没有改这一既有统计口径。

## 正确性

- 迁移前后内排 BAM：完整文件 `cmp` 返回 0。
- 迁移后内排与 `128M` 外排 BAM：完整文件 `cmp` 返回 0。
- 两种模式均得到 records=received=2,668,351，QNAME groups=1,328,196，
  BGZF data blocks=12,102。
- `run_all` 将迁移后内排 BAM 转为 SAM（8302696），成功解析并输出 2,668,351 条。
- `collate_verify.sh`（8302702）完整输出：records=2,668,351，
  groups=unique_qnames=1,328,196，`continuity_errors=0`、`read_order_errors=0`；
  `record_multiset_md5=38cf435adf7c60b4135660003b82ef42`。
- Sunway `RabbitBAM-MPI` 构建及 `git diff --check` 通过。

注意：验证脚本虽然打印了上述指标和 `collate_verify: done`，调度器仍报告任务
exit_code=1。本轮未定位该退出状态，不能将脚本任务记为正常退出；数据检查结果与
两次独立 `cmp`（均返回 0）分别记录。collate 迁移前、迁移后、外排和 BAM2SAM
四个计算任务本身均正常结束。

未测试全部内存档位、空 rank、多级 consolidation、压缩级别或 POSIX/MPI-IO 输入。
`MpiCollateMemoryToMemory` 使用同一核心，但本轮未重跑完整 dedup-pipeline。

## 复现命令

在 `build_sunway_sduhpc` 执行，迁移前后分别用 `before` / `after` 区分文件名：

```bash
bsub.py swls -q q_share -I -b -J collate_pipe_after \
  -N 1 -np 6 -cgsp 64 -share_size 14000 -cache_size 32 \
  -o collate_pipe_after_20260926.log \
  -x LD_LIBRARY_PATH=../ext/htslib-1.20:../ext/libdeflate-1.20/build \
  ./RabbitBAM-MPI@online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc \
  collate -i ../performance_test/WES_0.25G.bam \
  -o ../performance_test/collate_pipe_after_20260926.bam \
  --compress-level 1 -m 8G -n 66
```

外排使用同样资源参数，改为 `-m 128M -T ../performance_test/tmp/collate_pipe_20260926`，
输出/日志为 `collate_pipe_external_20260926.bam` / `.log`。
转换验证使用同样资源参数：

```bash
run_all -i ../performance_test/collate_pipe_after_20260926.bam \
  -o ../performance_test/collate_pipe_after_20260926.sam
```

QNAME 验证只需一个 rank：

```bash
bsub.py swls -q q_share -I -b -J collate_pipe_verify \
  -N 1 -np 1 -cgsp 64 -share_size 14000 -cache_size 32 \
  -o collate_pipe_verify_20260926.log -x SORT_MEMORY=8G \
  ../tools/collate_verify.sh@online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc \
  ../performance_test/collate_pipe_after_20260926.sam ../performance_test/tmp
```

在 SWLS `performance_test` 目录比较：

```bash
cmp collate_pipe_before_20260926.bam collate_pipe_after_20260926.bam
cmp collate_pipe_after_20260926.bam collate_pipe_external_20260926.bam
```

日志保存在 SWLS 的 `online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc/`：
`collate_pipe_{before,after,external,to_sam,verify}_20260926.log`。
