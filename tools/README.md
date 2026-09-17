# Tools

本目录集中放置性能微基准、资源探针和结果验证工具，不包含 RabbitBAM 的正式命令
实现。

## 资源与性能探针

### `swbam-cgs-memory-probe`

验证神威大共享模式下单个主核控制 6 个核组，并通过 cross segment 使用超过单核组
容量的内存。工具使用 `_sw_xmalloc()` 分配，384 个 CPE 分区触页，主核校验；指定
`--file` 后还会把文件直接 `pread` 到 cross buffer，并比较主核和从核校验和。

```bash
bsub.py swls -q q_share -I -b -J cgs_memory_probe \
  -N 1 -np 1 -mpecg 6 -cgsp 64 \
  -share_size 4000 -cross_size 24000 -cache_size 32 \
  -o cgs_memory_probe.log \
  ./swbam-cgs-memory-probe@online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc \
  --allocate 20G --file ../performance_test/WES_0.25G.bam
```

这里的 `share_size` 与 `cross_size` 是不同地址空间的配额。要聚合使用多个核组的
主存，应把大缓冲放到 cross segment。探针显式使用 `_sw_xmalloc()`；正式 APP 可在
作业参数中增加 `-xmalloc`，让普通 `malloc/realloc/free`（包括 `std::vector`）使用
cross segment。正式运行前应为程序代码、栈、控制结构、算法 workspace 和输出缓冲
预留空间，不要把物理内存全部配满。

探针中的 `address_span_rate` 仅表示 CPE 分页触及的地址跨度除以耗时，不是实际内存
带宽；该测试用于验证容量、地址可见性和 384 CPE 访问正确性。

### 现有 `RabbitBAM-MPI` 使用六核组内存

不需要改造 `RabbitBAM-CGS`，也不需要给 APP 增加专用 allocator。采用六节点、每节点
一个 MPI rank，并为每个 rank 预留六个核组和 cross segment：

```bash
bsub.py swls -q q_share -I -b -J mpi_large_memory \
  -N 6 -np 1 -mpecg 6 -cgsp 64 \
  -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc \
  -o mpi_large_memory.log \
  -x LD_LIBRARY_PATH=../ext/htslib-1.20:../ext/libdeflate-1.20/build \
  ./RabbitBAM-MPI@online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc \
  run_all -i input.bam -o output.bam \
  --io-backend memory --io-output-backend memory \
  --rank-body-backend memory --compress-level 1
```

该布局仍然是 6 个 MPI rank，每个 rank 的现有算法仍启动 64 个 CPE；另外五个核组
主要提供 cross memory，并未用于算法计算。`-xmalloc` 使完整输入、rank-local body、
模拟输出和大部分 STL workspace 自动进入 cross segment，因此 APP 的处理和计时逻辑
不变。代价是每个 rank 独占一个节点的六个核组，计算资源利用率较低。

2026-09-17 的实测结果：

- `cross_size=24000` 时，现有 BAM2BAM memory input/output 路径完整运行并写出 BAM；
  重新执行 flagstat 得到 2,668,351 条记录，统计与输入一致。
- 最大成功配置为 `share_size=100, cross_size=88000`：裸探针实际逐页提交并由 CPE
  校验 85 GiB/rank，现有 flagstat 和 BAM2BAM memory input/output 路径也运行成功。
- `share_size=100, cross_size=88500/89000` 均在加载器 `space_part` 阶段失败。因此当前
  节点和运行时下，可靠的单一连续 cross heap 上限约为 85 GiB，而不是物理理论值
  96 GiB。
- 增大 `share_size` 不会扩大 APP 可用堆。`share_size` 按核组预留，反而会挤占可分给
  cross segment 的总预算。例如 `share_size=14000, cross_size=60000` 在程序启动前即
  分段失败。若使用 `-xmalloc`，大部分 APP 动态分配都来自 cross，应压低 share 并把
  配额留给 cross。

裸探针通过 85 GiB 不代表 APP 可以处理 85 GiB 输入。memory 模式的峰值至少包含完整
输入、当前算法 workspace、各 rank 压缩 body，以及 rank0 聚合/验证输出。以
BAM-to-BAM 且输出大小约等于输入为例，rank0 的主体峰值约为
`input + final_output + local_body + simulated_local_body`，即约 `2.33 * input`，尚未
计入算法 workspace。纯 BAM-to-BAM 的理论输入上限约 36 GiB，实际建议先用
24--30 GiB 输入逐级测试；sort、collate 和 markdup 还应根据各自 workspace 进一步
降低。若超出 cross segment，现有程序会分配失败；此时应减小输入或改用
POSIX/MPI-IO、spool 等有界后端。

### `swbam-mpi-memory-probe`

用于验证 MPI rank 的节点映射、每个进程的神威运行时内存段配置、实际内存提交能力
以及共享文件系统读取吞吐。

```bash
./swbam-mpi-memory-probe --allocate 1G
./swbam-mpi-memory-probe --allocate 1G \
  --file input.bam --file-mode replicated
./swbam-mpi-memory-probe --allocate 1G \
  --file input.bam --file-mode partitioned
```

- `replicated`：每个 rank 读取完整文件，对应当前 memory input backend 的读取行为。
- `partitioned`：所有 rank 分段读取一次文件，对应 POSIX/MPI-IO 流式 backend 的目标
  行为。
- 文件读取结果可能受页缓存影响，应使用相同节点状态重复测试并报告中位数。
- 默认只逐页提交 512 MiB/rank。容量测试建议逐步增加，例如 4G、8G、12G，避免
  直接触发节点 OOM。

作业手册规定 `-np` 是每节点主核数。因此六节点、每节点一个 MPI rank 应使用：

```bash
bsub.py swls -q q_share -I -b -J memory_probe \
  -N 6 -np 1 -cgsp 64 \
  -share_size 14000 -cache_size 32 \
  -o memory_probe.log \
  ./swbam-mpi-memory-probe@online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc \
  --allocate 1G
```

`-N 6 -np 6` 表示 6 个节点、每节点 6 个 rank，即总计 36 个 rank，并不是每节点
一个 rank。

`share_size` 是每个进程所绑定核组的连续共享段配置；`cache_size` 是从核 LDM 中的
硬件 cache 配置。一个进程预留 1/2/3/6 个核组的大共享资源需要
`-mpecg 1/2/3/6`。`athread_spawn_cgs()` 用于同时启动这些核组的全部 CPE，但扩大
APP 堆内存并不要求把算法改为 384 CPE：实测在 `-mpecg 6 -xmalloc` 下，现有普通
`athread_spawn()` 的 64 CPE kernel 可以直接访问 cross buffer。

探针输出中的 `linux_private_total/available` 来自 Linux `sysinfo` 和
`/proc/meminfo`，在当前神威运行时中只反映约 2 GiB 的进程私有段视图；判断可用的
核组共享内存应查看 `share_size`，并以逐页提交测试为准。

## 2026-09-16 小规模验证

使用 `share_size=14000`、`cache_size=32` 和 6 个 MPI rank 得到：

| 布局 | 节点映射 | 每 rank 提交 | 总提交 | 结果 |
| --- | --- | ---: | ---: | --- |
| `-N 1 -np 6` | 1 节点 × 6 rank | 512 MiB | 3 GiB | 成功 |
| `-N 6 -np 1` | 6 节点 × 1 rank | 4 GiB | 24 GiB | 成功 |

两种布局下每个 rank 都报告 `share_size=14000`、`cross_size=0`、
`cache_size=32`。因此六节点布局不会给单个 rank 提供 84/96 GiB 连续地址空间，也
不会比单节点六 rank 增加这 6 个 rank 的理论总内存；它主要用于减少节点内争用并
利用多节点 I/O。

在同一 0.239 GiB BAM 上的一次小规模读取结果如下。数据受共享文件系统和页缓存
影响，只能作为可行性证据，正式性能结论需多次运行取中位数。

| 读取方式 | `-N 1 -np 6` | `-N 6 -np 1` | 物理读取量 |
| --- | ---: | ---: | ---: |
| replicated | 2.125 s | 2.857 s | 1.435 GiB（6 份输入） |
| partitioned | 1.872 s | 0.386 s | 0.239 GiB（1 份输入） |

这说明多节点对分段 POSIX/MPI-IO 流式输入有明显潜力；当前 memory backend 若仍由
每个 rank 加载完整 BAM，则文件仍须能装入单个 rank 的共享段，并不会因为增加节点
而支持大于单 rank 内存的输入。

### `swbam-mpi-memcpy-probe`

比较 `memcpy`、逐字节循环和 64 位循环在多个 MPI rank 上的内存复制吞吐。参数依次
为每个缓冲区 MiB、迭代次数和预热次数。

```bash
./swbam-mpi-memcpy-probe 256 20 3
```

## 正确性与内核微基准

- `collate_verify.sh`：检查 QNAME 连续性、read1/read2 次序和记录 multiset。
- `markdup_verify_fast`：快速比较统计量、mandatory fields 和 duplicate key 指纹。
- `RabbitBAM-CRC32-Bench` / `bench_crc32.cpp`：比较神威从核 CRC32 实现。

这些文件原位于 `bench/`，现统一归档到 `tools/`；可执行目标名称保持不变。
