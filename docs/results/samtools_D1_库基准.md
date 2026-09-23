# D1 库基准结果：Samtools 与 SWBAM

## 1. 测试目的

使用 Samtools 和 SWBAM 对 D1 执行只读扫描及 BAM read-write 回环。Samtools 按正常
方式直接读写 bigssd；SWBAM 使用 memory backend，并只统计排除输入加载和最终 dump
后的内存驻留核心时间。两组结果用于分别观察各平台的扩展性，不能直接作为严格的
跨平台端到端加速比。

## 2. 测试环境

| 项目 | 配置 |
|---|---|
| 日期 | 2026-09-21 |
| 主机 | `hpcl-fat` |
| CPU | 2 x Intel Xeon Platinum 8260 @ 2.40 GHz |
| 物理核 | 48 |
| Samtools | 1.20 |
| Samtools 参数 | `-@ 32` |
| 文件系统 | `/home/bigssd`，`/dev/md0` |
| 输入 | `D1-WES-9G.coordinate.bam` |
| 输入大小 | 9,891,133,304 bytes（9.891 GB） |
| 记录数 | 177,798,571 |
| 排序状态 | coordinate sorted |

输入文件：

```text
/home/bigssd/wzs_data/D1-WES-9G/initial/D1-WES-9G.coordinate.bam
```

测试没有使用 `/dev/shm` 或 `/dev/null`，也没有主动清理 Linux page cache。GNU time
记录的只读 `fs_inputs=0` 表明测试时输入主要命中了系统页缓存，这是本次正常运行环境
的实际状态。

## 3. 测试命令

公共环境：

```bash
export SAMTOOLS=/home/user_home/wzs/software/samtools/bin/samtools
export LD_LIBRARY_PATH=/home/user_home/wzs/software/libdeflate_dist/lib:\
/home/user_home/wzs/software/libdeflate_dist/lib64:${LD_LIBRARY_PATH:-}

IN=/home/bigssd/wzs_data/D1-WES-9G/initial/D1-WES-9G.coordinate.bam
```

### 3.1 只读扫描

```bash
/usr/bin/time \
  -f 'wall_s=%e user_s=%U sys_s=%S cpu=%P maxrss_kb=%M fs_inputs=%I fs_outputs=%O exit=%x' \
  $SAMTOOLS view --no-PG -@ 32 -c "$IN"
```

该命令完成 BGZF 解压、BAM record 读取及计数，不产生 BAM 输出。

### 3.2 Read-write 回环

```bash
OUT=/home/bigssd/wzs_data/samtools_bench/D1_library/D1.roundtrip.bam

/usr/bin/time \
  -f 'wall_s=%e user_s=%U sys_s=%S cpu=%P maxrss_kb=%M fs_inputs=%I fs_outputs=%O exit=%x' \
  $SAMTOOLS view --no-PG -@ 32 \
  -O BAM,level=1 -o "$OUT" "$IN"
```

这里使用 `-O BAM,level=1` 明确指定 BAM 压缩级别 1。Samtools 的 `view -l` 表示
按 read-group library 过滤，不能用于指定压缩级别。

## 4. 固定 `-@ 32` 初测结果

### 4.1 只读

| 轮次 | Wall/s | User/s | System/s | CPU | Max RSS/KiB | 记录数 |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 19.14 | 80.48 | 9.96 | 472% | 21,580 | 177,798,571 |
| 2 | 18.86 | 79.91 | 9.33 | 473% | 20,800 | 177,798,571 |
| 3 | 18.38 | 77.78 | 9.31 | 473% | 20,796 | 177,798,571 |

### 4.2 Read-write

| 轮次 | Wall/s | User/s | System/s | CPU | Max RSS/KiB | 输出大小/bytes |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 40.72 | 448.87 | 33.70 | 1185% | 39,624 | 11,134,791,931 |
| 2 | 39.09 | 431.91 | 32.37 | 1187% | 39,812 | 11,134,791,931 |
| 3 | 47.16 | 511.93 | 34.95 | 1159% | 39,588 | 11,134,791,931 |

## 5. 固定配置汇总

吞吐量采用十进制单位，输入 GB/s 只按输入 BAM 大小计算。

| 测试 | Wall 中位数/s | Records/s | 输入吞吐/GB/s | 输出吞吐/GB/s |
|---|---:|---:|---:|---:|
| Samtools read | 18.86 | 9.43 M | 0.524 | - |
| Samtools read-write | 40.72 | 4.37 M | 0.243 | 0.273 |

read-write 的输入与输出合计数据率约为 `0.516 GB/s`。三次输出大小完全一致，最终
输出通过 `samtools quickcheck`，输入与输出记录数均为 `177,798,571`。

该组为独立短测；线程扩展性结论以后续完整矩阵的三次中位数为准。长时间连续写回会
带来额外系统波动，因此完整矩阵中的 `-@ 32` read-write 中位数与本组略有差异。

## 6. 线程扩展性

### 6.1 方法

补测 `-@ 0、2、4、8、16、24、32、48`。每个线程点运行 3 次，三轮顺序分别为
升序、降序、升序，以减小持续写回和系统负载随时间变化造成的单向偏差。表中时间均为
三次中位数；加速比以 `-@ 0` 为基线。后续另行补测 read-write
`-@ 47`，用于观察 48 物理核主机上避免辅助线程过量配置的效果。

`-@` 表示 Samtools/HTSlib 配置的附加压缩或解压线程数，不应直接解释为程序实际占用
的总 CPU 核数。

### 6.2 只读扩展性

| `-@` | 三次 Wall/s | 中位数/s | 加速比 | Records/s | 输入 GB/s |
|---:|---|---:|---:|---:|---:|
| 0 | 56.72 / 56.41 / 56.80 | 56.72 | 1.00x | 3.13 M | 0.174 |
| 2 | 25.18 / 25.29 / 25.05 | 25.18 | 2.25x | 7.06 M | 0.393 |
| 4 | 17.31 / 17.32 / 17.77 | 17.32 | 3.27x | 10.27 M | 0.571 |
| 8 | 17.29 / 17.16 / 16.49 | **17.16** | **3.31x** | **10.36 M** | **0.576** |
| 16 | 17.84 / 17.22 / 17.73 | 17.73 | 3.20x | 10.03 M | 0.558 |
| 24 | 17.42 / 17.33 / 17.82 | 17.42 | 3.26x | 10.21 M | 0.568 |
| 32 | 18.23 / 17.76 / 18.91 | 18.23 | 3.11x | 9.75 M | 0.543 |
| 48 | 17.24 / 17.83 / 17.75 | 17.75 | 3.20x | 10.02 M | 0.557 |

只读在 `-@ 4～8` 后基本饱和。最佳中位数为 `-@ 8` 的 `17.16 s`，相比串行配置
加速 `3.31x`；继续增加线程没有稳定收益。

### 6.3 Read-write 扩展性

| `-@` | 三次 Wall/s | 中位数/s | 加速比 | Records/s | 输入 GB/s | 输入+输出 GB/s |
|---:|---|---:|---:|---:|---:|---:|
| 0 | 293.53 / 293.66 / 294.73 | 293.66 | 1.00x | 0.61 M | 0.034 | 0.072 |
| 2 | 154.62 / 165.72 / 155.14 | 155.14 | 1.89x | 1.15 M | 0.064 | 0.136 |
| 4 | 83.29 / 118.74 / 82.70 | 83.29 | 3.53x | 2.13 M | 0.119 | 0.252 |
| 8 | 59.28 / 67.27 / 53.58 | 59.28 | 4.95x | 3.00 M | 0.167 | 0.355 |
| 16 | 49.94 / 55.55 / 48.14 | 49.94 | 5.88x | 3.56 M | 0.198 | 0.421 |
| 24 | 47.28 / 51.05 / 48.25 | 48.25 | 6.09x | 3.68 M | 0.205 | 0.436 |
| 32 | 43.10 / 46.26 / 45.19 | 45.19 | 6.50x | 3.93 M | 0.219 | 0.465 |
| 47 | 39.41 / 37.49 / 38.79 | **38.79** | **7.57x** | **4.58 M** | **0.255** | **0.542** |
| 48 | 42.94 / 44.12 / 43.53 | 43.53 | 6.75x | 4.08 M | 0.227 | 0.483 |

read-write 随线程数增加持续加速，但从 `-@ 16` 开始收益递减。补测的 `-@ 47`
中位数为 `38.79 s`，相对串行配置加速 `7.57x`，也比 `-@ 48` 快约
10.9%。由于 `-@` 指定的是附加 HTSlib 工作线程，`-@ 47` 更接近这台 48 物理核
主机的满核配置；不过该点是后续独立补测，正式论文实验仍应用统一轮次顺序复测。
连续向
bigssd 写出约 11.13 GB 时存在明显环境波动，例如 `-@ 4` 第二轮为 `118.74 s`，
因此正式结果采用中位数而不是最小值。

所有只读测试均得到 `177,798,571` 条记录；所有 read-write 输出大小均为
`11,134,791,931` bytes，并逐次通过 `samtools quickcheck`。

## 7. 原始日志

fat 上保留的原始计时文件：

```text
/home/bigssd/wzs_data/samtools_bench/D1_library/raw_read_20260921_153710
/home/bigssd/wzs_data/samtools_bench/D1_library/raw_rw_20260921_153944
/home/bigssd/wzs_data/samtools_bench/D1_scaling/raw_20260921_154839
/home/bigssd/wzs_data/samtools_bench/D1_D2_final/raw_20260921_165619
```

全部临时回环 BAM 在完成 `quickcheck` 和记录数核对后已删除。

## 8. SWBAM 测试环境

| 项目 | 配置 |
|---|---|
| 日期 | 2026-09-22 |
| 平台 | `swls` 神威集群，`q_share` 队列 |
| 程序 | `RabbitBAM-MPI`，优化路径 |
| 节点数 | `1、2、4、8、16、24、32、48` |
| MPI 布局 | 每节点 1 rank；每 rank 使用 1 MPE + 64 CPE 计算 |
| 大共享配置 | `-mpecg 6 -share_size 100 -cross_size 88000 -xmalloc` |
| 输入 backend | `memory`，每个 rank 将完整输入加载到 cross memory |
| 输入 | `D1-WES-9G.coordinate.bam` |
| 输入大小 | 9,891,133,304 bytes（9.891 GB） |
| 记录数 | 177,798,571 |
| 读写输出大小 | 11,467,682,880 bytes（11.468 GB） |

大共享配置为每个 rank 预留一个节点的 6 个核组及约 85 GiB 可可靠使用的 cross heap，
其中算法仍只启动 64 个 CPE；其余 5 个核组主要提供内存容量。因此，下面的节点数等于
MPI rank 数，但不等价于 Samtools 的 `-@` 参数，也不表示每个节点的 384 个 CPE 都参与
计算。例如 48 节点配置实际执行 48 个 MPE 和 3,072 个 CPE。

## 9. SWBAM 测试命令

公共资源配置：

```bash
bsub.py swls -q q_share -I -b -J JOB_NAME \
  -N ${N} -np 1 -mpecg 6 -cgsp 64 \
  -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc \
  -o OUTPUT_LOG \
  -x LD_LIBRARY_PATH=../ext/htslib-1.20:../ext/libdeflate-1.20/build \
  ./RabbitBAM-MPI@online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc \
  COMMAND
```

只读使用融合 CPE flagstat 路径：

```bash
RabbitBAM-MPI flagstat \
  -i D1-WES-9G.coordinate.bam \
  --io-backend memory
```

Read-write 使用融合 BAM-to-BAM 路径，压缩级别为 1：

```bash
RabbitBAM-MPI run_all \
  -i D1-WES-9G.coordinate.bam \
  -o D1.roundtrip.bam \
  --io-backend memory \
  --io-output-backend memory \
  --rank-body-backend memory \
  --compress-level 1
```

## 10. SWBAM 计时口径

每个节点数运行 3 次，顺序为升序、降序、升序，结果取中位数。加速比以 1 节点核心
时间为基线，并行效率为 `speedup / N`。

- 只读主指标取 `444Complete the total body cost`，覆盖 4.1～4.4，不包含
  `222Complete the input open cost`。
- Read-write 主指标取 `Complete the total (4.1~4.6) cost`，不包含输入加载、4.5b
  模拟输出内存分配、4.7 统计、验证 gather 和最终文件 dump。
- `4.3` 单独列出，用于观察融合 CPE 算子的扩展性；主结果仍使用包含扫描、布局和模拟
  内存写出的核心总时间。
- 吞吐量使用十进制单位。Read-write 的合计吞吐按输入 9.891 GB 与输出 11.468 GB
  之和计算。

测试期间有一次 8 节点作业停在输入加载阶段，换节点后恢复正常；32 节点第三轮首次
提交因计算节点不足失败，资源可用后重试成功。失败尝试均未进入统计。

## 11. SWBAM 只读扩展性

| 节点/rank | 三次核心时间/s | 中位数/s | `4.3` 中位数/s | 加速比 | 并行效率 | Records/s | 输入 GB/s |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | 15.494172 / 15.494665 / 15.499825 | 15.494665 | 15.258401 | 1.00x | 100.0% | 11.47 M | 0.638 |
| 2 | 7.894712 / 7.884782 / 7.883362 | 7.884782 | 7.642904 | 1.97x | 98.3% | 22.55 M | 1.254 |
| 4 | 4.073929 / 4.074237 / 4.079054 | 4.074237 | 3.832719 | 3.80x | 95.1% | 43.64 M | 2.428 |
| 8 | 2.161533 / 2.163953 / 2.162430 | 2.162430 | 1.920965 | 7.17x | 89.6% | 82.22 M | 4.574 |
| 16 | 1.208032 / 1.207071 / 1.206850 | 1.207071 | 0.964886 | 12.84x | 80.2% | 147.30 M | 8.194 |
| 24 | 0.897962 / 0.889834 / 0.889033 | 0.889834 | 0.645927 | 17.41x | 72.6% | 199.81 M | 11.116 |
| 32 | 0.729295 / 0.731101 / 0.732116 | 0.731101 | 0.486637 | 21.19x | 66.2% | 243.19 M | 13.529 |
| 48 | 0.569906 / 0.576956 / 0.569457 | **0.569906** | **0.325972** | **27.19x** | **56.6%** | **311.98 M** | **17.356** |

只读融合算子从 1 节点的 `15.258 s` 降到 48 节点的 `0.326 s`，接近按 rank 数划分
工作量。与此同时，4.1 block 扫描、广播和分区约保持在 `0.24 s`，在 48 节点核心总
时间中已占约 42.5%，成为并行效率从 98.3% 下降到 56.6% 的主要原因。

## 12. SWBAM Read-write 扩展性

| 节点/rank | 三次核心时间/s | 中位数/s | `4.3` 中位数/s | 加速比 | 并行效率 | Records/s | 输入 GB/s | 输入+输出 GB/s |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 53.989948 / 53.978821 / 53.995399 | 53.989948 | 50.801560 | 1.00x | 100.0% | 3.29 M | 0.183 | 0.396 |
| 2 | 27.240919 / 27.233061 / 27.243198 | 27.240919 | 25.512551 | 1.98x | 99.1% | 6.53 M | 0.363 | 0.784 |
| 4 | 13.759391 / 13.760879 / 13.755252 | 13.759391 | 12.775471 | 3.92x | 98.1% | 12.92 M | 0.719 | 1.552 |
| 8 | 7.002306 / 7.012681 / 7.010257 | 7.010257 | 6.399941 | 7.70x | 96.3% | 25.36 M | 1.411 | 3.047 |
| 16 | 3.643280 / 3.642519 / 3.644592 | 3.643280 | 3.215383 | 14.82x | 92.6% | 48.80 M | 2.715 | 5.863 |
| 24 | 2.520041 / 2.521441 / 2.522046 | 2.521441 | 2.153714 | 21.41x | 89.2% | 70.51 M | 3.923 | 8.471 |
| 32 | 1.957067 / 1.954816 / 1.959419 | 1.957067 | 1.617354 | 27.59x | 86.2% | 90.85 M | 5.054 | 10.914 |
| 48 | 1.397325 / 1.397511 / 1.397086 | **1.397325** | **1.086101** | **38.64x** | **80.5%** | **127.24 M** | **7.079** | **15.286** |

Read-write 的核心计算量更大，固定开销占比更低，因此扩展性优于只读。48 节点时，
`4.3` 融合解压、解析、pack 和压缩约为 `1.086 s`，4.1 扫描/广播约为 `0.243 s`，
模拟分布式内存写约为 `0.067 s`；总体获得 `38.64x` 加速和 80.5% 并行效率。

## 13. SWBAM 正确性与结果说明

全部 24 次只读运行均统计到 `177,798,571` 条记录；全部 24 次 read-write 运行均报告：

```text
total_records=177798571
kept_records=177798571
dropped_records=0
```

保留 1 节点和 48 节点第一轮输出并分别运行 `stats --basic`。两个输出文件大小均为
`11,467,682,880 bytes`，全部 `SN` 统计逐行一致，关键结果如下：

```text
raw total sequences:       176997814
supplementary alignments:     800757
is sorted:                         1
total records including supplementary: 177798571
```

因此，1～48 节点的并行划分没有改变记录数量或基础统计语义。原始日志和两份验证输出
保存在：

```text
/home/export/base/shisuan/swls-CFD/online/guoshi/wzs/wzs_data/
  swbam_bench/D1_library/swbam_D1_library_20260922_103731/
```

需要特别注意，Samtools 表格记录的是 x86 主机从 bigssd 读写的 wall time，而 SWBAM
表格记录的是神威 memory backend 的核心时间，并且 48 节点远多于一台 x86 主机的资源。
因此当前结果可以分别证明两套程序的线程/节点扩展趋势，但不应直接用二者的最短时间
相除并宣称跨平台端到端加速比。正式论文若需要该加速比，应进一步统一 I/O 介质、计时
边界和计算资源口径。

## 14. Composable SDK 初步性能探测

为评估公共 SDK 路径的当前性能，新增 `swbam-sdk-mpi-benchmark` 示例。该示例只调用
公开的 `swbam::swbam` API，将 MPI input plan 分配给各 rank，并在一次 memory backend
加载后依次运行 Raw read、Raw read-write、`bam1_t` read 和 `bam1_t` read-write。
Read-write 只写入 rank-local `MemoryBamOutput`，不执行最终磁盘 dump。

本轮先完成 1 节点、1 rank、64 CPE 的单次探测：

| SDK 路径 | 核心时间/s | Records/s | 输入 GB/s | 输入+输出 GB/s |
|---|---:|---:|---:|---:|
| Raw read | 62.831 | 2.830 M | 0.157 | - |
| Raw read-write | 134.336 | 1.324 M | 0.074 | 0.159 |
| `bam1_t` read | 158.383 | 1.123 M | 0.062 | - |
| `bam1_t` read-write | 445.479 | 0.399 M | 0.022 | 0.048 |

四条路径均处理 `177,798,571` 条记录；两条 read-write 路径均产生
`11,467,688,536 bytes` 的 rank-local 内存输出。关键时间构成为：

| 路径 | 主要分项 |
|---|---|
| Raw read | CPE `13.282 s`；MPE pipeline post-processing `49.472 s` |
| Raw read-write | pack `26.219 s`；compress `29.048 s`；writer 总计 `56.604 s` |
| `bam1_t` read | CPE `13.285 s`；`bam1_t` materialize `101.066 s`；业务回调 `5.832 s` |
| `bam1_t` read-write | materialize `123.890 s`；encode `186.718 s`；pack `26.705 s`；compress `29.044 s` |

结果表明，Composable 路径的主要瓶颈不是 CPE BGZF 解压，而是 MPE 侧的 Raw record
view 构造/遍历、`bam1_t` payload 物化和反向编码。相对同一数据上的 Optimized path，
Raw read 和 Raw read-write 分别慢约 `4.06x` 和 `2.49x`；`bam1_t` read 和 read-write
分别慢约 `10.22x` 和 `8.25x`。

若仅与本文件中的 Samtools 最优 wall time 作数值参考，四条 SDK 路径分别慢约
`3.66x、3.46x、9.23x、11.48x`。但这里的 Samtools 使用 48 核 x86 配置，而 SWBAM
仅使用 1 rank/64 CPE，并且二者的 I/O 计时边界不同，因此这些比值不是正式跨平台
加速比。即便采用对 SWBAM 更有利的纯核心时间，当前单 rank SDK 仍明显偏慢，已经足以
说明需要优化 Composable 路径。

曾尝试直接运行 48 节点上界探测，但 48 个 rank 同时复制读取完整 9.9 GB 输入后，作业
始终未越过 memory input open，最终收到 `bkill`，没有产生任何核心计时。该失败反映的
是 replicated memory backend 对共享文件系统的启动压力，不是四条 SDK 核心路径的
执行失败；因此未继续提交扩展性作业，也未将其计入结果。

下一步应先优化后再重测：Raw 路径优先减少 MPE record-view 扫描成本；`bam1_t` 路径
优先减少对象池重置和 payload 复制，并考虑让 CPE 生成 mandatory-field descriptor。
对于未修改记录的 read-write，还可保留原 Raw view 并直接透传，避免无意义的
`bam1_t -> raw BAM` 全量重编码。完整 `bam1_t` 对象所有权和动态内存管理仍应留在 MPE，
不建议整体搬入容量受限的 CPE LDM。
