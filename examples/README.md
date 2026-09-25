# SWBAM SDK 示例

示例只包含业务逻辑，不依赖 `CmdInfo`、RabbitBAM 命令分发或内部 MPI helper。
这些示例都通过 `swbam/swbam.h` 和 CMake target `swbam::swbam` 使用公共库。

这里是理解新版 SDK 的推荐起点。先看只读示例，再看带 writer 的过滤示例；不建议
先从 `apps/mpi/algorithms/` 的大型实现反推公共接口。

## 构建

```bash
cmake -S . -B build_sunway_sduhpc -DPLATFORM=sunway
cmake --build build_sunway_sduhpc \
  --target swbam-sdk-record-count swbam-sdk-bam1-stats \
           swbam-sdk-filter-bam swbam-sdk-bam1-filter \
           swbam-sdk-mpi-benchmark swbam-sdk-optimized-count \
           swbam-sdk-optimized-bam2bam -j 8
```

## record count

`sdk_record_count.cpp` 展示 `PosixBamInput + RunComposableRawBamPipeline +
RawBamBatchPostProcessor`，统计 records、mapped 和 duplicate 数量：

```bash
./swbam-sdk-record-count input.bam
```

使用内存 backend 并打印不含文件加载的核心时间：

```bash
./swbam-sdk-record-count input.bam --memory-io
```

## bam1_t stats

`sdk_bam1_stats.cpp` 展示 `RunComposableBam1Pipeline + Bam1BatchPostProcessor`。记录由
SDK 按 batch 物化为标准 `bam1_t`，示例通过 HTSlib accessor 访问 mandatory fields、
QNAME/CIGAR/SEQ/QUAL 和 `NM` aux tag：

```bash
./swbam-sdk-bam1-stats input.bam
./swbam-sdk-bam1-stats input.bam --memory-io
```

输出中的 `timing_materialize` 是 Raw record 转换为 `bam1_t` 的额外成本；它与
`timing_bam1_post_process` 分开，便于比较易用路径和零拷贝 Raw 路径。batch
post-processor 不能保存
回调中的记录指针，确需保存时使用 `bam_dup1()`。

## bam1_t filter + writer

`sdk_bam1_filter.cpp` 展示完整的 HTSlib 对象路径：Composable read pipeline 批量物化
`bam1_t`，batch post-processor 通过 `record->core.qual` 过滤，再由 `Bam1Writer`
批量编码并
复用 `RawBamWriter + CPE BGZF compress` 写出 BAM：

```bash
./swbam-sdk-bam1-filter input.bam output.bam 30
./swbam-sdk-bam1-filter input.bam output.bam 30 --memory-io
```

输出分别列出 materialize、filter、encode、pack 和 compress 时间。该示例适合需要
HTSlib accessor 或修改 `bam1_t` 的业务；只访问少量 mandatory fields 时，下面的
Raw filter 路径更轻。

## filter BAM

`sdk_filter_bam.cpp` 展示 `RawBamFilterBatchPostProcessor + RawBamWriter +
PosixBamOutput`，保留 MAPQ 不低于阈值的记录：

```bash
./swbam-sdk-filter-bam input.bam output.bam 20
```

内存基准模式会先把输入加载到内存、将压缩输出保存在内存中，`timing_core`
不包含初始加载和最终 dump：

```bash
./swbam-sdk-filter-bam input.bam output.bam 20 --memory-io
```

Sunway CPE object 需要直接交给 hybrid linker，不能先封装进普通主核静态库。
该约束已由 `swbam::cpe_kernels` 隐藏，实例本身无需列出从核 object files。

## 示例与生产命令的区别

- 示例使用 Composable `RawBamBatchPostProcessor`，代码短，适合开发新统计和简单过滤。
- 需要直接调用 HTSlib record accessor 时，使用 `Bam1BatchPostProcessor`；它比 Raw
  路径易用，但需要物化并复制 record payload。
- `RabbitBAM-MPI` 的热点命令可使用专用融合 kernel，以减少中间记录和主从核搬运。
- 示例仍需显式管理 `athread_init()`/`athread_halt()`；当前 SDK 尚未提供运行时 RAII。
- `RawBamRecordView` 只在当前 batch post-processor 回调期间有效，不能跨 batch 保存指针。
- `Bam1BatchPostProcessor` 收到的 `bam1_t` 同样是 batch-local；需要保留时调用
  `bam_dup1()`，并由调用者负责 `bam_destroy1()`。

## MPI 扩展性微基准

`sdk_mpi_benchmark.cpp` 只使用公开 SDK 接口，将 MPI input plan 分配给各 rank，并在
一次 memory backend 加载后依次测试 Raw read、Raw read-write、`bam1_t` read 和
`bam1_t` read-write：

```bash
./swbam-sdk-mpi-benchmark input.bam
```

四条路径的 `core` 都不包含输入文件加载、block plan 和磁盘输出。Read-write 只写入
rank-local `MemoryBamOutput`，不 dump 文件；该工具用于比较 Composable SDK 路径和
`bam1_t` 物化开销，不代替生成单个分布式 BAM 文件的正式应用。

## Optimized 路径核心基准

下面两个独立示例使用与 MPI app 相同的 memory input backend、block plan 和
CPE batch 流水线；`core` 是各 rank 核心处理时间的最大值，不含初始加载、扫描
分区和最终磁盘 dump：

```bash
./swbam-sdk-optimized-count input.bam
./swbam-sdk-optimized-bam2bam input.bam
```

`sdk_optimized_count.cpp` 在 CPE 解压 BGZF 后逐条解析到复用的临时
`bam1_t` 再计数，不计算 flagstat 指标；可与 `samtools view -c` 的记录数核对。
它采用项目现有的 record 不跨 BGZF block 约束。此前只扫描 BAM record 长度
的 `core=0.047015 s` 是不同工作量的历史基线，不可作为此版本的性能结果。

`sdk_optimized_bam2bam.cpp` 复用 `RunBamTransformPipeline`，不设置过滤条件，
用 level 1 重新压缩，将各 rank 的 BGZF body 保存在内存中。它报告记录数和
压缩 body 总字节数；不拼接 header/EOF，也不生成最终输出文件，因此是库的
读写核心基准，不是完整的 BAM 转换命令。

WES_0.25G、6 ranks 的历史 smoke：count 得到 `2,668,351` 条；
bam2bam 同样处理 `2,668,351` 条，
`body_bytes=275,206,630`、两次 `core=0.196671/0.196398 s`。这些数字仅说明功能与
测试口径正常，正式性能比较仍需同配置重复运行取中位数。

改用 `bam1_t` 解析后，同一输入的 count 两次 smoke 均为 `2,668,351` 条，
`core=0.050270/0.050454 s`。与旧的只扫长度路径相比，工作量已不同；也不能把
这个核心时间与包含磁盘 I/O 的 Samtools wall time 直接等同比较。
