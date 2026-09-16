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
           swbam-sdk-filter-bam swbam-sdk-bam1-filter -j 8
```

## record count

`sdk_record_count.cpp` 展示 `PosixBamInput + RunGenericRawBamPipeline +
RawBamRecordConsumer`，统计 records、mapped 和 duplicate 数量：

```bash
./swbam-sdk-record-count input.bam
```

使用内存 backend 并打印不含文件加载的核心时间：

```bash
./swbam-sdk-record-count input.bam --memory-io
```

## bam1_t stats

`sdk_bam1_stats.cpp` 展示 `RunGenericBam1Pipeline + Bam1RecordConsumer`。记录由
SDK 按 batch 物化为标准 `bam1_t`，示例通过 HTSlib accessor 访问 mandatory fields、
QNAME/CIGAR/SEQ/QUAL 和 `NM` aux tag：

```bash
./swbam-sdk-bam1-stats input.bam
./swbam-sdk-bam1-stats input.bam --memory-io
```

输出中的 `timing_materialize` 是 Raw record 转换为 `bam1_t` 的额外成本；它与
`timing_bam1_consume` 分开，便于比较易用路径和零拷贝 Raw 路径。consumer 不能保存
回调中的记录指针，确需保存时使用 `bam_dup1()`。

## bam1_t filter + writer

`sdk_bam1_filter.cpp` 展示完整的 HTSlib 对象路径：通用读流水线批量物化
`bam1_t`，consumer 通过 `record->core.qual` 过滤，再由 `Bam1Writer` 批量编码并
复用 `RawBamWriter + CPE BGZF compress` 写出 BAM：

```bash
./swbam-sdk-bam1-filter input.bam output.bam 30
./swbam-sdk-bam1-filter input.bam output.bam 30 --memory-io
```

输出分别列出 materialize、filter、encode、pack 和 compress 时间。该示例适合需要
HTSlib accessor 或修改 `bam1_t` 的业务；只访问少量 mandatory fields 时，下面的
Raw filter 路径更轻。

## filter BAM

`sdk_filter_bam.cpp` 展示 `RawBamFilterConsumer + RawBamWriter +
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

- 示例使用通用 `RawBamRecordConsumer`，代码短，适合开发新统计和简单过滤。
- 需要直接调用 HTSlib record accessor 时，使用 `Bam1RecordConsumer`；它比 Raw
  路径易用，但需要物化并复制 record payload。
- `RabbitBAM-MPI` 的热点命令可使用专用融合 kernel，以减少中间记录和主从核搬运。
- 示例仍需显式管理 `athread_init()`/`athread_halt()`；当前 SDK 尚未提供运行时 RAII。
- `RawBamRecordView` 只在当前 consumer 回调期间有效，不能跨 batch 保存指针。
- `Bam1RecordConsumer` 收到的 `bam1_t` 同样是 batch-local；需要保留时调用
  `bam_dup1()`，并由调用者负责 `bam_destroy1()`。
