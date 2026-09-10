# SWBAM SDK 示例

示例只包含业务逻辑，不依赖 `CmdInfo`、RabbitBAM 命令分发或内部 MPI helper。
两者都通过 `swbam/swbam.h` 和 CMake target `swbam::swbam` 使用公共库。

这里是理解新版 SDK 的推荐起点。先看只读示例，再看带 writer 的过滤示例；不建议
先从 `apps/mpi/algorithms/` 的大型实现反推公共接口。

## 构建

```bash
cmake -S . -B build_sunway_sduhpc -DPLATFORM=sunway
cmake --build build_sunway_sduhpc \
  --target swbam-sdk-record-count swbam-sdk-filter-bam -j 8
```

## record count

`sdk_record_count.cpp` 展示 `PosixBamInput + RunGenericRawBamPipeline +
RawBamRecordConsumer`，统计 records、mapped 和 duplicate 数量：

```bash
./swbam-sdk-record-count input.bam
```

## filter BAM

`sdk_filter_bam.cpp` 展示 `RawBamFilterConsumer + RawBamWriter +
PosixBamOutput`，保留 MAPQ 不低于阈值的记录：

```bash
./swbam-sdk-filter-bam input.bam output.bam 20
```

Sunway CPE object 需要直接交给 hybrid linker，不能先封装进普通主核静态库。
该约束已由 `swbam::cpe_kernels` 隐藏，实例本身无需列出从核 object files。

## 示例与生产命令的区别

- 示例使用通用 `RawBamRecordConsumer`，代码短，适合开发新统计和简单过滤。
- `RabbitBAM-MPI` 的热点命令可使用专用融合 kernel，以减少中间记录和主从核搬运。
- 示例仍需显式管理 `athread_init()`/`athread_halt()`；当前 SDK 尚未提供运行时 RAII。
- `RawBamRecordView` 只在当前 consumer 回调期间有效，不能跨 batch 保存指针。
