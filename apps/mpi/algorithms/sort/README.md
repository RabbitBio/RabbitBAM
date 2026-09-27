# MPI Sort

本目录实现坐标排序，是“公共库 + 全局重排算法 + 专用 CPE kernel”的完整案例。

## 文件分工

- `swbam_sort_mpi.cpp`：命令级编排、输入 backend、内存策略、header 和分布式输出；
  入口为 `ProcessSortMPI()`，pipeline 还会使用 memory-to-memory helper。
- `swbam_sort_optimized_mpi.cpp`：排序核心及公共流水线适配，包括元数据提取、local
  sort、采样选取 splitter、MPI bucket exchange、k-way merge、外排 run/segment。
- 对应从核实现：`slave/algorithms/sort_mpi.cpp`。

## 两类临时存储

- 算法外排：工作集超过 `-m` 后保存有序 run，由 `-T` 指定位置。
- rank body spool：输出 body 超过内存策略后暂存压缩结果，由
  `--rank-body-*` 控制。

两者用途不同，阅读和统计性能时不要混为一类 I/O。

## 建议追踪

```text
ProcessSortMPI
-> RunCpeReadPipeline + MpiSortExtractOperator
-> app-owned raw/metadata or bounded run arena
-> local sort + sample/splitter
-> exchange buckets
-> merge
-> RunCpeWritePipeline + MpiSortCompressOperator
-> RankBodySink
-> DistributedBamOutput
```

## 读写两端如何复用库

内排和严格外排的读取都通过公共 BAM read 流水线。`MpiSortExtractOperator` 继承
`BamReadBatchOperator`，沿用原 `slave_mpi_sort_extract_raw`；内排的批后处理把 raw
和 metadata 追加到 app 存储，外排的批后处理负责追加 run arena、满后 spill。
raw view 不跨回调保存，稳定排序使用的全局 block/record 标识保持原语义。

写端把按最终顺序归并、打包的逻辑实现为 `UncompressedBgzfSource`。公共
`RunCpeWritePipeline` 使用 `MpiSortCompressOperator` 适配原压缩 kernel，并将
结果交给 `MpiSortCompressedOutput` 追加到 `RankBodySink`。
`overlap_source=true` 使 CPE 压缩当前批时，MPE flush 上一批并归并下一批；
两套 payload 与两套压缩输出缓冲互不覆盖。无需把 raw records 转成 `bam1_t`。

读写调度、批缓冲归库；键提取适配、排序、MPI 交换、run 文件和归并状态仍在 app。
memory-to-memory helper 通过相同核心复用以上路径。

严格外排的 `-m` 预算计入公共读流水线的四套 block batch 及 metadata，较原两套
block workspace 增加约 8 MiB/rank；run arena 会按可用预算缩减，不绕过内存限制。
外排仍区分真实临时 I/O 与模拟计时，合并/压缩的重叠时间不能直接相加。

验证命令与结果见 [sort 流水线迁移验证](../../../../docs/results/sort_流水线迁移验证.md)。
