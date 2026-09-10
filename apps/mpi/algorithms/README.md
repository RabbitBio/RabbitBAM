# MPI 全局算法实例

本目录存放需要全局数据关系、MPI exchange 或复杂状态的数据处理算法。它们作为
SWBAM 库的高性能应用案例，而不是通用 I/O 库本身。

## 子目录

- `sort/`：坐标排序、采样分区、交换、内排和外排。
- `collate/`：QNAME hash 分桶、交换、组内排序和外排。
- `markdup/`：候选提取、跨 rank 去重判断、重写和流式窗口。

这些实现复用 `BamInputBackend`、`BgzfBlockBatch`、`RankBodySink` 和
`DistributedBamOutput`，同时通过 `slave/algorithms/` 中的专用 kernel 保留性能。

阅读时先找每个文件末尾的 `Process*MPI()`，从阶段编排反向追踪内部函数；不要从
顶部所有结构体开始逐行阅读。
