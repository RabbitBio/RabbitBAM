# MPI Sort

本目录实现坐标排序，是“公共库 + 全局重排算法 + 专用 CPE kernel”的完整案例。

## 文件分工

- `swbam_sort_mpi.cpp`：命令级编排、输入 backend、内存策略、header 和分布式输出；
  入口为 `ProcessSortMPI()`，pipeline 还会使用 memory-to-memory helper。
- `swbam_sort_fused_mpi.cpp`：排序核心，包括 CPE 元数据提取、local sort、采样选取
  splitter、MPI bucket exchange、k-way merge、压缩，以及外排 run/segment。
- 对应从核实现：`slave/algorithms/sort_mpi.cpp`。

## 两类临时存储

- 算法外排：工作集超过 `-m` 后保存有序 run，由 `-T` 指定位置。
- rank body spool：输出 body 超过内存策略后暂存压缩结果，由
  `--rank-body-*` 控制。

两者用途不同，阅读和统计性能时不要混为一类 I/O。

## 建议追踪

```text
ProcessSortMPI
-> extract keys/raw
-> local sort + sample/splitter
-> exchange buckets
-> merge
-> compress to RankBodySink
-> DistributedBamOutput
```
