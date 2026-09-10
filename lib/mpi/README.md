# MPI 运行时

本目录构成 `swbam_mpi_runtime` 静态库，为公共 I/O 数据面增加 MPI 语义。

## 文件

- `swbam_mpi_runtime.cpp`：实现 `MpiBamInput` backend factory、MPI-IO BAM/SAM
  reader、全局 BGZF block plan、rank 分区、错误/计时规约和分布式输出。

## 输入路径

```text
rank0 ScanBlocks
-> broadcast spans
-> 按 block 数切分 rank_begin/rank_end
-> 每 rank 通过同一 backend ReadBatch
```

`MpiBamInput` 支持 `memory|posix|mpiio|auto`。`auto` 可因内存预算回退 POSIX，
MPI-IO 目前需要显式选择。

## 输出路径

`DistributedBamOutput` 收集每 rank 的压缩 body size，计算 header/body/EOF 的最终
offset，再选择 gather 到 rank0 或 MPI-IO 直写。它只负责字节布局，不修改 SAM
header，也不理解 sort/markdup 等算法语义。
