# SWBAM 库实现

本目录实现 `include/swbam/` 声明的主核侧公共能力。这里不放命令行解析，也不放
sort/collate/markdup 的完整算法语义。

## 子目录

- `io/`：`swbam_io`，输入输出和 rank body 存储。
- `mpi/`：`swbam_mpi_runtime`，MPI backend、分区和分布式输出。
- `cpe/`：`swbam_cpe_runtime` 的 MPE 侧流水线和通用算子包装。
- `bam/`：`swbam_cpe_runtime` 中的 Raw BAM adapter、filter 和 writer。

## 设计边界

```text
backend -> BGZF batch -> CPE pipeline -> raw/operator -> output backend
```

库的扩展调用发生在 batch 或 consumer 层。性能敏感命令可以保留专用融合 kernel，
但应继续复用这里的 backend、批缓冲、调度和输出接口。

阅读顺序建议为 `io -> cpe read -> bam raw -> cpe write -> mpi`。
