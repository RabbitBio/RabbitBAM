# CPE 算法专用 Kernel

本目录包含复杂 MPI 算法的从核热点，不作为普通 SDK 示例默认链接内容。它们通过
`swbam::algorithm_kernels` 加入 `RabbitBAM-MPI`。

## 文件

- `sort_mpi.cpp`：排序 key/raw 提取、payload 压缩和 bucket/range pack。
- `collate_mpi.cpp`：解压、QNAME/hash/meta 提取和 codec cache 管理。
- `fixmate_mpi.cpp`：mate plan、字段/tag 长度规划和 record rewrite。
- `markdup_mpi.cpp`：pair/single candidate 提取、key 构造和 duplicate rewrite。

MPE 侧对应实现位于 `apps/mpi/algorithms/` 或 `apps/mpi/commands/`。阅读时应把参数
结构定义、MPE `athread_spawn/join` 调用和这里的 `extern "C" slave_*` 入口对照起来。

这些 kernel 可以使用 `slave/core/`，但算法状态、MPI owner 和全局归并仍由 MPE
应用层负责。
