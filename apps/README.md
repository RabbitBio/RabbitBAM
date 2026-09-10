# 应用实例

本目录存放基于 SWBAM 公共库构建的可执行程序。这里负责命令语义和算法编排，
不负责实现通用 BGZF I/O、CPE 批调度或 MPI 输出布局。

## 目录

- `mpi/`：当前新版主线，生成 `RabbitBAM-MPI`。

旧版单机和 CGS 实现仍位于项目其他目录，不属于这里的阅读重点。

## 阅读建议

先看 `mpi/main.cpp` 理解命令入口，再按 `commands -> algorithms -> pipelines`
阅读。遇到 `swbam::` 类型时跳转到 `include/swbam/` 和 `lib/`，遇到
`slave_*` 入口时跳转到 `slave/`。
