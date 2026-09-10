# MPI Markdup

本目录实现坐标有序 BAM 的重复标记/删除，是当前状态和边界处理最复杂的实例。

## 文件分工

- `swbam_markdup_mpi.cpp`：命令编排、CPE candidate 提取、owner 判断、稀疏 remote
  exchange、flat hash 去重、bitmap 结果、记录重写、压缩与分布式输出。
- `swbam_markdup_stream_mpi.cpp`：坐标滑动窗口和有界状态的 duplicate 查找实现，由
  主文件的 streaming 路径调用。
- 对应从核实现：`slave/algorithms/markdup_mpi.cpp`。

## 两条路径

- 默认高性能路径：批量收集 candidate 后执行分布式去重，适合内存充足的基线。
- `--stream`：按坐标淘汰不再可能匹配的状态，面向大文件和有界候选内存。

`-r` 删除重复记录，`-c` 先清理已有 duplicate 标记和相关标签。输入必须是经过
fixmate `-m` 的坐标有序 BAM。

## 建议追踪

从 `ProcessMarkdupMPI()` 开始，依次理解 input plan、candidate pass、owner/group、
rewrite pass 和 output；边界 candidate 与 rank 间交换是正确性的重点。
