# MPI Collate

本目录实现按 QNAME 聚集记录，但不保证不同 QNAME group 之间的全局字典序。

## 主要逻辑

`swbam_collate_mpi.cpp` 同时包含：

- CPE 解压和 QNAME/hash/meta 提取；
- 按逻辑 bin 确定 owner rank；
- local sort、remote exchange 和 self-owned 快路径；
- memory 模式的 raw arena、segment 和 loser tree merge；
- external 模式的 run、临时 extent、consolidation 和文件 cursor；
- 压缩、`RankBodySink` 与分布式输出；
- `ProcessCollateMPI()` 和 pipeline memory helper。

对应从核实现位于 `slave/algorithms/collate_mpi.cpp`。

## 阅读建议

先看 `ProcessCollateMPI()` 如何选择 memory/external，再分别追踪
`FusedBamMemoryCollateMPI` 和 `FusedBamExternalCollateMPI`。最后再看 exchange 与
loser tree，避免同时展开两套模式。

验证重点是 record multiset、QNAME continuity 和 read1/read2 顺序，不要求输出与
SAMtools 具有完全相同的 group 排列。
