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

## 公共流水线接入

```text
MemReader / BamInputBackend + spans
  -> RunCpeReadPipeline + CollateExtractOperator
  -> app 持久 raw/meta 或外排 run arena
  -> app bin/hash 排序、MPI exchange、loser tree 归并
  -> CollatePayloadSource
  -> RunCpeWritePipeline + CollateCompressOperator
  -> RankBodySink -> 分布式输出
```

- `CollateMemoryExtractOperator` 让 CPE 直接解压到 app 的持久 raw arena，
  批后处理只追加 metadata；`CollateRunExtractOperator` 把当前 decoded batch
  保存到有界 run arena，满后按原规则 spill。二者共用 kernel 参数绑定和校验。
- 读流水线管理 compressed/decoded 双缓冲；本轮显式关闭输入覆盖，保留原有
  `extract_read` 与 CPE 时间不重叠的口径。算法持久数据不借用会被复用的 batch。
- 写流水线管理两套 payload 和两套压缩输出。开启 `overlap_source`，压缩当前批时
  MPE flush 上一批、归并打包下一批。内排仍直接遍历缓存 cursor 的 loser tree，
  不增加逐记录虚调用；外排保留一条记录的 lookahead，避免跨 payload 边界丢记录。
- hash/bin、owner、self/remote exchange、raw arena、consolidation、header/EOF
  均未改动。专用 key 提取和压缩 kernel 保留，本轮迁移的是批次调度两端。
- `MpiCollateMemoryToMemory` 间接复用同一核心；不改变 dedup-pipeline 的衔接方式。

读 workspace 比旧单批调度多约 8 MiB/rank，已计入外排工作区预算；实际 run arena
按读、交换、压缩三个阶段的最大 workspace 计算。`tracked_peak` 是算法工作集估计，
不是进程 RSS 或完整内存 I/O 总占用。

小数据结果见 [collate 流水线迁移验证](../../../../docs/results/collate_流水线迁移验证.md)。

## 阅读建议

先看 `ProcessCollateMPI()` 如何选择 memory/external，再分别追踪
`OptimizedBamMemoryCollateMPI` 和 `OptimizedBamExternalCollateMPI`。最后再看 exchange 与
loser tree，避免同时展开两套模式。

验证重点是 record multiset、QNAME continuity 和 read1/read2 顺序，不要求输出与
SAMtools 具有完全相同的 group 排列。
