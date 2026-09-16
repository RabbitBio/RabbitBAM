# MPI 组合流水线

本目录只负责把已经存在的阶段组合起来，不重新实现 collate、fixmate、sort 或
markdup 算法。

## 当前文件

- `swbam_dedup_pipeline_mpi.cpp`：实现
  `collate -> fixmate -m -> coordinate sort -> markdup`。
- `swbam_dedup_workflow_mpi.cpp`：直接顺序调用四个 standalone MPI 命令，以三个
  中间 BAM 连接阶段；用于大文件、逐阶段验证和端到端基线。

当前 pipeline 使用完整内存 BAM 作为阶段间格式：rank0 组装阶段输出并广播给下一
阶段。这样可以复用已验证的 memory-to-memory helper，但峰值内存和调试复杂度高于
依次执行独立命令，也不是大文件有界处理的主路径。

阅读时重点看阶段生命周期、广播、旧中间态释放、CPE cache reset、失败输入 dump
和计时口径。各阶段算法细节应回到对应 `commands/` 或 `algorithms/` 阅读。

## file-backed dedup-workflow

```bash
./RabbitBAM-MPI dedup-workflow \
  -i input.bam -o markdup.bam \
  -T ../performance_test/tmp/dedup-run \
  -m 8G -n 66 --compress-level 1
```

阶段文件依次为 `<prefix>.collate.bam`、`<prefix>.fixmate.bam` 和
`<prefix>.sort.bam`。默认在下一阶段成功后删除上一个中间文件；失败时保留已有文件，
便于从失败阶段继续排查。`--keep-intermediates` 可保留全部阶段结果。

最后一步复用全内存 pipeline 已验证的 `coordinate-stream-two-pass` 路径；新的单遍
`markdup --stream` 仍由 standalone 命令显式选择和单独评价。

该命令的 stage wall 包含当前服务器上的实际中间文件 load/dump；算法核心性能仍看
各 standalone 命令内部打印的 `4.1~4.6`。它与全内存 `dedup-pipeline` 的目的不同：
前者优先可扩展性和可调试性，后者用于估计消除中间存储开销后的性能上界。
