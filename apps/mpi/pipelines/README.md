# MPI 组合流水线

本目录只负责把已经存在的阶段组合起来，不重新实现 collate、fixmate、sort 或
markdup 算法。

## 当前文件

- `swbam_dedup_pipeline_mpi.cpp`：实现
  `collate -> fixmate -m -> coordinate sort -> markdup`。

当前 pipeline 使用完整内存 BAM 作为阶段间格式：rank0 组装阶段输出并广播给下一
阶段。这样可以复用已验证的 memory-to-memory helper，但峰值内存和调试复杂度高于
依次执行独立命令，也不是大文件有界处理的主路径。

阅读时重点看阶段生命周期、广播、旧中间态释放、CPE cache reset、失败输入 dump
和计时口径。各阶段算法细节应回到对应 `commands/` 或 `algorithms/` 阅读。
