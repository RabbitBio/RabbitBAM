# MPI 命令实例

本目录实现相对线性的 BAM/SAM 命令，以及复杂算法共用的应用层辅助逻辑。

## 文件

- `swbam_bam2bam_mpi.cpp`：BAM 过滤/重写，保留专用 decode+filter 路径。
- `swbam_bam2sam_mpi.cpp`：BAM 解压、解析并并行格式化为 SAM。
- `swbam_sam2bam_mpi.cpp`：切分 SAM 文本、解析、pack 并压缩为 BAM。
- `swbam_flagstat_mpi.cpp`：调用公共 `RunFlagstatPipeline`，负责 MPI 规划、归约和打印。
- `swbam_stats_mpi.cpp`：调用公共 `RunStatsBasicPipeline`，负责 MPI 归约、全局
  有序性判断、SN/diagnostic 输出。
- `swbam_fixmate_mpi.cpp`：公共 BAM read/write 两端之间，按 QNAME group 更新 mate
  字段与 `MC/MQ/ms` 标签；配对和 MPI 边界算法保留在 app。
- `swbam_io_check_mpi.cpp`：公共库微基准及 decode/raw/write 回环校验。
- `swbam_mpi_common.cpp`：旧调用点的兼容包装和应用层公共辅助函数。

## 阅读顺序

推荐 `io-check -> flagstat -> stats -> bam2bam -> bam2sam/sam2bam -> fixmate`。
前三个最容易看清公共 runtime 与专用 operator 的分工。

`swbam_mpi_common.cpp` 不是新 SDK 的首选入口；新增通用代码应优先调用
`swbam::io` 和 `swbam::mpi` 公共接口。

## Fixmate 的读写组合

```text
MemReader / BamInputBackend + spans
  -> RunCpeReadPipeline + FmReadOperator
  -> app 分组、plan/rewrite、记录打包计划
  -> CpeRecordWriteSession -> 公共 BAM write -> RankBodySink
```

- `FmReadSource` 适配内存输入和 backend batch；`FmReadOperator` 复用公共
  `BamReadBatchOperator`、`BindBamReadBatch` 及 passthrough 解压解析 kernel。
- `PostProcessBatch()` 只在 decode join 后进入 app 算法，再依次调用 plan、rewrite、
  写流水线。不是 decode 与压缩两个 CPE kernel 同时运行。
- 跨 batch 的 pending group 和跨 rank 的 leading/trailing group 独立保存记录字节，
  不长期引用即将复用的解析记录池。完整的中间 group 可直接处理本批记录指针。
- `FmCompressPlans()` 沿用原装块次序，将 `BamRecordPackBlock` 逐批提交给已有写适配器。
  每次 Submit 返回前已完成当前压缩，Finish 排空最后输出；之后才释放 rewrite arena。
- 配对、标签语义、边界交换及 `prefix -> middle -> suffix` 输出布局不变。
  standalone 与 `MpiFixmateMemoryToMemory` 共用上述核心。
- rewrite 的每个 CPE 输出区按 64 字节对齐，避免相邻从核非对齐拷贝覆盖记录尾部。
  间隔只存在于 workspace，序列化仍按 `l_data` 写出，不向 BAM 插入 padding。

读 workspace 比原单批读端多约 8 MiB/rank；写端仍是两套 scratch 加两套压缩输出。
默认 memory 输入及 body 缓存行为不变。`fixmate -m` 是生成 mate score 标签的开关，
不是 sort/collate 的算法内存限制。

结果见 [fixmate 流水线迁移验证](../../../docs/results/fixmate_流水线迁移验证.md)。
