# MPI 命令实例

本目录实现相对线性的 BAM/SAM 命令，以及复杂算法共用的应用层辅助逻辑。

## 文件

- `swbam_bam2bam_mpi.cpp`：BAM 过滤/重写，保留专用 decode+filter 路径。
- `swbam_bam2sam_mpi.cpp`：BAM 解压、解析并并行格式化为 SAM。
- `swbam_sam2bam_mpi.cpp`：切分 SAM 文本、解析、pack 并压缩为 BAM。
- `swbam_flagstat_mpi.cpp`：融合 decode+parse+flag count 的统计实例。
- `swbam_stats_mpi.cpp`：`stats --basic`、有序性和统计归并。
- `swbam_fixmate_mpi.cpp`：按 QNAME group 更新 mate 字段与 `MC/MQ/ms` 标签。
- `swbam_io_check_mpi.cpp`：公共库微基准及 decode/raw/write 回环校验。
- `swbam_mpi_common.cpp`：旧调用点的兼容包装和应用层公共辅助函数。

## 阅读顺序

推荐 `io-check -> flagstat -> stats -> bam2bam -> bam2sam/sam2bam -> fixmate`。
前三个最容易看清公共 runtime 与专用 operator 的分工。

`swbam_mpi_common.cpp` 不是新 SDK 的首选入口；新增通用代码应优先调用
`swbam::io` 和 `swbam::mpi` 公共接口。
