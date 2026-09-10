# CPE 可复用算子

本目录放通用 operator 和较简单的命令融合 operator。它们调用 `slave/core/`，由
`lib/cpe/` 或相应 MPI command 在 MPE 侧准备参数和启动。

## 文件

- `slave_decode_operator.cpp`：通用 BGZF decode。
- `slave_compress_operator.cpp`：通用 BGZF compress。
- `slave_bam_filter_operator.cpp`：BAM2BAM decode+parse+filter/passthrough。
- `slave_flagstat_operator.cpp`：decode+parse+flagstat count。
- `slave_stats_basic_operator.cpp`：decode+parse+basic stats 和有序性信息。

通用 operator 适合 SDK 组合；融合 operator 减少 decoded record 中间态和主从核传输，
用于内置命令性能路径。新增 operator 时应把通用 codec/parser 留在 `core/`，不要复制
解压和缓存实现。
