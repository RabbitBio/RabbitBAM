# I/O 与 BGZF 数据面

本目录构成 `swbam_io` 静态库，负责把存储介质转换成统一 BGZF batch，并为顺序或
分布式输出提供可组合的数据源/数据槽。

## 文件

- `swbam_io.cpp`：
  - `MemoryBamInput`、`PosixBamInput`；
  - BAM header/body offset；
  - 8 MiB 窗口 BGZF scanner；
  - `preadv` batch read；
  - SAM range 切分；
  - memory/POSIX output。
- `swbam_rank_body_sink.cpp`：
  - memory、spool、auto spill rank body；
  - segmented source；
  - byte-size 参数解析和 BGZF block append。

## 重要概念

- `BgzfBlockSpan` 只描述压缩块在输入中的位置和长度。
- `BgzfBlockBatch` 管理供 MPE/CPE 使用的对齐 block slots。
- `RankBodySink` 保存某 rank 已压缩的输出 body，不负责 header、EOF 或全局 offset。

算法外排临时文件不属于本目录；sort/collate 的 run 文件由各自算法管理。
