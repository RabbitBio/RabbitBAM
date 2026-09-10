# CPE 公共核心

这里放不依赖具体命令语义的从核基础能力，由多个融合算子共同调用。

## 文件

- `slave_bgzf_codec.cpp`：BGZF inflate/deflate、CRC/ISIZE、每 CPE 的
  libdeflate compressor/decompressor cache 及统一释放入口。
- `slave_bam_parser.cpp`：从 decoded BGZF block 顺序读取 BAM record 的 fast parser。

对应声明位于 `include/swbam/cpe_codec.h` 和
`include/swbam/cpe_bam_parser.h`。

本层只理解 BGZF 和 BAM 编码，不执行过滤、统计、排序或去重。修改这里会影响几乎
所有 CPE 路径，应优先进行 codec 正确性、CRC 和小/空 block 回归。
