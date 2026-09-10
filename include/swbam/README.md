# SWBAM 公共头文件

本目录定义新版 SDK 的公开编译接口。普通使用者优先包含 `swbam/swbam.h`；需要控制
依赖或开发底层算子时再包含细分头文件。

## 高层接口

- `swbam.h`：公共总入口。
- `io.h`：BAM/SAM backend、BGZF batch、输出 backend 和 RankBodySink。
- `raw_bam.h`：零拷贝 `RawBamRecordView`、consumer 和通用 raw pipeline。
- `raw_bam_filter.h`：通用过滤 consumer。
- `raw_bam_writer.h`：Raw BAM pack、压缩和完整 BAM writer。
- `mpi_runtime.h`：MPI input plan、MPI-IO 和分布式 BAM 输出。

## 流水线扩展接口

- `cpe_pipeline.h`、`generic_decode.h`：CPE 读取、解压和批级 operator。
- `cpe_write_pipeline.h`、`generic_compress.h`：CPE 压缩和输出流水。

## 底层/高级接口

- `cpe_codec.h`：CPE BGZF codec 与每核 codec cache。
- `cpe_bam_parser.h`：decoded block 内的 BAM record fast parser。

公共接口目前仍暴露 `bam_block`、`sam_hdr_t`、MPI 和部分项目类型，因此适合作为
当前神威工程 SDK，但尚未形成平台无关、ABI 稳定的外部软件包。
