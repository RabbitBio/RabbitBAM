# Raw BAM 适配层

本目录连接 decoded BGZF batch 与应用可消费的 Raw BAM records，同时提供通用过滤和
重新写出能力。它被编入 `swbam_cpe_runtime`。

## 文件

- `swbam_raw_bam.cpp`：扫描 decoded block 中的记录边界，构造零拷贝
  `RawBamRecordView` 并调用 consumer。
- `swbam_raw_bam_filter.cpp`：读取 mandatory fields，执行常用过滤并把保留视图传给
  downstream。
- `swbam_raw_bam_writer.cpp`：序列化 header、按记录边界填充 BGZF payload、调用 CPE
  压缩、写 output backend 并追加 EOF。

## 生命周期与约束

`RawBamRecordView` 指向当前 decoded batch，不拥有 payload，回调结束后即失效。
当前 parser 假设单条 BAM record 不跨 BGZF block。复杂命令若需要长期保存记录，必须
复制必要字段或 payload，不能保存 view 指针。
