# Raw BAM 适配层

本目录连接 decoded BGZF batch 与应用可消费的 Raw BAM records，同时提供通用过滤和
重新写出能力。它被编入 `swbam_cpe_runtime`。

## 文件

- `swbam_raw_bam.cpp`：扫描 decoded block 中的记录边界，构造零拷贝
  `RawBamRecordView` 并调用 consumer。
- `swbam_bam1.cpp`：将 raw record batch 物化为库管理的 `bam1_t` 对象池，供
  HTSlib 兼容 consumer 使用，并单独统计 materialize/consume 时间。
- `swbam_bam1_writer.cpp`：将只读或用户修改后的 `bam1_t` 批量编码回 raw BAM
  record，复用 arena、视图数组和 `RawBamWriter`，不为每条记录单独分配。
- `swbam_raw_bam_filter.cpp`：读取 mandatory fields，执行常用过滤并把保留视图传给
  downstream。
- `swbam_raw_bam_writer.cpp`：序列化 header、按记录边界填充 BGZF payload、调用 CPE
  压缩、写 output backend 并追加 EOF。

## 生命周期与约束

`RawBamRecordView` 指向当前 decoded batch，不拥有 payload，回调结束后即失效。
当前 parser 假设单条 BAM record 不跨 BGZF block。复杂命令若需要长期保存记录，必须
复制必要字段或 payload，不能保存 view 指针。

`RunGenericBam1Pipeline` 在 Raw adapter 之后增加一次 MPE 物化。`bam1_t` 和 data
由可复用对象池持有，同样只在当前 `ConsumeBam1()` 回调期间有效；需要跨回调保存时
调用 `bam_dup1()`。该路径支持标准 HTSlib accessor，但会产生 payload 复制和解析
成本，因此不会替代零拷贝 Raw 路径或性能敏感命令的融合 CPE kernel。

`Bam1Writer` 是该路径的对称写端。输入记录只在 `ConsumeBam1()` 调用期间读取，
writer 不取得其所有权；编码器移除 HTSlib 的 QNAME 内存 padding，校验数据布局与
CIGAR/query length，并把整批记录一次交给 `RawBamWriter`。当前仍要求单条编码记录
不超过一个 BGZF payload，且不写回 `n_cigar > 65535` 的长 CIGAR 记录。
