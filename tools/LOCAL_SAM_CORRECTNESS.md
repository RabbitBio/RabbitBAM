# SWLS 原地 SAM 内容校验

`swbam-bam-correctness-digest` 是独立 HTSlib 校验工具，不调用 SWBAM
解析器或 CPE kernel，也不修改 SDK/app 的执行路径。

```bash
cmake -S . -B build-check -DBUILD_SWBAM_TOOLS=ON
cmake --build build-check --target swbam-bam-correctness-digest

# 在神威计算节点运行，stdout 为小型 JSON，stderr 为进度或错误。
./swbam-bam-correctness-digest --bam output.bam --progress-every 5000000 > output.json
```

工具使用 HTSlib 1.20 的 `sam_read1()` / `sam_format1()` 将 BAM 转成
SAM 记录后立即计算摘要，不落地完整 SAM。仅排除头部 `@PG`；保留全部
11 个字段和全部 AUX，AUX 按字段排序以忽略合法 tag 排列差异。

- `header_without_pg_sha256`：其余头部的有序 SHA-256。
- `record_multiset_sha256_sum`：每条规范化 SAM 记录的 SHA-256 按模
  `2^256` 相加，保留重复次数，允许 sort/collate 的合法记录重排。
- `ordered_sha256`：记录有序摘要，另作严格行序诊断；sam2bam 验收仍要求一致。
- DUP 数量、清除 DUP 后的记录摘要和 DUP 成员摘要：辅助定位标记差异。
- 坐标逆序、QNAME 组内 READ1/READ2 顺序：算法输出约束检查。

Collate 还需要确认同一 QNAME 不会再次出现：

```bash
set -euo pipefail
export LC_ALL=C
./swbam-bam-correctness-digest --bam collate.bam --groups qnames.txt > collate.json
sort -u -S 1G --parallel=1 -T ./tmp qnames.txt | wc -l
```

唯一 QNAME 数应等于 JSON 中的 `qname_groups`。这个手工方案需要 GNU
sort；计算节点的 BusyBox sort 不支持这些资源选项。

本轮自动校验不需要再次排序组名：A4 参考已精确得到唯一 QNAME 数。
只有在完整记录多重集一致的前提下，输出的 QNAME 集合必与参考一致；
因此输出观测组数等于参考唯一名数，即可证明每个名字只出现一个连续
组。该推理依赖完整记录摘要通过，不能单独用组数相等代替内容校验。

自动控制器还可按应用日志中六段输出 body 的已知 BGZF 边界并行校验，
只占一个节点的六个 MPE。工具的 `--range-begin` / `--range-end` 为压缩
字节偏移，必须是已知记录对齐的 BGZF 边界，不能随意按文件字节均分。
控制器合并记录摘要，并修正分段接缝上的坐标、QNAME 组数及 READ 顺序。
全局有序 SHA-256 不能从各段 SHA-256 合并，分段模式将此诊断置为 null；
要求严格行序的 sam2bam 仍采用整文件校验。

两端必须使用相同摘要口径。与 A4 已有参考 JSON 比较即可，不回传整个
BAM。摘要匹配是强概率性一致性证据，并非压缩 BAM 逐字节一致的证明。
读取、格式化、截断输入或写入失败时工具返回非零，不输出成功 JSON。

原 SAM stdin + OpenSSL 构建仍保留；Sunway 构建使用带 MIT 许可的
PicoSHA2，便于在没有目标 OpenSSL 库时生成同一 SHA-256。该工具只用于
正确性验证，其耗时不计入论文性能成绩。
