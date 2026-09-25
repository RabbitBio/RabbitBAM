# MPE 侧可复用算子

本目录实现建立在公共 batch 流水线之上的 MPE 侧算子封装。算子负责准备 CPE
参数、选择 kernel、校验状态和归并 rank-local 结果；输入读取、双缓冲和
launch/join 由 `lib/cpe/` 的 runtime 负责。

- `swbam_flagstat_operator.cpp`：封装 decode + parse + flagstat count 融合
  CPE kernel，并公开 `FlagstatOperator` 与 `RunFlagstatPipeline`。
- `swbam_stats_basic_operator.cpp`：封装 decode + parse + basic stats 融合
  CPE kernel，并维护跨 batch 的计数直方图和坐标排序边界。
- `swbam_bam_transform_operator.cpp`：封装 BAM 直通或过滤 CPE kernel，
  将保留记录按 BGZF 容量打包，再逐批提交给公共 BAM record 写流水线。
- `swbam_sam_to_bam_operator.cpp`：接收 SAM read 流水线产生的 parsed batch，
  生成记录打包计划，并逐批提交给公共 BAM write 流水线。
- `swbam_bam_to_sam_operator.cpp`：复用 BAM read 调度与 passthrough 解析 kernel，
  在 MPE 批次后处理时按原始块顺序提交给 SAM write 流水线。

对应的从核 kernel 位于 `slave/operators/`。这里的“operator”是批级执行单元，
不会在逐记录热路径引入虚函数调用。
