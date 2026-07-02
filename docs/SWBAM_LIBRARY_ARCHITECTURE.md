# SWBAM 库架构说明

## 目标

SWBAM 将 BAM 数据搬运、异构流水调度与具体命令语义分离。库提供可复用的
批处理流水线，同时允许性能敏感的命令保留融合 CPE kernel，避免为了通用接口
而强制生成中间 record batch。

核心原则是：扩展接口的调用成本只发生在每个 BGZF batch，而不进入逐条 record
热路径。

## 分层结构

### `swbam_io`

- `BamInputBackend` 定义输入后端接口。
- `MemoryBamInput` 是当前实现，保持项目现有的“完整文件加载到内存”模拟 I/O
  方式。
- `BgzfBlockSpan`、`BgzfBlockBatch` 和 `BgzfSpanBatchReader` 提供 BGZF block
  索引、64 字节对齐批缓冲和按 span 批量读取。
- 后续可以增加 POSIX 或 MPI-IO backend，而无需修改上层算子。

### `swbam_mpi_runtime`

- `MpiBamInputPlan` 保存全局 BGZF 索引和当前 rank 的 block 范围。
- `PrepareMpiBamInputPlan` 完成 rank0 扫描、索引广播和 block 分区。
- 公共的 rank 状态检查和计时规约也放在这一层。

### `swbam_cpe_runtime`

- `RunCpeReadPipeline` 管理两套 compressed batch 和两套 decoded batch。
- 默认将 MPE 读取下一批与 CPE 处理当前批重叠。该行为由显式 option 控制，
  不再依赖命令头文件里的隐藏宏。
- runtime 统一负责 launch/join、算子生命周期、状态校验顺序、批计时和缓冲轮换。
- `CpeBatchOperator` 负责具体命令的参数准备、融合 CPE 入口、结果校验和 MPE
  侧归并。

`flagstat` 和 `stats --basic` 是首批两个算子实例。它们原有的
decode+parse+count 融合 CPE kernel 没有拆开，因此不会增加中间数据写回。

## 通用路径与融合路径

两条路径需要同时保留：

- 通用路径：后续由 `GenericDecodeOperator` 输出解压后的 BGZF 数据或 raw BAM
  record batch，供简单工具、库使用者和功能原型直接消费。
- 融合路径：内置性能敏感命令继续使用 decode+flagstat、decode+stats、
  decode+filter、decode+collate-extract 等专用融合算子。

二者共用输入 backend、批缓冲、双缓冲调度、CPE 生命周期和错误处理，仅算子
不同。markdup candidate ownership、collate QNAME exchange 等算法专用状态继续留在
具体命令中，不进入 BAM I/O 核心库。

## `flagstat` 实例

`FlagstatCpeOperator` 实现 `CpeBatchOperator`：

1. `Initialize` 分配每个 CPE 的 scratch 和计数 slice。
2. `Prepare` 将当前 compressed/decoded block、scratch 和 slice 绑定到
   `MpiFlagstatCountPara[64]`。
3. `kernel_entry` 返回现有 `slave_mpi_flagstat_count`。
4. CPE 在一次 kernel 中完成 BGZF 解压、BAM record 解析和 flagstat 计数。
5. `Validate` 检查每个 CPE 的状态，`Consume` 将 64 个计数 slice 合并到 rank
   结果。

命令本身不再维护读取、双缓冲和 launch/join 循环。

## `stats --basic` 实例

`StatsBasicCpeOperator` 使用同一个 runtime，但提供
`slave_mpi_stats_basic_count` 及其参数：

1. CPE 融合执行 BGZF 解压、BAM 解析、SN 统计、insert-size/orientation 统计和
   block 内坐标有序性检查。
2. `Consume` 合并 block sort state 和 record 数。
3. `Finish` 在全部 batch 完成后合并 64 个持续累积的统计 slice。

这说明同一套 pipeline 可以承载不同参数结构、结果结构和融合 kernel。

## 为什么从核代码暂时没有拆文件

当前已经抽出的是主核控制面和公开扩展接口。从核 kernel 仍位于 `slave.cpp`，原因是
它们依赖同文件中的 BAM fast parser、libdeflate cache、BGZF codec 和若干静态 helper。
直接移动两个入口会复制这些 helper，或者引入跨文件符号和初始化风险。

后续应先把真正通用的从核能力拆为独立 CPE object target：

- BGZF decode/encode；
- 每 CPE 的 libdeflate compressor/decompressor cache；
- raw BAM record iterator/fast parser；
- 公共状态码、容量检查和计时字段。

然后再把 `slave_mpi_flagstat_count`、`slave_mpi_stats_basic_count` 移入
`slave/operators/`。这一顺序可以避免为了整理目录而破坏已经验证的热路径。

## 如何增加新算子

新算子继承 `swbam::cpe::CpeBatchOperator`，实现：

1. `Initialize/Shutdown`：管理参数和结果 workspace。
2. `Prepare`：绑定当前 compressed/decoded batch 与 CPE 参数。
3. `kernel_entry/kernel_arguments`：指定每批启动的 CPE kernel。
4. `Validate/Consume`：处理状态并归并结果。
5. `Finish`：执行可选的跨 batch 最终规约。

命令只需把 rank-local spans 和算子交给 `RunCpeReadPipeline`，不再复制读取、
双缓冲和 CPE 调度循环。

## 对称 writer/compression pipeline

读流水线方向是：

`compressed BGZF -> CPE decode/parse/operator -> MPE consume`

写流水线方向则是：

`record/raw batch -> CPE pack/compress -> compressed BGZF batch -> output backend`

计划中的 `CpeWritePipeline` 将统一管理：

- 两套未压缩输入 batch 和两套压缩输出 batch；
- CPE pack/compress 与 MPE 写出上一批的重叠；
- BGZF block 大小、CRC、EOF block 和错误处理；
- compression level、codec cache 和统一计时；
- memory、POSIX、MPI-IO output backend。

通用 writer 接受 raw record batch；sort、markdup 等命令也可以提供专用
rewrite+pack+compress 融合算子，继续避免中间拷贝。

## 当前完成边界

已经完成：

- memory input backend、BGZF 索引和批缓冲；
- MPI 输入计划；
- 带显式 overlap 策略的 CPE read pipeline；
- `flagstat` 和 `stats --basic` 两个融合算子实例。

后续阶段：

- 将通用 slave-side BGZF codec/cache/helper 拆为独立 CPE object target；
- 增加通用 decoded/raw-record operator；
- 增加对称 writer/compression pipeline；
- 增加 POSIX/MPI-IO 输入输出 backend；
- 优先迁移 filter 和转换命令，再迁移 sort/collate/markdup 等算法型命令。

该路线使每个阶段都可以独立验证，同时在公共库逐步完善的过程中保留现有专用
优化路径。
