# SWBAM 库架构说明

## 目标

SWBAM 将 BAM 数据搬运、异构流水调度与具体命令语义分离。库提供可复用的
批处理流水线，同时允许性能敏感的命令保留融合 CPE kernel，避免为了通用接口
而强制生成中间 record batch。

核心原则是：扩展接口的调用成本只发生在每个 BGZF batch，而不进入逐条 record
热路径。

## 物理目录与构建目标

公共库和实现实例现在按职责物理分离：

```text
include/swbam/                  公共 SDK 头文件
lib/io/                         输入后端、BGZF batch、RankBodySink
lib/mpi/                        MPI input plan 与分布式输出
lib/cpe/                        MPE 侧 CPE 读写流水线
lib/bam/                        raw BAM adapter、filter、writer
slave/core/                     CPE codec/parser
slave/operators/                通用和统计类 CPE 算子
slave/algorithms/               sort/collate/markdup/fixmate 专用 CPE kernel
apps/mpi/commands/              转换、统计、fixmate 等命令实例
apps/mpi/algorithms/            sort、collate、markdup 全局算法实例
apps/mpi/pipelines/             dedup-pipeline 编排
examples/                       不依赖命令内部类型的 SDK 示例
```

正式静态库仍为 `swbam_io`、`swbam_mpi_runtime` 和 `swbam_cpe_runtime`。
`swbam_cpe_kernels` 是传播通用 CPE object files 的 INTERFACE target：Sunway hybrid
linker 要求从核对象直接参与最终链接，不能把它们先放进普通主核静态 archive。
sort/collate/markdup/fixmate 专用对象由独立的 `swbam_algorithm_kernels` 提供，
不会进入普通 SDK 示例。

应用只需链接聚合目标 `swbam::swbam`。`RabbitBAM-MPI` 自身也使用该目标，并额外
链接 `swbam::algorithm_kernels`，避免主程序和 SDK 示例采用两套不同链接方式。
公共总头文件是 `swbam/swbam.h`。

## 分层结构

### `swbam_io`

- `BamInputBackend` 定义输入后端接口。
- `MemoryBamInput` 保持项目现有的“完整文件加载到内存”模拟 I/O 方式，也是默认
  性能基线。
- `PosixBamInput` 使用 htslib 读取 BAM header，扫描阶段按大窗口预读并解析
  BGZF header，处理阶段
  通过 `preadv` 将一个连续 span batch 直接填入最终 block slots，不要求完整文件
  常驻内存。
- `BgzfBlockSpan`、`BgzfBlockBatch` 和 `BgzfSpanBatchReader` 提供 BGZF block
  索引、64 字节对齐批缓冲和按 span 批量读取。
- `SamInputBackend` 定义 SAM 文本输入的 header/body offset、按换行切分
  rank range 与 range read 接口。`PosixSamInput` 只读取每 rank 自己的
  文本范围，不要求每个 rank 常驻完整 SAM。
- MPI-IO input 放在 MPI runtime 层，而不让 `swbam_io` 依赖 MPI。

磁盘 BGZF scanner 使用 8 MiB 固定窗口顺序预读。窗口内只解析
18-byte header 中的 `BSIZE`，按 block length 跳转，只在剩余窗口不足
header/EOF 检查时再发起读取。对 WES_0.25G 约将 `12100` 次小
`pread` 降为几十次连续读。该策略会顺序预读大部分压缩文件，
而不是仅读取 header 字节；它以额外顺序带宽换取大量小 I/O 的消除，
并能预热后续 body read 所需的文件页。

### `swbam_mpi_runtime`

- `MpiBamInputPlan` 保存全局 BGZF 索引和当前 rank 的 block 范围。
- `PrepareMpiBamInputPlan` 完成 rank0 扫描、索引广播和 block 分区。
- `MpiIoBamInput` 复用 POSIX metadata/header 语义，body 通过 MPI-IO 读取。
  对连续 spans 创建 hindexed memory datatype，一次将文件连续字节直接散布到
  64 个最终 block slots，不使用中间 scratch buffer。
- `MpiIoSamInput` 实现同一 `SamInputBackend`，复用 POSIX header/range
  边界语义，使用 chunked `MPI_File_read_at` 读取当前 rank 的 SAM 范围。
- `MpiBamInput` 是统一 backend factory/lifetime：支持
  `memory|posix|mpiio|auto`，由
  rank0 选择后广播，保证所有 rank 使用同一种 backend。
- 公共的 rank 状态检查和计时规约也放在这一层。

`auto` 使用 `--io-memory-limit` 和当前可用内存的保守预算判断；若 memory 打开在
任一 rank 失败，所有 rank 一起回退 POSIX。MPI-IO 需要显式选择，避免在未知
并行文件系统特性时自动改变 I/O 策略。默认仍是显式 `memory`，论文性能测试
不会因机器内存状态不同而悄悄切换路径。

### `swbam_cpe_runtime`

- `RunCpeReadPipeline` 管理两套 compressed batch 和两套 decoded batch。
- 默认将 MPE 读取下一批与 CPE 处理当前批重叠。该行为由显式 option 控制，
  不再依赖命令头文件里的隐藏宏。
- runtime 统一负责 launch/join、算子生命周期、状态校验顺序、批计时和缓冲轮换。
- `CpeBatchOperator` 负责具体命令的参数准备、融合 CPE 入口、结果校验和 MPE
  侧归并。

`flagstat` 和 `stats --basic` 是首批两个算子实例。它们原有的
decode+parse+count 融合 CPE kernel 没有拆开，因此不会增加中间数据写回。
两个命令现已支持 `--io-backend memory|posix|mpiio|auto` 和
`--io-memory-limit SIZE`；backend 只改变 batch 来源，不改变算子实现。

## 通用路径与融合路径

两条路径需要同时保留：

- 通用路径：`RunGenericDecodePipeline` 通过 `DecodedBgzfConsumer` 输出解压后的
  BGZF batch；`RunGenericRawBamPipeline` 进一步生成零拷贝 `RawBamRecordView`
  batch，供简单工具、库使用者和功能原型消费。
- 融合路径：内置性能敏感命令继续使用 decode+flagstat、decode+stats、
  decode+filter、decode+collate-extract 等专用融合算子。

二者共用输入 backend、批缓冲、双缓冲调度、CPE 生命周期和错误处理，仅算子
不同。markdup candidate ownership、collate QNAME exchange 等算法专用状态继续留在
具体命令中，不进入 BAM I/O 核心库。

当前提供 `RabbitBAM-MPI io-check -i input.bam` 作为通用路径实例：它不执行任何
业务统计，只校验全局 BGZF block 数、decoded block 状态、解压字节数、raw record
边界和记录数，并打印 pipeline/kernel/consumer 计时。

`RawBamRecordView` 直接指向 decoded batch 中包含 4 字节 `block_size` 的原始编码，
不构造 `bam1_t`，也不复制 record payload。view 只在 consumer 回调期间有效。当前
与项目其他路径一致，要求单条 BAM record 不跨 BGZF block。

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

## 从核公共层与专用 kernel

通用从核能力目前已经拆入 `swbam_cpe_core_objects`：

- `slave/core/slave_bgzf_codec.cpp`：BGZF decode/encode、CRC、每 CPE 的
  libdeflate compressor/decompressor cache 和 cache 释放入口。
- `slave/core/slave_bam_parser.cpp`：从 decoded BGZF block 顺序读取 raw BAM
  record 的 fast parser。
- `include/swbam/cpe_codec.h` 与 `cpe_bam_parser.h`：稳定的 CPE core 接口。

原有 BAM2BAM、flagstat 和 stats 融合 kernel 已改为调用这些公共接口，不再各自
持有 codec/cache/parser 实现。

专用质检 kernel 也已拆入 `swbam_cpe_operator_objects`：

- `slave/operators/slave_bam_filter_operator.cpp`：BAM2BAM 的 MPI
  decode+parse+filter 与 passthrough 融合入口；
- `slave/operators/slave_flagstat_operator.cpp`：flagstat 逐记录计数和融合
  decode+parse+count 入口。
- `slave/operators/slave_stats_basic_operator.cpp`：stats SN、insert size、orientation、
  NM 和坐标有序性统计入口。

`slave.cpp` 不再定义这两个入口。parser 仍调用少量旧文件提供的 CIGAR/aux 修复
helper；这些 helper 需要在构造更多边界样例后再决定是否并入 core。

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

当前 `CpeWritePipeline` 已统一管理：

- 两套未压缩输入 batch 和两套压缩输出 batch；
- CPE compress 当前批与 MPE 消费上一批的重叠；
- BGZF block 大小、CRC、compression level、codec cache、错误处理和统一计时；
- `UncompressedBgzfSource` 与 `CompressedBgzfConsumer` 两侧扩展接口；
- `GenericCompressOperator` 通用 BGZF 压缩实例。

`RawBamWriter` 已补齐通用 raw record 写路径：它将 `RawBamRecordView` 按记录
边界装入 BGZF payload，以固定 256-block chunk 控制临时内存，再交给
`CpeWritePipeline` 压缩并追加到 `BamOutputBackend`。writer 还负责序列化标准 BAM
header、让 body 从独立 BGZF block 开始并追加标准 EOF block，因此输出是可重新打开
的完整 BAM，而不是只有 body 的内部片段。

`MemoryBamInput::OpenMemoryCopy()` 可把完整内存 BAM 重新作为下一阶段输入，形成
`memory input -> generic decode/raw view -> memory writer -> memory input` 闭环。当前
chunk source 到 write pipeline input 尚有一次 block copy；MPI-IO output backend
尚未实现。sort、markdup 等命令仍可提供专用 rewrite+pack+compress 融合算子，
避免为了通用接口强制落地中间 record batch。

## 输出 backend

`BamOutputBackend` 只暴露批外 `Write/Flush/bytes_written`，不会进入逐记录热路径。
目前有两个实现：

- `MemoryBamOutput`：默认模拟高速存储，支持 reserve、连续 append 和完整内存 BAM；
- `PosixBamOutput`：真实顺序文件输出，为后续高速 SSD/并行文件系统保留直接路径。

`RawBamWriter`、BGZF block consumer、BAM2BAM rank-local body append 和最终 memory
dump 已经通过该接口。BAM2BAM 使用一个只包装既有 `MemWriter` 的 adapter，因此
rank gather/layout 完全不变；实测 backend 化前后 `4.3` 只相差约 0.15%。

## 通用过滤与专用过滤

`RawBamFilterConsumer` 是由公共库构建的第一个通用 BAM 应用算子。它直接从 raw
BAM mandatory fields 读取 `tid/MAPQ/FLAG/l_qseq`，不构造 `bam1_t`；通过的记录
仍以 `RawBamRecordView` 交给下游，因此 aux、序列、质量值等 payload 不发生重建。
它可连接 `RawBamWriter`，并按配置使用 memory 或 POSIX output backend。

默认 BAM2BAM 命令没有因此降级。其
`slave_mpi_decompress_filterfunc` 和 passthrough 入口作为专用融合算子保留在
`swbam_cpe_operator_objects`，CPE 一次完成 decode+parse+filter，主机继续使用成熟的
flat pack 与 record serialize+compress 路径。通用算子用于快速构建实例和功能验证，
专用算子用于论文性能结果；二者共享 codec、parser、backend 和数据语义。

`io-check` 分两层验证写流水线：确定性合成 payload 用来检查 BSIZE、ISIZE、CRC
和逐字节内容；真实 BAM 回环则重新序列化 header、打包全部 raw records、CPE 压缩
并追加 EOF，再由 `MemoryBamInput` 重新打开，校验 header、body 索引、CRC、记录数和
原始记录字节数。

## BAM2BAM 命令输入迁移

`run_all` 的 BAM→BAM 路径已接入 `MpiBamInput` factory，支持
`--io-backend memory|posix|mpiio|auto` 和 `--io-memory-limit`。各 backend 共用
`MpiBamInputPlan` 的全局 BGZF 索引与 rank 分区，但热路径有意保留两种实现：

- memory backend 直接将该 rank 的连续 BGZF 窗口交给原
  `FusedBamToBamMPI`，不新增逐 block 拷贝；
- POSIX backend 在同一个融合核心内通过 `ReadBatch`/`preadv` 将每批 BGZF
  block 直接读入既有 CPE input slots，工作内存不随输入文件增长。
- MPI-IO backend 使用派生 memory datatype 在一次 `MPI_File_read_at` 中将连续
  文件范围散布到同样的 CPE input slots。

BAM header 也按 backend 分流：memory 保留原文件 header 字节，POSIX 使用
htslib 解析结果重建合法压缩 header。CPE decode/filter、MPE flat pack、CPE
compress 和双缓冲 overlap 没有替换。默认 memory 输出仍按既有逻辑在
rank0 组装完整输出内存并 dump；显式选择 MPI-IO 时则使用下述分布式路径。

## MPI-IO 分布式输出

`swbam_mpi_runtime` 提供 `MpiFileOutput`，封装 collective open/resize/sync/close
和支持超过 `INT_MAX` 字节的 chunked `WriteAt`。它不继承顺序
`BamOutputBackend`：前者表示多 rank 全局 offset 布局，后者表示单一顺序字节流，
将两种语义强行放进同一虚接口会隐藏 collective 约束。

在此之上，`DistributedBamOutput` 统一处理最终 BAM 的公共编排：

- `PrepareBody()` 收集每个 rank 的 body size，计算 body prefix 和总长度；
- `SetEnvelope()` 广播序列化 header 与 BGZF EOF 的长度，完成全局文件布局；
- `GatherToRoot()` 以固定8 MiB消息把任意 `RankBodySource` 汇总成完整内存 BAM；
- `WriteMpiIo()` 使用相同布局将 header、各 rank body 和 EOF 分布式写出。

该对象不负责修改 SAM header，也不参与压缩热路径。命令实例仍决定 `SO/GO/PG`
语义并生成序列化 header，公共运行时只负责可靠放置这些字节。这样既避免把算法
语义塞入 I/O 库，也消除了 sort/collate 中重复的 size/prefix、chunk gather 和
MPI-IO 代码。

`run_all --io-output-backend mpiio` 的布局是：

- rank0 写 `[0, output_body_start)` 的 SAM/BAM header；
- 各 rank 在 `output_body_start + prefix(rank)` 直接写自己的压缩 body；
- BAM 输出由 rank0 在末尾写标准 BGZF EOF；
- 文件在打开时 collective resize，避免覆盖旧文件时留下尾部数据。

该路径不再将所有 body 发往 rank0，也不分配 rank0 完整输出 buffer。转换命令、
standalone fixmate、markdup、sort 和 collate 已接入下述 `RankBodySink`，本地
body 可以独立选择 memory 或 spool。默认 `memory` 路径保留不变，继续用于当前
核心性能模拟。

## RankBodySink：rank 本地输出抽象

`swbam_io` 提供只读的 `RankBodySource` 和可追加的 `RankBodySink`。它们只出现
在每批压缩或格式化完成之后，不进入逐 record 热循环，接口包括顺序
`Append()`、最终 `size()`、分块 `ReadAt()` 和生命周期管理。当前实现包括：

- `MemoryRankBodySink`：兼容包装原 `MemWriter`，保留原容量增长和连续内存布局；
- `SpoolRankBodySink`：压缩块顺序写入 rank-local 临时文件，最终按 8 MiB
  分块读出并交给 `MpiFileOutput::WriteAt()`；
- `AdaptiveRankBodySink`：支持 `memory|spool|auto`。`auto` 先写内存，body 超过
  `--rank-body-memory-limit` 时将已有内容一次转入 spool，后续继续顺序追加；
- `SegmentedRankBodySource`：逻辑拼接多个 source，`ReadAt()` 可跨 segment 读取，
  用于 fixmate 的 prefix/middle/suffix，避免整段 `memmove`。

spool/auto 要求显式提供 `--rank-body-temp-dir`，例如
`../performance_test/tmp`。临时文件通过 `mkstemp` 创建并立即 unlink，不搜索或
自动迁移到其他目录；创建失败、空间不足或写入失败会直接报错并提示检查该目录。
这使性能实验的存储位置和计时口径可控。`Flush()` 不调用 `fsync`：同一文件
描述符上的后续 `pread` 已能看到此前写入内容，而强制落盘只会增加额外等待。

本地 body 与最终输出已经解耦：`--rank-body-backend` 控制本地存储，
`--io-output-backend memory|mpiio` 控制最终文件。memory 最终输出仍会在 rank0
分配完整输出内存，因此只有 `spool/auto + mpiio` 才是输出侧有界内存路径。
默认 `memory + memory` 保持原性能基线。

## 文件转换实例

`run_all` 的三种转换现在共用库的输入分区和输出布局：

- BAM→BAM：使用 `BamInputBackend + MpiBamInputPlan`，保留专用 CPE
  decode/filter/compress 融合算子；
- BAM→SAM：`FusedBamToSamMPI` 提供 `MemReader` 和
  `BamInputBackend + spans` 两个入口，memory 继续走原路径，POSIX/MPI-IO
  按 batch 读取，格式化 CPE kernel 不变；
- SAM→BAM：memory 保留原完整内存模拟，显式 `posix/mpiio` 时使用
  `SamInputBackend` 仅读本 rank 范围，然后交给原 CPE copy/count/parse/compress
  融合路径。SAM 的 `auto` 保守选择 POSIX streaming，避免大文本在每 rank
  上完整复制；要求 memory 模拟时需显式使用默认 `memory`。

三种转换的最终 memory gather 与 MPI-IO offset write 由同一布局逻辑处理。
三种融合核心现在直接接受 `RankBodySink`，因此 BAM→SAM 的文本膨胀也可在
阈值处动态 spill。公共库负责存储、分区、batch 与生命周期，命令仍负责 SAM
语义和高性能融合算子。

## Sort 与 Collate：算法型实例

standalone `sort` 和 `collate` 已完整接入公共 I/O 闭环，同时保留算法专用核心：

- 输入统一由 `MpiBamInput` 打开，由 `PrepareMpiBamInputPlan` 建立全局 BGZF
  block 索引并切分 rank 范围；
- memory backend 继续从完整输入的连续 rank 窗口构造 `MemReader`，保持原热路径；
- POSIX/MPI-IO backend 使用 `BgzfSpanBatchReader` 直接填充算法已有的64个
  `bam_block` slots，不聚合整段 rank 压缩输入；
- sort 的内排 extract、严格外排 run generation，以及 collate 的 memory extract、
  external run generation 都消费相同 block batch；全局 block index 仍来自 input
  plan，稳定排序和 QNAME 元数据语义不变；
- 压缩端直接向 `RankBodySink` 追加 BGZF block。memory-to-memory pipeline 通过
  `MemoryRankBodySink` 包装旧 `MemWriter`，因此四阶段融合流程的存储策略不变；
- standalone 与 pipeline helper 的最终 body 组装均使用
  `DistributedBamOutput`；同一 `RankBodySource` 可服务 memory gather 或
  MPI-IO 分布式写入。

这里公共库只管理“数据如何进入批次、body 如何保存和最终如何输出”。sort 的
sample/partition/exchange/k-way merge、严格外排，collate 的 hash/bin、ring
exchange、raw arena、loser tree 和外排归并仍在命令实例中。它们属于算法专用
融合算子，不应塞进通用 BAM I/O 接口。

两类临时文件必须区分：`-T/--temp-prefix` 管理 sort/collate 算法外排的
run/segment；`--rank-body-temp-dir` 只管理最终压缩 body 的 spool。两者可以指向
同一目录，但分别配置、分别计时；库不会自动寻找其他磁盘。

## 当前完成边界

已经完成：

- memory input backend、BGZF 索引和批缓冲；
- POSIX streaming input backend，支持 BGZF span scan 和 batch `preadv`；
- MPI-IO streaming input backend，支持派生 memory datatype 直接 scatter read；
- SAM text input backend，支持 POSIX/MPI-IO rank-range streaming；
- MPI backend factory 与 flagstat/stats 的 `memory|posix|mpiio|auto` 选择；
- MPI 输入计划；
- 带显式 overlap 策略的 CPE read pipeline；
- 独立的 CPE BGZF codec/cache 和 raw BAM fast parser object target；
- 独立的 CPE operator object target，以及 `flagstat`、`stats --basic` 两个融合
  算子实例。
- `GenericDecodeOperator`、`DecodedBgzfConsumer` 和可运行的 `io-check` 通用解码
  实例。
- `RawBamRecordView`、`RawBamRecordConsumer` 和通用 raw BAM batch adapter。
- `CpeWritePipeline`、`GenericCompressOperator` 和通用 source/consumer 接口。
- `RawBamWriter`、标准 BAM header/EOF 序列化、`MemoryBamOutput` 与
  `PosixBamOutput`；完整内存
  BAM 可以通过 `OpenMemoryCopy()` 再次进入读流水线。
- `RawBamFilterConsumer` 通用二进制过滤算子，以及独立的 BAM2BAM CPE 融合
  filter/passthrough operator。
- BAM2BAM 命令的 `memory|posix|mpiio|auto` 输入选择；memory 继续使用
  连续窗口快速路径，POSIX/MPI-IO 使用有界 batch 读取。
- `MpiFileOutput` 与 `run_all --io-output-backend mpiio`，各 rank 可按全局
  offset 直接写 header/body/EOF，不再在 rank0 汇总完整输出。
- BAM→BAM、BAM→SAM、SAM→BAM 三种转换均已接入公共 input/output
  backend，且保留各自专用 CPE 融合核心。
- `fixmate -m` 已接入同一套 `MpiBamInput`、`MpiBamInputPlan` 和
  `MpiFileOutput`。默认 memory backend 仍从 rank 的连续压缩窗口构造
  `MemReader`；POSIX/MPI-IO backend 则通过 `BgzfSpanBatchReader` 每批直接
  填充 64 个 CPE input slots。两条路径共用原 fixmate 解压、QNAME 分组、
  边界 group 交换、mate tag 重写和压缩核心。
- standalone `markdup --stream` 已迁移到公共 backend。内部使用可重置的
  block-pass 适配器，对同一 rank spans 完成首坐标探测、流式候选扫描和最终
  rewrite 三遍读取；每遍都只保留双槽、每槽64个 BGZF block。memory backend
  继续走连续 `MemReader`，POSIX/MPI-IO 直接填充原 markdup input slots，
  duplicate window、跨 rank coordinate ownership 和 CPE rewrite/compress
  算子保持不变。MPI-IO 输出沿用 header/body-prefix/EOF 分布式布局。
- `RankBodySource`、`AdaptiveRankBodySink` 和 `SegmentedRankBodySource` 已进入
  `swbam_io`。三种转换、standalone fixmate 和 markdup 均支持独立的
  `--rank-body-backend memory|spool|auto`；fixmate 通过三段 source 保持边界
  group 输出顺序，pipeline 的 memory-to-memory 接口保持连续内存实现。
- standalone `sort` 和 `collate` 已成为算法型公共 I/O 实例，支持
  `memory|posix|mpiio|auto` 输入、`memory|spool|auto` rank body 和
  `memory|mpiio` 最终输出；内排、外排和 CPE 融合核心保持原实现。
- `DistributedBamOutput` 已用于 sort/collate 的 standalone 和
  memory-to-memory helper，统一 body size/prefix、header/EOF 布局、rank0
  memory gather 与 MPI-IO 写出。

后续阶段：

- 评估非流式 markdup 是否值得接入 streaming backend，或继续作为 memory
  flat-hash 极致性能实例；
- 减少 packer chunk 到 write pipeline input 的一次 block copy；
- 将转换、fixmate 和 markdup 中仍然保留的最终输出编排逐步切换到
  `DistributedBamOutput`；
- dedup pipeline 暂时仍采用完整内存 BAM 广播，是否进一步改成分布式阶段交接
  取决于论文是否需要证明超大文件端到端有界处理。

该路线使每个阶段都可以独立验证，同时在公共库逐步完善的过程中保留现有专用
优化路径。
