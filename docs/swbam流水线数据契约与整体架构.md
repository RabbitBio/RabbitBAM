# SWBAM 流水线数据契约与整体架构

## 目标与边界

SWBAM 采用一套 batch 级运行时承载不同 BAM/SAM 操作，而不是为每个命令重复实现
I/O、双缓冲和 CPE 调度。流水线在 **batch 边界**组合，逐记录热点由一个 CPE
operator 静态融合执行，因此库化不应引入逐记录虚调用或强制中间对象物化。

当前项目约束是：一条 BAM record 完整位于一个 BGZF block 内。违反该约束的输入
必须显式失败，不能返回半条记录或静默丢失数据。

## 整体结构

```text
Input backend
    -> rank/block plan
    -> typed batch pipeline (buffer ownership + double buffering)
    -> CPE operator (codec/parser/business logic)
    -> MPE batch post-processing
    -> output sink/backend
```

规划中的四条类型化流水线为：

| 流水线 | 输入 | 主要 operator | 输出 |
| --- | --- | --- | --- |
| BAM read | BGZF spans | decode / decode+parse+业务 | 结果 slice、raw view 或算法元数据 |
| BAM write | raw/record batch | pack+compress | compressed BGZF batch |
| SAM read | 文本 range/chunk | split+parse | BAM record 或算法元数据 |
| SAM write | BAM record batch | format | SAM text batch |

应用负责 CLI、MPI 全局划分/归约和结果展示；库负责 backend、batch、调度、operator
及输出能力。sort、collate、fixmate、markdup 的全局算法状态仍属于应用层，但其
输入和输出边界复用公共流水线。

## Batch 数据契约

### 所有权与生命周期

- backend 拥有文件或内存输入；pipeline 只在运行期间借用 backend。
- pipeline 拥有 compressed/decoded 双缓冲，operator 不得释放或长期保存其地址。
- batch view 只在当前回调期间有效；跨 batch 保存记录必须显式复制。
- operator 拥有参数数组、CPE scratch 和 rank-local result slice。
- output sink 拥有写出状态；pipeline 只按顺序提交完整 batch。

### 顺序与索引

- `BgzfBlockSpan` 的顺序等于全局 BGZF 顺序，rank 只处理计划分配的连续区间。
- batch 内 block 顺序不得改变；需要全局重排的算法必须显式产生稳定 key/index。
- BAM read operator 返回的 record 数必须等于该批所有成功解析记录之和。
- writer 只接收完整 record；BGZF EOF 由最终输出层追加一次。

### 容量与错误

- `active_blocks <= operator.batch_capacity()`，当前 CPE 标准批容量为 64 blocks。
- 所有 arena 必须在 launch 前确定容量，CPE 不进行不可控动态扩容。
- backend read、kernel status、容量上限、解析和 post-processing 任一失败均终止
  当前 rank；MPI command 再通过 `AllRanksOk` 保证所有 rank 一致退出。
- 失败 batch 不得向后续处理暴露部分有效结果。

### 并发与同步

- 默认重叠 MPE 读取下一批和 CPE 执行当前批。
- kernel 完成并通过 `Validate` 后，MPE 才能读取结果 slice。
- `PostProcessBatch` 是 MPE 侧批后处理，不是 CPE 逐记录 consumer。
- operator 的公开多态调用只发生一次/batch，不进入 record loop。

### 计时口径

- pipeline 统一记录 `read`、`kernel`、`post_process` 和 `total`。
- operator 可细分 kernel 内的 decode/parse/business 时间，但不得重复计入 total。
- MPI scan/partition、跨 rank reduce 和最终打印由 command 单独计时。
- memory backend 的文件预加载和最终 dump 不计入算法核心时间。

## Operator 组合原则

operator 可以复用 codec、parser、counter、filter 等源代码组件，但生产路径在一次
CPE kernel 中静态融合，例如 `decode+parse+flagstat`。开发与验证路径可以输出
`RawBamRecordView` 或 `bam1_t` 供 MPE 自定义后处理，但这不是性能关键命令必须经过
的中间层。

原则概括为：**batch 层动态组合，record 层静态融合**。

当前 BAM read 的从核基础步骤为 `swbam_cpe_decode_bam_block`（BGZF 解压）和
`swbam_cpe_walk_bam_records`（BAM record 解析与遍历）。遍历接收编译期动作：
filter、passthrough、flagstat、stats-basic、count 复用同一骨架，仍由一次 CPE
launch 完成。主核 `BindBamReadBatch` 统一绑定 BGZF 输入/输出块；结果 arena、
状态 slice 和后处理仍由各算子拥有。bam2sam 也使用 passthrough 解析动作，解析得到
的记录仅在当前 batch 内交给 SAM write 流水线。

## Flagstat 首个迁移实例

```text
MpiBamInput / MpiBamInputPlan
    -> RunFlagstatPipeline
       -> RunCpeReadPipeline
          -> FlagstatOperator
             -> slave_mpi_flagstat_count
                (BGZF decode + BAM parse + count)
          -> merge 64 CPE count slices
    -> MPI_Reduce
    -> samtools-style output
```

`FlagstatOperator`、计数结果和 rank-local metrics 属于公共库；命令只负责参数、MPI
规划、跨 rank 归约与打印。迁移不改变 CPE kernel、参数布局、批容量、双缓冲或计时
含义，因而可作为后续 stats、filter 和格式转换库化的模板。

WES_0.25G、单节点 6 ranks 的迁移后 smoke 得到 `12,100` blocks、`2,668,351`
records，全部 flagstat 计数与既有基线一致；`4.3=0.051421 s`，相对历史库化基线
`0.051994 s` 无性能回退。

## Stats --basic 第二个迁移实例

`StatsBasicOperator` 与 flagstat 复用相同 BAM read runtime，但保持两种额外状态：

- 每个 CPE slot 跨 batch 持续累积 value、insert-size 和 orientation slices，全部
  batch 完成后由 `Finish()` 合并一次；
- 每个 block 返回首尾坐标和 block-local sorted 状态，`PostProcessBatch()` 按
  block 顺序合并为 rank-local `StatsSortState`。

MPI command 对 values 求和、对 maximum-length 字段取最大值、归约直方图，并收集
各 rank 的首尾坐标判断全局有序性。该边界保证 operator 负责 rank-local 计算，
command 负责跨 rank 语义和用户可见格式。

WES_0.25G、单节点 6 ranks 的迁移后 smoke 得到 `2,668,351` records、
NM mismatches=`878,007`，全部 SN 字段和 `is sorted=0` 与既有输入基线一致；
`4.3=0.072359 s`，低于历史 `0.099174 s`，确认本次迁移没有性能回退。

## BAM2BAM / Filter 的读写流水线组合

```text
BamInputBackend + spans / MemReader
    -> RunCpeReadPipeline (64-block double buffer)
       -> BamTransformOperator
          -> decode + parse + passthrough/filter CPE kernel
          -> MPE 按记录长度生成最多 64 个 BGZF 打包计划
          -> CpeRecordWriteSession::Submit
             -> CpeWritePipelineSession (output double buffer)
             -> record-serialize + compress CPE operator
             -> RankBodySink
```

读流水线支持 backend span 和借用 `MemReader` 的 batch 输入。后者直接读取已有内存
中的压缩块，不再复制完整 BAM。算子保留直通时的原块边界和过滤时的重新装块规则；
`CpeWritePipelineSession` 是唯一的 BAM 写调度器：已打包 BGZF 块使用
`ComposableCompressOperator`，记录打包计划使用 `RecordCompressOperator`；
`RunCpeWritePipeline` 是拉取式 source 的适配入口，`CpeRecordWriteSession` 是
`RankBodySink` 和记录计划的推送式适配入口。两者共享压缩输出双缓冲及同一套
`Prepare -> spawn -> flush(previous) -> join -> validate -> swap` 状态机。
两类写算子由 `CpeWriteKernelSpec` 描述入口、batch 容量与是否需要未压缩块 scratch；
从核共用 `swbam_cpe_compress_bam_block`。记录输入在同一次 CPE launch 中先序列化
再压缩，已打包块输入直接压缩；不会先在主核物化另一份完整 BAM。
写 session 在 `Submit()` 返回前完成当前批的压缩，上一批已压缩块可与压缩并行写出。
因此下一个 read batch 复用记录 arena 时，不会留下悬空的 `bam1_t` 指针。
这不表示 BAM read 与 BAM write 的 CPE kernel 同时执行。

命令层的 `OptimizedBamToBamMPI` 现在只适配原有四种输入/输出签名，并继续负责
header、BGZF EOF、MPI 输出布局、归约和打印。`t_read`、`t_decomp_filter`、
`t_pack`、`t_compress`、`t_write` 保留原字段；read 与 decode、write 与 compress
存在重叠，分项相加不能作为 `4.3` wall time。

WES_0.25G、单节点 6 ranks 的迁移后 smoke：无过滤 `2,668,351` 条，
`4.3=0.182463 s`；MAPQ >= 30 保留 `2,481,669` 条，`4.3=0.178360 s`。
两份输出都经 `stats --basic` 重新读取；分别与同输入的旧版输出
`io_refactor.bam`、`io_refactor_mapq30.bam` 做 `cmp`，均逐字节一致。
已有同输入旧版无过滤日志的 `4.3` 为 `0.181363 s`，单次差约 `0.6%`；
这只是一组节点运行结果，尚未作为稳定性能结论。

共享 BAM write 调度器改造后，同数据、同过滤参数的 6-rank `run_all` 保留
`2,481,669` 条，`4.3=0.179002 s`（改造前单次 `0.178360 s`），输出与
`io_refactor_mapq30.bam` 逐字节一致；`stats --basic` 可复读。
`io-check` 的合成块压缩、完整 BAM 回环和 MAPQ 过滤回环也全部通过。

## SAM2BAM 的流水线组合

```text
MemReader (rank-local SAM body)
    -> RunCpeSamReadPipeline
       -> SAM 行边界分块 -> CPE copy/count -> CPE parse
       -> SamToBamPostProcessor (MPE 生成 BGZF 记录打包计划)
       -> CpeRecordWriteSession -> 公共 BAM write 流水线 -> RankBodySink
```

SAM read 流水线拥有文本 chunk、解析后的记录池和数据 arena；解析结果只在同步
`PostProcessParsedBatch()` 回调期间有效。后处理只保存本批记录指针，提交写流水线
后等待 CPE 完成，再允许下一波解析复用记录池。写出压缩仍与上一批的 MPE flush
重叠；进入下一组 SAM chunk 的 CPE copy/count 时，也可 flush 最后待写的压缩批。
`SamReadKernelSpec` 描述已有的 copy/count 与 parse 两阶段 CPE 入口；新应用通常只
实现 `SamParsedBatchPostProcessor`，不改解析 kernel。自定义 kernel 必须遵守
`MpiSamParseBatch` 参数和 status/arena 契约，两阶段不能任意交换顺序。
命令层只保留 `MemWriter`/`RankBodySink` 适配、header 和最终 MPI 输出布局。

WES_0.25G、6 ranks、level 1 的迁移前后均处理 `2,668,351` 条，输出 BAM
逐字节一致；单次 `4.3` 从 `0.235816 s` 到 `0.235334 s`。这里只说明未观察到
明显回退，不将单次运行作为加速结论。

## BAM2SAM 的流水线组合

```text
BamInputBackend / MemReader
    -> RunCpeReadPipeline (BGZF 双缓冲 + passthrough decode/parse kernel)
    -> BamToSamOperator::PostProcessBatch (MPE 按 block 顺序收集记录指针)
    -> CpeSamWriteSession (CPE format + 双文本缓冲 + RankBodySink)
```

record arena 由 BAM read 算子持有，`Submit()` 完成格式化后才可在下一批解码时
复用。SAM write session 持有格式化文本，下一批 BAM read CPE kernel 运行时由
MPE flush 上批；最后一批由 `Finish()` 写出。命令层仅适配输入、输出及 MPI 布局。
这是 batch 级组合，不在逐记录格式化热路径加入虚调用。
该算子选择读流水线的可选调度：在解码期间 flush 上批 SAM，在本批 CPE
formatter 期间预读下一批 BGZF；其他 BAM read 算子继续使用默认预读调度。
`SamWriteKernelSpec` 描述格式化入口与记录容量。新应用可使用
`CpeSamWriteSession::SubmitRecords()` 批量提交指针；bam2sam 继续直接填写
`RecordSlots()`，避免额外指针复制。两种入口都要求记录在本次 `Submit()` 返回前
有效，文本缓冲由 session 持有直到下一批 flush 或 `Finish()`。

WES_0.25G、6 ranks、memory backend：旧版与迁移版均输出 `2,668,351` 条，
最终 SAM 逐字节一致；`4.3` 单次为旧版 `0.164583 s`、迁移版 `0.162327 s`。
该结果仅用于确认本轮未出现明显性能回退，不代表稳定加速比。

同一输入的 POSIX backend smoke 得到 `2,668,351` 条，`stats --basic` 复读正常。
由于命令行产生的输出 header 长度不同，完整 BAM 与 memory backend 输出不能
直接逐字节比较；跳过各自 header 后，`275,206,630` 字节的压缩 body 逐字节一致。
