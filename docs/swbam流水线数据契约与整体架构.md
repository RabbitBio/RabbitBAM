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

当前四条类型化流水线为：

| 流水线 | 输入 | 主要 operator | 输出 |
| --- | --- | --- | --- |
| BAM read | BGZF spans | decode / decode+parse+业务 | 结果 slice、raw view 或算法元数据 |
| BAM write | raw/record batch | pack+compress | compressed BGZF batch |
| SAM read | 文本 range/chunk | split+parse | BAM record 或算法元数据 |
| SAM write | BAM record batch | format | SAM text batch |

应用负责 CLI、MPI 全局划分/归约和结果展示；库负责 backend、batch、调度、operator
及输出能力。复杂算法的全局状态属于应用层。sort、collate 和 fixmate 的读写两端
已接入公共 BAM 流水线；markdup 的单遍 `--stream` 主路径也已接入，两遍 fallback
的第二遍 rewrite/write 仍保留原调度。

## Batch 数据契约

### 所有权与生命周期

- backend 拥有文件或内存输入；pipeline 只在运行期间借用 backend。
- pipeline 拥有 compressed/decoded 双缓冲，operator 不得释放或长期保存其地址。
- batch view 只在当前回调期间有效；跨 batch 保存记录必须显式复制。
- 算法也可在 launch 前绑定自己拥有的持久 arena，例如 collate 内排直接解压到
  最终 raw 位置；此时持久的是 app arena，不是借用的 pipeline decoded buffer。
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
CPE kernel 中静态融合，例如 `decode+parse+flagstat`。算子可输出 batch-local
`bam1_t`、解压块或算法元数据供 MPE 后处理；旧 `RawBamRecordView`/MPE 物化的
高层 adapter 已退役，不再作为另一套 SDK 实现。

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
`BgzfCompressOperator`，记录打包计划使用 `RecordCompressOperator`；
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

## Sort 的读写端点

```text
MemReader / BamInputBackend + spans
  -> RunCpeReadPipeline + MpiSortExtractOperator
     -> 原 CPE raw/key 提取 kernel
     -> app 持久 raw/metadata，或严格外排 run arena
  -> app 排序、采样、MPI 交换、run/segment 管理、归并
  -> UncompressedBgzfSource（app 按最终顺序打包完整 raw records）
  -> RunCpeWritePipeline + MpiSortCompressOperator
     -> 原 CPE payload 压缩 kernel
     -> RankBodySink -> DistributedBamOutput
```

这是两个阶段的组合：读取阶段完成后，app 的持久数据仍然有效；排序/交换完成后
才启动最终写出。内排和严格外排使用同一批次提取适配，区别在批后处理如何保存数据。
排序键提取保留轻量 raw 路径，不强制物化 `bam1_t`，也不在逐记录热路径增加虚调用。

BAM write 新增 `CpeWritePipelineOptions::overlap_source`，默认 false。
sort 显式开启后，调度顺序为：

```text
Prepare -> spawn compress(current)
        -> flush(previous) -> Fill(next) [MPE]
        -> join -> validate -> swap
```

推送式调用可使用 `SubmitWithPrefetch(input, callback, context)`。回调只使用下一批
独立存储，不能启动 CPE kernel；当前 input 在 Submit 返回前始终有效。回调返回失败
时也要完成 join 后再释放输入。最终 `Finish()` 排空最后一批压缩输出。
这保留了 sort 原有“压缩 + flush + 归并”的覆盖关系，不代表两个 CPE kernel 并行。

外排预算现在覆盖两套 compressed、两套 decoded batch 和提取 metadata；约增加
8 MiB/rank 的读 workspace 会从 run arena 可用预算扣除。写端仍为两套 payload 加
两套压缩结果。外排模拟时间沿用实际临时 I/O 扣除与压缩周期估算模型；`source`、
`kernel`、`post_process` 有重叠，不能相加当作核心 wall time。

迁移验证见 [sort 流水线迁移验证](results/sort_流水线迁移验证.md)。

## Collate 的读写端点

```text
MemReader / BamInputBackend + spans
  -> RunCpeReadPipeline + CollateExtractOperator
     -> 原 CPE 解压、QNAME/hash/meta 提取 kernel
     -> 内排：app raw arena 直写 + metadata append
        外排：decoded batch -> app run arena -> spill / resident run
  -> app 按 bin/hash 排序、MPI ring exchange、segments、loser tree
  -> CollatePayloadSource（按归并次序打包完整记录）
  -> RunCpeWritePipeline + CollateCompressOperator
  -> RankBodySink -> 分布式输出
```

公共库负责 batch 生命周期、spawn/join、校验后回调、buffer 轮换和写端覆盖。
collate app 保留 owner/bin、QNAME 分组次序、持久 raw/meta、外排文件及归并状态。
`CollateExtractOperator` 继承公共 `BamReadBatchOperator`，复用
`BindBamReadBatch`；读端的两种后处理按算法模式选择，不另写一套 I/O 调度。

内排 operator 使用自己的输出描述符绑定最终 raw arena，不改动库拥有的 decoded
arena 指针，避免一次解压后 raw 拷贝；remote receive 直接写最终 segment 的优化
也保留。外排 decoded 数据必须在批回调返回前追加到 run arena，不能长期保存
pipeline buffer 的地址。全局稳定序号仍由 `global_block_begin + local_block + b`
和 block 内记录号生成。

写端复用 sort 已接入的 `overlap_source=true`，保持 `compress(current)` 与
`flush(previous) + Fill(next)` 的覆盖。内排 fill 直接调用 loser tree，外排 fill
保留原 cursor 的延迟 advance 和一条记录 lookahead；batch 切换不改变记录次序。
这仍是“读端 -> app 全局算法 -> 写端”的两阶段组合，不是边读边输出已分组结果。

本轮读端 `overlap_input=false`，不把新增输入覆盖混入架构迁移；两套 compressed、
两套 decoded buffer 已进入外排 workspace 预算。内排虽直接使用持久 raw arena，
公共读端仍分配 decoded 双缓冲，本轮未扩展自定义 arena 分配接口。

核心 `4.3`、`4.1~4.6` 计时边界保持不变；`compress_pipeline` 含被覆盖的 MPE
flush/fill，不能把它与 `compress_fill`、`write` 相加。`extract_alloc/free` 现在只
统计算子 metadata 池，公共 batch 的分配释放仍包含在核心总时间内。
外排保留 `actual - temp_read_actual - temp_write_actual` 的模拟时间模型，模拟
memcpy 已在 actual 内，不能重复相加。旧 `accounted/unaccounted` 分项未统一扣除
真实临时 I/O，外排会出现负 unaccounted，不能将其解释为漏计的 CPU 时间。

迁移验证见 [collate 流水线迁移验证](results/collate_流水线迁移验证.md)。

## Fixmate 的读写端点

```text
MemReader / BamInputBackend + spans
  -> RunCpeReadPipeline + FmReadOperator
     -> passthrough CPE 解压、解析为批量 bam1_t
     -> PostProcessBatch：app QNAME 分组、补齐 pending group
        -> 原 fixmate plan CPE kernel
        -> MPE 计算输出 offsets / 分配 rewrite arena
        -> fixmate rewrite CPE kernel（允许 workspace 对齐间隔）
        -> app 生成 BamRecordPackBlock
        -> CpeRecordWriteSession -> 公共 BAM write -> middle sink
  -> app MPI 边界交换 + 同一组处理/写出逻辑 -> prefix/suffix sink
  -> prefix + middle + suffix -> 分布式输出
```

fixmate 是状态型实例：公共读流水线管理 BGZF batch 和 decode 调度，app 的
`FmReadOperator` 管理解析记录池和未完成的 QNAME group。跨 batch/rank 的 group
复制 core/data 到 `FmStoredGroup`；其余记录只在当前批回调内使用，不缓存完整
rank 的解压结果。极端情况下单个 QNAME group 仍可能增长，本轮没有改变该边界。

`Validate()` 先检查全部 active block，才进入 app 批后处理。read kernel 已 join，
此时可以依次启动配对 plan、rewrite 和 BAM write kernel；本轮不引入 CPE kernel
并发。读端设置 `overlap_input=false`，保留原 read/decompress 不重叠的计时口径。
内存 batch 适配器将正常 EOF 与读取失败分开处理，坏块不再作为正常 EOF 返回。

写端不新建 fixmate 专用压缩 operator，直接复用已有 `CpeRecordWriteSession`、
`RecordCompressOperator` 和 `CpeWritePipelineSession`。app 仍负责装块计划，
同一次 CPE launch 完成序列化和压缩；MPE flush 上一压缩批与本批压缩重叠。
每次 `FmCompressPlans()` 都完成 Finish，再释放被借用的 rewrite 记录和 data arena。
这保留原 group 边界上的装块方式，没有先把整个中间 BAM 复制回 MPE 打包。

prefix/middle/suffix 独立写入，最后仍由 app 按顺序组合；库不处理 mate 关系或决定
跨 rank group 的 owner。standalone 和 pipeline memory helper 复用同一实现，
但完整 dedup-pipeline 的跨阶段衔接不属于本轮迁移。

公共读端增加约 8 MiB/rank 固定 buffer；记录池、写端双缓冲和 I/O backend 策略
不变。`4.3` 和 `4.1~4.6` 核心范围保持原样，压缩和写出分项仍有重叠，不能重复
加总。默认 memory backend 的整文件 load、最终 gather/dump 继续在核心范围外。

本次逐字节回归发现一条跨批 QNAME group 中的 `ms` 标签尾字节不同。原 rewrite
紧密拼接不同 CPE 的输出记录，非对齐拷贝存在相邻字覆盖风险。现将各非空 CPE
输出区起点按 64 字节对齐，边界容量校验允许间隔，但真实记录长度仍取 plan 的
`output_data_len`。最多增加 63 个不足 64 字节的间隔，不改变输出顺序或 BAM 编码。
修正后完整 BAM 与旧版基线 `cmp` 一致；小数据单次 4.3 为 0.507465 -> 0.518867 s，
约 +2.25%，未进行多次统计或所有 backend 回归。

迁移验证见 [fixmate 流水线迁移验证](results/fixmate_流水线迁移验证.md)。

## 流式 Markdup 的读写端点

```text
MdInputPass (MemReader / backend spans)
  -> RunCpeReadPipeline + MdCandidateReadOperator
     -> passthrough decode/parse
     -> app candidate extract CPE / candidate merge
  -> app owner exchange / sliding window / boundary halo
  -> app decoded-block pending ring / final decisions / flag patch
  -> CpeWritePipelineSession + MdOnePassCompressOperator
     -> 原 payload compress CPE -> MdOnePassOutput -> RankBodySink
```

不是把一个完整 markdup 功能塞入库：公共库管理 batch 调度、校验和写缓冲，app
管理候选、坐标窗口、MPI 决策及块何时能够退休。candidate 提取仍是 app 的独立
CPE 阶段，不是 MPE 逐记录虚调用，也不在本轮合并或改写判重 kernel。

### 保持原有覆盖

读端保存两槽记录和候选池：当前槽 decode 运行时，`DuringKernel()` 只归并上一槽
已完成提取的 candidate。当前 decode join/validate 后，才处理上一槽的窗口和
输出，再启动当前槽的 candidate extract。最后一槽在 `Finish()` 归并、后处理。
全局进度、ordinal 和各 rank 补齐到 `target_batches` 的空回调保持原语义。

写端新增受控 split 提交接口，复用原写调度器而非新建第五条流水线：

```text
Start(ready payload)   -> CPE 压缩已确定记录
app MPE/MPI 工作      -> 当前候选交换、窗口决策
Complete()            -> join / validate
Flush()               -> 写入 sink
app retire/recycle    -> 回收已完成块
```

Start 到 Complete 之间不能启动其他 CPE kernel，也不能修改/释放所提交的 payload、
描述符或参数。`Submit()` 和 `SubmitWithPrefetch()` 仍是同步接口，内部复用此逻辑。
失败清理时 session 先 join，再释放 pending arena；析构不代替成功路径 Finish/Flush。
markdup 本轮仍在 Complete 后立即 Flush，保留原来的“压缩覆盖窗口计算”，没有声称
额外实现 flush 与后续压缩的覆盖；其他 writer 的双缓冲 flush 策略不变。

### 所有权与计时

pending 环保存 app 自己的 decoded arena，批回调通过交换 arena 所有权留住原始字节，
不是复制整个 raw BAM，也不是窃取库的 decoded buffer。原有宽块 fallback 仍按完整
记录拆块，不扩大支持范围到跨 BGZF record。

公共读端另外分配了两套 decoded buffer，本轮算子实际用自有可转移 arena；加上公共
写端多一套输出 buffer，固定 workspace 增加约 12 MiB/rank，并纳入原 `-m` 检查。
后续可以在公共接口支持算子管理 decoded 存储后消除这部分冗余，无需修改窗口算法。

核心 4.3、4.1~4.6 的边界不变，memory load/dump 仍排除。`candidate_decomp` 现在取
读 runtime kernel wall 减掉被覆盖的 candidate merge，含少量 launch/调度开销；
`compress` 仍取 Complete 的暴露等待和校验时间，不是完整压缩 wall。
公共 write 的 `timing.kernel` 是 Start 到 Complete 的区间，包含被覆盖的窗口/通信，
不能与 app 分项简单相加。`tracked_peak` 是工作区估计，不是 RSS。

此轮未迁移第二遍 rewrite/write、首坐标预探测或 pipeline 跨阶段编排，也未改
`swbam_markdup_stream_mpi.cpp`。共享 candidate 入口同时影响两遍和非流式路径，
验证覆盖默认单遍、两遍 `-r` 及公共同步 writer 的 BAM2BAM smoke。

结果见 [markdup 流水线迁移验证](results/markdup_流水线迁移验证.md)。

2026-09-27 边界回归补充：candidate 批回调返回后，extract 清空候选及 QNAME 容器。
需跨回调保存的 halo probe 数据必须在回调内同时转移这两个容器的所有权，不能只
保留坐标后假定返回参数仍有内容。此次修复用 swap 留存并校验候选数量；非 probe
范围的批次跳过预判查表。完整串联输入的 np=6 输出已与同输入 np=1 BAM 逐字节
一致，详见 [去重串联与边界回归](results/去重串联_边界回归验证.md)。
