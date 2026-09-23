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
