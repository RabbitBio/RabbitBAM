# SWBAM SDK 构建与使用说明

## 1. 当前定位

SWBAM 当前既是一个面向神威异构平台的 BAM 命令行工具集合，也是一个可在项目内
复用的 C++11 SDK。SDK 已经把输入输出、MPI 分区、MPE/CPE 批处理流水线、通用
BGZF 编解码和 Raw BAM 记录处理拆成独立 CMake 目标，并用两个示例证明可以在不依赖
RabbitBAM 命令分发代码的情况下构建新工具。

目前它属于“项目内 SDK”，还不是可独立安装的系统软件包：尚未提供 `install()`、
导出的 `SWBAMConfig.cmake`、版本化 ABI 或 `find_package(SWBAM)`。因此这里的“库”指
CMake 构建目标和静态库，而不是 DEB/RPM 包或已经安装到系统目录的开发包。

## 2. 编译后得到什么

Sunway 配置下的主要产物如下：

| 类型 | CMake 目标 | 典型产物 | 作用 |
| --- | --- | --- | --- |
| 主库 | `swbam_io` / `swbam::io` | `libswbam_io.a` | BAM/SAM 输入后端、BGZF block 扫描与批读取、顺序输出、rank body 内存或 spool 存储 |
| 主库 | `swbam_mpi_runtime` / `swbam::mpi_runtime` | `libswbam_mpi_runtime.a` | MPI-IO 输入、全局 block 分区、rank 状态规约、分布式 BAM 布局和写出 |
| 主库 | `swbam_cpe_runtime` / `swbam::cpe_runtime` | `libswbam_cpe_runtime.a` | MPE 侧 CPE 读写流水线、双缓冲、通用解压/压缩、Raw BAM 扫描、过滤和 writer |
| 聚合目标 | `swbam_sdk` / `swbam::swbam` | 无独立文件 | 聚合三个主库、通用 CPE objects、MPI/HTSlib/libdeflate 等依赖和神威 hybrid 链接选项 |
| 通用 CPE 目标 | `swbam_cpe_kernels` / `swbam::cpe_kernels` | INTERFACE，无独立 archive | 把通用 codec、parser、统计和过滤从核对象直接交给 hybrid linker |
| 算法 CPE 目标 | `swbam_algorithm_kernels` / `swbam::algorithm_kernels` | INTERFACE，无独立 archive | 提供 sort、collate、fixmate、markdup 的专用从核 kernel，仅主程序需要 |

神威 hybrid linker 要求 CPE object files 直接参与最终链接，所以
`swbam_cpe_kernels` 和 `swbam_algorithm_kernels` 不能按普通主核静态库处理。SDK 使用者
只需链接 `swbam::swbam`，不需要手工列出这些从核对象。

主要可执行程序包括：

| 产物 | 用途 |
| --- | --- |
| `RabbitBAM-MPI` | MPI 主程序，承载转换、统计、排序、collate、fixmate、markdup 和去重流程 |
| `RabbitBAM-X` | 原有主程序 |
| `RabbitBAM-CGS` | 原有 CGS 转换实现 |
| `swbam-sdk-record-count` | SDK 只读示例 |
| `swbam-sdk-filter-bam` | SDK 读、过滤、重写 BAM 示例 |
| `RabbitBAM-CRC32-Bench`、`markdup_verify_fast` | 性能微基准与结果验证工具 |

构建 SDK 示例：

```bash
cmake -S . -B build_sunway_sduhpc -DPLATFORM=sunway
cmake --build build_sunway_sduhpc \
  --target swbam-sdk-record-count swbam-sdk-filter-bam -j 8
```

## 3. 代码与库的对应关系

```text
include/swbam/          对外头文件
lib/io/                 swbam_io 实现
lib/mpi/                swbam_mpi_runtime 实现
lib/cpe/                CPE 流水线的 MPE 侧实现
lib/bam/                Raw BAM adapter、filter 和 writer
slave/core/             通用 CPE BGZF codec 与 BAM parser
slave/operators/        通用 decode/compress/filter/statistics kernel
slave/algorithms/       sort/collate/fixmate/markdup 专用 kernel
apps/mpi/               命令、复杂算法和 pipeline 实例
examples/               独立于命令内部类型的 SDK 示例
```

这一组织表达了两条边界：通用 BAM 数据面进入 SDK；全局排序、候选交换、重复判断等
算法语义留在应用实例中。专用 kernel 可以复用库的 I/O 和运行时，但不会被强行改写为
较慢的逐记录通用回调。

## 4. 对外公开的接口

总入口为：

```cpp
#include <swbam/swbam.h>
```

该头文件汇总下面四组接口。也可以只包含需要的细分头文件。

### 4.1 输入输出与 BGZF batch

定义在 `swbam/io.h`：

- `BamInputBackend`：BAM 输入抽象，暴露 header、body offset、block scan 和 batch read。
- `MemoryBamInput`：完整 BAM 常驻内存，作为当前性能模拟路径。
- `PosixBamInput`：按 BGZF span 批量 `preadv`，不要求完整输入常驻内存。
- `SamInputBackend`、`PosixSamInput`：SAM header、按换行切分 range 和局部读取。
- `BgzfBlockSpan`：一个 BGZF block 的文件偏移与压缩长度。
- `BgzfBlockBatch`、`BgzfSpanBatchReader`：对齐批缓冲及顺序批读取。
- `BamOutputBackend`：顺序输出抽象。
- `MemoryBamOutput`、`PosixBamOutput`：内存输出和真实 POSIX 文件输出。
- `RankBodySink`：保存每个 MPI rank 生成的压缩 body。
- `MemoryRankBodySink`、`SpoolRankBodySink`、`AdaptiveRankBodySink`：内存、临时文件及自动 spill 策略。

### 4.2 通用 CPE 流水线

定义在 `swbam/cpe_pipeline.h`、`generic_decode.h`、
`cpe_write_pipeline.h` 和 `generic_compress.h`：

- `RunCpeReadPipeline`：读取下一批与 CPE 处理当前批重叠。
- `CpeBatchOperator`：批级专用算子扩展点，不在逐记录热路径增加虚调用。
- `RunGenericDecodePipeline`：BGZF 解压后把 decoded block batch 交给 consumer。
- `RunCpeWritePipeline`：未压缩 block 生产、CPE 压缩和输出消费流水。
- `RunGenericCompressPipeline`：通用 BGZF 压缩实现。
- 对应 timing/metrics 结构用于拆分 read、kernel、consume、inflate、CRC、deflate 等时间。

### 4.3 Raw BAM 记录接口

定义在 `swbam/raw_bam.h`、`raw_bam_filter.h` 和 `raw_bam_writer.h`：

- `RawBamRecordView`：指向 decoded batch 中原始 BAM record 的零拷贝只读视图。
- `RawBamRecordConsumer`：应用实现的批量记录回调。
- `RunGenericRawBamPipeline`：通用解压、记录边界解析和 consumer 调用。
- `RawBamFilterConsumer`：读取 mandatory fields 并完成常用过滤。
- `RawBamWriter`：将 Raw BAM records 重新分块、CPE 压缩并写入 output backend。

`RawBamRecordView` 的内存只在当前 `ConsumeRaw()` 回调期间有效，不能保存指针供后续
batch 使用。当前实现还沿用项目的输入约束，即单条 BAM record 不跨 BGZF block。

### 4.4 MPI 运行时

定义在 `swbam/mpi_runtime.h`：

- `MpiBamInput`：统一选择 `memory|posix|mpiio|auto` 输入后端。
- `MpiBamInputPlan`、`PrepareMpiBamInputPlan`：rank0 扫描并广播 block plan，再划分 rank 范围。
- `MpiIoBamInput`、`MpiIoSamInput`：MPI-IO 输入。
- `DistributedBamOutput`：汇总 rank body size，计算全局 offset，并选择 gather-to-root 或 MPI-IO 写出。
- `MpiFileOutput`：大于 `INT_MAX` 时仍可分块执行 MPI offset 写入。
- `AllRanksOk`、`ReduceMaxCost`：统一错误传播和最大 wall time 规约。

## 5. 如何使用

### 5.1 只读统计工具

典型流程是：

```text
PosixBamInput::Open
  -> ScanBlocks
  -> 自定义 RawBamRecordConsumer
  -> RunGenericRawBamPipeline
  -> Close
```

`examples/sdk_record_count.cpp` 展示了完整实现。使用者只需要在
`ConsumeRaw()` 中读取所需 BAM mandatory fields；公共库负责 BGZF 扫描、批读取、
CPE 解压、记录边界识别和双缓冲。

### 5.2 过滤并生成 BAM

典型流程是：

```text
PosixBamInput
  -> RunGenericRawBamPipeline
  -> RawBamFilterConsumer
  -> RawBamWriter
  -> PosixBamOutput
```

对应代码见 `examples/sdk_filter_bam.cpp`。`RawBamWriter::InitializeBam()` 负责序列化
header，`Finish()` 刷新最后一个 payload 并追加 BGZF EOF。业务代码主要负责设置
`BamFilterOptions` 和选择 downstream consumer。

### 5.3 链接方式

项目内新增程序只需：

```cmake
add_executable(my_swbam_tool my_swbam_tool.cpp)
target_compile_options(my_swbam_tool PRIVATE ${HOST_FLAGS})
target_link_libraries(my_swbam_tool PRIVATE swbam::swbam)
```

当前示例仍显式调用 `athread_init()`/`athread_halt()`。MPI 工具还需负责
`MPI_Init()`/`MPI_Finalize()`，再使用 `swbam::mpi` 接口。由于尚未导出 CMake package，
独立仓库目前不能直接写 `find_package(SWBAM)`；最稳妥的使用方式仍是在 RabbitBAM
构建树中增加 target。

## 6. 已有实例

### 6.1 SDK 示例

- `sdk_record_count`：通用 Raw BAM 扫描，统计 records、mapped 和 duplicates。
- `sdk_filter_bam`：通用 Raw BAM 过滤及完整 BAM 重写。

### 6.2 RabbitBAM-MPI 命令实例

- 线性扫描：`flagstat`、`stats --basic`、BAM/SAM 转换和过滤。
- 全局重排：`sort`、`collate`。
- 有状态处理：`fixmate -m`、默认或流式 `markdup`。
- 组合流程：`dedup-pipeline`。
- 库微基准及完整性检查：`io-check`。

这些命令并不都通过 `RawBamRecordConsumer` 实现。`flagstat`、`stats`、转换、sort、
collate、fixmate 和 markdup 保留了针对数据流融合的专用 CPE kernel，同时复用公共
backend、BGZF batch、双缓冲调度、rank body sink 和分布式输出。这是通用性与性能
并存的设计，而不是两套互不相关的实现。

## 7. 与 SAMtools 的对齐程度

### 7.1 已对齐的层面

- 命令概念对齐：提供转换/过滤、`flagstat`、`stats`、`sort`、`collate`、
  `fixmate` 和 `markdup` 等典型 BAM 工作流。
- 核心输出语义对齐：已经针对记录统计、name continuity、coordinate sorted、
  mandatory fields 和 duplicate key/flag 集合做过与 SAMtools 的验证。
- 大文件策略方向对齐：线性工具支持流式输入；sort/collate 提供算法外排；rank body
  可通过 spool 控制内存；最终结果可通过 MPI-IO 分布式写出。

### 7.2 尚未对齐的层面

- **不是 SAMtools CLI 的直接替代品**：当前通常要求 `-i/-o`，部分选项名称不同，
  压缩等级只接受 `0/1/6`，`stats` 目前只实现 `--basic`。
- **不是 HTSlib API 兼容层**：SWBAM 提供批级 C++11 接口，而不是
  `sam_open/sam_read1/sam_write1` 的同名替代。
- 格式覆盖不完整：重点是 BAM；CRAM、索引、region iterator、远程 URL、参考序列
  管理和 SAMtools 的大量边缘选项没有完整实现。
- `fixmate`、`markdup`、collate 等只覆盖当前项目验证过的主要选项和输入约束，不能
  宣称与 SAMtools 全参数等价。
- `dedup-pipeline` 是 SWBAM 自有组合命令，不是 SAMtools 单命令接口。

因此更准确的表述是：**SWBAM 在论文涉及的核心 BAM 工作流和主要结果语义上与
SAMtools 对齐，但 CLI、API、格式和选项覆盖仍是一个面向神威的性能子集。**

## 8. 易用性评价

当前易用性可评价为“项目内较好，外部 SDK 尚可但仍需收口”：

### 已经易用的部分

- 普通应用只链接一个目标 `swbam::swbam`，无需了解从核 object 的链接规则。
- `swbam/swbam.h` 提供单一总入口。
- 输入、输出、过滤、Raw BAM writer 可以组合，两个示例分别覆盖只读和读写工具。
- backend 和算法解耦，同一 consumer 可以使用 memory 或 POSIX 输入。
- 扩展发生在 batch/consumer 层，开发新统计或简单过滤无需复制 BGZF/CPE 调度代码。

### 当前使用门槛

- 没有 `install()` 和 `find_package()`，外部项目接入仍不够自然。
- 使用者需要手动管理 `athread`，MPI 应用还要管理 MPI 生命周期。
- `RawBamRecordView` 是低层编码视图，复杂字段访问缺少完善的类型安全 accessor。
- 公共头文件仍暴露 `BamTools.h`、HTSlib `sam_hdr_t`、MPI 类型和 `bam_block`，接口
  与内部依赖及神威平台耦合较强。
- 错误接口主要是整数返回值和 stderr，尚无统一错误对象或诊断上下文。
- 缺少独立 SDK API reference、版本策略和外部工程示例。

## 9. 推荐的下一步

若目标是论文中的“可复用高性能 BAM 库”，当前结构已经足以支撑架构与应用案例。
若目标是让第三方真正方便地使用，建议按以下顺序收尾：

1. 增加 `install()`、目标 export 和 `SWBAMConfig.cmake`，支持
   `find_package(SWBAM CONFIG REQUIRED)`。
2. 增加一个真正位于独立目录、只通过安装包编译的 consumer 示例，验证依赖没有
   从源码树泄漏。
3. 提供 `RuntimeGuard`/`MpiRuntimeGuard`，隐藏 `athread` 和 MPI 生命周期样板代码。
4. 为 FLAG、MAPQ、坐标、QNAME 和 aux tag 增加只读 accessor，避免 SDK 用户手工按
   byte offset 解析 mandatory fields。
5. 将低层 CPE codec/parser 头文件标注为 advanced API，并为高层 API 建立稳定性和
   版本策略。

完成前两项后，SWBAM 才适合称为“可安装、可被外部 CMake 项目直接消费的 SDK”；
后面三项主要提升开发体验，不影响现有专用融合命令的性能路径。
