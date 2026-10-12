# SWBAM SDK 使用说明

当前 SDK 以四条流水线和可组合 operator 为入口，独立项目位于 `swbam-sdk/`，
使用 `find_package(SWBAM CONFIG REQUIRED)` 后链接统一目标 `swbam::swbam`。
旧 Raw view/MPE bam1_t adapter 与旧五个示例已退役，不再维护两套
示例架构；CPE BAM parser 对 `bam1_t` 的解析能力保留。

## 构建产物

| 目标 | 用途 |
| --- | --- |
| `swbam_io` / `libswbam_io.a` | 输入/输出后端、BGZF batch、rank body storage |
| `swbam_mpi_runtime` / `libswbam_mpi_runtime.a` | MPI input plan、MPI-IO、分布式输出 |
| `swbam_cpe_runtime` / `libswbam_cpe_runtime.a` | MPE 侧四条流水线与算子适配 |
| `swbam_cpe_kernels` | 共享 CPE object 的接口目标，不是普通 `.a` |
| `swbam::host` | MPE 编译选项、公共头文件和布局配置 |
| `swbam::cpe` | 自定义 CPE kernel 编译选项、公共头文件和 HTSlib 类型 |
| `swbam::swbam` | SDK 总链接目标，传递头文件、库、MPI 依赖和 hybrid 链接选项 |

应用算法 object 位于 `swbam-app/slave/`，由 app 自己构建，不随 SDK 安装。
三个静态库有对应的 `swbam::io`、`swbam::mpi_runtime`、`swbam::cpe_runtime`
导出目标；通常直接使用 `swbam::swbam` 即可。

## 独立安装与使用

从仓库根目录构建并安装库：

```bash
cmake -S swbam-sdk -B build_sunway_sduhpc/sdk \
  -DBUILD_SWBAM_EXAMPLES=OFF \
  -DCMAKE_INSTALL_LIBDIR=lib \
  -DCMAKE_INSTALL_PREFIX="$PWD/build_sunway_sduhpc/sdk-install"
cmake --build build_sunway_sduhpc/sdk -j 8
cmake --install build_sunway_sduhpc/sdk
```

SDK 源码可单独复制出来。此时配置 `-DSWBAM_DEPENDENCY_ROOT=/path/to/ext`，其中
应已有神威构建的 `htslib-1.20`、`libdeflate-1.20`；也可分别指定
`SWBAM_HTSLIB_ROOT` 和 `SWBAM_LIBDEFLATE_ROOT`。库构建不访问 app 或 tools。

安装目录内容为：

```text
include/swbam/         SDK 公共接口、MPE/CPE 共享类型和 CPE 编解码声明
include/htslib/        同版本 HTSlib 公共头文件
include/libdeflate.h   主核 libdeflate 头文件
lib/libswbam_*.a       三个 MPE 静态库
lib/libhts.a           适配的主核 HTSlib
lib/libdeflate.a       适配的主核 libdeflate
lib/swbam/slave/         通用 CPE object，按原方式参与 hybrid 链接
lib/cmake/SWBAM/       Config、版本文件、导出 targets 和 Sunway 工具链
share/SWBAM/licenses/  随包依赖的许可证
```

MPI、athread、zlib 和系统库由 Sunway 环境提供。`SUNWAY_MPI_ROOT` 默认沿用构建时
的位置，可在使用安装包时重新指定。安装前缀可移动；CMake 导出文件不引用 SDK
源码或原构建目录。二进制对象仍要求兼容的神威工具链、硬件和主从核 ABI。

外部程序的 CMake 文件可以写为：

```cmake
cmake_minimum_required(VERSION 3.13)
project(MyBamTool LANGUAGES C CXX)
find_package(SWBAM 0.1 CONFIG REQUIRED)
add_executable(my_tool my_tool.cpp)
target_link_libraries(my_tool PRIVATE swbam::swbam)
```

配置时指定安装位置和工具链：

```bash
cmake -S /path/to/my-tool -B /path/to/my-build \
  -DCMAKE_PREFIX_PATH=/path/to/sdk-install \
  -DCMAKE_TOOLCHAIN_FILE=/path/to/sdk-install/lib/cmake/SWBAM/SunwayToolchain.cmake
cmake --build /path/to/my-build -j 8
```

`swbam-app/` 和 `swbam-sdk/examples/` 均可用同样方式单独构建。根目录的联合构建
入口继续可用，便于日常开发；与独立构建共享同一份源码和目标定义。

程序负责 MPI 初始化和 athread 生命周期。导出目标自动带入 CPE object 和
`-mhybrid`，不能自行换成普通主核 archive 链接。安装包版本目前为 `0.1.0`。

## 应用自定义 CPE 算子

app 中间算法需要专用 kernel 时，可链接 SDK 的编译接口：

```cmake
add_library(my_kernel OBJECT my_kernel.cpp)
target_link_libraries(my_kernel PRIVATE swbam::cpe)
target_sources(my_tool PRIVATE $<TARGET_OBJECTS:my_kernel>)
```

这样 MPE 调度仍链接 `swbam::swbam`，CPE 编译使用 `swbam::cpe`；应用自己的
MPE/CPE 参数结构由应用头文件定义。参数字段、对齐与类型必须在两端一致。
库侧保留 `cpe_bam_read_steps.h`、`cpe_bam_write_steps.h` 等可组合步骤，用户无需
包含 app 的头文件，也不需要访问 SDK 的私有实现目录。

## 两个示例

### 只读计数

`swbam-sdk/examples/sdk_optimized_count.cpp`：

```text
MpiBamInput(posix or borrowed RAM) -> PrepareMpiBamInputPlan
    -> RunRecordCountPipeline -> RunCpeReadPipeline
    -> CPE decode + parse bam1_t + count
    -> MPI sum(records), max(core time)
```

命令体：`./swbam-sdk-optimized-count input.bam`。只计数，不运行完整 flagstat。
输出 `sdk_count ranks=... records=... blocks=... input_bytes=... core=...`。

### BAM 读写回环

`swbam-sdk/examples/sdk_optimized_bam2bam.cpp`：

```text
MpiBamInput(posix or borrowed RAM) -> PrepareMpiBamInputPlan
    -> RunBamTransformPipeline (no-op filter, level 1)
    -> BAM read -> pack plan -> BAM write -> rank body storage
    -> MPI sum(records/body bytes), max(core time)
```

命令体：`./swbam-sdk-optimized-bam2bam input.bam -o output.bam`。输出包含完整
header、body 和 EOF。默认走真实磁盘输入和 MPI-IO 输出；RAM 模式的实际保存发生在
MPI 结束后。两示例的 `core` 只作为流水线诊断，不代表新端到端成绩。

## RAM 存储与端到端计时

`benchmark-support/` 是 app 和示例共用的私有使用端模块，不安装进 SDK，也不被
SDK 库反向依赖。独立构建可指定
`-DSWBAM_BENCHMARK_SUPPORT_DIR=/path/to/benchmark-support`；设为空或该目录不存在时
仍可构建、运行磁盘模式，RAM 参数会明确报错。

```bash
./swbam-sdk-optimized-count input.bam \
  --storage ram --ram-capacity 8G --ram-ranks 6 \
  --ram-session /absolute/unique-session/count
./swbam-sdk-optimized-bam2bam input.bam -o output.bam \
  --storage ram --ram-capacity 8G --ram-ranks 6 \
  --ram-session /absolute/unique-session/copy
```

以上命令体必须由匹配的 MPI 启动器启动。例如六节点配置为
`-N 6 -np 1 -mpecg 6 -cgsp 64 -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc`。
`--ram-ranks` 必须与实际 rank 数一致；会话目录的父目录须存在，本次目录必须不存在。
默认就绪超时为 7200 秒，可用 `--ram-ready-timeout` 修改。普通 `malloc` 配合
`-xmalloc`，禁止直接在 MPI 初始化前调用 `_sw_xmalloc`。

初始化前创建并预触碰固定容量池、装入本 rank 的完整输入分片及原始头/扫描信息、
用文件标记同步。MPI 内重新解析 header、验证扫描信息并分区，运行原有流水线。
存储池按 4 MiB 页管理；临时文件实际保存数据，支持位置读写、覆盖、追加、EOF
和关闭回收。`--ram-capacity` 为每 rank 存储总预算，包含输入、临时文件和完整输出，
**不包含**算法 `-m` 工作区。容量不足明确报错，绝不落盘或循环覆盖输出；外排中
容量失败会终止 MPI 作业，避免其他 rank 留在通信中等待。

主成绩取所有 `[swbam-e2e rank=...] seconds=...` 的最大值，边界为 MPI_Init 前
到 MPI_Finalize 后，包含 athread、正常清理、实际 RAM I/O 及通信。预加载和
Finalize 后保存不计入；不扣减真实临时 I/O、不累加阶段最大值。
默认各 rank 按已确定偏移保存完整输出，不再使用 MPI。`--no-dump` 仅允许 RAM 模式，
仍保留完整输出和 header/EOF，只跳过最终物理保存。普通四步串联的前三步不能使用它。

公共库新增非拥有式 `MpiBamInput::Attach(backend, memory_resident)`：后端需已打开，
寿命覆盖所有调用；`Close()` 只解除借用，不关闭或删除它，不能 Attach 自己拥有的后端。
`ResidentRange(offset, size)` 可选地返回文件全局偏移对应的连续驻留地址；不存在时
返回空指针，使用端不得修改输入，地址在后端关闭/存储释放前有效。
四条流水线和 CPE kernel 不依赖这套私有 RAM 存储。

可用 `-DSWBAM_RAM_IO_AUDIT=ON` 单独构建诊断版。它在 MPI 区间包装链接进程序的
POSIX 文件调用，检查输入/输出及内部临时数据路径；`physical_data_io_calls=0`
与实际 RAM 临时读写计数一起检查。该审计不是全系统 syscall trace，正式性能构建关闭它。

## 扩展入口

| 流水线 | 接口 | 定制位置 |
| --- | --- | --- |
| BAM read | `RunCpeReadPipeline` | `CpeBatchOperator` / `BamReadBatchOperator`，结果 slice 与批后处理 |
| BAM write | `CpeWritePipelineSession` / `RunCpeWritePipeline` | `CpeWriteBatchOperator`，块压缩或记录序列化+压缩 |
| SAM read | `RunCpeSamReadPipeline` | `SamReadKernelSpec` + `SamParsedBatchPostProcessor` |
| SAM write | `CpeSamWriteSession` | `SamWriteKernelSpec` + record batch + `RankBodySink` |

常用封装为 `RunFlagstatPipeline`、`RunStatsBasicPipeline`、`RunRecordCountPipeline`、
`RunBamTransformPipeline`、`RunBamToSamPipeline`、`RunSamToBamPipeline`；准确签名见
`swbam-sdk/include/swbam/operators/`。`RunBgzfCompressPipeline` 将已打包块交给公共 BAM writer。

record 层动作通过 `cpe_bam_read_steps.h` / `cpe_bam_write_steps.h` 在 CPE kernel 中
静态组合；不是把每条 record 送回 MPE 再逐条虚调用。用户可在 MPE 批后处理实现
自己的算法，跨批保留记录时必须自己持有内存。

## 验证与限制

`RabbitBAM-MPI io-check -i small.bam` 使用当前 count、transform 和 BGZF 压缩算子，
检查 memory/POSIX/MPI-IO 输入、CRC/ISIZE、完整 BAM 回环和 MAPQ>=30 过滤。
RAM 模式仅检查已预加载的 RAM 后端，避免诊断时重新读物理输入。
诊断包含独立 MPE 校验，耗时不能代表库核心性能。

输入仍要求 record 不跨 BGZF block。内存容量、batch 指针生命周期及 CPE launch
互斥约束见[流水线契约](swbam流水线数据契约与整体架构.md)。构建与 SWLS 提交示例
统一放在[根 README](../README.md)。

本轮构建、小数据全字段回归、RAM 外排审计和 D1 六核组探测见
[RAM 存储与端到端重构验证](results/RAM_存储与端到端重构验证.md)。复测脚本为
`tools/ram_storage_regression.sh`，D1 单次探测脚本为 `tools/ram_d1_probe.sh`；
正式矩阵需另行固定版本、容量和轮次，不直接把单次探测写成中位数。
