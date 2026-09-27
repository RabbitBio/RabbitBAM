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
MpiBamInput(memory) -> PrepareMpiBamInputPlan
    -> RunRecordCountPipeline -> RunCpeReadPipeline
    -> CPE decode + parse bam1_t + count
    -> MPI sum(records), max(core time)
```

命令体：`./swbam-sdk-optimized-count input.bam`。只计数，不运行完整 flagstat。
输出 `sdk_count ranks=... records=... blocks=... input_bytes=... core=...`。

### BAM 读写回环

`swbam-sdk/examples/sdk_optimized_bam2bam.cpp`：

```text
MpiBamInput(memory) -> PrepareMpiBamInputPlan
    -> RunBamTransformPipeline (no-op filter, level 1)
    -> BAM read -> pack plan -> BAM write -> AdaptiveRankBodySink(memory)
    -> MPI sum(records/body bytes), max(core time)
```

命令体：`./swbam-sdk-optimized-bam2bam input.bam`。它是库核心微基准，输出 rank body
到内存，不组装最终 header/EOF、不 dump 文件。需要可用的磁盘 BAM 时用 `run_all`。
两示例的 `core` 不含输入 open/load、全局 scan/plan、最终归约和磁盘 dump；不要直接
等同于应用的 `4.1~4.6` 或端到端 wall time。

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
诊断包含独立 MPE 校验，耗时不能代表库核心性能。

输入仍要求 record 不跨 BGZF block。内存容量、batch 指针生命周期及 CPE launch
互斥约束见[流水线契约](swbam流水线数据契约与整体架构.md)。构建与 SWLS 提交示例
统一放在[根 README](../README.md)。
