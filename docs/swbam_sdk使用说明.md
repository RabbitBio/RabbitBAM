# SWBAM SDK 使用说明

当前 SDK 以四条流水线和可组合 operator 为入口，源码树内链接统一目标
`swbam::swbam`。旧 Raw view/MPE bam1_t adapter 与旧五个示例已退役，不再维护两套
示例架构；CPE BAM parser 对 `bam1_t` 的解析能力保留。

## 构建产物

| 目标 | 用途 |
| --- | --- |
| `swbam_io` / `libswbam_io.a` | 输入/输出后端、BGZF batch、rank body storage |
| `swbam_mpi_runtime` / `libswbam_mpi_runtime.a` | MPI input plan、MPI-IO、分布式输出 |
| `swbam_cpe_runtime` / `libswbam_cpe_runtime.a` | MPE 侧四条流水线与算子适配 |
| `swbam_cpe_kernels` | 共享 CPE object 的接口目标，不是普通 `.a` |
| `swbam_algorithm_kernels` | 应用算法 CPE object，SDK 示例不需要链接 |
| `swbam::swbam` | SDK 总链接目标，传递头文件、库、MPI 依赖和 hybrid 链接选项 |

尚未提供 install/find_package 包。新示例在本仓库 CMake 中按已有示例添加：

```cmake
add_executable(my_tool my_tool.cpp)
target_compile_options(my_tool PRIVATE ${HOST_FLAGS})
target_link_libraries(my_tool PRIVATE swbam::swbam)
```

程序负责 MPI 初始化和 athread 生命周期；不能将神威从核对象当成普通主核 archive。

## 两个示例

### 只读计数

`examples/sdk_optimized_count.cpp`：

```text
MpiBamInput(memory) -> PrepareMpiBamInputPlan
    -> RunRecordCountPipeline -> RunCpeReadPipeline
    -> CPE decode + parse bam1_t + count
    -> MPI sum(records), max(core time)
```

命令体：`./swbam-sdk-optimized-count input.bam`。只计数，不运行完整 flagstat。
输出 `sdk_count ranks=... records=... blocks=... input_bytes=... core=...`。

### BAM 读写回环

`examples/sdk_optimized_bam2bam.cpp`：

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
`include/swbam/operators/`。`RunBgzfCompressPipeline` 将已打包块交给公共 BAM writer。

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
