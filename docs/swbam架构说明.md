# SWBAM 架构说明

当前主线是 **四条类型化流水线 + 可组合 operator + 应用算法**。详细数据契约与迁移
记录统一见[流水线数据契约与整体架构](swbam流水线数据契约与整体架构.md)，不再维护
另一套 Raw/MPE bam1_t adapter 架构。

```text
swbam-sdk/lib/io + lib/mpi
    backend、block plan、rank body sink、distributed output
                       |
swbam-sdk/lib/cpe       四条 MPE 调度流水线
    BAM read / BAM write / SAM read / SAM write
                       |
swbam-sdk/lib/operators batch 参数、校验、归并和流水线连接
                       |
swbam-sdk/slave         CPE 编解码、解析与编译期业务动作组合
                       |
swbam-app              CLI、MPI 协同、复杂算法状态和流程编排
swbam-app/slave         应用专用 CPE kernel
```

`slave/` 专门存放实际从核源码；`lib/cpe/` 保留原名，存放在 MPE 上执行的
流水线调度。`swbam::cpe` 等接口名仍表示从核编译环境，不是源码目录名。

## 职责边界

- 库管理输入输出与批缓冲，不控制 sort/collate 的全局重排或 markdup 判重策略。
- app 的 sort、collate、fixmate、单遍流式 markdup 在读写两端复用库，中间状态留在
  app；仍在使用的两遍 markdup 回退保留，不能当成死代码删除。
- BAM read 的解压、解析和 count/filter 等动作在一次 CPE launch 内静态融合。
- BAM write 统一调度已打包块压缩和记录序列化+压缩；record adapter 只负责连接。
- 双缓冲覆盖 MPE 预取/flush 与 CPE 工作，不允许同核组未 join 就启动另一 CPE kernel。
- `dedup-workflow` 与 `dedup-pipeline` 保留；前者阶段文件串联，后者全内存中间态，
  两者不是相同的内存规模和计时模型。

## 当前产物

主程序 `RabbitBAM-MPI`；SDK 示例只有 `sdk_optimized_count.cpp` 和
`sdk_optimized_bam2bam.cpp`。三个主核静态库与 CPE object 通过 `swbam::swbam`
统一链接。SDK 和 app 是两个独立 CMake 项目；app 可通过安装包
`find_package(SWBAM CONFIG REQUIRED)` 构建。见[SDK 使用说明](swbam_sdk使用说明.md)。

SDK 公共头文件仅位于 `swbam-sdk/include/swbam/`。原 `BamTools.h` 的库类型迁入
`bam_types.h`；sort/collate/fixmate/markdup 的参数原样迁入
`swbam-app/include/algorithm_types.h`，由 app 主从核共同包含。
`CmdInfo.h`、`swbam_mpi.h` 和 CLI11 为应用私有依赖，不随 SDK 安装。

安装包保留 CPE object 的直接 hybrid 链接，未增加运行时分发。构建选项通过
`swbam::host`、`swbam::cpe` 传递；源码树不再使用全局 include/link 搜索目录。
应用自定义 CPE kernel 用 `swbam::cpe` 编译，仍可使用库的解压/解析/压缩步骤。

## I/O 与实验边界

memory backend 预加载整个文件后模拟内存读写；POSIX/MPI-IO 是可选真实输入后端，
rank body 可使用 memory/spool/auto。算法外排临时文件与输出 body spool 独立。
核心模拟时间排除真实 load/dump，不能标成真实磁盘端到端耗时。

单条记录不跨 BGZF block 是当前显式约束。全局 block index、某些算法状态与内存
输出仍可能随数据规模增长，不承诺任意大输入均有界。默认测试为 level 1；大规模、
其他参数及外排需要各自的功能和资源验证。
