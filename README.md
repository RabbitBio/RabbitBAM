# SWBAM

**面向神威异构众核平台的分布式 BAM/SAM 处理库与应用。**

SWBAM 在 MPI rank 间划分任务，在每个核组内由 MPE 调度、CPE 执行编解码与记录
处理。应用通过四条批处理流水线和可组合 operator 构建，复用输入后端、双缓冲、
BGZF 编解码、解析、批次生命周期管理与分布式输出。

当前维护 Sunway 版本，命令行程序名保留为 **`RabbitBAM-MPI`**。
旧单机/x86、旧 CGS 应用和旧 SDK 示例已退役；历史实现可从 Git 历史查阅。

## 功能

| 命令 | 功能 |
| --- | --- |
| `run_all` | SAM→BAM、BAM→SAM、BAM→BAM；BAM→BAM 可附加 MAPQ、FLAG、参考序列和长度过滤 |
| `flagstat` | 比对标志分类统计 |
| `stats --basic` | SN 摘要统计及排序状态检查，不等同于 SAMtools 完整 stats |
| `sort` | 坐标排序，支持内排与外排 |
| `collate` | 按 QNAME 聚合，支持内排与外排，不保证与其他实现具有相同组间顺序 |
| `fixmate` | mate 信息和标签修复；`-m` 生成后续 markdup 使用的 mate score |
| `markdup --stream` | 坐标有序输入上的流式重复标记；保留已有两遍回退与删除模式 |
| `dedup-workflow` | 使用阶段文件串联 collate、fixmate、sort、markdup |
| `dedup-pipeline` | 使用内存 BAM 中间态连接四阶段；不等同于全程有界内存处理 |
| `io-check` | 新 BAM read/write 流水线、count/filter、BGZF 校验及输入后端一致性检查 |

具体选项以各子命令 `--help` 为准。SWBAM 并不实现 SAMtools 的全部命令和选项。

## 架构与目录

```text
input backend → MPI input plan → BAM/SAM read pipeline
             → CPE operator → MPE batch post-processing / app algorithm
             → BAM/SAM write pipeline → rank body sink → distributed output

include/swbam/     库接口、四条流水线及 operator 契约
lib/io/           memory/POSIX 输入输出、rank body memory/spool
lib/mpi/          MPI-IO、输入分区与分布式输出
lib/cpe/          MPE 侧四条流水线和写出适配器
lib/operators/    批次参数准备、状态校验、结果归并及读写连接
slave/core/       CPE 共享 BGZF codec、BAM parser 与辅助函数
slave/operators/  CPE 计数、过滤、解析、格式化、压缩 kernel
slave/algorithms/ sort/collate/fixmate/markdup 的算法专用 kernel
apps/mpi/         CLI、命令、复杂算法和去重编排
examples/         两个与应用复用同一流水线的 SDK 示例
tools/            正确性验证、内存与 CRC 探测工具
docs/             架构契约、实验方案、结果与研究记录
ext/              HTSlib、libdeflate 等依赖及神威适配
overleaf/         论文材料
```

**batch 层组合，record 层静态融合。** 例如解压、解析、count/filter 可以在一次 CPE
kernel 中完成；MPE 回调只发生在批次边界。双缓冲重叠输入预取与计算、压缩与上一批
输出，不表示同一核组的读写 CPE kernel 同时运行。

`include/BamTools.h` 保存仍在使用的 MPE/CPE 共享参数布局；`apps/mpi/CmdInfo.h`
和 `apps/mpi/swbam_mpi.h` 是应用内部接口，不属于 SDK。

## 构建

需要 Sunway `swgcc`/`swg++`、athread、MPI、CMake ≥ 3.13，以及本仓库神威适配的
HTSlib 1.20 / libdeflate 1.20。不是直接在 x86 上运行的程序。

在仓库根目录执行：

```bash
# 依赖未构建时执行；需要本地可用的神威工具链。
bash build_ext_sw.sh

cmake -S . -B build_sunway_sduhpc \
  -DBUILD_SWBAM_EXAMPLES=ON \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build_sunway_sduhpc -j 8
```

MPI 默认路径为 `/usr/sw/mpi/mpi_20220608_SEA`，可通过 `-DSUNWAY_MPI_ROOT=...`
指定。保留 `-DPLATFORM=sunway` 旧脚本兼容；不再提供 x86 平台分支。

默认产物：

- `RabbitBAM-MPI`：应用入口。
- `swbam-sdk-optimized-count`、`swbam-sdk-optimized-bam2bam`：两个 SDK 示例。
- `libswbam_io.a`、`libswbam_mpi_runtime.a`、`libswbam_cpe_runtime.a`：三个主核静态库。
- CPE kernel 以 object target 参与神威 hybrid 链接，不作为普通主核静态库打包。

验证/探测工具通过 `-DBUILD_SWBAM_TOOLS=ON` 开启。对外统一 CMake 链接目标为
`swbam::swbam`；当前采用源码树内集成，尚未提供 install/find_package 分发包。

## 运行示例

以下命令体需要通过集群的 MPI/作业启动器执行，不应在登录节点直接运行大数据任务。
`run_all` 根据输入格式和输出后缀选择转换方向。

```bash
./RabbitBAM-MPI flagstat -i input.bam --io-backend memory
./RabbitBAM-MPI run_all -i input.bam -o output.bam --compress-level 1
./RabbitBAM-MPI run_all -i input.bam -o filtered.bam --min-mapq 30 --compress-level 1
./RabbitBAM-MPI run_all -i input.bam -o output.sam
./RabbitBAM-MPI run_all -i input.sam -o output.bam --compress-level 1
./RabbitBAM-MPI io-check -i small.block-aligned.bam
```

SWLS 示例（在构建目录，路径按实际部署位置调整）：

```bash
bsub.py swls -q q_share -I -b -J swbam_check \
  -N 1 -np 6 -cgsp 64 -share_size 14000 -cache_size 32 \
  -o io_check.log \
  -x LD_LIBRARY_PATH=../ext/htslib-1.20:../ext/libdeflate-1.20/build \
  ./RabbitBAM-MPI@online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc \
  io-check -i ../performance_test/WES_0.25G.bam
```

普通四步去重的命令体：

```bash
./RabbitBAM-MPI collate -i input.bam -o name.bam -m 8G --compress-level 1
./RabbitBAM-MPI fixmate -i name.bam -o fixmate.bam -m --compress-level 1
./RabbitBAM-MPI sort -i fixmate.bam -o sorted.bam -m 8G --compress-level 1
./RabbitBAM-MPI markdup -i sorted.bam -o markdup.bam --stream -l 300 -m 8G --compress-level 1
```

这里 `fixmate -m` 表示生成 mate score，其他三个命令的 `-m 8G` 为算法内存配置；
`markdup -l` 必须按输入的实际读长/比对特点设置，不应对任意数据照搬示例。

## SDK 示例

两个示例都默认使用内存输入，由 MPI 分区后调用与应用相同的库入口：

| 示例 | 库入口与行为 |
| --- | --- |
| [sdk_optimized_count.cpp](examples/sdk_optimized_count.cpp) | `RunRecordCountPipeline`：CPE 解压、解析为可复用 `bam1_t` 并计数，不执行完整 flagstat |
| [sdk_optimized_bam2bam.cpp](examples/sdk_optimized_bam2bam.cpp) | `RunBamTransformPipeline`：BAM read → 打包 → BAM write，压缩 body 写入 memory sink |

命令体分别是 `./swbam-sdk-optimized-count input.bam` 和
`./swbam-sdk-optimized-bam2bam input.bam`。后者是库读写微基准，**不生成带 header/EOF
的最终磁盘 BAM**；需要文件输出使用应用 `run_all`。

## I/O、计时与限制

- 默认 memory 模式先把完整输入加载到各 rank 内存，以内存读写模拟高速 I/O。
  核心时间不包含预加载和最终文件 dump，不能描述为真实磁盘端到端时间。
- Standalone 命令提供 `--io-backend memory|posix|mpiio|auto`，输出可选
  `--io-output-backend memory|mpiio`；具体支持范围以命令帮助为准。
- `--rank-body-backend memory|spool|auto` 控制 rank 本地输出缓存；其临时目录
  由 `--rank-body-temp-dir` 指定，与 sort/collate 的 `-T` 外排临时文件独立。
- **要求单条 BAM record 不跨 BGZF block**。一般来源的 BAM 不保证满足此条件；
  `io-check` 可检查小输入是否满足当前解析契约，不负责自动重分块。
- BAM record 长度、批容量和 CPE scratch 有边界限制；不能将流式输入等同于支持
  任意大小单条记录。全局 BGZF 索引、内存输入/输出及算法状态仍需计入容量规划。
- 默认 BAM 压缩 level 1。正式性能结果应记录数据、资源、压缩配置、版本和计时范围。
- `io-check` 是诊断，不是性能基准：它重复读取并用 MPE 校验 CRC 和结果，建议用小数据。

## 文档

- [流水线数据契约与整体架构](docs/swbam流水线数据契约与整体架构.md)
- [SDK 使用说明](docs/swbam_sdk使用说明.md)
- [实验方案](docs/swbam实验方案.md)与[结果目录](docs/results/)
- [大共享内存及工具说明](docs/swbam工具使用说明.md)

实验结果与论文草稿包含历史版本信息，应以对应测试记录中的版本、命令和计时口径为准。
