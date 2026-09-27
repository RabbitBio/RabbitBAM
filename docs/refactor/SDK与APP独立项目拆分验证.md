# SDK 与 app 独立项目拆分验证

日期：2026-09-27。拆分前版本：`064b166`（Streamline the project）。

## 拆分内容

- `swbam-sdk/`：公共头文件、三个主核静态库、四条流水线、算子、通用 CPE 实现、
  CMake 安装包和两个 SDK 示例；没有 app 头文件或应用算法依赖。
- `swbam-app/`：CLI、命令、复杂算法、去重编排和四个应用专用 CPE kernel。
  独立构建时使用 `find_package(SWBAM 0.1 CONFIG REQUIRED)`。
- `BamTools.h` 的库内容迁入 `swbam/bam_types.h`，应用参数迁入
  `swbam-app/include/algorithm_types.h`；字段定义和辅助函数逐字核对一致。
- `Globals.h` 迁入 `swbam/platform.h`；CLI11 归入 app；CRC benchmark 头文件和
  从核代码归入 tools。算法源文件只移动或替换 include，没有修改算法体。
- 根目录继续支持一键构建，产物仍位于 `build_sunway_sduhpc/`；两个去重命令保留。

## 安装包与链接

SDK 版本暂定 `0.1.0`。导出统一入口 `swbam::swbam`，另有主核/CPE 编译接口
`swbam::host`、`swbam::cpe`。库只传递 target 级 include、编译和链接选项。

安装包携带通用 CPE object，继续直接进入神威 hybrid 链接；主核 HTSlib 和
libdeflate 采用静态依赖。安装的 CMake 文件不含原仓库路径、源码路径或构建路径。
MPI 与神威系统运行时使用 `SUNWAY_MPI_ROOT` 和系统工具链。

## 构建验证

1. 根目录联合构建：app、两个 SDK 示例、验证/内存/CRC 工具全部通过。
2. `cmake -S swbam-sdk ... -DBUILD_SWBAM_EXAMPLES=OFF` 独立构建 SDK，通过。
3. 安装 SDK 后，将 `sdk-install` 移为 `sdk-relocated`，验证前缀可迁移。
4. 将 app 和 examples 分别复制为独立源码目录，仅通过新前缀与安装后的
   Sunway toolchain 配置；app 与两个示例均构建、链接通过。
5. 检查独立构建的 compile/link 输入，不引用原 SDK 的源码/构建目录或仓库 ext。
6. `git diff --check HEAD` 通过。

验证产物位于 `build_sunway_sduhpc/standalone/`：`sdk-build`、`sdk-relocated`、
`app-source`、`app-build`、`examples-source` 和 `examples-build`。
这些是忽略的本地构建产物，不作为源码版本提交。

独立构建与安装的日常命令见 [SDK 使用说明](../swbam_sdk使用说明.md)。

## SWLS 小数据验证

测试的是 **从安装包构建的 app**，复制为 `RabbitBAM-MPI-sdk-split` 后提交。
单节点、6 ranks、每 rank 64 CPE；输入 `WES_0.25G.bam`，12,100 BGZF blocks，
2,668,351 条记录，解压 payload 788,045,119 bytes。

```bash
bsub.py swls -q q_share -I -b -J sdk_split_io_check \
  -N 1 -np 6 -cgsp 64 -share_size 14000 -cache_size 32 \
  -o sdk_split_io_check.log \
  ./RabbitBAM-MPI-sdk-split@online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc \
  io-check -i ../performance_test/WES_0.25G.bam
```

作业 `8318149`，返回 0。日志：本地构建目录 `sdk_split_io_check.submit.log`，
SWLS 对应构建目录 `sdk_split_io_check.log`。

| 检查 | 结果 |
| --- | --- |
| memory / POSIX / MPI-IO read/count | 均 PASS，records=2,668,351 |
| BGZF 压缩与校验 | PASS，每 rank 128 个合成块 |
| BAM read → BAM write 回环 | PASS，total=kept=2,668,351，body bytes=275,253,532 |
| MAPQ≥30 过滤回环 | PASS，total=2,668,351，kept=2,481,669，body bytes=252,880,296 |

本轮未做完整算法回归或重复性能测量。io-check 包含额外独立校验，不能用于证明
所有命令性能无变化；本次证据覆盖项目边界、可迁移安装包、混合链接及基本读写语义。

## 从核目录命名统一

同日将 SDK 和 app 的从核源码目录统一为 `swbam-sdk/slave/` 与
`swbam-app/slave/`，从核专用头文件目录同步为 `include/swbam/slave/`。
`swbam-sdk/lib/cpe/` 仍是 MPE 侧调度实现，不更名；未修改算法和 kernel 逻辑。

重命名后联合构建通过；重新安装 SDK，并通过该安装包独立构建两个 SDK 示例通过。
从核 object 的安装位置同步为 `lib/swbam/slave/`。验证产物在
`build_sunway_sduhpc/slave-layout/`，`git diff --check HEAD` 通过。
本次命名调整未另行提交 SWLS 作业，上述 io-check 为重命名前的项目拆分验证。
