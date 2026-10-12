# RAM 存储与端到端重构验证

日期：2026-10-11。本记录是单次功能回归和 D1 探测，不是三次正式性能矩阵。
历史 D1/D2/D3 核心成绩保持原口径，不替换或改名。

## 实现与计时边界

- 默认 `--storage posix`，磁盘输入、真实临时文件和最终磁盘输出均在 MPI 区间内。
- `--storage ram` 使用私有 `benchmark-support`，固定 4 MiB 分页池实际保存完整数据。
- MPI 前：普通 `malloc` 配合 `-xmalloc`，预触碰、按 BGZF 块/完整 SAM 行装入本 rank
  分片、原始头和扫描信息，并通过唯一会话目录同步。MPI 内重新解析、扫描、分区。
- 主时间是每 rank 的 `MPI_Init` 前到 `MPI_Finalize` 后实测区间，取最大值。
  包含 athread、CLI、工作区清理、通信、流水线和实际 RAM 临时/输出 I/O，不扣时。
- 输出复用 body 页并保留 header/EOF；MPI 结束后各 rank 按全局偏移保存。
  `--no-dump` 仅跳过物理保存，不省略完整输出。普通四步串联前三步必须保存。
- `--ram-capacity` 为每 rank 输入、临时和输出存储预算，向下取整到 4 MiB 页；
  算法 `-m` 工作区另计。容量不足明确失败，绝不自动落盘或循环覆盖数据。
- `dedup-workflow` / `dedup-pipeline` 本轮保持磁盘模式；RAM 参数明确拒绝。

四条公共流水线和从核 kernel 未改。SDK 公共部分仅增加借用式
`MpiBamInput::Attach` 和可选 `ResidentRange`，不安装或依赖私有 RAM 模块。

## 环境与构建

Sunway 本地交叉编译器 `/usr/sw/bin/swg++`，C++11，现有 SDK、MPI、HTSlib 1.20
和 libdeflate 构建配置。SWLS `q_share`；每节点 1 rank、64 个实际计算 CPE，
另外 5 个核组只提供跨核组内存。

```bash
cmake -S . -B build_ram_e2e -DBUILD_SWBAM_TOOLS=OFF -DBUILD_SWBAM_EXAMPLES=ON
cmake --build build_ram_e2e -j8

# 独立构建：先安装 SDK，再构建使用端。
cmake -S swbam-sdk -B build_ram_e2e/sdk_standalone \
  -DCMAKE_TOOLCHAIN_FILE="$PWD/swbam-sdk/cmake/SunwayToolchain.cmake" \
  -DCMAKE_INSTALL_PREFIX="$PWD/build_ram_e2e/sdk-install"
cmake --build build_ram_e2e/sdk_standalone -j8
cmake --install build_ram_e2e/sdk_standalone

PREFIX="$PWD/build_ram_e2e/sdk-install"
cmake -S swbam-app -B build_ram_e2e/app_standalone \
  -DCMAKE_TOOLCHAIN_FILE="$PREFIX/lib/cmake/SWBAM/SunwayToolchain.cmake" \
  -DSWBAM_DIR="$PREFIX/lib/cmake/SWBAM" \
  -DSWBAM_BENCHMARK_SUPPORT_DIR="$PWD/benchmark-support"
cmake --build build_ram_e2e/app_standalone -j8
cmake -S swbam-sdk/examples -B build_ram_e2e/examples_standalone \
  -DCMAKE_TOOLCHAIN_FILE="$PREFIX/lib/cmake/SWBAM/SunwayToolchain.cmake" \
  -DSWBAM_DIR="$PREFIX/lib/cmake/SWBAM" \
  -DSWBAM_BENCHMARK_SUPPORT_DIR="$PWD/benchmark-support"
cmake --build build_ram_e2e/examples_standalone -j8
```

工作区、独立 SDK/app/示例、没有私有模块的磁盘 app/示例均已完成构建。
磁盘独立构建设置 `SWBAM_BENCHMARK_SUPPORT_DIR=/nonexistent`。
SDK 安装头和 CMake 导出中没有私有模块引用。旧二进制保存在
`build_sunway_sduhpc/{RabbitBAM-MPI,sdk-count,sdk-bam2bam}-before-ram-e2e`。

## 小数据回归

输入 `performance_test/WES_0.25G.bam`：256,845,825 字节、12,100 BGZF 数据块、
2,668,351 条记录。两节点、每节点 1 rank。基本命令为：

```bash
bsub -N 2 -np 1 -health 1 -mpecg 6 -cgsp 64 \
  -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc \
  -a mpi -b -K -q q_share -J ram_test -o "$UNIQUE_LOG" \
  "$APP" collate -i "$INPUT" -o "$OUTPUT" -n 66 -m 128M \
  --storage ram --ram-capacity 8G --ram-ranks 2 \
  --ram-session "$UNIQUE_ABSOLUTE_SESSION" --ram-ready-timeout 300
```

完整命令及串行提交实现见 `tools/ram_storage_regression.sh`；它依次运行两个 SDK、
统计/过滤、BAM→SAM→BAM、四步串联和两个磁盘 dedup 命令，分别保留日志。
重新执行须使用新根目录或独立 `SESSION_TAG`，不复用旧日志/会话目录。

实际小数据输入绝对路径为
`/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/code/RabbitBAM/performance_test/WES_0.25G.bam`。
小数据根目录为下述 `ram_e2e_20261011_small2`，最后串联复测日志带 `_v8` 后缀。

| 检查 | 结果 |
|---|---|
| SDK count | 2,668,351；仍由 CPE 构建 `bam1_t` 后计数 |
| SDK 完整 BAM | 275,215,131 字节；记录全集与输入一致 |
| `--min-mapq 30` | 2,481,669，与独立 Pysam 筛选一致 |
| BAM→SAM→BAM | 2,668,351，记录全集字段一致 |
| collate，强制外排 `-m 128M -n 66` | 记录全集保持，1,328,196 个 QNAME 组连续 |
| fixmate `-m` | mate 坐标、MC 和 ms 独立逐对检查 |
| sort，强制外排 `-m 128M` | 坐标排序、记录全集及标签保持检查 |
| markdup `--stream` | 2,668,351，15,695 个重复标记 |
| 磁盘 dedup-workflow / dedup-pipeline | 完整可解析；记录数、排序和重复标记检查 |

独立验证使用 Pysam 0.24.1，见 `tools/verify_ram_outputs.py`。记录文本的 128-bit
多集合指纹同时覆盖顺序无关的全集字段，不只核对计数；另检查排序、QNAME 连续性、
mate 坐标和标签。SHA/多集合指纹不能代替所有算法语义测试，本次是输入上的回归。

独立全记录复读已通过：fixmate 和 sort 的指纹均为
`40286709157607695747077692683583223282`；最终 RAM markdup、磁盘 workflow
和融合 pipeline 均为 `295976077210537395667207712327857750245`。
逐对检查了 1,323,599 组已映射 primary pair，`ms` 与 mate 的 Q>=15 质量和一致。

最终版普通串联单轮结果（每次 MPI 区间独立取 rank 最大值）：

| 阶段/作业 | MPI 区间 max | 4.3 诊断 |
|---|---:|---:|
| collate 8499578，`-m 128M -n 66` | 9.801075 s | 2.643227 s |
| fixmate 8499582，`-m` | 8.701562 s | 1.553360 s |
| sort 8499590，`-m 128M` | 9.691323 s | 2.583234 s |
| markdup 8499598，`--stream` | 9.076337 s | 1.963667 s |
| 逐阶段 MPI 区间合计 | **37.270297 s** | 不作为正式中位数 |

输出在前三条命令的 Finalize 后实际保存，下一条重新预加载；这些准备、保存、
排队时间不计入上述区间合计。最终验证作业 8499606 / 8499609 已成功。
外排的新 `4.3` 包含真实 RAM 临时 I/O，不能与历史扣时后的模型 `4.3` 直接比较。

最后两个磁盘模式回归也已结束：`dedup-workflow` 8499620 的 MPI 区间 max 为
43.133725 s，`dedup-pipeline` 8499623 为 20.412566 s；随后复读作业
8499625 / 8499627 均得到 2,668,351 条、`is sorted: 1`。最后输出的 SHA-256
与已完成独立全字段验证的本地副本一致：

- workflow：`e0ced2bd1fa060eadf480c76028994ed7d5b96d05f6396e752c9bcc922c2ecf9`。
- pipeline：`0d1b1c9ab1e860525b95c5833892e1c521254943f9b92fa826a8b76df1bdb587`。

最终小数据回归驱动于 2026-10-11 11:27:11 成功结束；这些单次回归耗时不是正式成绩。

回归额外修复了两个潜在问题，算法规则和 kernel 保持不变：

1. 磁盘融合 pipeline 中 collate 的 HTSlib 内存文件关闭会释放借用输入，造成后续
   stage 悬空；改用已有非拥有式头解析 helper。
2. 外排 sort 的 CPE 交换打包在非对齐边界偶发损坏记录尾部。按照内排路径的做法，
   join 后由 MPE 补写最多每 CPE 两条边界记录，不增加全批复制。

## 外排与边界检查

RAM 文件单元测试通过 ASan/UBSan：分页读写、覆盖、空洞清零、EOF、缩短/增长、
关闭回收、碎片化页、增长不移动已有地址及容量失败回滚。
SWLS 额外检查空 BAM、单条 BAM/SAM 的空 rank 分片、默认保存、Finalize 后数据存活、
`--no-dump`、预加载/输出/临时容量不足、准备超时及非法 CLI 的退出。

可选 `-DSWBAM_RAM_IO_AUDIT=ON` 诊断构建包装链接进程序的 POSIX 文件调用。
强制外排实测如下，所有 rank 的 `physical_data_io_calls=0`：

| 作业 | RAM 临时读字节合计 | RAM 临时写字节合计 | 存储池峰值/rank | MPI 区间 max |
|---|---:|---:|---:|---:|
| collate 8499210，`-m 128M` | 1,645,372,344 | 1,644,929,604 | 1,107,296,256 | 9.715908 s |
| sort 8499211，`-m 128M` | 1,777,419,983 | 1,776,932,582 | 1,170,210,816 | 9.675086 s |

审计不是全系统 syscall trace；结论针对被包装调用，并结合输入/输出/临时适配器
的源代码检查。控制目录的就绪标记不是用户数据 I/O；正式性能构建关闭审计。
512 MiB 强制外排容量失败测试明确触发 `ENOSPC` 和 `MPI_Abort(28)`；SWLS
调度器未及时回收此失败作业，已人工取消该自建测试，不作为成功性能结果。

## D1 六核组单次探测

输入 9,891,133,304 字节、814,215 数据块、177,798,571 条记录。
命令完整保存在 `tools/ram_d1_probe.sh`，所有提交资源为：

```bash
-N 6 -np 1 -health 1 -mpecg 6 -cgsp 64 \
-share_size 100 -cross_size 88000 -cache_size 32 -xmalloc \
-a mpi -b -K -q q_share
```

SDK count：`--storage ram --ram-capacity 2G --ram-ranks 6`。
SDK 回环：`INPUT -o OUTPUT --storage ram --ram-capacity 8G --ram-ranks 6`，
level 1，默认保存完整输出并通过 `stats --basic` 复读。
每条命令还有唯一 `--ram-session` 和 `--ram-ready-timeout 900`。

| 命令/作业 | 新核心诊断 | 旧二进制同日核心 | 回退 | 新 MPI 区间 max |
|---|---:|---:|---:|---:|
| SDK count 8499090 / 8499168 | 2.470854 s | 2.470581 s | +0.011% | 9.714538 s |
| SDK BAM 回环 8499229 / 8499316 | 8.237332 s | 8.577784 s | -3.967% | 15.508824 s |

count 存储池峰值/rank 为 1,694,498,816 字节；BAM 回环最大峰值/rank 为
3,619,684,352 字节，最终完整 BAM 为 11,467,690,272 字节。
MPI 区间新增约 7 s 的初始化、退出和正常控制成本，不能据此判定流水线回退。
两项单次核心比较均在本轮约 5% 回退目标内，不代表正式三次中位数。

| app/作业 | 4.3 诊断 | MPI 区间 max | 完整输出字节 |
|---|---:|---:|---:|
| flagstat 8499528 | 3.380914 s | 10.678473 s | 无输出文件 |
| MAPQ>=30 过滤 8499547 | 8.448405 s | 15.721746 s | 10,899,180,596 |
| collate 8499566，memory，`-n 66 -m 24G` | 57.248594 s | 64.526982 s | 18,755,398,492 |
| sort 8499579，memory，`-m 24G` | 47.282966 s | 54.562130 s | 11,467,679,314 |

过滤保留 170,359,591 条；collate 的 records/received 均为 177,798,571，
groups=88,498,907；sort 的 local/received 均为 177,798,571。
这三项写 app 使用 `--no-dump`，完整输出驻留；只另保存 SDK D1 输出用于全量复读。
SDK 复读 8499381 得到 177,798,571 条、`is sorted: 1`，MPI 区间 max 为 10.512741 s。
D1 没做全量逐记录指纹比较，不能把上述计数检查称为 D1 全字段一致性验证。

相同 D1 原始输入的历史 `4.3` 中位数：flagstat 3.373142 s、过滤 8.793768 s、
collate 58.188765 s。此次分别 +0.230%、-3.927%、-1.616%；这是历史参考，
不替代同日 SDK 旧二进制配对。sort 本次输入是原始 coordinate BAM，历史流程
使用 fixmate 输出，不能拿两者的时间直接证明回退或提升。

本轮没有启动完整 D1/D2 矩阵，也没有取消或覆盖旧实验任务。

本次 9 个串行 D1 作业（含两项旧版对照和复读）于 10:54:58 至 11:22:49 完成。
RAM 输入准备实测每条约 106～148 s，不计成绩但影响等待；SDK 保存也在区间外。
后续仅六核组两项 SDK 各三次，按当前 I/O 约需 15～25 分钟；若按 6 个资源点、
两项、各三次共 36 个命令顺序跑 D1，先估计 1.5～2.5 小时，排队和额外验证另计。
D2 应先容量/准备耗时探测，再给出矩阵估计，不能直接复用 D1 的每 rank 容量。
扩展性矩阵也须随 rank 数调整容量；上述 count 2G / 回环 8G 只针对六 rank，
不能直接用于单 rank。D1 单 rank 的完整读写存储预算可先按约 24G 探测。

## 日志与版本

SWLS 独立目录均位于
`/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data/swbam_bench/`：

- `ram_e2e_20261011_small2`：小数据各版本回归、完整输出、失败诊断与修复后的重跑。
- `ram_e2e_20261011_final`：冻结版本的 D1 单次探测、旧 SDK 二进制对照及容量边界测试。
- `ram_e2e_20261011_audit`：强制 RAM 外排的 I/O 审计构建和日志。

D1 冻结版与最后小数据版的区别是清理边界、磁盘融合输入所有权和外排边界保护，
SDK 的 CPE 计算路径没有改变。应随正式矩阵重新固定最终二进制 SHA-256。

工作区基于 `9d6f72251bcbe7182e1e1636af80d84474cbf5bd`，本轮修改未自动 commit。

| 二进制 | SHA-256 |
|---|---|
| D1 冻结 app | `7fb10077a794d8481b63d9c385dadef54a855548b12c29df5ca5001c20cde47c` |
| D1 冻结 count | `4ceb087dc9fd8c875214537088df4b68e0695d637935564b77effe4b0c139563` |
| D1 冻结回环 | `aaa67aeeac7f6282d8bc35166fdb7191de6af9cd33d24dd4fc55b0bbd62ddd0e` |
| 最终小数据 app | `5fcc7b5e6758b8c503d75019c2a14e03db49db96bdf8c031687200f1bd9a8399` |
| 最终 count | `925af1c5e20d79cb57bddf56ab0a5c570056b27feaab467291d54420d1b80b0a` |
| 最终回环 | `8c6ee35ca1cb7b8ce34881384de5bfbb24fcf89303b730bd1759b0968096bd8f` |
| 旧 count 对照 | `32e16d1f6d39e19d921b0cff45776ad899ee9b7f579a3a76153a61e7f5352a2c` |
| 旧回环对照 | `426660818471b928c992f1b840a2384a75ff09513a73449306c67404f1910f36` |

本地原始日志归档为 `build_ram_e2e/benchmarks/{small,D1,audit}/`，完整输出留在
SWLS 各独立目录的 `out/`。原有论文文件及历史结果的用户修改未回退。
