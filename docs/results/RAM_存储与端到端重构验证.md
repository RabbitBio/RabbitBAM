# RAM 存储与端到端重构验证

## 2026-10-13 完整处理计时补充

为 app 与两个 SDK 增加从 `MPI_Init` 返回后到 `MPI_Finalize` 调用前的
`[swbam-processing rank=...]`。每 rank 使用同一单调时钟取两个时间戳，
Finalize 后打印；不增加 MPI 测量通信。处理时间包含正常初始化、扫描分区、
算法、完整最终 RAM 输出和清理，与 `core` 和 MPI 区间分列。
原算法、四条流水线、从核 kernel、RAM 容量和算法预算没有因计时改变。

工作区 RAM 构建及独立无 RAM app/示例构建通过，30 项脚本测试通过，
覆盖字段缺失/重复、处理时长大于同 rank MPI 区间，以及两个 max 来自不同 rank。
`git diff --check` 通过。

两节点 WES_0.25G 小数据的五条作业完成，均为 2,668,351 records：

| 项目 | 作业 | 处理 max/s | MPI 区间 max/s |
|---|---:|---:|---:|
| RAM SDK count | 8532055 | 0.146074 | 7.276649 |
| RAM SDK BAM 回环 | 8532057 | 0.511552 | 7.621294 |
| RAM flagstat | 8532062 | 0.209146 | 7.284005 |
| 回环输出完整 stats 复读 | 8532071 | 0.209470 | 7.292564 |
| 无私有 RAM 模块的磁盘 SDK count | 8532123 | 2.753961 | 9.791406 |

RAM SDK 回环默认在 Finalize 后保存完整 BAM，随后完整复读。
回环处理 max 来自 rank 0，MPI 区间 max 来自 rank 1；这是独立取时汇总，
并非将两种全局最大值相减。小数据用于验证字段和正确性，不是 D1 正式性能点。
验证作业 DONE 后日志短暂延迟可见，验证脚本等待完整日志并复用同一作业，
没有重新计算或将无日志结果记为通过。

远端验证目录：

```text
/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data/swbam_bench/processing_timer_smoke_20261013_043000
```

完整处理计时版冻结 SHA-256：

| 二进制 | SHA-256 |
|---|---|
| app | `81b68788a9d897b80f494fd1d84945f15b2fe8c7791a28f3559779c77c683cb1` |
| count | `abdca8692b6b3fa6cb888da4da8ade0c7a268403628807e095a5a5c83a67a7f7` |
| bam2bam | `121d764161e8f6737aad9ada8f3ce5c85195f46456fa3daaba8f682b404df9e4` |

04:43 启动对应 D1 正式矩阵和 D2 容量探测，分别见其结果文档。
下文为此前重构版本的验证和成绩记录。

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

## 2026-10-12 六节点 I/O 探测

独立作业 **8525736**，21:59 提交，22:02:10 确认三轮完成。
固定六个健康空闲节点 **43987-43992**，每节点一个 MPI rank；配置与正式矩阵相同：
`-N 6 -np 1 -health 1 -mpecg 6 -cgsp 64 -share_size 100 -cross_size 88000
-cache_size 32 -xmalloc`。探测只在 MPE 搬运数据，CPE 不参与，不是六核组的
总 DRAM 峰值测试。没有修改 SDK/app 或冻结程序，也没有取消 D2 任务。
D2 矩阵同期在其他节点运行，共享存储负载未隔离；下表是当时的快照，不是系统空载上限。

工具为 `tools/mpi_io_bandwidth_probe.cpp`，提交入口为
`tools/run_mpi_io_bandwidth_probe.sh`。每 rank 的 source/destination 各 512 MiB，
另外一个 8 MiB scratch，通过普通 malloc 分配并预触碰；复制单项重复 16 次，
每项每 rank 搬运 8 GiB。内存复制、两次复制的组装和文件 I/O 各测三轮。
每次文件测试为每 rank 512 MiB，六节点共 3 GiB；传输块为 8 MiB。
私有文件写计入 fsync，共享文件写使用 MPI_File_write_at 并计入 MPI_File_sync。
分配、初始化、打开/布局、内容校验和清理在带宽计时外；共享写前的完成检查在计时内。

单位全部为十进制 GB/s；合计取 `6 * bytes_per_rank / max(rank_elapsed)`，
再对三轮取中位数，不将六个独立 rank 的带宽简单相加。
内存项按逻辑输出字节数计，不把源读取、目标写入或两次 memcpy 的流量重复计入分子。

| 项目 | 单 rank 范围，三轮所有 ranks / GB/s | 六节点合计中位数 / GB/s |
|---|---:|---:|
| 64 KiB 分块 memcpy | 2.409-2.417 | **14.461387** |
| 8 MiB 分块 memcpy | 4.794-4.814 | **28.793163** |
| 8 MiB source -> scratch -> destination，两次 memcpy | 2.368-2.376 | **14.215216** |
| 独立文件顺序写 + fsync | 0.0795-0.0904 | **0.480160** |
| 独立文件读，先请求 FADV_DONTNEED | 0.1077-0.1658 | **0.737345** |
| 同文件立即再次读取 | 0.3975-0.5373 | **2.587589** |
| 共享文件分布式位置写 + MPI_File_sync | 0.0575-0.0605 | **0.353735** |

测试目录位于当前共享文件系统，不是本地 SSD。FADV_DONTNEED 只是建议，
无法保证清空客户端和服务端所有缓存；立即复读也不能当作真实冷盘带宽。
O_DIRECT 打开后的读取没有形成有效样本，三轮均跳过；原日志没有记录读取失败
的具体 errno，不据此断言整个文件系统或硬件不支持直读。没有运行 drop_caches。
上述读数是共享文件接口的有效带宽，不是后端磁盘阵列的纯硬件指标。

独立 POSIX 文件和 MPI 共享文件均全量读回，对比每 rank 的随机源数据通过；
全部自建测试数据已清理，`remaining_files.txt` 为空。
本地两 rank 小数据验证、Sunway 交叉构建、三轮六节点验证、bash -n 和
git diff --check 通过。单文件测量规模不能保证大数据持续写的带宽相同。

D1 fixmate 输出约 19.29 GB；若最终组装使用本次两次复制的有效带宽，
纯搬运约 `19.29 / 14.215 = 1.36 s`。这仅为数量级估算，不包含页分配、
生命周期管理、MPI 控制及应用并发干扰，不是新的 fixmate 实测成绩。
相比耗时，更应先解决同时保留旧 body 和新输出造成的峰值容量。
RAM 池是固定分配且预触碰的，预算要算完整池容量而非仅文件页 peak；
输入、临时文件和输出已经包含在池里，不能再次重复加算，算法工作区另计。

远端根目录：

```text
/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data/swbam_bench/io_n6_20261012_215000
```

原始日志、三轮摘要、命令、节点、作业和校验和已归档到本地
`build_ram_e2e/benchmarks/io_n6_20261012_215000/`。
探测二进制 SHA-256：
`a89ecabb64a1a235c0e65e8a09e5c4ae5d765084e1c36e1bc8c695b229a0c13d`。
提交前曾因登录目录在不允许提交的 base 文件系统而被拒绝，没有产生计算作业；
失败信息保留在远端 `submit-cwd-failure.log`，切换到测试目录后才提交上述成功作业。

## 2026-10-12 A4 单线程 NVMe 与热缓存 I/O

通过 `wzs@a4.rabbitbio.org` 实测，主机为 `a4.localdomain`，
`/home/bigssd` 是 `/dev/nvme1n1` 上的 ext4，设备型号
`WUS4BB076D7P3E3`，容量约 7.0 TiB。启动时空闲约 1.4 TiB，
可用内存约 181 GiB，load average 为 0.35/0.31/0.28；瞬时 iostat
显示该 NVMe 没有读写。其他用户仍有轻量进程，不称为严格独占机器。
没有修改 sysctl、清理全局 page cache、取消旧实验或更改 app/SDK。

使用系统 `fio 3.19`，`psync`、一个工作线程、`iodepth=1`，顺序位置读写
32 GiB 普通文件；64 KiB、8 MiB 两种请求大小，每种模式三次。
结果为十进制 GB/s，按有效数据字节数除以 fio 的 I/O runtime，再取中位数。
打开、预分配和校验不计入带宽时间；带同步的写入使用 `end_fsync=1`，
末尾同步包含在 I/O runtime 内。不主动 fsync 的补测每次创建新输出文件，
`end_fsync=0`，完成后在计时外同步该文件，再删除；没有用 `/dev/null`、
立即删除未同步的文件或只保留循环覆盖窗口来替代完整输出。

| A4 单工作线程顺序 I/O | 64 KiB 请求 / GB/s | 8 MiB 请求 / GB/s |
|---|---:|---:|
| `O_DIRECT` 读，绕过 OS page cache | **1.031824** | **2.320663** |
| `O_DIRECT` 写 + 末尾 fsync | **0.875340** | **1.564722** |
| 缓冲写 + 末尾 fsync | **1.223768** | **1.239887** |
| 缓冲读，先请求失效本文件缓存 | **1.824346** | **1.812796** |
| 随即缓存复读，不请求失效 | **6.906480** | **6.535997** |
| 新输出缓冲写，不主动 fsync | **1.609582** | **1.630123** |

64 KiB 缓存复读三轮为 6.898161/6.913428/6.906480 GB/s，
8 MiB 为 6.535997/6.573510/6.468324；热读比较稳定。
不主动 fsync 写入后，每个 32 GiB 文件额外同步耗时 6.77--7.75 s，
不计入该模式带宽。这是文件系统接受写入的速度，不是完整持久化速度。
缓存失效是针对测试文件的建议，不能声称彻底排除所有缓存；
`O_DIRECT` 也不排除设备内部缓存。这是单线程接口带宽，不是 NVMe 的
多队列、多线程硬件峰值。相关语义见
[Linux page cache](https://docs.kernel.org/mm/page_cache.html)、
[write(2)](https://man7.org/linux/man-pages/man2/write.2.html)。

### 与单核组 RAM 的关系

此前 SWLS 六节点探测中，单 rank/MPE 的 64 KiB 单次 memcpy 约
2.41 GB/s、8 MiB 约 4.80 GB/s。这是内存搬运参考值，不能直接称为
分别测得的 RAM 文件读速和写速，也不包含页管理和接口开销。

另外核对已归档的 D1 `app_fixmate_n6_r3.log`：本地输入 body 合计
18,755,367,916 bytes，批读最大 1.383491 s；输出合计
19,289,393,734 bytes，批写最大 1.447287 s。采用
`总字节数 / 最大 rank 分项时间 / 6`，归一化到每节点，
应用 RAM 批读约 **2.259 GB/s**、批写约 **2.221 GB/s**。
这是一次真实应用的阶段计时参照，不是对各 rank 完全等量负载的独立测速；
输出字节数包含微小的 header/EOF，不能把这两个分项当成端到端时间。
写入与 CPE 压缩部分重叠，写分项也不是全部额外关键路径耗时。

因此，以当前 BGZF 批读写而言，可用约 **2.2--2.4 GB/s/核组**作为
RAM 存储接口的数量级参照。与 A4 的 64 KiB 热缓存读相比，A4 单线程
约快 3 倍；与正常新文件、不主动 fsync 的写入相比，RAM 批写约快 1.4 倍。
前者结合应用阶段计时与 I/O 微测试，只用于判断数量级，不是两边应用的
严格 I/O 阶段消融；大块 memcpy 的 4.80 GB/s 不能替代所有 app 的 I/O 速度。
所有数值按逻辑 payload 计，不把同一次复制的读与写重复计入字节数。

### 旧成绩的缓存证据与限制

旧 A4 库基准中，D1 最佳只读 `-@8` 三轮、D2 最佳只读 `-@32` 三轮
GNU time 日志都为 `fs_inputs=0`，支持这两组只读成绩主要来自热 page cache。
D2 应用补测最初 flagstat `-@0` 的 coarse 日志为 `fs_inputs=68858320`，
后续 `-@4` 为 0。不能把“重复运行、不主动清缓存”统一描述成所有命令
完全热缓存；排序工作区、较大 SAM 输出和四步中间文件会改变缓存占用。
输出每轮新建，写入也不是“热读同样的缓存命中”问题。

本次结果说明：在已缓存输入的场景，SWLS RAM 单核组 I/O 不一定快于
A4 文件接口；但仍不能据此认定两个完整存储流程已等价。
当前 RAM 最终输出复用已有 body 页，磁盘模式还需 body 读出与最终位置写入，
该差异没有通过本次微测试消除。历史应用和库成绩不改名、不覆盖。

主测试在 A4 日志时间 22:34:53--22:47:04 完成，补测为
22:48:08--22:52:24；各模式字节数和返回码检查通过，两套 32 GiB
独立全文件 CRC32C 写后复读通过，最后两个测试数据文件均已删除。
只新增 `tools/run_a4_io_bandwidth_probe.sh`，支持 `full` 和
`buffered-write` 两个模式；先完成 512 MiB 小规模验证，再跑正式测试。
`bash -n` 和 `git diff --check` 通过，没有启动额外应用矩阵。

远端两目录为 `/home/bigssd/wzs_data/samtools_bench/` 下的
`io_probe_20261012_223024`、`io_probe_warm_write_20261012_223024`。
原始 JSON、完整命令、环境、iostat、同步分项和旧基准缓存证据归档于
`build_ram_e2e/benchmarks/a4_io_20261012_223024/`。
两份实际执行脚本也保留在本地；主测试脚本 SHA-256 为
`2fcd68d0535fbcd92de27faf39c166e6348ab66252b6ceb27a4423d26514c17c`，
加入不主动 fsync 模式后的脚本为
`163515835cda88811fcfbd573369bcd3cecf12bffab83247ed626641cfe3a2af`。

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

## 2026-10-13 自动预算与统一最终输出

本节是新的存储实现验证，不覆盖上述页引用版成绩；不启动 D1/D2 正式矩阵，
也不取消正在执行的旧实验。实现顺序为公共预算/有界暂存、最终输出、fixmate、
markdup 延后数据缓存，普通 CPE 流水线保持原有处理方式。

### 实现边界

- `ResolveWorkMemory` 按节点额度、同节点 rank 数和 RAM 池占用给出工作预算，
  保留 10%（至少 256 MiB）的运行时余量。跨核组启动使用
  `SWBAM_NODE_MEMORY=88000M`；不是根据标称物理内存假定整个节点均可用。
- 边界/延后数据缓存自动取 `min(256 MiB, max(8 MiB, work / 16))`，
  fixmate 为组缓存和交换副本再划分额度。只在阈值触发时创建临时文件。
  RAM 临时文件仍保存在显式 RAM 池，禁止溢出到物理磁盘。
- `DistributedBamOutput::Write` 是 RAM/磁盘共同的最终输出循环，
  使用至多 8 MiB 工作缓冲执行完整 body 读取与按偏移写入，再包含 header/EOF。
  源页读完即可回收，但最终文件的每个字节都实际写入，不能仅登记页引用。
  磁盘默认 body 为 spool，临时目录默认 `TMPDIR` 或 `/tmp`。
- fixmate 普通组仍交给原 CPE planner/rewrite；大型保留组与边界交换改用有界缓存，
  超限后记录和处理计划分批暂存、回放。主核异常路径与从核共用标志/TLEN/tag helper。
- markdup 的普通边界探测保持不变；未决解码块超限后暂存，输出时按原顺序恢复，
  暂存文件约 64 MiB 轮换，避免已处理数据一直占用 RAM 池。
  必要判重索引仍在内存，超预算明确失败，不宣称可以处理任意规模的极端输入。

这些是应用级预算，不是操作系统 RSS 硬上限；分配器、页缓存、MPI 仍有额外占用。
本轮遵照性能优先要求，不开展极端数据实验矩阵。

### 构建与小数据

工作区构建、独立 SDK 构建/安装、无私有 RAM 模块的独立 app/示例构建、
审计 app 构建均通过。`git diff --check`、24 项 Python 脚本测试、
`ram_file_test.cpp` 和 `bounded_cache_test.cpp` 通过。
RAM 文件测试另经 AddressSanitizer/UndefinedBehaviorSanitizer 检查通过。
新增页回收测试用已占满的 32 MiB 池，完整读取源文件并写入新文件，逐字节验证，
证明回收消费页可以避免第二份完整 body 的存储峰值，不是跳过复制。

SWLS 两节点小数据为 WES_0.25G，2,668,351 条记录。两个 SDK、BAM/SAM 回环、
过滤、flagstat、io-check、普通四步串联和磁盘 dedup-workflow/dedup-pipeline
均运行并通过从核复读。sort/collate 使用 `-m 128M` 强制外排。
最终 app 的独立 Pysam 检查已确认：collate 全字段多重集与输入一致、
1,328,196 个 QNAME 组连续；1,323,599 对主比对的 mate 坐标/MC/ms 校验通过；
sort 与 fixmate 全字段多重集一致且坐标有序；markdup 保留全部记录，重复标记 15,695 条。
SDK RAM/磁盘回环与 SAM/BAM 转换的全字段多重集均与输入一致；MAPQ>=30
过滤的 2,481,669 条记录及全字段多重集均匹配。磁盘 dedup-workflow 和
dedup-pipeline 的有序输出与新 RAM 四步末端的全字段多重集、记录数和标记数一致。
独立验证使用 Pysam 0.24.1，三份日志为 `four-step.log`、`streaming.log`、
`native-filter.log`，不仅检查退出码或 BAM 记录数。

审计作业 `8529383`：最终 app 的 RAM sort（`-m 128M`，每 rank RAM 池 8G）
产生 14 个 runs、28 个 segments，实际 RAM 临时读/写累计分别为
1,777,419,983 / 1,776,932,582 bytes。两个 rank 均报告
`physical_data_io_calls=0`；最终 body 分别完整读写 120,421,598 和
121,432,067 bytes。该审计覆盖程序链接的 POSIX 数据调用，不是系统级全调用跟踪。

小数据目录为 `bounded_io_20261013` 与 `bounded_io_20261013_final`，
位于上文 SWLS `swbam_bench/` 下；后一目录复测最终版四步及磁盘 dedup 命令。
本地完整输出和独立验证日志位于 `build_bounded_io/verification/`。

### D1 六核组单次对照

使用 `tools/bounded_io_d1_probe.sh`，六节点各一个核组，真实完整 RAM 输出，
`--no-dump` 只省略 MPI 结束后的物理保存。输入、压缩级别及算法参数保持一致，
旧二进制来自 `ram_D1_matrix_20261012_010800`，新作业和日志保存在
`bounded_io_D1_20261013`。这是单次回归探测，不是三轮中位数。

实际探测参数（脚本由 SWLS 登录节点执行）：

```bash
BASE=/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data
export SWBAM_NODE_MEMORY=88000M
bash tools/bounded_io_d1_probe.sh \
  "$BASE/swbam_bench/bounded_io_D1_20261013" \
  "$BASE/swbam_bench/ram_D1_matrix_20261012_010800" \
  "$BASE/D1-WES-9G/initial/D1-WES-9G.coordinate.bam" \
  "$BASE/swbam_bench/D1_apps_current/apps_D1_n6_20261001_031800/out/name.bam" \
  "$BASE/swbam_bench/D1_apps_current/apps_D1_n6_20261001_031800/out/fixmate.sorted.bam"
```

RAM 池分别为 count 2G、SDK 回环 8G、fixmate 12G、markdup 8G；
markdup 两版均为 `--stream -m 8G`，fixmate 两版均带 mate-score 标志 `-m`。
脚本拒绝复用已有日志，重跑需新目录；此次脚本完成后不再自动提交后续矩阵。

| 项目 | 旧核心 s | 新核心 s | 旧 MPI 区间 s | 新 MPI 区间 s | 新最终装配 s |
|---|---:|---:|---:|---:|---:|
| SDK BAM 读写 | 8.242165 | 8.236799 | 15.476310 | 16.277601 | 0.816692 |
| fixmate | 34.485228 | 34.461240 | 41.755768 | 43.268465 | 1.522206 |
| markdup --stream | 47.891899 | 47.618762 | 55.147587 | 55.781325 | 0.865098 |

SDK 两版均处理 177,798,571 条记录，body 均为 11,467,657,614 bytes。
核心时间基本一致，新增端到端耗时与必要的最终读写相符；两版各 rank RAM 池峰值
逐项相同，最大 3,619,684,352 bytes。新 count 的 MPI 区间为 9.696166 s。
fixmate 两版的记录数 177,798,571、配对组数 88,498,907 相同，核心时间也基本一致。
markdup 两版记录数 177,798,571、重复标记数 50,084,605 一致，RAM 池最大峰值
同为 4,102,029,312 bytes。新版 fixmate/markdup 在 D1 均未触发延后数据暂存，
普通路径继续使用原先的从核处理方式。

三项单次对照均未观察到核心时间退化，新增最终装配为 0.82--1.52 s，且没有
出现第二份完整输出导致的 RAM 池峰值增长。这只能作为重构回归证据，不能将
小幅时间下降表述为统计显著加速。新端到端时间已包含完整最终读写；后续正式
D1/D2 矩阵需固定此路径，不能混用历史页引用版成绩。

本轮七个作业在 01:36--02:01 完成，较长墙钟主要是计时区间外的预加载；
最后作业为 `8529718`，全部 rank 返回 status=0。原始日志、作业号和汇总归档于
`build_bounded_io/verification/D1-logs/`、`jobs.tsv`、`e2e.tsv`、`driver.retry.log`。
普通数据回归已完成；没有专门构造极端 QNAME 组或极端判重窗口开展数据集级验证。

冻结新二进制 SHA-256：

| 二进制 | SHA-256 |
|---|---|
| app | `43bda2270a14a913221cae425f1d9d04ec49315ff9915c3fcc3ab652ecfda852` |
| count | `912551514173b255611d12e9cffa69f67f0f7d9418cdfa4afc5f2a3cc07fcd84` |
| bam2bam | `b62c7778622fb25b53307da546ad88fda1dbbb495caa673111e460ced84da4ce` |
