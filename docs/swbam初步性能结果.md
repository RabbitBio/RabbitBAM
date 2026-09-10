# SWBAM 初步性能结果

## 测试基线

- 数据集：`WES_0.25G`，2,668,351 条记录，12,100 个 BGZF 数据块。
- 平台配置：6 个 MPE、384 个 CPE（`np=6`）。
- 以下结果来自已有测试记录，尚未完成统一环境下的正式重复实验。
- 加速比计算：`SAMtools time / SWBAM time`。

## 公共库微基准

| 测试 | SWBAM 时间/s | 记录吞吐量 | 解压数据吞吐量 | 说明 |
|---|---:|---:|---:|---|
| 通用 BGZF 解压 | 0.0448 |  | 17.60 GB/s | 788.0 MB 解压数据，不含逐记录遍历 |
| 通用 Raw BAM 扫描 | 0.1255 | 21.27 M records/s | 6.28 GB/s | 零拷贝 `RawBamRecordView` |
| BAM read-write 回环 | 0.3193 | 8.36 M records/s | 2.47 GB/s | 解压、解析、pack、压缩 |
| MAPQ >= 30 过滤回环 | 0.3899 | 6.84 M records/s | 2.02 GB/s | 保留 2,481,669 条记录 |
| MPI-IO BGZF 批读取 |  |  |  | 待测 |
| BGZF 压缩吞吐量 |  |  |  | 待测 |
| `np=1/2/4/6` 扩展性 |  |  |  | 待测 |

解压数据吞吐量统一按输入 BAM body 解压后的 `788,045,119` 字节计算，使用十进制
单位（`1 GB = 10^9 bytes`）。read-write 和过滤回环按输入字节计一次，不重复累加
输出字节；这些数值表示 6 个 rank 的聚合核心处理吞吐量，不包含初始磁盘加载。

## 库设计评估

| 对比项目 | 基线/s | 优化后/s | 结果 |
|---|---:|---:|---|
| I/O 库重构：flagstat | 0.052073 | 0.051994 | 性能持平，约快 0.15% |
| I/O 库重构：stats --basic | 0.098962 | 0.099174 | 回退约 0.21% |
| markdup 默认路径 -> 流式路径 | 0.815248 | 0.786508 | 约快 3.6% |
| markdup workspace 峰值 | 333.12 MB | 106.63 MB | 降低约 68% |
| collate 外排优化 | 1.550668 | 1.365146 | 1.136x，耗时下降约 12.0% |
| overlap 关闭 -> 开启 |  |  | 待正式 A/B 测试 |
| generic -> fused operator |  |  | 待同工作量 A/B 测试 |

## 应用性能

| 功能 | SWBAM/s | SAMtools/s | 初步加速比 |
|---|---:|---:|---:|
| collate `-n 66` | 0.708 | 4.538 | **6.41x** |
| fixmate `-m` | 0.585 | 1.581 | **2.70x** |
| coordinate sort | 0.675 | 2.131 | **3.16x** |
| markdup | 0.831 | 1.985 | **2.39x** |
| BAM -> SAM | 0.234 | 1.241 | **5.30x** |
| SAM -> BAM |  |  |  |
| BAM -> BAM |  |  |  |
| filter |  |  |  |
| flagstat | 0.052 |  |  |
| stats --basic | 0.099 |  |  |

## 完整去重流程

流程为 `collate -> fixmate -m -> sort -> markdup`。

| 方案 | collate/s | fixmate/s | sort/s | markdup/s | 合计/s | 初步加速比 |
|---|---:|---:|---:|---:|---:|---:|
| SWBAM 四阶段 | 0.708 | 0.585 | 0.675 | 0.831 | **2.799** | **3.66x** |
| SAMtools 四阶段 | 4.538 | 1.581 | 2.131 | 1.985 | **10.235** | 1.00x |
| SWBAM 内存 pipeline | 0.701 | 0.585 | 0.689 | 0.780 | **2.755** |  |

## 正确性与资源结果

| 项目 | 当前结果 |
|---|---|
| collate QNAME 连续错误 | 0 |
| collate READ1/READ2 顺序错误 | 0 |
| sort 坐标有序 | 通过 |
| markdup 重复记录数 | 15,695 |
| markdup mandatory fields / duplicate key set | 与基线一致 |
| flagstat 全部计数 | 与 SAMtools 一致 |
| stats --basic SN 项 | 与 SAMtools 一致 |

## 计时口径说明

当前 SWBAM 部分数字采用核心处理时间，而已有 SAMtools 数字可能包含文件 I/O，
因此表中加速比仅作为初步观察。正式论文实验应在相同高速存储或 `tmpfs` 上比较
双方完整 wall time，并至少重复 3 次报告中位数；SWBAM 核心时间另用于分析加速来源。
