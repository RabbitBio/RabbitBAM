# MPI Markdup

本目录实现坐标有序 BAM 的重复标记/删除，是当前状态和边界处理最复杂的实例。

## 文件分工

- `swbam_markdup_mpi.cpp`：命令编排、CPE candidate 提取、owner 判断、稀疏 remote
  exchange、flat hash 去重、bitmap 结果、记录重写、压缩与分布式输出。
- `swbam_markdup_stream_mpi.cpp`：坐标滑动窗口、winner 淘汰以及最终
  KEEP/DUP 决策生成，由主文件的 streaming 路径调用。
- 对应从核实现：`slave/algorithms/markdup_mpi.cpp`。

## 两条路径

- 默认高性能路径：批量收集 candidate 后执行分布式去重，适合内存充足的基线。
- `--stream` 默认标记模式：单遍解压和解析输入。每批 decoded block 进入待定环，
  owner rank 的滑动窗口在 winner 不再可能被后续记录替换时返回 KEEP/DUP 决策；
  本 rank 按顺序压缩并写出已经全部确定的 block。CPE 压缩与下一批 MPE
  candidate 交换、窗口判重重叠执行。记录采用 optimistic KEEP：本地 candidate
  直接更新状态，只有跨 rank candidate 才交换最终 decision；只有仍存活的 carry
  winner 进入 UNKNOWN 状态，避免为全部记录构造和排序 KEEP decision。窗口还使用
  batch-local 坐标位图预过滤 paired marker；位图只排除必然不可能命中的 marker，
  其余项仍执行精确 hash 查询，因此不改变判重语义。remote candidate 的 UNKNOWN
  设置与 owner 分类合并为一次扫描，并为本 rank 的首坐标区间提供直接命中快路径。
  为避免 rank 边界状态长期阻塞 decoded block 释放，主扫描前用每个 rank 首尾各
  32 个 BGZF block 构造 boundary halo，沿用正式 candidate/hash 逻辑预先求出
  边界记录的 KEEP/DUP 状态。若 halo 未覆盖 `max_read_length` 所需坐标范围，路径会
  自动回退到原保守边界策略，不以错误结果换取内存下降。
- `--stream -r`、`--stream -c`：暂时保留已验证的两遍路径。设置
  `RABBITBAM_MARKDUP_TWO_PASS=1` 也可让默认标记回退到两遍实现，便于正确性和
  性能对照。

单遍路径不保存完整 duplicate bitmap，也不重新读取、解压和解析输入；内存主要由
candidate 双槽、尚未完成决策的 decoded block 环、滑动窗口状态和固定输出压缩槽
组成，并统一受 `-m` 限制。若单个合法 BGZF block 解压后超过 CPE 安全压缩大小，
仅对该块按完整 BAM record 边界重新分包。

当前 6 rank、2,668,351 records 的速度优先单次结果中，优化后单遍路径的 4.3 为
0.673 s，较初版单遍路径 1.098 s 降低约 38.7%；`result` 从约 0.135 s
降至 0.005 s，`group_scan` 从约 0.186 s 降至 0.001 s，`stream_buffer`
从约 0.107 s 降至 0.001 s。输出通过
`markdup_verify_fast`，mandatory fields 与 duplicate key 指纹均与正确结果一致。
该结果比已验证的两遍 `--stream` 4.3 基线 0.768 s 快约 12.4%。

加入 boundary halo 后的低内存单次结果中，4.3 为 0.709 s；相对 0.673 s
速度优先基线慢约 5.3%，但仍比两遍基线快约 7.7%。`tracked_peak_max` 从
257,029,337 B（245.1 MiB/rank）降至 126,461,103 B（120.6 MiB/rank），下降
约 50.8%；`stream_pending_peak_max` 从 143,173,139 B（136.5 MiB）降至
12,619,908 B（12.0 MiB），下降约 91.2%。输出包含完整的 2,668,351 条记录和
12,675 个 BGZF block，并再次通过 `markdup_verify_fast`。增量输出路径同时增加
block/record 守恒检查，防止部分 batch 已准备但尚未启动压缩时被提前释放。

“BAM record 不跨 BGZF block”假设只简化单块解析，并不能消除相邻 MPI rank 在
坐标窗口内的判重依赖；boundary halo 处理的是后者。因此两项约束分别服务于解析
正确性和分布式判重正确性，不能相互替代。

`-r` 删除重复记录，`-c` 先清理已有 duplicate 标记和相关标签。输入必须是经过
fixmate `-m` 的坐标有序 BAM。

## 建议追踪

从 `ProcessMarkdupMPI()` 开始，依次理解 input plan、candidate pass、owner/group、
rewrite pass 和 output；边界 candidate 与 rank 间交换是正确性的重点。
