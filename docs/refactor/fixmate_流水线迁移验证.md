# Fixmate 公共 BAM 读写流水线迁移验证

日期：2026-09-26。WES 小数据 smoke，验证端点迁移、记录生命周期和明显性能回退，
不代表全输入范围的正确性证明或稳定性能基线。

## 迁移范围

```text
公共 BAM read + passthrough decode/parse
  -> app QNAME group / pending / plan / rewrite
  -> BamRecordPackBlock -> 公共 BAM write -> RankBodySink
  -> app prefix + middle + suffix -> 分布式输出
```

- `FmReadSource` 连接 MemReader 或 backend spans；`FmReadOperator` 继承公共
  `BamReadBatchOperator`，由 `RunCpeReadPipeline` 管理缓冲、spawn/join 和批回调。
- QNAME 分组、跨批 pending、配对与标签计划、MPI 边界交换仍属于 app。
  普通记录仅在批回调内使用；未完成组和 rank 边界组保存独立的 core/data。
- 写端复用 `CpeRecordWriteSession`，由已有记录压缩 operator 接入公共 BAM write。
  保留从核序列化和压缩，以及压缩当前批时 MPE flush 上一批的双缓冲覆盖。
- `Finish()` 排空压缩结果后再释放 rewrite arena。read kernel 已 join 后才进入
  app plan/rewrite/write，不同时启动两个 CPE kernel。
- CLI、配对规则、header、EOF、输出顺序及 memory/spool 策略不变。
  standalone 和 `MpiFixmateMemoryToMemory` 共用迁移后的核心。

## 配置与性能

SWLS 单节点、6 MPI ranks、每 rank 64 CPE；memory 输入/输出、压缩 level 1。
输入为上轮已验证的 `collate_pipe_after_20260926.bam`，来自 WES_0.25G：
12,102 个 BGZF 数据块、2,668,351 条记录、1,328,196 个 QNAME group。
`fixmate -m` 表示生成 mate score，不是算法内存上限。

| 版本 | 作业号 | 4.3 / s | 4.1~4.6 / s | 结果 |
| --- | ---: | ---: | ---: | --- |
| 迁移前 | 8302786 | 0.507465 | 0.525128 | 基线 |
| 首次迁移，尚未修正 rewrite 对齐 | 8302851 | 0.505338 | 0.522896 | 一条记录的 ms 值不同，不合格 |
| 迁移及对齐修正后 | 8303040 | 0.518867 | 0.536536 | 完整 BAM 与基线逐字节一致 |

最终版本 4.3 单次变化 **+2.25%**，核心总时间 **+2.17%**，处于本轮约 5% 的目标内。
这里只做单次比较，不解释为稳定回退或加速。文件完整 load、最终 gather/dump 不计入
核心时间；4.1~4.6 还排除模拟输出内存分配 4.5b。

解压、分组及压缩 max 分别约 0.070、0.121、0.192 s，基本不变；输出索引 max 从
约 0.029 变为 0.041 s，是本次对齐布局实现的主要新增开销。各 rank max 不可直接
相加，compress 已包含被覆盖的 flush，write 不是额外独立阶段。

公共读端比旧单批读端多约 8 MiB/rank 固定 buffer；写端和记录池规模不变。
每次 rewrite 最多增加 63 个不足 64 字节的 CPE 间隔，总计不足 4 KiB。
这不包括 backend 的整文件输入、rank body 和最终输出驻留内存。

## 回归发现与修正

不能只看统计计数：首次迁移所有记录数和标签更新数一致，`stats --basic` 也能读取，
但完整 `cmp` 失败。解压后逐条比较发现仅第 1,298,912 条记录不同：

- QNAME：`SRR7890919.10801832`，`ms:i` 从 4693 变为 85，差异为一个字节。
- 该 QNAME 的四条记录跨两个读 batch，作为 stored group 一起处理；rewrite 将
  它们分给不同 CPE，输出 payload 之前却紧密相邻。
- 相邻 CPE 的非对齐拷贝可能以整字读改写形式覆盖前一条记录的尾字节。现将每个
  非空 CPE 的输出区起点对齐到 64 字节，隔离并行写区；kernel 容量检查从精确长度
  改为允许对齐间隔，`l_data` 和序列化长度仍是真实记录长度。
- 修正后整个压缩 BAM 与迁移前文件 `cmp` 返回 0，既非仅计数相等，也非忽略 AUX。

这属于迁移验证暴露的已有紧密布局风险，配对规则和 tag 计算公式没有修改。
`fixmate_pipe_after_20260926.bam` 是不合格的首次产物，后续使用 `aligned` 版本。

## 正确性与范围

- 完整输出 `cmp` 返回 0。
- 最终输出复读 `stats --basic`（8303056）：total_records=2,668,351，
  primary=2,656,392，supplementary=11,959，mismatches=878,007。
- `is sorted: 0` 符合预期：输入为 collate 输出，fixmate 不进行坐标排序。
- 配对组 1,328,196；MQ/MC/ms 更新数分别为 2,651,742 / 2,656,286 / 2,656,392；
  rank 边界合并组 3、边界字节 991；输出 BGZF 数据块 13,015，均与旧版相同。
- Sunway `RabbitBAM-MPI` 构建和 `git diff --check` 通过。

未重跑完整 dedup-pipeline，未测试所有 POSIX/MPI-IO/spool 组合、空 rank、异常输入、
压缩级别或超大 QNAME group；本轮不修改这些算法边界或宣称覆盖完整。

## 复现命令

在 `build_sunway_sduhpc` 执行：

```bash
bsub.py swls -q q_share -I -b -J fixmate_pipe_aligned \
  -N 1 -np 6 -cgsp 64 -share_size 14000 -cache_size 32 \
  -o fixmate_pipe_aligned_20260926.log \
  -x LD_LIBRARY_PATH=../ext/htslib-1.20:../ext/libdeflate-1.20/build \
  ./RabbitBAM-MPI@online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc \
  fixmate -i ../performance_test/collate_pipe_after_20260926.bam \
  -o ../performance_test/fixmate_pipe_aligned_20260926.bam \
  --compress-level 1 -m
```

基线使用迁移前二进制，同参数并将输出/日志后缀改为 `before_20260926`。
最终复读使用相同资源参数，将命令部分替换为：

```bash
stats -i ../performance_test/fixmate_pipe_aligned_20260926.bam --basic
```

在 SWLS `performance_test` 目录执行：

```bash
cmp fixmate_pipe_before_20260926.bam fixmate_pipe_aligned_20260926.bam
```

原始日志位于 SWLS `online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc/`：
`fixmate_pipe_{before,after,aligned,final_stats}_20260926.log`。
