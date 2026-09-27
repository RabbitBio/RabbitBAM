# SDK 与普通四步去重：最终版小数据回归

日期：2026-09-27。测试代码为 `a47c66850ac0c0141c92c8c85038994fad41d5cb`。
本轮只更新实验文档，未修改算法、SDK 或计时实现。

## 1. 配置与结论

- SWLS，单节点 6 MPI ranks，每 rank 64 CPE；`q_share`，share_size=14000，cache_size=32。
- 输入 `WES_0.25G.bam`：256,845,825 bytes，12,100 BGZF blocks，2,668,351 条记录。
- memory 输入/输出、压缩 level 1；每个性能项仅运行一次，不是正式中位数基准。
- 两个 SDK 成功；四条独立命令重新串联成功；io-check 与最终 stats 检查通过。
- collate、fixmate、sort 输出分别与先前已验证阶段输出完整 BAM `cmp` 一致。
- 最终 markdup 输出与先前同一 sorted 输入的 np=1 基线完整 BAM `cmp` 一致，
  包括 header、全部记录及 AUX、重复标记、压缩字节和 EOF。

参照基线与此前跨 rank 漏标修复见
[普通四步去重与边界回归](../refactor/去重串联_边界回归验证.md)。
不要误用其中修复前的 `workflow_smoke_20260927.markdup.bam`。

## 2. SDK 结果

| 项目 | 作业号 | core / s | 记录数 | M records/s | 输入 GB/s |
|---|---:|---:|---:|---:|---:|
| CPE 解析 bam1_t 后 count | 8319004 | 0.052791 | 2,668,351 | 50.55 | 4.87 |
| BAM read/write 回环 | 8319028 | 0.200435 | 2,668,351 | 13.31 | 1.28 |

读写输出 rank body 合计 **275,206,630 bytes**；输入+body 吞吐量为 **2.65 GB/s**。
GB 为十进制。两个示例的 core 为各 rank 流水线墙钟时间的最大值，不包含整文件
load、scan/plan、初始 body reserve、最终统计归约。读写示例只将压缩 body 写入内存，
没有组装最终完整文件，也没有 dump；不能把其 core 当成 app 的 4.1~4.6。

SDK 读写例内部要求 kept=total，但不逐字节验证记录。另运行作业 **8319081** 的
`io-check`，复用相同公共 count/transform 流水线进行独立 MPE 校验：

| 检查 | 结果 |
|---|---|
| memory / POSIX / MPI-IO read/count | 全部 PASS，records=2,668,351，decoded_bytes=788,045,119 |
| BGZF 写出、解压与校验 | PASS，每 rank 128 个合成块 |
| 无过滤 BAM 回环 | PASS，total=kept=2,668,351，output_bytes=275,253,532 |
| MAPQ >= 30 过滤回环 | PASS，kept=2,481,669，output_bytes=252,880,296 |

io-check 的 output_bytes 包含每 rank 独立验证 BAM 的 header/EOF，因此不等于 SDK
body_bytes。它校验 header、CRC/ISIZE、记录边界、计数、总字节数和 MAPQ，不逐字节
比较转换前后的所有字段；本轮没有增加 SDK 写出文件接口或完整字段对照测试。
io-check 含额外诊断及多后端操作，不作为微基准性能成绩。

## 3. 普通四步串联结果

每一步均读取本轮上一步生成的文件，未复用历史中间 BAM，也未调用内存融合命令。

| 阶段 | 作业号 | 4.3 / s | 4.1~4.6 / s | 验证 |
|---|---:|---:|---:|---|
| collate，-n 66 -m 8G | 8319046 | 0.695453 | 0.713140 | memory mode；与旧 collate BAM 一致 |
| fixmate -m | 8319050 | 0.518072 | 0.534761 | 与旧 fixmate BAM 一致 |
| sort -m 8G | 8319052 | 0.652880 | 0.671407 | memory mode；与旧 sort BAM 一致 |
| markdup --stream -l 300 -m 8G | 8319062 | 0.687252 | 0.704129 | 与旧 np=1 正确 BAM 一致 |
| **合计** | | **2.553657** | **2.623437** | 四步核心时间求和 |

4.1~4.6 排除输入 load、输出内存分配 4.5b 和最终真实 gather/dump；不是完整作业 wall。
此前对应核心合计为 2.602898 s，本次约增加 **0.79%**。单次回归未见明显性能退化，
但不能据此证明所有规模和配置的性能无变化。

markdup 日志确认 `algorithm=coordinate-stream-one-pass`，boundary halo 开启；没有
退回两遍模式。marked=15,695（pair=15,550，single=145），removed=0。
算法 tracked_peak_max=139,051,288 bytes，pending_peak_max=12,619,909 bytes；这些
不是包含完整输入/输出驻留的进程 RSS。

最终 `stats --basic` 作业 **8319068**：

- total_records=2,668,351；`is sorted: 1`；`reads duplicated: 15695`。
- `raw total sequences: 2656392` 排除了 supplementary/secondary，不应与总记录数
  混淆；其差额 11,959 为 supplementary。
- 最终 BAM 与 np=1 参照的 SHA-256 均为
  `b0ea844de76e8afac783b4790a5e0e4a6b92b93dd256102eb792fbbfdd19ddc2`。

完整 BAM 相同说明本次保留了已验证基线的全部输出，而不只是 duplicate 数相同。
本轮没有新增独立 SAMtools 四阶段重算，结论是对既有基线的正确性回归。

## 4. 复现命令与日志

本地在 `build_sunway_sduhpc/` 执行。为固定二进制，复制构建产物为本轮专用名称：

```bash
cmake --build . --target RabbitBAM-MPI swbam-sdk-optimized-count swbam-sdk-optimized-bam2bam -j 8
cp -p RabbitBAM-MPI RabbitBAM-MPI-final-a47c668
cp -p swbam-sdk-optimized-count sdk-count-final-a47c668
cp -p swbam-sdk-optimized-bam2bam sdk-bam2bam-final-a47c668

submit() {
    label=$1
    exe=$2
    shift 2
    bsub.py swls -q q_share -I -b -J "final_$label" \
      -N 1 -np 6 -cgsp 64 -share_size 14000 -cache_size 32 \
      -o "final_a47c668_${label}.log" \
      "./${exe}@online/guoshi/wzs/code/RabbitBAM/build_sunway_sduhpc" \
      "$@" > "final_a47c668_${label}.submit.log" 2>&1
}

set -e
I=../performance_test/WES_0.25G.bam
P=../performance_test/final_a47c668
APP=RabbitBAM-MPI-final-a47c668
submit count sdk-count-final-a47c668 "$I"
submit bam2bam sdk-bam2bam-final-a47c668 "$I"
submit collate "$APP" collate -i "$I" -o "$P.collate.bam" --compress-level 1 -m 8G -n 66
submit fixmate "$APP" fixmate -i "$P.collate.bam" -o "$P.fixmate.bam" --compress-level 1 -m
submit sort "$APP" sort -i "$P.fixmate.bam" -o "$P.sort.bam" --compress-level 1 -m 8G
submit markdup "$APP" markdup -i "$P.sort.bam" -o "$P.markdup.bam" --compress-level 1 -m 8G --stream -l 300
submit stats "$APP" stats -i "$P.markdup.bam" --basic
submit io_check "$APP" io-check -i "$I"
```

这段函数汇总实际逐项执行的命令；复跑前换新输出前缀，避免覆盖本轮证据。
二进制静态链接已配置的 HTSlib/libdeflate，本轮无需设置 LD_LIBRARY_PATH。

SWLS 的 `online/guoshi/wzs/code/RabbitBAM/performance_test/` 下执行过以下比较，均退出 0：

```bash
cmp final_a47c668.collate.bam workflow_smoke_20260927.collate.bam
cmp final_a47c668.fixmate.bam workflow_smoke_20260927.fixmate.bam
cmp final_a47c668.sort.bam workflow_smoke_20260927.sort.bam
cmp final_a47c668.markdup.bam workflow_smoke_20260927.markdup.np1.bam
```

本地日志：`build_sunway_sduhpc/final_a47c668_*.submit.log`。
远端日志：SWLS 对应构建目录中的 `final_a47c668_*.log`。
阶段文件保留在上面的 performance_test 目录，供后续检查。

| 固定二进制 | SHA-256 |
|---|---|
| RabbitBAM-MPI-final-a47c668 | eada8f056411f3a3e1043f6bba85b3033fb241c0310b0af3572278471e5a4120 |
| sdk-count-final-a47c668 | 788c5e60661adb831eeee5844325b76920c659e2e0570aca4038b80272c1a1df |
| sdk-bam2bam-final-a47c668 | c130b9a57f90ac83b3e4c62f13763fcae1fa87216a3f52ba593b3a5ab55c721d |

## 5. 下一步

可以开始 D1 正式 SDK 实验，再推进应用与 D2。消融开关尚未透传到 SDK 示例，
需单独补充受控配置后再测试，不能将历史两路径结果当作当前架构消融。
本轮未覆盖外排、大数据容量、跨块记录、删除模式和所有边界输入。
