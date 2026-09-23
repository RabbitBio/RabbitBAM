# Samtools D2 库基准结果

## 1. 测试目的

在 D2 WGS 数据上测量 Samtools 只读扫描与 BAM read-write 回环的线程扩展性，作为
SWBAM Raw BAM 和 `bam1_t` 库路径的对照。Samtools 按正常方式直接读写 bigssd；
本文档仅记录实测口径，不将其与 SWBAM 排除文件 load/dump 的核心时间混为同一口径。

## 2. 测试环境

| 项目 | 配置 |
|---|---|
| 日期 | 2026-09-21 |
| 主机 | `hpcl-fat` |
| CPU | 2 x Intel Xeon Platinum 8260 @ 2.40 GHz |
| 物理核 | 48，每核 1 线程 |
| Samtools / HTSlib | 1.20 / 1.20 |
| 文件系统 | `/home/bigssd`，`/dev/md0` |
| 输入 | `D2-WGS-25G.coordinate.bam` |
| 输入大小 | 35,255,807,945 bytes（35.256 GB） |
| 记录数 | 365,079,492 |
| 排序状态 | coordinate sorted |

输入文件：

```text
/home/bigssd/wzs_data/D2-WGS-25G/initial/D2-WGS-25G.coordinate.bam
```

测试没有使用 `/dev/shm` 或 `/dev/null`，也没有主动清理 Linux page cache。本轮
GNU time 的 `fs_inputs=0` 表明输入主要命中系统页缓存；read-write 仍将约
38.37 GB 压缩 BAM 直接写入 bigssd。

## 3. 测试方法

公共环境：

```bash
export SAMTOOLS=/home/user_home/wzs/software/samtools/bin/samtools
export LD_LIBRARY_PATH=/home/user_home/wzs/software/libdeflate_dist/lib:\
/home/user_home/wzs/software/libdeflate_dist/lib64:${LD_LIBRARY_PATH:-}

IN=/home/bigssd/wzs_data/D2-WGS-25G/initial/D2-WGS-25G.coordinate.bam
```

测试 `-@ 0、2、4、8、16、24、32、47、48`，每个点运行 3 次，表中报告中位数。
`-@ 8` 只读复用紧邻扩展性测试前完成的三轮数据；其余点使用下述相同命令完成。

只读扫描：

```bash
/usr/bin/time -f \
  'wall_s=%e user_s=%U sys_s=%S cpu=%P maxrss_kb=%M fs_inputs=%I fs_outputs=%O exit=%x' \
  $SAMTOOLS view -@ "$threads" -c "$IN"
```

Read-write 回环：

```bash
OUT=/home/bigssd/wzs_data/samtools_bench/D2_scaling/D2.roundtrip.bam

/usr/bin/time -f \
  'wall_s=%e user_s=%U sys_s=%S cpu=%P maxrss_kb=%M fs_inputs=%I fs_outputs=%O exit=%x' \
  $SAMTOOLS view -@ "$threads" -O BAM,level=1 -o "$OUT" "$IN"

$SAMTOOLS quickcheck -v "$OUT"
```

`-@` 表示 Samtools/HTSlib 的附加压缩或解压线程数，不等同于程序最终实际使用的
CPU 核数。`-@ 47` 是针对 48 物理核主机增加的近似满核配置。

## 4. 只读扩展性

| `-@` | 三次 Wall/s | 中位数/s | 加速比 | Records/s | 输入 GB/s |
|---:|---|---:|---:|---:|---:|
| 0 | 181.01 / 180.72 / 180.62 | 180.72 | 1.00x | 2.02 M | 0.195 |
| 2 | 82.10 / 82.69 / 82.54 | 82.54 | 2.19x | 4.42 M | 0.427 |
| 4 | 42.38 / 42.05 / 42.40 | 42.38 | 4.26x | 8.61 M | 0.832 |
| 8 | 43.04 / 43.95 / 43.30 | 43.30 | 4.17x | 8.43 M | 0.814 |
| 16 | 40.95 / 41.93 / 43.54 | **41.93** | **4.31x** | **8.71 M** | **0.841** |
| 24 | 41.66 / 42.08 / 42.10 | 42.08 | 4.29x | 8.68 M | 0.838 |
| 32 | 43.58 / 42.98 / 43.89 | 43.58 | 4.15x | 8.38 M | 0.809 |
| 47 | 42.91 / 42.15 / 41.20 | 42.15 | 4.29x | 8.66 M | 0.836 |
| 48 | 43.98 / 42.56 / 43.52 | 43.52 | 4.15x | 8.39 M | 0.810 |

只读从 `-@ 4` 开始已基本饱和，`-@ 4～48` 的中位数均在 `41.93～43.58 s`。
最佳观测点为 `-@ 16`，相对 `-@ 0` 加速 `4.31x`；但它仅比 `-@ 4` 快约
1.1%，因此实用推荐值是 `-@ 4`，继续增加线程无稳定收益。

## 5. Read-write 扩展性

| `-@` | 三次 Wall/s | 中位数/s | 加速比 | Records/s | 输入 GB/s | 输入+输出 GB/s |
|---:|---|---:|---:|---:|---:|---:|
| 0 | 863.46 / 861.87 / 858.31 | 861.87 | 1.00x | 0.42 M | 0.041 | 0.085 |
| 2 | 465.04 / 470.85 / 473.29 | 470.85 | 1.83x | 0.78 M | 0.075 | 0.156 |
| 4 | 269.54 / 269.47 / 274.29 | 269.54 | 3.20x | 1.35 M | 0.131 | 0.273 |
| 8 | 165.54 / 169.97 / 163.55 | 165.54 | 5.21x | 2.21 M | 0.213 | 0.445 |
| 16 | 120.15 / 115.58 / 115.97 | 115.97 | 7.43x | 3.15 M | 0.304 | 0.635 |
| 24 | 116.56 / 117.21 / 117.94 | 117.21 | 7.35x | 3.11 M | 0.301 | 0.628 |
| 32 | 107.66 / 108.03 / 114.14 | 108.03 | 7.98x | 3.38 M | 0.326 | 0.682 |
| 47 | 95.79 / 98.27 / 100.06 | **98.27** | **8.77x** | **3.72 M** | **0.359** | **0.749** |
| 48 | 96.24 / 99.14 / 99.24 | 99.14 | 8.69x | 3.68 M | 0.356 | 0.743 |

read-write 随线程增加持续加速，到 `-@ 16` 后收益明显递减。`-@ 47` 与
`-@ 48` 的中位数仅相差 0.9%，可视为同一性能水平；考虑主程序仍需要 CPU，
在这台 48 核机器上建议使用 `-@ 47` 作为满核 Samtools 对照。若希望减少资源占用，
`-@ 16` 的时间只比最佳点高约 17.0%。

## 6. 正确性与输出

- 所有只读轮次都得到 `365,079,492` 条记录。
- 所有 read-write 输出均通过 `samtools quickcheck -v`。
- `-@ 47` 回环输出复核得到 `365,079,492` 条记录，与输入一致。
- 统一命令口径下输出大小为 `38,367,990,867～38,367,990,871` bytes；所有输出的
  BAM 结构和记录数均通过校验。

## 7. 原始日志

fat 上保留的原始计时文件：

```text
/home/bigssd/wzs_data/samtools_bench/D1_D2_final/raw_20260921_165619
/home/bigssd/wzs_data/samtools_bench/D2_scaling/raw_20260921_170704
```

所有临时回环 BAM 在完成 `quickcheck` 后已删除。
