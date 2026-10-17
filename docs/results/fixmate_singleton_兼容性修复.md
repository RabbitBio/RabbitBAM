# fixmate singleton 兼容性修复

更新时间：2026-10-16 09:12 CST（Asia/Shanghai）。

## 最小修改

孤立 QNAME 组只有一个 primary 时，除已有的 mate 坐标/TLEN 清理外，还应像 samtools 1.20 一样清除 `BAM_FPAIRED`。本次只改两个已有位掩码：

- `swbam-app/slave/fixmate_mpi.cpp`：普通 CPE plan 的 singleton 分支。
- `swbam-app/commands/swbam_fixmate_mpi.cpp`：有界持久组主核 replay 的 singleton 分支。

两处均从 `~(BAM_FMREVERSE | BAM_FPROPER_PAIR)` 改为 `~(BAM_FPAIRED | BAM_FMREVERSE | BAM_FPROPER_PAIR)`。不改正常配对分支、secondary/supplementary 处理、operator、批次大小、流水线调度、缓存预算、记录布局或 I/O；不增加循环、分配、复制和 kernel 启动。

## D1/D2 不重跑矩阵的依据

已直接复核 swls 原始 fixmate 日志，19 份全部满足全局 `singleton=0`，且各份六 rank 的 singleton 也均为 0：

| 原始目录 | fixmate 日志数 | 覆盖 |
| --- | ---: | --- |
| `ram_D1_processing_matrix_20261013_043000` | 6 | 独立 app 三次、四步流程三轮 |
| `ram_D1_markdup_retest_20261013_215500` | 3 | 新版 markdup 的四步重测三轮 |
| `ram_D2_processing_matrix_20261013_114500` | 6 | 独立 app 三次、四步流程三轮 |
| `ram_D2_external_20261014_165000` | 4 | 一次探测、外排四步流程三轮 |

D1 每份为 177,798,571 条记录，D2 每份为 365,079,492 条。因此这些已完成样本没有执行本次改动的 singleton 分支。无需因这个修复重跑 D1/D2 的 SDK、八项 app、外排矩阵或完整串联。

旧冻结程序、命令、日志、输出和成绩不改写，也不把旧成绩标成新版本实测。论文/最终版本说明中应记录：历史 D1/D2 使用原冻结版本，本次兼容性修复仅改变未被这些样本触发的 singleton 分支。

## 性能证据与验证

使用与已有新版 markdup app 相同的 SDK 安装、编译器和选项，在 `build_bounded_io/app_fixmate_singleton/` 独立构建，未覆盖旧构建或冻结二进制。

对比现有 `app_markdup_ready` 的两份 fixmate 对象反汇编，忽略对象文件名后：

- CPE 对象仅一条指令由 `bic r47,0x22,r49` 变为 `bic r47,0x23,r49`。
- MPE 对象也仅 singleton 位置一条指令编码变化，指令数量与相邻地址不变。
- 正常配对路径以及这两份对象中的其他指令不变。此证据说明没有引入正常配对路径的额外执行工作，不是新增的端到端性能测量或对任意平台抖动的保证。

新增 `tools/fixmate_singleton_test.cpp`，在主机上直接编译实际 CPE planner、令 `_PEN=0` 做轻量功能回归，覆盖孤立 primary 的 flag/坐标/TLEN、原记录不被修改、secondary/supplementary 不误改以及正常配对和 mate score 计划。回归通过；完整 Sunway app 构建通过。未专门构造超大持久组数据矩阵。

```bash
g++ -std=c++11 -O2 -D_PEN=0 -Iswbam-app/include -Iswbam-sdk/include \
    -Iext/htslib-1.20 tools/fixmate_singleton_test.cpp \
    swbam-app/slave/fixmate_mpi.cpp -o /tmp/swbam-fixmate-singleton-test
/tmp/swbam-fixmate-singleton-test
```

历史日志核查、两份反汇编差异和两行代码补丁归档在：

```text
build_bounded_io/verification/fixmate-singleton-20261016/
  historical-singletons.json
  cpe-disassembly.diff
  mpe-disassembly.diff
  code.patch
```

## D3 版本隔离与回归

新版 SHA-256：`fc7f37cbc87a6fc1a65b01cb87037c012da3df7851e2e28bf082f7692ea64602`。

仅通过探测脚本的 `--fixmate-app` 指向独立的 `bin/RabbitBAM-MPI-fixmate-singleton`；原基础 app、markdup 和两个 SDK 的四个冻结指纹不变。另存 `fixmate.sha256`；依赖 fixmate 的五个小数据步骤使用含新指纹的标签和独立输出目录，旧输出/日志不覆盖。

只重测 D3 的 fixmate、sort、sort 复读、markdup、markdup 复读；其余小数据关卡复用已通过结果。字段/AUX 对照仍使用未经归一化的原输出，必须消除先前 7,878 条 singleton 差异。通过后才恢复原定的全量容量关卡，不自动提交正式矩阵。

09:11 已完成本轮硬件回归与独立字段/AUX 对照：

| 小数据关卡 | 作业号 | 结果 |
| --- | --- | --- |
| 新版 fixmate | 8599980 | 六 rank 成功，200,000 条 |
| 下游 sort | 8599996 | 内排，无 RAM 临时读写 |
| sort 复读 | 8600002 | 200,000 条，sorted=1 |
| 下游 markdup | 8600017 | duplicate bit 22,259，与 samtools 一致 |
| markdup 复读 | 8600025 | 200,000 条，sorted=1 |

09:11:00 接收 a4 的 `small-verified.json`，`passed=true`：先前 7,878 条 fixmate singleton 差异全部消除，完整字段/AUX 与 samtools 一致；SDK/转换、过滤、collate、sort 内容与顺序，以及 markdup 的记录数/顺序/重复总数检查均通过。原失败证据和旧输出保留，未放宽验收规则。该 JSON 与本次小数据结果、fixmate 指纹已归档到上文验证目录。

09:11:03 已恢复原 D3 全量容量关卡，首项 count `8600047` 为 RUN，六节点 `21972-21977`，每 rank RAM 池 8G。D1/D2 不重新提交，D3 正式三轮矩阵也没有自动启动。后续实时状态见 [D3 容量与正确性探测](D3_小鼠_RAM容量与正确性探测.md)。
