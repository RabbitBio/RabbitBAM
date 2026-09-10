# RabbitBAM-MPI 应用

本目录是新版 MPI 应用层。它把公共 I/O/CPE/MPI 库与具体 BAM 功能组合成
`RabbitBAM-MPI`，应用层可以使用专用融合 kernel，但不重复实现公共 backend。

## 结构

- `main.cpp`：CLI11 参数定义、MPI/CPE 生命周期和命令分发。
- `commands/`：转换、统计、fixmate 和公共命令辅助代码。
- `algorithms/`：sort、collate、markdup 等全局算法。
- `pipelines/`：多个阶段的组合流程。

## 总体数据流

```text
CLI
-> input backend + MPI block plan
-> rank-local CPE/算法处理
-> RankBodySink
-> DistributedBamOutput
-> memory gather 或 MPI-IO 最终输出
```

应用层仍依赖 `CmdInfo` 和部分历史数据结构，因此不是 SDK 公共接口。新建简单工具
应优先参考 `examples/`；只有需要全局交换或专用融合时才参考这里。
