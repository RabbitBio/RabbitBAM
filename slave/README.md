# 神威 CPE 代码

本目录包含使用 `-mslave` 编译的从核代码。新版相关部分已经按“公共核心、可复用
算子、算法专用 kernel”分层。

## 新版目录

- `core/`：通用 BGZF codec、codec cache 和 BAM parser。
- `operators/`：通用 decode/compress 及 flagstat/stats/filter 融合算子。
- `algorithms/`：sort/collate/fixmate/markdup 的专用 kernel。

## 其他内容

- `slave_zlib/`：底层压缩依赖源码，不建议作为理解新版架构的阅读入口。
- 根目录的 `slave.cpp`、`sam_parse.*` 和 CGS 文件仍承载兼容 helper 或旧版路径；
  阅读新版时只在链接符号追踪需要时查看。

从核 object 需要直接交给神威 hybrid linker。CMake 通过 INTERFACE target 隐藏了
这一限制，普通 SDK 使用者不需要手工处理对象列表。
