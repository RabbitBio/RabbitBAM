# CPE 流水线主核运行时

本目录实现 `swbam_cpe_runtime` 的 MPE 侧调度。真正执行的从核函数位于
`slave/core/` 和 `slave/operators/`。

## 文件

- `swbam_cpe_read_pipeline.cpp`：两套 compressed/decoded batch，重叠 MPE 读取下一批与
  CPE 处理当前批。
- `swbam_composable_decode.cpp`：把 BGZF decode kernel 包装成 Composable
  `CpeBatchOperator`。
- `swbam_cpe_write_pipeline.cpp`：两套 uncompressed/compressed batch，重叠 CPE
  压缩与 MPE 输出 batch post-processing。
- `swbam_composable_compress.cpp`：Composable BGZF compress operator。

## 扩展方式

简单工具使用 `RunComposableDecodePipeline` 或 `RunComposableRawBamPipeline`。需要融合
decode+parse+业务逻辑时，实现 `CpeBatchOperator`，提供参数准备、kernel entry、
状态校验和 MPE 归并即可。

虚接口按 batch 调用，不进入逐条记录热循环。
