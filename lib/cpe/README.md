# CPE 流水线主核运行时

本目录实现 `swbam_cpe_runtime` 的 MPE 侧调度。真正执行的从核函数位于
`slave/core/` 和 `slave/operators/`。

## 文件

- `swbam_cpe_read_pipeline.cpp`：两套 compressed/decoded batch，重叠 MPE 读取下一批与
  CPE 处理当前批；BAM2SAM 可选择在批后格式化期间预读。
- `swbam_cpe_sam_read_pipeline.cpp`：按完整 SAM 行分块，调度 CPE 计数和解析，
  以同步的 parsed batch 回调交给 MPE 后处理。
- `swbam_cpe_sam_write_pipeline.cpp`：提交 BAM 记录批次给 CPE formatter，
  保留上批文本，在下一批读取 kernel 运行时写入 `RankBodySink`。
- `swbam_composable_decode.cpp`：把 BGZF decode kernel 包装成 Composable
  `CpeBatchOperator`。
- `swbam_cpe_write_pipeline.cpp`：统一的 BAM write 调度器，使用压缩输出双缓冲，
  在 CPE 压缩本批时让 MPE 写出上批；支持拉取式 BGZF source 和推送式 batch。
- `swbam_cpe_record_write_adapter.cpp`：记录计划 operator 与 `RankBodySink`
  适配器；调用统一写调度器，CPE 仍直接序列化 `bam1_t` 并压缩。
- `swbam_composable_compress.cpp`：Composable BGZF compress operator。

写算子使用 `CpeWriteKernelSpec` 选择已打包块压缩或记录序列化+压缩，
共享 `CpeWritePipelineSession` 的压缩输出双缓冲。SAM read 使用
`SamReadKernelSpec` 描述 copy/count 与 parse 两阶段入口；SAM write 使用
`SamWriteKernelSpec` 描述 format 入口，应用既可直接填写 `RecordSlots()`，
也可用 `SubmitRecords()` 提交当前批记录指针。上述可替换入口都必须遵守现有
batch 参数布局；它们不是逐记录动态插件。

## 扩展方式

简单工具使用 `RunComposableDecodePipeline` 或 `RunComposableRawBamPipeline`。需要融合
decode+parse+业务逻辑时，实现 `CpeBatchOperator`，提供参数准备、kernel entry、
状态校验和 MPE 归并即可。

虚接口按 batch 调用，不进入逐条记录热循环。
