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
  拉取式 source 可设置 `overlap_source=true`，同时预读/归并下一批输入；推送式
  session 可用 `SubmitWithPrefetch()` 提供同类 MPE 回调。
  状态型 app 可用 `Start()` / `Complete()` 分开启动与等待，在中间执行 MPE/MPI
  工作；输入和参数必须保持有效且不变，期间不能启动另一个 CPE kernel。
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

写端预取回调只能处理下一批的独立缓冲区，不能启动其他 CPE kernel；即使回调失败，
session 也先 join 当前 kernel 再报告失败。`Submit()` 的同步输入借用语义保持不变。
`timing.source`、`timing.post_process` 可被 `timing.kernel` 覆盖，总时间以 `total`
为准。新选项默认关闭，sort 的归并输出显式开启。

`Submit()` / `SubmitWithPrefetch()` 现在复用同一组 Start/Complete 内部逻辑，同步
接口仍在返回前 join。split 接口 Complete 之后才能复用输入或调用 Flush/Finish；
空 Start 不启动任务，也不需要 Complete。session 析构会等待未完成的 kernel，
但不隐式 flush，错误路径必须先销毁 session，再释放它借用的 app arena。
单遍 markdup 使用此接口保留“压缩已确定块 + MPE 窗口判重”的覆盖。
