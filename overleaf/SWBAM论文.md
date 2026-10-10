# SWBAM论文

[TOC]

## 2.提示词

```bash
#codex卡住的问题
ps -fu wzs | grep -E "codex"
kill xxx
```



```bash
其他行文要求如下：请参考我的项目SWBAM的具体代码还有Rabbitbam论文以及SWBWA论文的行文风格，行文要保持科技学术论文风格，使用中短句，先用一句话总结，然后直观明了阐述过程，不要有生涩词汇，用词保持高性能专业领域常用词汇，保持直观明了简洁易懂，同时一定要保持严谨不能与我的代码严重不符    
```



## 3.论文

### 0. Abstract



### 0. Introduction

```tex
1.测序数据增长 + BAM/SAM 的重要性
2.SAMtools/HTSlib、Sambamba、quickBAM、RabbitBAM 已有工作
3.引入 Sunway 异构架构
4.明确三个新的 challenge
5.引出 SWBAM 的 library + applications
6.突出核心思想：batch-level composition + record-level fusion
7.复杂 distributed apps + contributions
```

```tex
\section{Introduction}
The rapid growth of sequencing data continues to increase the computational and storage requirements of genomic data analysis~\cite{amarasinghe2020opportunities,goodwin2016coming,schmidt2017ngs}.After sequence alignment, results are commonly stored in the SAM (Sequence Alignment/Map) or BAM (Binary Alignment/Map) format~\cite{li2009samtools}.These files are subsequently processed by a series of operations, including statistics, filtering, format conversion, sorting, mate information updating, and duplicate marking, before being used by downstream analyses such as variant calling.Many analysis workflows therefore read, reorganize, and rewrite large BAM files multiple times.As sequencing datasets continue to grow, efficient BAM and SAM processing becomes increasingly important for the overall throughput of large-scale genomic analysis.
%随着测序数据规模持续增长，基因组数据分析对计算和存储资源的需求也不断提高。序列比对完成后，结果通常以 SAM（Sequence Alignment/Map）或 BAM（Binary Alignment/Map）格式进行存储。这些文件在进入下游分析之前，通常还需要经过统计、过滤、格式转换、排序、mate 信息更新以及重复标记等一系列处理。许多分析流程因此需要多次读取、重组和重新写出大规模 BAM 文件。随着测序数据集不断增大，高效的 BAM 和 SAM 数据处理对于大规模基因组分析流程的整体吞吐量变得越来越重要。

SAMtools and HTSlib provide the widely used software infrastructure for BAM and SAM processing and support operations such as format conversion, sorting, indexing, and duplicate marking~\cite{danecek2021samtools,bonfield2021htslib}.Several tools have further improved the performance of specific BAM operations.Sambamba uses multi-threaded processing for common alignment-file operations~\cite{tarasov2015sambamba}, while quickBAM provides parallel access to sorted and indexed BAM files~\cite{pitman2023quickbam}.Our previous work, RabbitBAM, focuses on general BAM I/O on multi-core CPU platforms and improves parallel efficiency through pre-parsing, parallel parsing, memory pools, and lock-free queues~\cite{yan2025rabbitbam}.These approaches demonstrate that BAM processing can benefit substantially from architecture-aware optimization.However, their parallel models are mainly designed around conventional CPU threads and do not directly address heterogeneous many-core execution and distributed BAM processing.
%SAMtools 和 HTSlib 提供了目前广泛使用的 BAM/SAM 处理基础软件，并支持格式转换、排序、索引和重复标记等操作。此后也有一些工具针对特定 BAM 操作进行了性能优化。例如，Sambamba 使用多线程方式加速常见的比对文件处理，而 quickBAM 针对已经排序并建立索引的 BAM 文件提供并行访问能力。我们此前提出的 RabbitBAM 主要面向传统多核 CPU 平台上的通用 BAM I/O，通过 pre-parsing、parallel parsing、memory pool 和 lock-free queue 提高并行效率。这些工作表明，针对硬件架构特征进行优化能够显著提升 BAM 处理性能。然而，这些方法的并行模型主要围绕传统 CPU 线程设计，并没有直接解决异构众核执行以及分布式 BAM 处理中的问题。

The new Sunway architecture provides strong computing capability but also introduces a different execution model.The SW26010-Pro processor consists of multiple Core Groups (CGs), each containing one Management Processing Element (MPE) and 64 Compute Processing Elements (CPEs), while MPI is typically used for communication among processes on different CGs~\cite{zhu2021sunway}.Previous work such as SWQC and SWBWA has shown that bioinformatics applications can achieve high performance on this architecture by moving computation-intensive operations to CPEs and redesigning the original parallel workflow for MPE-CPE execution~\cite{yan2025swqc,yan2026swbwa}.However, BAM processing introduces additional challenges because it combines compressed file I/O, record parsing, format conversion, and application-specific data dependencies.
%新一代神威架构具有较强的计算能力，同时也采用了与传统多核 CPU 不同的执行模式。SW26010-Pro 处理器由多个 Core Group（CG）组成，每个 CG 包含一个 Management Processing Element（MPE）和 64 个 Compute Processing Elements（CPE），不同 CG 上的进程通常通过 MPI 进行通信。此前的 SWQC 和 SWBWA 已经表明，通过将计算密集型操作下放到 CPE，并针对 MPE-CPE 执行模式重新设计原有并行流程，可以在该架构上获得较好的生物信息计算性能。然而，BAM 数据处理还包含压缩文件 I/O、record 解析、格式转换以及不同应用特有的数据依赖，因此面临更加复杂的并行设计问题。

Efficient BAM processing on Sunway faces three main challenges.First, BAM and SAM libraries naturally expose record-oriented operations, whereas launching and synchronizing CPE tasks for individual records would introduce excessive fine-grained overhead.A coarser execution unit is therefore needed to efficiently utilize the CPE cores.Second, distributed BAM processing requires scalable data access and output organization.Multiple MPI ranks should directly process their assigned input regions and generate a valid global output without relying on a single rank to distribute or collect the complete file.Third, BAM applications have very different data dependencies.Simple operations such as statistics and filtering can be performed independently on individual batches, whereas sorting requires global data redistribution, fixmate requires QNAME-based grouping, and duplicate marking requires coordinate-ordered cross-record and cross-rank state.Implementing separate I/O and MPE-CPE execution frameworks for every operation would introduce substantial code duplication, while a fine-grained general-purpose abstraction could add overhead to the record-processing hot path.
%在神威平台上实现高效 BAM 处理主要面临三个方面的挑战。首先，BAM 和 SAM 库通常以 record 为基本处理对象，而如果直接针对每条 record 启动和同步 CPE 任务，会产生较大的细粒度调度和同步开销，因此需要更大的执行粒度以充分利用 CPE。其次，分布式 BAM 处理还需要高效的数据访问和输出组织机制。多个 MPI 进程应该能够直接处理各自负责的输入区域，并生成合法的全局输出，而不应依赖单个进程完成完整数据的分发和汇总。第三，不同 BAM 应用具有完全不同的数据依赖关系。统计和过滤等简单操作可以基于独立 batch 完成，而 sorting 需要全局数据重分布，fixmate 需要基于 QNAME 的分组，duplicate marking 则需要维护按坐标组织的跨 record 和跨 MPI 进程状态。如果针对每个操作单独实现一套 I/O 和 MPE-CPE 执行框架，会产生大量重复代码；而如果引入过于细粒度的通用抽象，又可能增加 record 处理热路径上的额外开销。

To address these challenges, we develop SWBAM, a high-performance BAM and SAM processing library and application framework for the Sunway many-core platform.At the library level, SWBAM provides four batch-based MPE-CPE pipelines for BAM reading, BAM writing, SAM reading, and SAM writing.The library separates storage access from computation and supports reusable I/O, batch management, CPE scheduling, and output components.Different processing operators can be attached to these pipelines according to the required operation.SWBAM also optimizes the major CPE kernels for BGZF compression and decompression, BAM parsing, SAM parsing and formatting, and record serialization.At the application level, these common components are reused to implement both streaming operations and distributed BAM algorithms, including statistics, filtering, BAM/SAM conversion, sorting, collation, fixmate, and duplicate marking.The distributed applications further manage global data redistribution, cross-rank communication, and algorithm-specific state.
%为了解决这些问题，我们提出 SWBAM，一个面向神威众核平台的高性能 BAM/SAM 处理库和应用框架。在 library 层，SWBAM 提供四条基于 batch 的 MPE-CPE 流水线，分别用于 BAM read、BAM write、SAM read 和 SAM write。该库将底层数据访问与计算过程分离，并提供可复用的 I/O、batch 管理、CPE 调度和输出组件。根据不同处理需求，可以为这些流水线选择不同的 operator。SWBAM 还针对 BGZF 压缩与解压、BAM parsing、SAM parsing 和 formatting 以及 record serialization 等主要 CPE kernel 进行了优化。在 application 层，这些公共组件被进一步复用，用于实现 streaming operations 和 distributed BAM algorithms，包括统计、过滤、BAM/SAM 格式转换、sorting、collation、fixmate 和 duplicate marking。对于分布式应用，应用层进一步负责全局数据重分布、跨 rank 通信以及算法相关状态的维护。

A key design principle of SWBAM is to provide flexibility at batch boundaries while keeping record-level processing efficient.Pipeline and operator selection is performed at the batch level, whereas performance-critical record operations are statically combined inside CPE kernels.For example, BGZF decompression, BAM parsing, and record-level statistics or filtering can be completed within one CPE launch.This avoids per-record runtime dispatch and unnecessary intermediate processing.Different pipelines can then be connected through MPE-side batch post-processing to construct complete operations such as BAM-to-BAM, BAM-to-SAM, and SAM-to-BAM conversion.In this way, SWBAM improves code reuse while preserving the optimized CPE execution path.
%SWBAM 的一个核心设计原则是在 batch 边界提供灵活组合，同时保持 record 级处理路径的高效执行。流水线和 operator 的选择在 batch 层完成，而性能关键的 record 操作则在 CPE kernel 内进行静态融合。例如，BGZF 解压、BAM 解析以及 record 级统计或过滤可以在一次 CPE launch 中完成。这样可以避免逐 record 的运行时分发以及不必要的中间处理。不同流水线还可以通过 MPE 侧的 batch 后处理进行连接，从而构建 BAM-to-BAM、BAM-to-SAM 和 SAM-to-BAM 等完整操作。通过这种方式，SWBAM 在提高代码复用性的同时，仍然保留了经过优化的 CPE 执行路径。

Based on this library, we further develop distributed implementations for more complex BAM operations.
Sorting and collation redistribute records among MPI ranks according to coordinate or QNAME information, respectively.Fixmate handles QNAME groups that may cross batch and rank boundaries.Duplicate marking uses a one-pass coordinate-streaming design with a sliding window, cross-rank boundary handling, and incremental output to control memory usage.These applications are further combined into an end-to-end deduplication workflow consisting of collation, fixmate, coordinate sorting, and duplicate marking.
%在该 library 基础上，我们进一步实现了多种更复杂的分布式 BAM 操作。Sorting 和 collation 分别根据坐标信息或 QNAME 信息在多个 MPI 进程之间重新分布 records。Fixmate 用于处理可能跨 batch 或跨 rank 的 QNAME groups。Duplicate marking 则采用 one-pass coordinate-streaming 设计，通过 sliding window、跨 rank 边界处理以及增量输出控制内存占用。这些应用进一步被组合成一条完整的 deduplication workflow，包括 collation、fixmate、coordinate sorting 和 duplicate marking。

The main contributions of this work are as follows:
\begin{itemize}
    \item We design SWBAM, a reusable BAM/SAM processing library for the heterogeneous Sunway many-core architecture. SWBAM provides four batch-based MPE-CPE pipelines together with common I/O, double buffering, distributed data access, and output support.
    \item We develop reusable and optimized CPE processing components for BGZF compression and decompression, BAM/SAM parsing and formatting, and record-level operations. SWBAM combines pipelines at batch boundaries and statically fuses record processing inside CPE kernels, reducing fine-grained scheduling and intermediate processing overhead.
    \item We implement distributed BAM applications on top of SWBAM, including sorting, collation, fixmate, and duplicate marking. These applications reuse the common library while introducing application-specific data redistribution, cross-rank state management, and a one-pass coordinate-streaming duplicate-marking algorithm.
    \item We integrate these applications into a complete distributed deduplication workflow and evaluate SWBAM at both the library and application levels on the Sunway platform.
\end{itemize}
%本文的主要贡献如下：
%- 我们设计了 SWBAM，一个面向神威异构众核架构的可复用 BAM/SAM 处理库。SWBAM 提供四条基于 batch 的 MPE-CPE 流水线，以及统一的 I/O、双缓冲、分布式数据访问和输出支持。
%- 我们实现了可复用并经过优化的 CPE 数据处理组件，用于 BGZF 压缩与解压、BAM/SAM parsing 和 formatting，以及 record 级操作。SWBAM 在 batch 边界组合不同流水线，并在 CPE kernel 内静态融合 record 处理，从而降低细粒度调度和中间处理开销。
%- 我们基于 SWBAM 实现了 sorting、collation、fixmate 和 duplicate marking 等分布式 BAM 应用。这些应用复用公共 library，同时结合应用相关的数据重分布、跨 rank 状态管理，以及 one-pass coordinate-streaming duplicate-marking 算法。
%- 我们进一步将这些应用整合为完整的分布式 deduplication workflow，并从 library 和 application 两个层面对 SWBAM 在神威平台上的性能进行了评估。
```

```tex
amarasinghe2020opportunities：Amarasinghe et al., Opportunities and challenges in long-read sequencing data analysis, Genome Biology, 2020。   a97b7858-44d8-4e0c-af4a-46e8764…
goodwin2016coming：Goodwin et al., Coming of age: Ten years of next-generation sequencing technologies, Nature Reviews Genetics, 2016。    a97b7858-44d8-4e0c-af4a-46e8764…
schmidt2017ngs：Schmidt and Hildebrandt, Next-generation sequencing: Big data meets high performance computing, Drug Discovery Today, 2017。    a97b7858-44d8-4e0c-af4a-46e8764…
li2009samtools：Li et al., The sequence alignment/map format and SAMtools, Bioinformatics, 2009。    a97b7858-44d8-4e0c-af4a-46e8764…
danecek2021samtools：Danecek et al., Twelve years of SAMtools and BCFtools, GigaScience, 2021。    a97b7858-44d8-4e0c-af4a-46e8764…
bonfield2021htslib：Bonfield et al., HTSlib: C library for reading/writing high-throughput sequencing data, GigaScience, 2021。    a97b7858-44d8-4e0c-af4a-46e8764…
tarasov2015sambamba：Tarasov et al., Sambamba: Fast processing of NGS alignment formats, Bioinformatics, 2015。    a97b7858-44d8-4e0c-af4a-46e8764…
pitman2023quickbam：Pitman et al., quickBAM: A parallelized BAM file access API for high-throughput sequence analysis informatics, Bioinformatics, 2023。    a97b7858-44d8-4e0c-af4a-46e8764…
yan2025rabbitbam：你的 RabbitBAM 论文，RabbitBAM: Accelerating BAM File Manipulation on Multi-Core Platforms, IEEE TCBB, 2025。    a97b7858-44d8-4e0c-af4a-46e8764…
zhu2021sunway：Zhu et al., Enabling and scaling the HPCG benchmark on the newest generation Sunway supercomputer with 42 million heterogeneous cores, SC 2021。    01c877c7-7454-4301-a26b-ef464f8…
yan2025swqc：Yan et al., SWQC: efficient sequencing data quality control on the next-generation Sunway platform, FGCS, 2025。    01c877c7-7454-4301-a26b-ef464f8…
yan2026swbwa：Yan et al., SWBWA: A Highly Efficient NGS Aligner on the New Sunway Architecture, Euro-Par 2025 / LNCS 15902, 2026。
```





### 3.1 Methods

#### 3.1.1 Overview

```tex
\section{Methods}
\subsection{Overview}
% 这里要画一个总的架构图主要是库的整体流程和库与app的联系等等
% 总体介绍 -> 采取的优化方法
As shown in Fig.~1, SWBAM consists of two main parts: the SWBAM library and a set of applications built on top of it.The library provides four batch-based processing pipelines for BAM reading, BAM writing, SAM reading, and SAM writing.These pipelines manage data input and output, batch organization, and MPE-CPE execution, providing reusable processing components for different BAM and SAM operations.
For BAM input, compressed BGZF blocks are organized into batches, while SAM input is prepared as text chunks with complete record boundaries.Each batch is processed by the CPEs according to the operation required by the current pipeline.The results are then returned to the MPE for batch-level post-processing and can be further passed to another processing pipeline when required.Based on these common pipelines, SWBAM applications implement complete functions such as format conversion, statistics and filtering, sorting, collating, fixmate, and duplicate marking.
%如图 1 所示，SWBAM 主要由两部分组成：SWBAM 库以及基于该库构建的一组应用程序。SWBAM 库提供四条基于 batch 的处理流水线，分别用于 BAM 读取、BAM 写入、SAM 读取和 SAM 写入。这些流水线负责数据输入输出、batch 组织以及 MPE-CPE 执行，为不同的 BAM 和 SAM 操作提供可复用的处理组件。对于 BAM 输入，压缩后的 BGZF block 被组织成 batch；对于 SAM 输入，则按照完整 record 边界组织为文本数据块。每个 batch 根据当前流水线所需的操作交由 CPE 进行处理。处理结果随后返回 MPE 进行 batch 级后处理，并在需要时继续传递给后续处理流水线。基于这些公共流水线，SWBAM 应用程序实现了格式转换、统计与过滤、排序、collate、fixmate 和重复标记等完整功能。

SWBAM applies performance optimizations at both the library and application levels.At the library level, batch-based MPE-CPE processing, double buffering, reusable buffers, and optimized BAM/SAM parsing and BGZF compression and decompression reduce scheduling, synchronization, and data movement overhead.The processing pipelines can also use different operators for different tasks, allowing record-level operations to be performed inside CPE kernels when appropriate.For BAM processing, SWBAM further supports multiple I/O backends and distributed data access across MPI ranks.At the application level, common pipelines and operators are reused to avoid repeated implementation of basic BAM/SAM processing.Applications that require global ordering, record grouping, or cross-record state further employ application-specific distributed algorithms and MPI communication.Together, these designs provide efficient and scalable BAM and SAM processing on the Sunway many-core platform.
%SWBAM 分别在库层和应用层进行性能优化。在库层，SWBAM 采用基于 batch 的 MPE-CPE 处理、双缓冲、可复用缓冲区，以及优化的 BAM/SAM 解析和 BGZF 压缩解压，以降低任务调度、同步和数据移动开销。不同流水线还可以根据具体任务使用不同的 operator，使适合的 record 级操作能够直接在 CPE kernel 内完成。对于 BAM 处理，SWBAM 进一步支持多种 I/O 后端以及多个 MPI 进程之间的分布式数据访问。在应用层，公共流水线和 operator 被不同应用复用，从而避免重复实现基础的 BAM/SAM 数据处理。对于需要全局排序、record 分组或跨 record 状态的应用，SWBAM 进一步设计了应用相关的分布式算法和 MPI 通信。通过这些设计，SWBAM 在神威众核平台上实现了高效且具有良好扩展性的 BAM 和 SAM 数据处理。
```



#### 3.1.2 SWBAM Library Design

##### 1. BAM/SAM input and output

```tex
\subsection{SWBAM Library Design}
\subsubsection{BAM/SAM input and output}
%这里要画一个输入输出的流程图，感觉可以跟后一节一起搞一个详细的数据流程图
%io这么设计的原因 -> BAM输入设计及优化 -> SAM输入设计及优化 -> 输出设计
SWBAM provides common input and output interfaces that separate data access from subsequent processing. Different BAM and SAM applications perform different operations, but they share similar requirements for file reading, data partitioning, and result output. Implementing these operations separately for each application would duplicate I/O logic and couple the application algorithm with a specific data access method. SWBAM therefore hides backend-specific I/O operations behind common interfaces and converts the input into batch-oriented data structures required by the processing pipelines. This design allows applications to reuse the same processing implementation when switching between different I/O methods and also provides a common basis for distributed processing across MPI ranks.
%SWBAM 提供统一的输入输出接口，将数据访问与后续处理过程分离。不同的 BAM 和 SAM 应用虽然执行不同的操作，但在文件读取、数据划分以及结果输出方面具有相似的需求。如果每个应用分别实现这些功能，不仅会产生重复的 I/O 代码，还会使应用算法与具体的数据访问方式紧密耦合。因此，SWBAM 通过统一接口隐藏不同后端的具体 I/O 操作，并将输入数据转换为后续处理流水线所需的基于 batch 的数据结构。该设计使应用在切换不同 I/O 方式时能够复用相同的处理实现，同时也为多个 MPI 进程之间的分布式处理提供统一基础。

For BAM input, SWBAM first scans the compressed data region to identify BGZF block boundaries and builds an ordered block-span index. In MPI execution, the global block index is generated by the root rank and broadcast to all ranks. The indexed BGZF blocks are then divided into contiguous ranges with approximately equal numbers of blocks, and each rank processes only its assigned range. SWBAM supports memory, POSIX, and MPI-IO BAM input backends. The memory backend loads the input file into memory and is mainly used when sufficient memory is available or when evaluating the processing pipeline with limited storage interference. For file-based access, SWBAM scans BGZF headers using an 8-MiB window, reducing the number of small file reads required to identify block boundaries. For contiguous batches, the POSIX backend uses \texttt{preadv} to read multiple BGZF blocks directly into their corresponding block buffers. The MPI-IO backend follows the same block organization and uses an MPI derived datatype to place a contiguous file region directly into multiple BGZF block buffers. These designs reduce small I/O requests and avoid an additional intermediate copy before batch processing.
%对于 BAM 输入，SWBAM 首先扫描压缩数据区域，识别 BGZF block 的边界，并建立有序的 block-span 索引。在 MPI 执行模式下，全局 block 索引由根进程建立，然后广播给所有进程。随后，索引后的 BGZF block 按照近似相同的 block 数量划分为连续区间，每个 MPI 进程只处理分配给自己的部分。SWBAM 支持 memory、POSIX 和 MPI-IO 三种 BAM 输入后端。memory 后端将完整输入文件加载到内存中，主要用于内存充足的情况，或者用于减少存储访问影响后评估核心处理流水线的性能。对于基于文件的输入，SWBAM 使用 8 MiB 的窗口扫描 BGZF block header，从而减少识别 block 边界时产生的大量小规模文件读取。对于连续的 batch，POSIX 后端使用 preadv 将多个 BGZF block 直接读取到对应的 block buffer 中。MPI-IO 后端采用相同的 block 组织方式，并使用 MPI derived datatype 将文件中的一段连续数据直接放入多个 BGZF block buffer。通过这些设计，SWBAM 减少了小规模 I/O 请求，并避免了 batch 处理之前额外的中间数据复制。

SAM input follows a similar separation between file access and processing, but uses text ranges instead of BGZF block spans. SWBAM identifies the SAM body after the header and divides it into approximately equal file ranges for distributed execution. Each range boundary is adjusted to a line boundary so that a SAM record is not split between two ranks. The assigned text range is then read through POSIX or MPI-IO and passed to the SAM read pipeline as rank-local input. In-memory input can also be used when the complete input is already available in memory. In this way, BAM and SAM use different data representations while exposing a consistent batch-oriented processing flow to the upper layers.
%SAM 输入采用类似的 I/O 与处理分离方式，但其划分单位是文本范围，而不是 BGZF block span。SWBAM 首先确定 SAM header 之后的数据区域，然后在分布式执行时按照文件位置将其划分为大小近似相同的范围。每个范围的边界会进一步调整到完整的文本行边界，从而避免一条 SAM record 被划分到两个 MPI 进程。随后，每个进程通过 POSIX 或 MPI-IO 读取自己负责的文本范围，并将其作为 rank-local 输入交给 SAM read 流水线。当完整输入已经位于内存中时，也可以采用内存输入方式。通过这种方式，BAM 和 SAM 虽然使用不同的数据表示形式，但在上层均采用统一的 batch 处理流程。

For output, SWBAM separates the generation of rank-local results from the construction of the final output file. Each rank first stores its generated output body in a rank-local sink. The body can be kept in memory, written to a temporary spool file, or selected automatically according to a configured memory limit. For distributed output, SWBAM collects the body sizes of all ranks and computes their offsets in the final file. The ranks can then write their local bodies directly to the corresponding file regions through MPI-IO, avoiding the need to gather the complete output on a single rank. For BAM output, the global layout contains the BAM header, the ordered rank-local compressed bodies, and a single BGZF EOF block. For SAM output, the text header is followed by the rank-local SAM bodies. This design supports both high-performance in-memory execution and bounded-memory processing of large files.
%对于输出，SWBAM 将各 MPI 进程本地结果的生成与最终输出文件的构建分离。每个进程首先将生成的输出数据保存到自己的 rank-local sink 中。rank-local body 可以保存在内存中、写入临时 spool 文件，或者根据配置的内存限制自动选择存储方式。在分布式输出阶段，SWBAM 收集所有进程的 body 大小，并计算每个进程在最终输出文件中的偏移位置。随后，各进程可以通过 MPI-IO 将本地 body 直接写入最终文件的对应区域，从而避免将完整输出集中到单个进程。对于 BAM 输出，最终文件依次由 BAM header、各进程有序的压缩 body 以及一个 BGZF EOF block 组成；对于 SAM 输出，则由文本 header 和各进程生成的 SAM body 组成。该设计既支持高性能的内存处理，也能够在处理大规模文件时通过有界内存方式控制内存占用。
```



##### 2. Batch-based MPE-CPE pipelines

```tex
\subsubsection{Batch-based MPE-CPE pipelines}
%这里可以跟上一节一起画一个图，有io、有流水线、有双缓冲
%为什么以batch为核心 -> 介绍4条流水线 -> 介绍每一条流水线 -> 介绍双缓冲优化
SWBAM adopts batch-based MPE-CPE pipelines to reduce fine-grained scheduling and synchronization overhead and to improve CPE utilization. Launching CPE tasks at the granularity of individual BAM or SAM records would introduce frequent MPE-CPE coordination and make it difficult to fully utilize the many-core architecture. SWBAM therefore increases the processing granularity from individual records to batches. For BAM processing, multiple BGZF blocks are grouped and processed together, while SAM processing uses batches of text chunks or BAM records according to the direction of conversion. In this way, CPE execution, synchronization, and buffer management are performed at batch boundaries rather than for each individual record.
%SWBAM 采用基于 batch 的 MPE-CPE 流水线，以减少细粒度的任务调度和同步开销，并提高 CPE 的利用率。如果在神威平台上以单条 BAM 或 SAM record 为粒度启动 CPE 任务，会引入频繁的 MPE-CPE 协调，从而难以充分利用众核架构。因此，SWBAM 将处理粒度从单条 record 提升到 batch。对于 BAM 数据，多个 BGZF block 被组织在一起进行处理；对于 SAM 数据，则根据转换方向，以文本数据块或 BAM record 批次作为处理单位。通过这种方式，CPE 执行、同步以及缓冲区管理都在 batch 边界进行，而不是针对每条 record 单独执行。

SWBAM provides four batch-based pipelines for BAM reading, BAM writing, SAM reading, and SAM writing. The BAM read pipeline accepts compressed BGZF batches and invokes a CPE operator that performs decompression and, when required, BAM parsing and application-specific processing. The BAM write pipeline accepts either uncompressed BGZF data or batches of BAM records and performs BGZF compression or record serialization followed by compression. The SAM read pipeline processes text chunks aligned to record boundaries and parses them into BAM records, while the SAM write pipeline formats batches of BAM records into SAM text. Although these pipelines operate on different data representations, they follow a common batch-oriented execution model and provide reusable processing components for upper-level applications.
%SWBAM 提供四条基于 batch 的流水线，分别用于 BAM 读取、BAM 写入、SAM 读取和 SAM 写入。BAM read 流水线接收压缩后的 BGZF batch，并调用相应的 CPE operator 完成解压；在需要时，还可以进一步执行 BAM 解析和应用相关处理。BAM write 流水线可以接收未压缩的 BGZF 数据，也可以接收 BAM record 批次，并分别完成 BGZF 压缩或者 record 序列化后再压缩。SAM read 流水线处理按照完整 record 边界组织的文本数据块，并将其解析为 BAM record；SAM write 流水线则将一批 BAM record 格式化为 SAM 文本。虽然四条流水线处理的数据表示不同，但它们均采用统一的 batch 级执行模式，并为上层应用提供可复用的处理组件。

For each batch, the MPE prepares the input buffers and operator parameters before launching the corresponding CPE kernel. After CPE execution is completed, SWBAM validates the processing status before exposing the results to MPE-side post-processing. The post-processing stage can aggregate batch-level results, submit records to another pipeline, or store application-specific intermediate data. The pipeline owns and reuses its working buffers, while data that need to remain valid across batches are stored in application-owned memory. This design defines clear data ownership and allows different BAM and SAM operations to share the same scheduling framework.
%对于每个 batch，MPE 首先准备输入缓冲区和 operator 参数，然后启动相应的 CPE kernel。CPE 执行结束后，SWBAM 首先检查处理状态，只有验证通过后才将结果交给 MPE 侧的后处理阶段。后处理阶段可以完成 batch 级结果汇总、将 record 提交给另一条流水线，或者保存应用相关的中间数据。流水线负责管理并复用自己的工作缓冲区，而需要跨 batch 保留的数据则存储在应用自身管理的内存空间中。该设计明确了不同数据的生命周期和所有权，同时使不同 BAM 和 SAM 操作能够复用相同的调度框架。

SWBAM further uses reusable alternating buffers to overlap MPE-side work with CPE execution. In the BAM read pipeline, while the CPEs process the current batch, the MPE reads the next compressed batch into another buffer. MPE-side post-processing starts after the current CPE kernel has completed and its results have been validated. In the BAM write pipeline, the MPE can write the previously compressed batch while the CPEs compress the current batch, and the next input batch can also be prepared during CPE execution when required. The SAM write pipeline similarly keeps current and pending output batches, while the SAM read pipeline uses separate CPE stages for text copy/count and record parsing and overlaps MPE-side output handling where possible. These pipeline-specific overlap strategies reduce idle time and reuse allocated buffers across batches without requiring multiple CPE kernels to execute concurrently within the same core group.
%SWBAM 进一步通过可复用的交替缓冲区，将 MPE 侧工作与 CPE 执行进行重叠。在 BAM read 流水线中，当 CPE 处理当前 batch 时，MPE 同时将下一批压缩数据读入另一个缓冲区；当前 CPE kernel 完成并通过结果验证后，MPE 才执行相应的后处理。在 BAM write 流水线中，当 CPE 压缩当前 batch 时，MPE 可以同时写出上一批已经完成压缩的数据；在需要时，还可以同时准备下一批输入。SAM write 流水线同样维护当前和待写出的输出 batch；SAM read 流水线则采用文本复制/计数和 record 解析两个 CPE 阶段，并在适合的阶段与 MPE 侧待处理输出进行重叠。通过这些针对不同流水线设计的重叠方式，SWBAM 减少了 MPE 和 CPE 的空闲时间，并在不同 batch 之间复用已分配的缓冲区，而不需要在同一核组中并发启动多个 CPE kernel。
```



##### 3. Pipeline composition and operators

```tex
\subsubsection{Pipeline composition and operators}
%这一节可以加两个SDK的伪代码展示算法流程
%怎么组合流水线 -> 从核内的operator如何组合 -> 介绍两个SDK示例 -> 一些优化
SWBAM supports flexible BAM and SAM processing by composing pipelines at batch boundaries and integrating record-level operations within CPE kernels. The four processing pipelines provide the common mechanisms for I/O, buffer management, CPE scheduling, and synchronization, while operators define the processing performed on each batch. For the BAM read pipeline, an operator prepares the batch parameters, selects the corresponding CPE kernel, validates the processing results, and performs batch-level post-processing on the MPE. These operator interfaces are invoked once for each batch rather than for each BAM record, avoiding fine-grained dispatch overhead in record processing.
%SWBAM 通过在 batch 边界组合不同流水线，并在 CPE kernel 内融合 record 级操作，实现灵活的 BAM 和 SAM 数据处理。四条处理流水线提供统一的 I/O、缓冲区管理、CPE 调度和同步机制，而 operator 则决定每个 batch 具体执行的处理操作。对于 BAM read 流水线，operator 负责准备 batch 参数、选择相应的 CPE kernel、检查处理结果，并在 MPE 上执行 batch 级后处理。这些 operator 接口每个 batch 只调用一次，而不是针对每条 BAM record 分别调用，从而避免在 record 处理过程中引入细粒度的调度开销。

Within a CPE kernel, SWBAM further combines common BAM processing steps with different record-level operations. The BAM read implementation reuses the same decompression and record traversal framework for operations such as record counting, flag statistics, basic statistics, filtering, and passthrough processing. The selected record operation is determined when the CPE kernel is built, so decompression, parsing, and the required record processing can be completed within a single CPE launch. This avoids per-record virtual calls and reduces intermediate processing between separate kernels. The operator can then return batch-level records, statistics, or application-specific metadata to the MPE for subsequent processing.
%在 CPE kernel 内部，SWBAM 进一步将通用 BAM 处理步骤与不同的 record 级操作进行组合。BAM read 的实现复用了相同的解压和 record 遍历框架，用于 record count、flag statistics、basic statistics、filter 以及 passthrough 等操作。具体的 record 操作在构建 CPE kernel 时确定，因此 BGZF 解压、BAM 解析和所需的 record 处理可以在一次 CPE launch 中完成。该设计避免了逐 record 的虚函数调用，同时减少了多个独立 kernel 之间的中间处理。随后，operator 可以将 batch 级的 records、统计结果或者应用相关的元数据返回给 MPE 进行后续处理。

Different pipelines can be further connected through MPE-side batch post-processing to construct complete BAM and SAM operations. For example, record counting only uses the BAM read pipeline with a counting operator and aggregates the batch results on the MPE. BAM-to-BAM processing connects a BAM read operator with the BAM write pipeline. The read operator performs decompression, parsing, and optional filtering, after which the MPE builds record packing plans and submits them to the write pipeline for serialization and BGZF compression. Similarly, BAM-to-SAM processing connects the BAM read pipeline with the SAM write pipeline, while SAM-to-BAM processing connects the SAM read pipeline with the BAM write pipeline. In these cases, records are exchanged between pipelines at batch boundaries rather than through a new record-level dispatch layer.
%不同流水线还可以通过 MPE 侧的 batch 后处理连接起来，从而构建完整的 BAM 和 SAM 操作。例如，record counting 仅使用 BAM read 流水线和 counting operator，并在 MPE 上汇总各 batch 的结果。BAM-to-BAM 处理则连接 BAM read 和 BAM write 两条流水线。read operator 完成解压、解析以及可选的过滤操作，随后 MPE 根据保留下来的 records 生成打包计划，并提交给 BAM write 流水线完成序列化和 BGZF 压缩。类似地，BAM-to-SAM 处理将 BAM read 流水线与 SAM write 流水线连接，而 SAM-to-BAM 则将 SAM read 流水线与 BAM write 流水线连接。在这些情况下，不同流水线之间的数据交换发生在 batch 边界，而不是引入新的逐 record 调度层。

SWBAM also defines clear ownership rules for data shared between pipelines. Pipeline buffers are reused across batches and are valid only during their corresponding processing stage. Records or intermediate data that must remain valid across batch boundaries are explicitly copied to application-owned memory. This design allows the same pipelines and operators to be reused by different applications while keeping the record-processing hot path on CPEs. As a result, SWBAM improves code reuse without introducing additional per-record abstraction overhead.
%SWBAM 还为流水线之间共享的数据定义了明确的生命周期和所有权。流水线内部的缓冲区会在不同 batch 之间重复使用，并且只在对应的处理阶段内有效。如果某些 records 或中间结果需要跨 batch 保存，则必须显式复制到由应用管理的持久内存中。通过这种设计，不同应用可以复用相同的流水线和 operator，同时将性能关键的 record 处理保留在 CPE 上。因此，SWBAM 在提高代码复用性的同时，不会在逐 record 热路径中引入额外的抽象开销。
```



##### 4. CPE kernel optimizations

```tex
\subsubsection{CPE kernel optimizations}
%讲从核侧的一些优化
%压缩解压侧的优化 -> 输入侧的融合优化 -> 输出侧的融合优化 -> sam的相关从核优化 
SWBAM further optimizes the major BAM and SAM processing kernels executed on CPEs to reduce computation, memory allocation, and intermediate data movement. BGZF compression and decompression are repeatedly invoked during continuous BAM processing, and repeatedly creating the corresponding codec objects introduces unnecessary overhead. SWBAM therefore maintains reusable libdeflate compressor and decompressor objects for individual CPEs. A decompressor is created when first used and is then reused by the same CPE across subsequent batches. Compressors are similarly cached according to the compression level and recreated only when the required level changes. For the performance-oriented level-1 compression path, SWBAM further adjusts the match-search strategy to reduce the work spent on compression while retaining valid BGZF output.
%SWBAM 进一步对运行在 CPE 上的主要 BAM 和 SAM 处理 kernel 进行优化，以减少计算开销、内存分配以及中间数据移动。BGZF 压缩和解压在连续的 BAM 处理过程中会被反复调用，如果每次都重新创建相应的编解码对象，会产生额外开销。因此，SWBAM 为各个 CPE 维护可复用的 libdeflate compressor 和 decompressor。decompressor 在第一次使用时创建，之后由同一个 CPE 在后续 batch 中重复使用。compressor 同样按照压缩级别进行缓存，只有在所需压缩级别发生变化时才重新创建。对于面向性能的 level-1 压缩路径，SWBAM 进一步调整匹配搜索策略，以减少压缩过程中需要执行的搜索工作，同时保持合法的 BGZF 输出。

BAM parsing is also performed directly on CPEs and is integrated with record-level processing. After a BGZF block is decompressed, SWBAM traverses the BAM records in the decoded block using a common parser. Frequently used record operations, including counting, flag statistics, basic statistics, filtering, and passthrough processing, reuse the same traversal framework. The selected operation is resolved when the CPE kernel is built and is executed directly during record traversal. Thus, BGZF decompression, BAM parsing, and the required record operation can be completed within one CPE launch without per-record runtime dispatch. SWBAM also reuses preallocated record buffers and scratch spaces across records and batches, reducing repeated memory allocation during parsing.
%BAM 解析同样直接在 CPE 上执行，并与 record 级处理进行融合。BGZF block 解压完成后，SWBAM 使用统一的解析器遍历解压 block 中的 BAM records。record count、flag statistics、basic statistics、filter 和 passthrough 等常用操作复用同一套 record 遍历框架。具体操作在 CPE kernel 构建时确定，并在 record 遍历过程中直接执行。因此，BGZF 解压、BAM 解析以及所需的 record 操作可以在一次 CPE launch 中完成，不需要在逐 record 热路径中进行运行时分发。SWBAM 还在不同 records 和 batches 之间复用预分配的 record buffer 和 scratch 空间，从而减少解析过程中的重复内存分配。

SWBAM applies a similar fusion strategy to BAM output. When the input is already organized as uncompressed BGZF blocks, the CPEs directly compress these blocks. When the input consists of BAM records, SWBAM first serializes the records into BGZF-ready blocks and then performs compression within the same CPE execution. This avoids constructing another complete BAM representation on the MPE between serialization and compression. The output buffers and temporary scratch spaces are reused across batches, further reducing memory management overhead.
%SWBAM 在 BAM 输出阶段采用类似的融合策略。当输入已经组织成未压缩的 BGZF block 时，CPE 直接完成这些 block 的压缩；当输入为 BAM records 时，SWBAM 首先将 records 序列化成可用于 BGZF 压缩的数据块，并在同一次 CPE 执行中继续完成压缩。这样可以避免在序列化和压缩之间由 MPE 再构建一份完整的中间 BAM 数据。输出 buffer 和临时 scratch 空间同样在不同 batch 之间重复使用，从而进一步降低内存管理开销。

The SAM pipelines are also implemented with CPE kernels for the major processing stages. SAM reading first partitions text into record-aligned chunks and uses CPE kernels for text copy/count and record parsing. SAM writing formats batches of BAM records into SAM text on CPEs. Together with the batch-based execution model described above, these kernel-level optimizations move the major parsing, formatting, compression, and decompression work to CPEs while reducing repeated allocation and intermediate processing on the MPE.
%SAM 流水线中的主要处理阶段也由 CPE kernel 完成。SAM 读取首先将文本划分为完整 record 边界的 chunks，并使用 CPE kernel 完成文本 copy/count 和 record 解析。SAM 写出则在 CPE 上将批量 BAM records 格式化为 SAM 文本。结合前面介绍的基于 batch 的执行模型，这些 kernel 级优化将主要的解析、格式化、压缩和解压工作交给 CPE，同时减少 MPE 上的重复内存分配和中间处理开销。
```



#### 3.1.3 Applications Built on SWBAM

##### 1. Streaming BAM operations

```tex
\subsection{Applications Built on SWBAM}
\subsubsection{Streaming BAM operations}
%主要介绍简单的这种流式处理
%可以画一张图介绍这些流程，不用伪代码
%什么是流式处理 -> 统计类操作 -> bam2bam -> bam2sam和sam2bam -> 简单总结
SWBAM implements several streaming BAM and SAM operations that process the input in a single pass without global record reordering. These operations include flag statistics, basic statistics, BAM filtering and rewriting, and format conversion between BAM and SAM. Since the processing of one batch does not depend on distant records, each MPI rank can independently process its assigned input range using the common SWBAM pipelines. Only small final results, such as statistics or output layout information, need to be exchanged among MPI ranks.
%SWBAM 实现了多种流式 BAM 和 SAM 操作，这些操作只需要顺序扫描一次输入数据，不需要进行全局 record 重排。这类操作包括 flag statistics、basic statistics、BAM 过滤与重写，以及 BAM 和 SAM 之间的格式转换。由于一个 batch 的处理通常不依赖较远位置的 records，因此每个 MPI 进程可以使用 SWBAM 的公共流水线独立处理分配给自己的输入范围。只有统计结果或最终输出布局等少量信息需要在 MPI 进程之间进行交换。

For statistical operations, SWBAM uses the BAM read pipeline with application-specific CPE operators. The flagstat and basic-statistics operators combine BGZF decompression, BAM parsing, and record-level statistics within the CPE kernel. Each batch produces rank-local counters, which are accumulated on the MPE after CPE processing. After all local batches are processed, the rank-local results are combined through MPI reduction to generate the final statistics. This design avoids materializing an additional BAM stream and keeps the main record processing on CPEs.
%对于统计类操作，SWBAM 使用 BAM read 流水线配合应用相关的 CPE operator。flagstat 和 basic-statistics operator 在 CPE kernel 内将 BGZF 解压、BAM 解析和 record 级统计融合执行。每个 batch 生成当前进程的局部统计结果，并在 CPE 处理完成后由 MPE 进行累加。当所有本地 batch 处理完成后，再通过 MPI reduction 合并各进程的局部结果，生成最终统计信息。该设计不需要额外生成中间 BAM 数据，并将主要的逐 record 处理保留在 CPE 上执行。

BAM filtering and rewriting are implemented by connecting the BAM read and BAM write pipelines. The read-side transform operator performs decompression, parsing, and optional record filtering. The retained records are then organized into record packing plans on the MPE and submitted to the BAM write pipeline, where record serialization and BGZF compression are performed on CPEs. For records that are removed by the filter, no output record is generated. The same pipeline can also operate in passthrough mode when no filtering is required. Each MPI rank produces its own ordered output body, which is later placed at the corresponding position in the final BAM file.
%BAM 过滤和重写通过连接 BAM read 和 BAM write 两条流水线实现。读端的 transform operator 负责执行解压、解析以及可选的 record 过滤。保留下来的 records 随后由 MPE 组织为 record packing plans，并提交给 BAM write 流水线，由 CPE 完成 record 序列化和 BGZF 压缩。对于被过滤掉的 records，则不会产生对应的输出记录。当不需要过滤时，同一套流水线也可以使用 passthrough 模式。每个 MPI 进程生成自己保持顺序的输出 body，之后再写入最终 BAM 文件中的相应位置。

Format conversion is implemented by composing the corresponding read and write pipelines. For BAM-to-SAM conversion, the BAM read pipeline decompresses and parses the input records using the passthrough operator, and the parsed records are passed in batches to the SAM write pipeline for parallel formatting on CPEs. For SAM-to-BAM conversion, the SAM read pipeline first identifies complete SAM records and parses them into BAM records, after which the BAM write pipeline serializes and compresses the records into BGZF blocks. The intermediate records are only kept for the current batch and are released or reused after the downstream pipeline completes its processing. This batch-level composition preserves record order while avoiding an additional record-level dispatch layer between the two pipelines.
%格式转换通过组合对应的读写流水线实现。对于 BAM-to-SAM 转换，BAM read 流水线使用 passthrough operator 对输入 records 进行解压和解析，然后以 batch 的方式将解析后的 records 交给 SAM write 流水线，在 CPE 上并行完成 SAM 格式化。对于 SAM-to-BAM 转换，SAM read 流水线首先识别完整的 SAM records，并将其解析为 BAM records，随后 BAM write 流水线将 records 序列化并压缩为 BGZF blocks。中间 records 仅在当前 batch 的处理过程中有效，当下游流水线完成后即可释放或复用。该 batch 级流水线组合在保持 record 顺序的同时，避免在两条流水线之间增加新的逐 record 调度层。

Compared with the distributed algorithms described below, these streaming operations require little cross-rank communication. Most computation and data movement are performed locally within each rank, while MPI is mainly used for final statistics aggregation and distributed output organization. Therefore, they directly benefit from the batch-based MPE-CPE pipelines and provide representative applications for evaluating the basic processing efficiency of SWBAM.
%与后面介绍的分布式算法相比，这些 streaming operations 只需要很少的跨 MPI 进程通信。大部分计算和数据移动都在各进程本地完成，MPI 主要用于最终统计结果汇总以及分布式输出组织。因此，这类操作能够直接利用 SWBAM 基于 batch 的 MPE-CPE 流水线，同时也适合作为评价 SWBAM 基础处理性能的代表性应用。
```



##### 2. Distributed BAM processing

```tex
\subsubsection{Distributed BAM processing}
%前面的算法可以总的有一张图，最后markdup可以用伪代码描述
%先总体介绍 -> sort/collate -> fixmate -> markdup
SWBAM further implements distributed BAM operations that require global data reorganization or cross-record state. Unlike streaming operations, these applications cannot process each MPI rank completely independently because records belonging to the same global order, QNAME group, or duplicate group may be distributed across different ranks. SWBAM therefore reuses the common BAM read and write pipelines for local data processing, while the application layer manages algorithm-specific metadata, MPI communication, and global state. This design keeps BGZF I/O and CPE processing reusable while allowing each application to adopt a distributed algorithm suited to its data dependency.
%SWBAM 进一步实现了需要全局数据重组或跨 record 状态的分布式 BAM 操作。与 streaming operations 不同，这些应用不能让每个 MPI 进程完全独立地处理数据，因为属于同一全局顺序、QNAME group 或 duplicate group 的 records 可能分布在不同的 MPI 进程中。因此，SWBAM 复用公共 BAM read 和 BAM write 流水线完成各进程上的本地数据处理，而应用层负责维护算法相关的 metadata、MPI 通信和全局状态。通过这种设计，BGZF I/O 和 CPE 处理可以在不同应用之间复用，同时每个应用又可以根据自身的数据依赖设计相应的分布式算法。

\paragraph{Sorting and collation}
Sorting and collation both require global redistribution of BAM records, but use different partitioning keys. For coordinate sorting, SWBAM first uses CPEs to extract the coordinate key and raw BAM data from each record. Each rank sorts its local metadata and selects samples from the ordered records. The samples are collected to determine global splitters, which divide the coordinate space among MPI ranks. Local records are then partitioned according to these splitters and exchanged with their destination ranks. Since the records sent from each source rank are already ordered, the receiving rank performs a k-way merge of the received segments instead of sorting all received records again. The final raw records are directly packed into BGZF payloads and compressed through the common BAM write pipeline. When the data exceed the configured memory limit, SWBAM generates sorted runs and performs bounded data exchange and streaming multi-way merging, allowing the same distributed sorting scheme to operate in external-memory mode.
%Sorting 和 collate 都需要在全局范围内重新分布 BAM records，但两者采用不同的划分键。对于 coordinate sorting，SWBAM 首先利用 CPE 从每条 record 中提取坐标排序键以及原始 BAM 数据。每个 MPI 进程首先对本地 metadata 进行排序，并从有序 records 中进行采样。随后汇总各进程的 samples，生成全局 splitters，并利用这些 splitters 将全局坐标空间划分给不同 MPI 进程。本地 records 根据 splitters 被划分并交换到对应的目标进程。由于来自每个源进程的数据本身已经有序，因此接收进程采用 k-way merge 合并不同的有序数据段，而不需要再次对所有接收到的 records 进行完整排序。最终的 raw records 被直接组织为 BGZF payload，并通过公共 BAM write 流水线完成压缩。当数据规模超过配置的内存限制时，SWBAM 生成多个有序 runs，并采用有界的数据交换和流式多路归并，使同一套分布式排序方法能够在 external-memory 模式下运行。

Collation follows a similar distributed organization but groups records according to QNAME instead of genomic coordinates. During BAM reading, CPEs extract the QNAME hash, QNAME-related metadata, and raw record data. Records are assigned to bins, and each bin is mapped to an owner rank so that records with the same QNAME are processed by the same rank. After local ordering, self-owned records are retained locally, while records owned by other ranks are exchanged in bounded chunks. The received segments are then merged according to their QNAME order, keeping records from the same QNAME together and preserving the required ordering within a group. Similar to sorting, SWBAM provides both in-memory and external-memory implementations. The external mode stores intermediate runs under the configured memory limit and uses streaming multi-way merging before the results are passed to the BAM write pipeline.
%Collate 采用类似的分布式组织方式，但它根据 QNAME 而不是基因组坐标对 records 进行分组。在 BAM 读取过程中，CPE 提取 QNAME hash、与 QNAME 相关的 metadata 以及原始 record 数据。Records 首先被划分到不同 bins，每个 bin 再映射到一个 owner rank，从而使具有相同 QNAME 的 records 由同一个 MPI 进程处理。完成本地排序后，属于当前进程的 records 直接保留，而属于其他进程的数据则以有界大小的数据块进行交换。接收到的数据段随后按照 QNAME 顺序进行归并，使相同 QNAME 的 records 保持连续，并维持 group 内所需的 record 顺序。与 sorting 类似，SWBAM 同时提供 in-memory 和 external-memory 两种实现。external 模式在给定内存限制下生成中间 runs，并在最终写出之前进行流式多路归并。

\paragraph{Fixmate}
Fixmate requires records with the same QNAME to be processed as a group and therefore introduces state across batch and rank boundaries. SWBAM reads the name-collated input through the common BAM read pipeline and groups consecutive records with the same QNAME on the MPE. Complete groups are processed immediately, while incomplete groups at batch boundaries are retained until the following batch is available. A QNAME group may also cross an MPI rank boundary. SWBAM identifies these boundary groups and transfers their fragments to a selected owner rank before mate information is updated. For each complete group, a CPE planning kernel determines the required mate-field and tag updates and the output size. The MPE then computes output offsets, and a second CPE kernel rewrites the records in parallel. The current implementation updates mate coordinates, flags, template length, and the \texttt{MQ}, \texttt{MC}, and \texttt{ms} tags required by the supported fixmate mode. The rewritten records are finally serialized and compressed through the common BAM write pipeline.
%Fixmate 需要以相同 QNAME 的 records 为单位进行处理，因此需要维护跨 batch 和跨 MPI 进程的状态。SWBAM 使用公共 BAM read 流水线读取经过 name-collated 的输入，并在 MPE 上将连续且具有相同 QNAME 的 records 组织为 group。完整的 group 可以立即进行处理，而 batch 边界上不完整的 group 会被暂时保留，直到下一批数据到达。同一个 QNAME group 还可能跨越 MPI rank 边界，因此 SWBAM 会识别这些边界 group，并在更新 mate 信息之前将不同 rank 上的 group fragments 传输到选定的 owner rank。对于每个完整 group，首先由 CPE planning kernel 计算需要修改的 mate 字段、tags 以及输出长度。随后 MPE 计算输出 offset，第二个 CPE kernel 再并行完成 record rewrite。当前实现更新所支持 fixmate 模式需要的 mate coordinates、flags、template length 以及 MQ、MC 和 ms tags。修改后的 records 最后通过公共 BAM write 流水线完成序列化和压缩。

\paragraph{Duplicate marking}
Duplicate marking is implemented for coordinate-ordered BAM data and requires both cross-record comparison and cross-rank coordination.SWBAM first uses the BAM read pipeline to decode and parse input batches, while application-specific CPE processing extracts only the information required for duplicate detection.Candidates are assigned to owner ranks according to their coordinate ranges, so that most candidates remain on the local rank and only candidates belonging to another coordinate range need to be exchanged.Each owner rank maintains separate states for paired and single-read candidates and compares records with the same duplicate key to determine the record to retain.
%Duplicate marking 面向 coordinate-ordered BAM 数据，需要同时处理跨 record 比较和跨 MPI 进程协调。SWBAM 首先通过 BAM read 流水线完成输入 batch 的解码和解析，同时利用应用相关的 CPE 处理提取 duplicate detection 所需的信息。Candidates 根据其坐标范围被分配给对应的 owner rank，因此大部分 candidates 可以保留在当前进程，只有属于其他坐标区间的 candidates 才需要进行 MPI 交换。每个 owner rank 分别维护 paired 和 single-read candidates 的状态，并比较具有相同 duplicate key 的 records，从中确定需要保留的 record。

To bound memory usage, SWBAM performs duplicate marking in a one-pass coordinate-streaming manner.Candidates are processed in coordinate order using a sliding window, and candidates that can no longer interact with future records are retired from the active state.Records close to MPI rank boundaries require additional handling because their duplicate partners may reside on neighboring ranks.SWBAM therefore builds a boundary halo around rank partitions and exchanges the required candidate information before final decisions are made.Decoded blocks whose duplicate decisions are not yet complete are temporarily retained in a bounded pending queue.Once all local and remote decisions for a block are available, SWBAM updates the duplicate flags and incrementally submits the ready block to the common BAM write pipeline.
The common library manages batch scheduling and output, while the application maintains the sliding window, MPI decisions, boundary state, and pending records.
%为了控制内存使用，SWBAM 采用 one-pass coordinate-streaming 方式进行 duplicate marking。Candidates 按照坐标顺序通过 sliding window 进行处理，当某个 candidate 不可能再与后续 records 形成 duplicate 关系时，就从当前活动状态中移除。位于 MPI rank 边界附近的 records 需要额外处理，因为其 duplicate partner 可能位于相邻 rank。因此，SWBAM 在 rank 分区边界构建 boundary halo，并在最终判定之前交换所需的 candidate 信息。对于 duplicate decision 尚未完成的 decoded blocks，SWBAM 将其暂时保存在一个有界的 pending queue 中。当一个 block 的本地和远程判定全部完成后，SWBAM 更新其中的 duplicate flags，并将已经可以输出的 block 增量提交给公共 BAM write 流水线。公共库负责 batch 调度和结果写出，而应用层负责维护 sliding window、MPI decisions、边界状态以及 pending records。
```



##### 3. Deduplication workflow

```tex
\subsubsection{Deduplication workflow}
%可以简单画一个流程图
%这四步为什么要按顺序执行 -> 怎么组织起来的
SWBAM provides an end-to-end deduplication workflow by combining collation, fixmate, coordinate sorting, and duplicate marking. These four stages are executed in sequence because duplicate marking requires mate information and coordinate-ordered input, while mate information is most efficiently updated after records with the same QNAME have been grouped together. The workflow therefore first applies collate to group records by QNAME, then runs fixmate with mate-score calculation enabled, reorders the updated records by genomic coordinate, and finally performs duplicate marking.
%SWBAM 通过组合 collate、fixmate、coordinate sorting 和 duplicate marking，提供了一条完整的去重处理流程。这四个阶段需要按顺序执行，因为 duplicate marking 需要完整的 mate 信息以及按基因组坐标排序的输入，而 mate 信息则适合在具有相同 QNAME 的 records 被组织到一起之后进行更新。因此，该流程首先使用 collate 按 QNAME 对 records 进行分组，然后执行开启 mate-score 计算的 fixmate，随后按照基因组坐标重新排序更新后的 records，最后进行 duplicate marking。

Each stage is implemented using the corresponding distributed SWBAM application described above and reuses the common BAM input, processing, and output interfaces. The output of one stage is written as an intermediate BAM file and used as the input of the next stage. Specifically, the workflow generates collated, fixmate-updated, and coordinate-sorted intermediate BAM files before producing the final duplicate-marked BAM file. These intermediate files can be automatically removed after they are no longer needed, or retained for debugging and validation. The memory limit and compression settings are propagated to the corresponding stages, while sort and collate can further use their external-memory implementations when required.
%每个阶段都由前面介绍的相应 SWBAM 分布式应用实现，并复用公共的 BAM 输入、处理和输出接口。一个阶段的输出会写成中间 BAM 文件，并作为下一阶段的输入。具体来说，workflow 依次生成 collate、fixmate 更新和 coordinate-sorted 的中间 BAM 文件，最终输出完成 duplicate marking 的 BAM 文件。当中间文件不再需要时可以自动删除，也可以选择保留用于调试和正确性验证。内存限制和压缩参数会传递给相应阶段，而 sort 和 collate 在需要时还可以使用各自的 external-memory 实现。
```



### 3.2 Results

#### 3.2.1 Overview

```tex
\section{Results}
\subsection{Overview}
We evaluate SWBAM at both the library and application levels to examine its processing performance, scalability, and effectiveness for complete BAM workflows.The experiments cover basic BAM input and output, MPE-CPE processing optimizations, streaming BAM operations, distributed BAM processing, and the complete deduplication workflow.
%我们从 library 和 application 两个层面对 SWBAM 进行评估，以考察其处理性能、可扩展性以及在完整 BAM workflow 中的性能表现。实验包括基础 BAM 输入输出、MPE-CPE 处理优化、streaming BAM operations、distributed BAM processing 以及完整的 deduplication workflow。

Table~\ref{tab:platforms} summarizes the experimental platforms.
SWBAM is evaluated on the SWLS Sunway platform, while SAMtools is evaluated on the A4 server equipped with two AMD EPYC 7T83 processors.
The SWLS platform provides six MPEs and 384 CPEs per node, whereas the A4 server provides 128 CPU cores.
SAMtools 1.20 is used as the main reference implementation.
Unless otherwise specified, the same logical BAM operation and comparable compression settings are used when comparing SWBAM with SAMtools.
%表~\ref{tab:platforms} 给出了实验所使用的硬件平台。SWBAM 在 SWLS 神威平台上进行测试，而 SAMtools 在配备两颗 AMD EPYC 7T83 处理器的 A4 服务器上进行测试。SWLS 每个节点包含 6 个 MPE 和 384 个 CPE，而 A4 服务器包含 128 个 CPU cores。实验采用 SAMtools 1.20 作为主要参考实现。除非特别说明，在 SWBAM 和 SAMtools 之间进行比较时，我们采用相同的 BAM 操作以及可比的压缩设置。

%表1-平台对照表
\begin{table}[!t]
\centering
\caption{Specification of the utilized hardware platforms}
\label{tab:platforms}
\begin{tabular}{lcc}
\toprule
 & SWLS (SWBAM) & A4 (SAMtools) \\
\midrule
CPU      & Sunway (rev. sw39000) & AMD EPYC 7T83 $\times$ 2 \\
Cores    & 6 MPE + 384 CPE & 128 \\
RAM      & 96 GB & 251 GiB \\
Storage  & Shared file system & 7.0 TiB NVMe SSD \\
OS       & Buildroot 2020.05 & Rocky Linux 8.9 \\
Compiler & SWGCC 7.1.0 & GCC 8.5.0 \\
\bottomrule
\end{tabular}
\end{table}

Table~\ref{tab:datasets} lists the three BAM datasets used in the experiments.
Their sizes range from 9.89~GB to 89.99~GB and contain approximately 178 million to 1.21 billion BAM records.
These datasets provide different input sizes for evaluating both basic BAM processing and more complex distributed applications.
%表~\ref{tab:datasets} 列出了实验使用的三个 BAM 数据集。数据规模从 9.89 GB 到 89.99 GB，共包含约 1.78 亿到 12.1 亿条 BAM records。这些不同规模的数据用于评估基础 BAM processing 以及更复杂的分布式应用。

%表2-实验数据介绍表
\begin{table}[!t]
\centering
\caption{Description of the utilized datasets}
\label{tab:datasets}
\begin{tabular}{llrrr}
\toprule
NCBI accession & Abbr. & \#Records & Average length & Data size \\
\midrule
SRR7890919  & D1 & 177,798,571   & 126 & 9.89G \\
SRR12898326 & D2 & 365,079,492   & 150 & 35.26G \\
SRR7890854  & D3 & 1,210,633,370 & 151 & 89.99G \\
\bottomrule
\end{tabular}
\end{table}
```



#### 3.2.2 Performance of the SWBAM Library

##### 1. BAM input and output performance

```tex
\subsection{Performance of the SWBAM Library}
\subsubsection{BAM input and output performance}
%对比库的只读和读写性能还有拓展性

```



##### 2. MPE-CPE processing performance

```tex
\subsubsection{MPE-CPE processing performance}
%针对库做类似消融实验，双缓冲等等的性能提升

```



#### 3.2.3 Performance of SWBAM Applications

##### 1. Streaming BAM operations

```tex
\subsection{Performance of SWBAM Applications}
\subsubsection{Streaming BAM operations}
%对比普通处理的性能
```



##### 2. Distributed BAM processing

```tex
\subsubsection{Distributed BAM processing}
%对比复杂算法的性能，包括内外排
```



##### 3. Deduplication Workflow Performance

```tex
\subsection{Deduplication Workflow Performance}
%对比完整去重流水线的性能

```



### 3.3 Conclusions



### 3.4 Acknowledgments



### 3.5 References









