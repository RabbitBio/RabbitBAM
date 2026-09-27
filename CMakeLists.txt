cmake_minimum_required(VERSION 3.13)

# Accept existing -DPLATFORM=sunway scripts, but do not offer an x86 build.
if(DEFINED PLATFORM AND NOT PLATFORM STREQUAL "sunway")
    message(FATAL_ERROR "SWBAM requires the Sunway host/CPE toolchain")
endif()
find_program(SWBAM_SUNWAY_C_COMPILER swgcc HINTS /usr/sw/bin)
find_program(SWBAM_SUNWAY_CXX_COMPILER swg++ HINTS /usr/sw/bin)
if(NOT SWBAM_SUNWAY_C_COMPILER OR NOT SWBAM_SUNWAY_CXX_COMPILER)
    message(FATAL_ERROR "Sunway compilers swgcc/swg++ were not found")
endif()
set(CMAKE_C_COMPILER "${SWBAM_SUNWAY_C_COMPILER}" CACHE FILEPATH "C compiler" FORCE)
set(CMAKE_CXX_COMPILER "${SWBAM_SUNWAY_CXX_COMPILER}" CACHE FILEPATH "C++ compiler" FORCE)
project(SWBAM LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(USE_SWLU "Enable SWLU profiling library" OFF)
option(BUILD_SWBAM_EXAMPLES "Build the two MPI SDK examples" ON)
option(BUILD_SWBAM_TOOLS "Build validation and memory/CRC probes" OFF)
set(SUNWAY_MPI_ROOT "/usr/sw/mpi/mpi_20220608_SEA" CACHE PATH "Sunway MPI installation root")
set(SUNWAY_MPI_LIB_DIR "${SUNWAY_MPI_ROOT}/lib/single_static")

# Keep the validated host/CPE code generation and hybrid linking flags.
set(COMMON_FLAGS "-g -O3 -ffast-math -DTIMING -funroll-loops -fstrict-aliasing")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_FLAGS}")
set(HOST_FLAGS -mhost -faddress_align=64)
set(SLAVE_FLAGS -mslave -msimd -faddress_align=64)
add_compile_definitions(PLATFORM_SUNWAY __STDC_LIMIT_MACROS)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib-1.20
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/libdeflate-1.20)
link_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib-1.20
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/libdeflate-1.20/build)
if(USE_SWLU)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ext/swlu/include)
    link_directories(${CMAKE_CURRENT_SOURCE_DIR}/ext/swlu)
    add_compile_definitions(USE_SWLU)
endif()

add_library(swbam_io STATIC
    lib/io/swbam_io.cpp
    lib/io/swbam_rank_body_sink.cpp)
target_compile_options(swbam_io PRIVATE ${HOST_FLAGS})
target_include_directories(swbam_io PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib-1.20)

add_library(swbam_mpi_runtime STATIC lib/mpi/swbam_mpi_runtime.cpp)
target_compile_options(swbam_mpi_runtime PRIVATE ${HOST_FLAGS})
target_include_directories(swbam_mpi_runtime PUBLIC ${SUNWAY_MPI_ROOT}/include)
target_link_libraries(swbam_mpi_runtime PUBLIC swbam_io)

add_library(swbam_cpe_runtime STATIC
    lib/cpe/swbam_cpe_bam_read_pipeline.cpp
    lib/cpe/swbam_cpe_bam_write_pipeline.cpp
    lib/cpe/swbam_cpe_record_write_adapter.cpp
    lib/cpe/swbam_cpe_sam_read_pipeline.cpp
    lib/cpe/swbam_cpe_sam_write_pipeline.cpp
    lib/operators/swbam_bgzf_compress_operator.cpp
    lib/operators/swbam_flagstat_operator.cpp
    lib/operators/swbam_record_count_operator.cpp
    lib/operators/swbam_stats_basic_operator.cpp
    lib/operators/swbam_bam_transform_operator.cpp
    lib/operators/swbam_bam_to_sam_operator.cpp
    lib/operators/swbam_sam_to_bam_operator.cpp)
target_compile_options(swbam_cpe_runtime PRIVATE ${HOST_FLAGS})
target_link_libraries(swbam_cpe_runtime PUBLIC swbam_io)

add_library(swbam_cpe_core_objects OBJECT
    slave/core/slave_bgzf_codec.cpp
    slave/core/slave_bam_parser.cpp
    slave/core/slave_bam_support.cpp)
add_library(swbam_cpe_operator_objects OBJECT
    slave/operators/slave_decode_operator.cpp
    slave/operators/slave_compress_operator.cpp
    slave/operators/slave_flagstat_operator.cpp
    slave/operators/slave_record_count_operator.cpp
    slave/operators/slave_stats_basic_operator.cpp
    slave/operators/slave_bam_filter_operator.cpp
    slave/operators/slave_record_compress_operator.cpp
    slave/operators/slave_sam_format_operator.cpp
    slave/operators/slave_sam_parse_operator.cpp)
add_library(swbam_cpe_support_objects OBJECT
    slave/sam_parse.c
    slave/slave_zlib/adler32.c
    slave/slave_zlib/crc32.c
    slave/slave_zlib/deflate_compress.c
    slave/slave_zlib/deflate_decompress.c
    slave/slave_zlib/gzip_compress.c
    slave/slave_zlib/gzip_decompress.c
    slave/slave_zlib/utils.c
    slave/slave_zlib/zlib_compress.c
    slave/slave_zlib/zlib_decompress.c)
add_library(swbam_cpe_algorithm_objects OBJECT
    slave/algorithms/sort_mpi.cpp
    slave/algorithms/collate_mpi.cpp
    slave/algorithms/markdup_mpi.cpp
    slave/algorithms/fixmate_mpi.cpp)
foreach(target swbam_cpe_core_objects swbam_cpe_operator_objects
               swbam_cpe_support_objects swbam_cpe_algorithm_objects)
    target_compile_options(${target} PRIVATE ${SLAVE_FLAGS})
    target_include_directories(${target} BEFORE PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/slave
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib-1.20
        ${CMAKE_CURRENT_SOURCE_DIR}/ext/libdeflate-1.20
        ${CMAKE_CURRENT_SOURCE_DIR}/slave/slave_zlib)
endforeach()

# CPE objects must reach the hybrid linker directly, not via a host archive.
add_library(swbam_cpe_kernels INTERFACE)
target_sources(swbam_cpe_kernels INTERFACE
    $<TARGET_OBJECTS:swbam_cpe_support_objects>
    $<TARGET_OBJECTS:swbam_cpe_core_objects>
    $<TARGET_OBJECTS:swbam_cpe_operator_objects>)
add_library(swbam_algorithm_kernels INTERFACE)
target_sources(swbam_algorithm_kernels INTERFACE
    $<TARGET_OBJECTS:swbam_cpe_algorithm_objects>)
add_library(swbam_sdk INTERFACE)
target_link_options(swbam_sdk INTERFACE -mhybrid -Wl,-zmuldefs)
target_link_directories(swbam_sdk INTERFACE ${SUNWAY_MPI_LIB_DIR} /usr/sw/lib)
target_link_libraries(swbam_sdk INTERFACE
    swbam_cpe_runtime swbam_mpi_runtime swbam_io swbam_cpe_kernels
    pthread rt deflate hts z mpicxx mpi swverbs_ft_perf swutils swtm m dl)
add_library(swbam::io ALIAS swbam_io)
add_library(swbam::mpi_runtime ALIAS swbam_mpi_runtime)
add_library(swbam::cpe_runtime ALIAS swbam_cpe_runtime)
add_library(swbam::cpe_kernels ALIAS swbam_cpe_kernels)
add_library(swbam::algorithm_kernels ALIAS swbam_algorithm_kernels)
add_library(swbam::swbam ALIAS swbam_sdk)

add_executable(RabbitBAM-MPI
    apps/mpi/main.cpp
    apps/mpi/CmdInfo.cpp
    apps/mpi/commands/swbam_bam2bam_mpi.cpp
    apps/mpi/commands/swbam_bam2sam_mpi.cpp
    apps/mpi/commands/swbam_fixmate_mpi.cpp
    apps/mpi/commands/swbam_flagstat_mpi.cpp
    apps/mpi/commands/swbam_io_check_mpi.cpp
    apps/mpi/commands/swbam_mpi_common.cpp
    apps/mpi/commands/swbam_sam2bam_mpi.cpp
    apps/mpi/commands/swbam_stats_mpi.cpp
    apps/mpi/algorithms/sort/swbam_sort_mpi.cpp
    apps/mpi/algorithms/sort/swbam_sort_optimized_mpi.cpp
    apps/mpi/algorithms/collate/swbam_collate_mpi.cpp
    apps/mpi/algorithms/markdup/swbam_markdup_mpi.cpp
    apps/mpi/algorithms/markdup/swbam_markdup_stream_mpi.cpp
    apps/mpi/pipelines/swbam_dedup_pipeline_mpi.cpp
    apps/mpi/pipelines/swbam_dedup_workflow_mpi.cpp)
target_compile_options(RabbitBAM-MPI PRIVATE ${HOST_FLAGS})
target_compile_definitions(RabbitBAM-MPI PRIVATE ENABLE_MPI_BAM2BAM)
target_include_directories(RabbitBAM-MPI PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/apps/mpi)
target_link_libraries(RabbitBAM-MPI PRIVATE swbam::swbam swbam::algorithm_kernels)
if(USE_SWLU)
    target_link_libraries(RabbitBAM-MPI PRIVATE swlu)
endif()

if(BUILD_SWBAM_EXAMPLES)
    add_executable(swbam-sdk-optimized-count examples/sdk_optimized_count.cpp)
    add_executable(swbam-sdk-optimized-bam2bam examples/sdk_optimized_bam2bam.cpp)
    foreach(target swbam-sdk-optimized-count swbam-sdk-optimized-bam2bam)
        target_compile_options(${target} PRIVATE ${HOST_FLAGS})
        target_link_libraries(${target} PRIVATE swbam::swbam)
    endforeach()
endif()

if(BUILD_SWBAM_TOOLS)
    add_executable(markdup_verify_fast tools/markdup_verify_fast.cpp)
    target_compile_options(markdup_verify_fast PRIVATE ${HOST_FLAGS})
    add_executable(swbam-mpi-memory-probe tools/mpi_memory_capacity_probe.cpp)
    add_executable(swbam-mpi-memcpy-probe tools/mpi_memcpy_probe.cpp)
    foreach(target swbam-mpi-memory-probe swbam-mpi-memcpy-probe)
        target_compile_options(${target} PRIVATE ${HOST_FLAGS})
        target_link_libraries(${target} PRIVATE swbam::swbam)
    endforeach()
    add_library(cgs_memory_probe_slave_object OBJECT tools/slave_cgs_memory_probe.cpp)
    target_compile_options(cgs_memory_probe_slave_object PRIVATE ${SLAVE_FLAGS})
    target_include_directories(cgs_memory_probe_slave_object PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tools)
    add_executable(swbam-cgs-memory-probe tools/cgs_memory_probe.cpp
        $<TARGET_OBJECTS:cgs_memory_probe_slave_object>)
    target_compile_options(swbam-cgs-memory-probe PRIVATE ${HOST_FLAGS})
    target_link_libraries(swbam-cgs-memory-probe PRIVATE pthread rt m dl)
    target_link_options(swbam-cgs-memory-probe PRIVATE -mhybrid)
    add_library(crc32_bench_slave_objects OBJECT slave/slave_crc32_bench.c)
    target_compile_options(crc32_bench_slave_objects PRIVATE ${SLAVE_FLAGS})
    target_include_directories(crc32_bench_slave_objects PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/slave ${CMAKE_CURRENT_SOURCE_DIR}/slave/slave_zlib)
    add_executable(RabbitBAM-CRC32-Bench tools/bench_crc32.cpp
        $<TARGET_OBJECTS:crc32_bench_slave_objects>)
    target_compile_options(RabbitBAM-CRC32-Bench PRIVATE ${HOST_FLAGS})
    target_link_libraries(RabbitBAM-CRC32-Bench PRIVATE pthread rt m dl)
    target_link_options(RabbitBAM-CRC32-Bench PRIVATE -mhybrid)
endif()
