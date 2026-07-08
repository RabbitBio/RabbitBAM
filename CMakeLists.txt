cmake_minimum_required(VERSION 3.13)

# Platform selection (must be before project() to set compiler)
if(NOT DEFINED PLATFORM)
    set(PLATFORM "x86" CACHE STRING "Platform: x86 or sunway")
endif()

# SWLU option (only available for Sunway platform)
option(USE_SWLU "Enable SWLU profiling library" OFF)
option(BUILD_SWBAM_EXAMPLES "Build SWBAM SDK examples" ON)

if(PLATFORM STREQUAL "sunway")
    # Set Sunway compiler before project()
    find_program(SWBAM_SUNWAY_C_COMPILER swgcc HINTS /usr/sw/bin)
    find_program(SWBAM_SUNWAY_CXX_COMPILER swg++ HINTS /usr/sw/bin)
    if(NOT SWBAM_SUNWAY_C_COMPILER OR NOT SWBAM_SUNWAY_CXX_COMPILER)
        message(FATAL_ERROR "Sunway compilers swgcc/swg++ were not found")
    endif()
    set(CMAKE_C_COMPILER "${SWBAM_SUNWAY_C_COMPILER}"
        CACHE FILEPATH "C compiler" FORCE)
    set(CMAKE_CXX_COMPILER "${SWBAM_SUNWAY_CXX_COMPILER}"
        CACHE FILEPATH "C++ compiler" FORCE)
    message(STATUS "Using Sunway compiler: C=${CMAKE_C_COMPILER}, CXX=${CMAKE_CXX_COMPILER}")
endif()

project(RabbitBAM-X)

# Platform configuration
if(PLATFORM STREQUAL "x86")
    set(PLATFORM_X86 TRUE)
    message(STATUS "Building for x86 platform")
elseif(PLATFORM STREQUAL "sunway")
    set(PLATFORM_SUNWAY TRUE)
    message(STATUS "Building for Sunway platform")
else()
    message(FATAL_ERROR "Invalid PLATFORM: ${PLATFORM}. Must be 'x86' or 'sunway'")
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags
# set(COMMON_FLAGS "-g -O0 -ffast-math -DTIMING")
set(COMMON_FLAGS "-g -O3 -ffast-math -DTIMING -funroll-loops -fstrict-aliasing")

if(PLATFORM_X86)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_FLAGS}")
    add_definitions(-DPLATFORM_X86)
elseif(PLATFORM_SUNWAY)
    # Common flags for all files (no -fopenmp for Sunway)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_FLAGS}")
    add_definitions(-DPLATFORM_SUNWAY)
    # Define SIZE_MAX if not available (needed for some Sunway toolchains)
    add_definitions(-D__STDC_LIMIT_MACROS)
    # Host core flags (for src files) - 64-byte address alignment
    set(HOST_FLAGS -mhost -faddress_align=64)
    # Slave core flags (for slave files) - 64-byte address alignment
    set(SLAVE_FLAGS -mslave -msimd -faddress_align=64)
endif()

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/ext/htslib-1.20
    ${CMAKE_SOURCE_DIR}/ext/libdeflate-1.20
)

# Add SWLU include directory if USE_SWLU is enabled
if(USE_SWLU)
    include_directories(${CMAKE_SOURCE_DIR}/ext/swlu/include)
    add_definitions(-DUSE_SWLU)
    message(STATUS "SWLU profiling enabled")
endif()

# Library directories
link_directories(
    ${CMAKE_SOURCE_DIR}/ext/htslib-1.20
    ${CMAKE_SOURCE_DIR}/ext/libdeflate-1.20/build
)

# Add SWLU library directory if USE_SWLU is enabled
if(USE_SWLU)
    link_directories(${CMAKE_SOURCE_DIR}/ext/swlu)
endif()

# Collect source files from src directory
file(GLOB_RECURSE SRC_FILES "src/*.cpp")
set(MAIN_SRC "${CMAKE_SOURCE_DIR}/src/main.cpp")
set(CGS_MAIN_SRC "${CMAKE_SOURCE_DIR}/src/main_cgs.cpp")
file(GLOB CGS_ONLY_SRC "src/*cgs*.cpp")
set(SINGLE_SRC_FILES ${SRC_FILES})
list(REMOVE_ITEM SINGLE_SRC_FILES ${CGS_ONLY_SRC})
set(CGS_SRC_FILES ${SRC_FILES})
list(REMOVE_ITEM CGS_SRC_FILES ${MAIN_SRC})
set(MPI_SRC_FILES ${SRC_FILES})
list(REMOVE_ITEM MPI_SRC_FILES ${MAIN_SRC})
list(REMOVE_ITEM MPI_SRC_FILES ${CGS_ONLY_SRC})
set(SWBAM_MPI_APP_FILES
    apps/mpi/main.cpp
    apps/mpi/commands/swbam_bam2bam_mpi.cpp
    apps/mpi/commands/swbam_bam2sam_mpi.cpp
    apps/mpi/commands/swbam_fixmate_mpi.cpp
    apps/mpi/commands/swbam_flagstat_mpi.cpp
    apps/mpi/commands/swbam_io_check_mpi.cpp
    apps/mpi/commands/swbam_mpi_common.cpp
    apps/mpi/commands/swbam_sam2bam_mpi.cpp
    apps/mpi/commands/swbam_stats_mpi.cpp
    apps/mpi/algorithms/sort/swbam_sort_mpi.cpp
    apps/mpi/algorithms/sort/swbam_sort_fused_mpi.cpp
    apps/mpi/algorithms/collate/swbam_collate_mpi.cpp
    apps/mpi/algorithms/markdup/swbam_markdup_mpi.cpp
    apps/mpi/algorithms/markdup/swbam_markdup_stream_mpi.cpp
    apps/mpi/pipelines/swbam_dedup_pipeline_mpi.cpp
)
list(APPEND MPI_SRC_FILES ${SWBAM_MPI_APP_FILES})

if(PLATFORM_SUNWAY)
    set(SUNWAY_MPI_ROOT "/usr/sw/mpi/mpi_20220608_SEA" CACHE PATH "Sunway MPI installation root")
    set(SUNWAY_MPI_LIB_DIR "${SUNWAY_MPI_ROOT}/lib/single_static")

    add_library(swbam_io STATIC
        lib/io/swbam_io.cpp
        lib/io/swbam_rank_body_sink.cpp
    )
    target_compile_options(swbam_io PRIVATE ${HOST_FLAGS})
    target_include_directories(swbam_io PUBLIC
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/htslib-1.20
    )

    add_library(swbam_mpi_runtime STATIC
        lib/mpi/swbam_mpi_runtime.cpp
    )
    target_compile_options(swbam_mpi_runtime PRIVATE ${HOST_FLAGS})
    target_include_directories(swbam_mpi_runtime PUBLIC
        ${CMAKE_SOURCE_DIR}/include
        ${SUNWAY_MPI_ROOT}/include
    )
    target_link_libraries(swbam_mpi_runtime PUBLIC swbam_io)

    add_library(swbam_cpe_runtime STATIC
        lib/cpe/swbam_cpe_pipeline.cpp
        lib/cpe/swbam_cpe_write_pipeline.cpp
        lib/cpe/swbam_generic_compress.cpp
        lib/cpe/swbam_generic_decode.cpp
        lib/bam/swbam_raw_bam.cpp
        lib/bam/swbam_raw_bam_filter.cpp
        lib/bam/swbam_raw_bam_writer.cpp
    )
    target_compile_options(swbam_cpe_runtime PRIVATE ${HOST_FLAGS})
    target_include_directories(swbam_cpe_runtime PUBLIC
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/htslib-1.20
    )
    target_link_libraries(swbam_cpe_runtime PUBLIC swbam_io)

    # For Sunway, also compile slave code
    file(GLOB SLAVE_FILES "slave/*.cpp")
    file(GLOB SLAVE_C_FILES   "slave/*.c")
    file(GLOB CGS_SLAVE_FILES "slave/*cgs*.cpp")
    set(SWBAM_ALGORITHM_SLAVE_FILES
        slave/algorithms/sort_mpi.cpp
        slave/algorithms/collate_mpi.cpp
        slave/algorithms/markdup_mpi.cpp
        slave/algorithms/fixmate_mpi.cpp
    )
    list(REMOVE_ITEM SLAVE_C_FILES "${CMAKE_SOURCE_DIR}/slave/slave_crc32_bench.c")
    set(SWBAM_CPE_SUPPORT_FILES ${SLAVE_FILES})
    list(REMOVE_ITEM SWBAM_CPE_SUPPORT_FILES ${CGS_SLAVE_FILES})
    # Include slave_zlib source files (compile with -mslave)
    # Only files in slave_zlib directory, not subdirectories
    file(GLOB SLAVE_ZLIB_FILES "slave/slave_zlib/*.c")

    set(SWBAM_CPE_CORE_FILES
        slave/core/slave_bgzf_codec.cpp
        slave/core/slave_bam_parser.cpp
    )
    add_library(swbam_cpe_core_objects OBJECT ${SWBAM_CPE_CORE_FILES})
    target_compile_options(swbam_cpe_core_objects PRIVATE ${SLAVE_FLAGS})
    target_include_directories(swbam_cpe_core_objects PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/htslib-1.20
        ${CMAKE_SOURCE_DIR}/ext/libdeflate-1.20
        ${CMAKE_SOURCE_DIR}/slave/slave_zlib
    )
    target_compile_definitions(swbam_cpe_core_objects PRIVATE PLATFORM_SUNWAY)

    set(SWBAM_CPE_OPERATOR_FILES
        slave/operators/slave_decode_operator.cpp
        slave/operators/slave_compress_operator.cpp
        slave/operators/slave_flagstat_operator.cpp
        slave/operators/slave_stats_basic_operator.cpp
    )
    add_library(swbam_cpe_operator_objects OBJECT ${SWBAM_CPE_OPERATOR_FILES})
    target_compile_options(swbam_cpe_operator_objects PRIVATE ${SLAVE_FLAGS})
    target_include_directories(swbam_cpe_operator_objects PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/htslib-1.20
        ${CMAKE_SOURCE_DIR}/ext/libdeflate-1.20
        ${CMAKE_SOURCE_DIR}/slave/slave_zlib
    )
    target_compile_definitions(swbam_cpe_operator_objects PRIVATE PLATFORM_SUNWAY)

    add_library(swbam_cpe_bam_transform_objects OBJECT
        slave/operators/slave_bam_filter_operator.cpp)
    target_compile_options(swbam_cpe_bam_transform_objects PRIVATE ${SLAVE_FLAGS})
    target_include_directories(swbam_cpe_bam_transform_objects PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/htslib-1.20
        ${CMAKE_SOURCE_DIR}/ext/libdeflate-1.20
        ${CMAKE_SOURCE_DIR}/slave/slave_zlib
    )
    target_compile_definitions(swbam_cpe_bam_transform_objects PRIVATE PLATFORM_SUNWAY)
    
    # Shared CPE support and algorithm-specific kernels are separate so SDK
    # clients do not pull sort/collate/markdup/fixmate implementations.
    add_library(swbam_cpe_support_objects OBJECT
        ${SWBAM_CPE_SUPPORT_FILES}
        ${SLAVE_C_FILES}
        ${SLAVE_ZLIB_FILES})
    target_compile_options(swbam_cpe_support_objects PRIVATE ${SLAVE_FLAGS})
    target_include_directories(swbam_cpe_support_objects BEFORE PRIVATE
        ${CMAKE_SOURCE_DIR}/slave
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/htslib-1.20
        ${CMAKE_SOURCE_DIR}/ext/libdeflate-1.20
        ${CMAKE_SOURCE_DIR}/slave/slave_zlib
    )
    target_compile_definitions(swbam_cpe_support_objects PRIVATE
        PLATFORM_SUNWAY
        SWBAM_CPE_OPERATOR_SPLIT)

    add_library(swbam_cpe_algorithm_objects OBJECT
        ${SWBAM_ALGORITHM_SLAVE_FILES})
    target_compile_options(swbam_cpe_algorithm_objects PRIVATE ${SLAVE_FLAGS})
    target_include_directories(swbam_cpe_algorithm_objects BEFORE PRIVATE
        ${CMAKE_SOURCE_DIR}/slave
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/htslib-1.20
        ${CMAKE_SOURCE_DIR}/ext/libdeflate-1.20
        ${CMAKE_SOURCE_DIR}/slave/slave_zlib
    )
    target_compile_definitions(swbam_cpe_algorithm_objects PRIVATE
        PLATFORM_SUNWAY
        SWBAM_CPE_OPERATOR_SPLIT)

    add_library(slave_cgs_objects OBJECT ${SLAVE_FILES} ${SLAVE_C_FILES} ${SLAVE_ZLIB_FILES})
    target_compile_options(slave_cgs_objects PRIVATE ${SLAVE_FLAGS})
    target_include_directories(slave_cgs_objects BEFORE PRIVATE
        ${CMAKE_SOURCE_DIR}/slave
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/htslib-1.20
        ${CMAKE_SOURCE_DIR}/ext/libdeflate-1.20
        ${CMAKE_SOURCE_DIR}/slave/slave_zlib
    )
    target_compile_definitions(slave_cgs_objects PRIVATE PLATFORM_SUNWAY ENABLE_CGS_BAM2BAM ENABLE_CGS_BAM2SAM ENABLE_CGS_SAM2BAM)

    # CPE objects must reach the Sunway hybrid linker directly; placing them
    # in a normal host archive triggers "hybrid linking using master ld".
    add_library(swbam_cpe_kernels INTERFACE)
    target_sources(swbam_cpe_kernels INTERFACE
        $<TARGET_OBJECTS:swbam_cpe_support_objects>
        $<TARGET_OBJECTS:swbam_cpe_bam_transform_objects>
        $<TARGET_OBJECTS:swbam_cpe_core_objects>
        $<TARGET_OBJECTS:swbam_cpe_operator_objects>
    )

    add_library(swbam_algorithm_kernels INTERFACE)
    target_sources(swbam_algorithm_kernels INTERFACE
        $<TARGET_OBJECTS:swbam_cpe_algorithm_objects>)

    add_library(swbam_sdk INTERFACE)
    target_link_options(swbam_sdk INTERFACE
        -mhybrid
        -Wl,-zmuldefs)
    target_link_directories(swbam_sdk INTERFACE
        ${SUNWAY_MPI_LIB_DIR}
        /usr/sw/lib)
    target_link_libraries(swbam_sdk INTERFACE
        swbam_cpe_runtime
        swbam_mpi_runtime
        swbam_io
        swbam_cpe_kernels
        pthread
        rt
        deflate
        hts
        z
        mpicxx
        mpi
        swverbs_ft_perf
        swutils
        swtm
        m
        dl
    )
    add_library(swbam::io ALIAS swbam_io)
    add_library(swbam::mpi_runtime ALIAS swbam_mpi_runtime)
    add_library(swbam::cpe_runtime ALIAS swbam_cpe_runtime)
    add_library(swbam::cpe_kernels ALIAS swbam_cpe_kernels)
    add_library(swbam::algorithm_kernels ALIAS swbam_algorithm_kernels)
    add_library(swbam::swbam ALIAS swbam_sdk)

    if(BUILD_SWBAM_EXAMPLES)
        add_executable(swbam-sdk-record-count
            examples/sdk_record_count.cpp)
        add_executable(swbam-sdk-filter-bam
            examples/sdk_filter_bam.cpp)
        target_compile_options(swbam-sdk-record-count PRIVATE ${HOST_FLAGS})
        target_compile_options(swbam-sdk-filter-bam PRIVATE ${HOST_FLAGS})
        target_link_libraries(swbam-sdk-record-count PRIVATE swbam::swbam)
        target_link_libraries(swbam-sdk-filter-bam PRIVATE swbam::swbam)
    endif()
    
    # Main executable with host code and slave objects
    add_executable(${PROJECT_NAME} ${SINGLE_SRC_FILES}
        $<TARGET_OBJECTS:swbam_cpe_support_objects>
        $<TARGET_OBJECTS:swbam_cpe_algorithm_objects>
        $<TARGET_OBJECTS:swbam_cpe_core_objects>)
    # Apply host flags to all source files in src directory
    target_compile_options(${PROJECT_NAME} PRIVATE ${HOST_FLAGS})
    add_executable(RabbitBAM-CGS ${CGS_SRC_FILES}
        $<TARGET_OBJECTS:slave_cgs_objects>
        $<TARGET_OBJECTS:swbam_cpe_core_objects>)
    target_compile_options(RabbitBAM-CGS PRIVATE ${HOST_FLAGS})
    target_compile_definitions(RabbitBAM-CGS PRIVATE PLATFORM_SUNWAY ENABLE_CGS_BAM2BAM ENABLE_CGS_BAM2SAM ENABLE_CGS_SAM2BAM)
    add_executable(RabbitBAM-MPI ${MPI_SRC_FILES})
    target_compile_options(RabbitBAM-MPI PRIVATE ${HOST_FLAGS})
    target_compile_definitions(RabbitBAM-MPI PRIVATE PLATFORM_SUNWAY ENABLE_MPI_BAM2BAM)
    target_include_directories(RabbitBAM-MPI PRIVATE ${SUNWAY_MPI_ROOT}/include)
    target_link_directories(RabbitBAM-MPI PRIVATE ${SUNWAY_MPI_LIB_DIR} /usr/sw/lib)

    add_library(crc32_bench_slave_objects OBJECT slave/slave_crc32_bench.c)
    target_compile_options(crc32_bench_slave_objects PRIVATE ${SLAVE_FLAGS})
    target_include_directories(crc32_bench_slave_objects PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/libdeflate-1.20
        ${CMAKE_SOURCE_DIR}/slave
        ${CMAKE_SOURCE_DIR}/slave/slave_zlib
    )
    target_compile_definitions(crc32_bench_slave_objects PRIVATE PLATFORM_SUNWAY)
    add_executable(RabbitBAM-CRC32-Bench
        bench/bench_crc32.cpp
        $<TARGET_OBJECTS:crc32_bench_slave_objects>
    )
    target_compile_options(RabbitBAM-CRC32-Bench PRIVATE ${HOST_FLAGS})
    target_include_directories(RabbitBAM-CRC32-Bench PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_link_directories(RabbitBAM-CRC32-Bench PRIVATE /usr/sw/lib)
    target_compile_definitions(RabbitBAM-CRC32-Bench PRIVATE PLATFORM_SUNWAY)
else()
    # For x86, only compile src files
    add_executable(${PROJECT_NAME} ${SINGLE_SRC_FILES})
endif()

add_executable(markdup_verify_fast
    bench/markdup_verify_fast.cpp
)
if(PLATFORM_SUNWAY)
    target_compile_options(markdup_verify_fast PRIVATE ${HOST_FLAGS})
endif()

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
    pthread
    rt
    deflate
    hts
    z
)
if(PLATFORM_SUNWAY)
    target_link_libraries(RabbitBAM-CGS PRIVATE
        pthread
        rt
        deflate
        hts
        z
    )
    target_link_libraries(RabbitBAM-MPI PRIVATE
        swbam::swbam
        swbam::algorithm_kernels)
    target_link_libraries(RabbitBAM-CRC32-Bench PRIVATE
        pthread
        rt
        m
        dl
    )
endif()

# Link SWLU library if enabled
if(USE_SWLU)
    target_link_libraries(${PROJECT_NAME} PRIVATE swlu)
    if(PLATFORM_SUNWAY)
        target_link_libraries(RabbitBAM-CGS PRIVATE swlu)
        target_link_libraries(RabbitBAM-MPI PRIVATE swlu)
    endif()
endif()

# Add static linking option
# target_link_options(${PROJECT_NAME} PRIVATE -static)

# Sunway platform specific link options
if(PLATFORM_SUNWAY)
    target_link_options(${PROJECT_NAME} PRIVATE -mhybrid)
    target_link_options(RabbitBAM-CGS PRIVATE -mhybrid)
    target_link_options(RabbitBAM-MPI PRIVATE -mhybrid -Wl,-zmuldefs)
    target_link_options(RabbitBAM-CRC32-Bench PRIVATE -mhybrid)
endif()

# OpenMP support (only for x86)
if(PLATFORM_X86)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

# Set executable name and output directory
set_target_properties(${PROJECT_NAME} PROPERTIES
    OUTPUT_NAME "RabbitBAM-X"
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

if(PLATFORM_SUNWAY)
    set_target_properties(RabbitBAM-CGS PROPERTIES
        OUTPUT_NAME "RabbitBAM-CGS"
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    set_target_properties(RabbitBAM-MPI PROPERTIES
        OUTPUT_NAME "RabbitBAM-MPI"
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    set_target_properties(RabbitBAM-CRC32-Bench PROPERTIES
        OUTPUT_NAME "RabbitBAM-CRC32-Bench"
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endif()
