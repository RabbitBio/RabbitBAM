cmake_minimum_required(VERSION 3.10)

# Platform selection (must be before project() to set compiler)
if(NOT DEFINED PLATFORM)
    set(PLATFORM "x86" CACHE STRING "Platform: x86 or sunway")
endif()

# SWLU option (only available for Sunway platform)
option(USE_SWLU "Enable SWLU profiling library" OFF)

if(PLATFORM STREQUAL "sunway")
    # Set Sunway compiler before project()
    set(CMAKE_C_COMPILER swgcc CACHE PATH "C compiler" FORCE)
    set(CMAKE_CXX_COMPILER swg++ CACHE PATH "C++ compiler" FORCE)
    message(STATUS "Using Sunway compiler: C=swgcc, CXX=swg++")
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
file(GLOB MPI_ONLY_SRC "src/*mpi*.cpp")
set(SINGLE_SRC_FILES ${SRC_FILES})
list(REMOVE_ITEM SINGLE_SRC_FILES ${CGS_ONLY_SRC})
list(REMOVE_ITEM SINGLE_SRC_FILES ${MPI_ONLY_SRC})
set(CGS_SRC_FILES ${SRC_FILES})
list(REMOVE_ITEM CGS_SRC_FILES ${MAIN_SRC})
list(REMOVE_ITEM CGS_SRC_FILES ${MPI_ONLY_SRC})
set(MPI_SRC_FILES ${SRC_FILES})
list(REMOVE_ITEM MPI_SRC_FILES ${MAIN_SRC})
list(REMOVE_ITEM MPI_SRC_FILES ${CGS_ONLY_SRC})
list(APPEND MPI_SRC_FILES
    "${CMAKE_SOURCE_DIR}/src/swbam_flagstat_mpi.cpp"
    "${CMAKE_SOURCE_DIR}/src/swbam_stats_mpi.cpp"
    "${CMAKE_SOURCE_DIR}/src/swbam_sort_mpi.cpp"
    "${CMAKE_SOURCE_DIR}/src/swbam_sort_fused_mpi.cpp"
    "${CMAKE_SOURCE_DIR}/src/swbam_markdup_mpi.cpp"
    "${CMAKE_SOURCE_DIR}/src/swbam_markdup_stream_mpi.cpp"
    "${CMAKE_SOURCE_DIR}/src/swbam_fixmate_mpi.cpp"
    "${CMAKE_SOURCE_DIR}/src/swbam_dedup_pipeline_mpi.cpp"
)
list(REMOVE_DUPLICATES MPI_SRC_FILES)

if(PLATFORM_SUNWAY)
    set(SUNWAY_MPI_ROOT "/usr/sw/mpi/mpi_20220608_SEA" CACHE PATH "Sunway MPI installation root")
    set(SUNWAY_MPI_LIB_DIR "${SUNWAY_MPI_ROOT}/lib/single_static")

    # For Sunway, also compile slave code
    file(GLOB SLAVE_FILES "slave/*.cpp")
    file(GLOB SLAVE_C_FILES   "slave/*.c")
    file(GLOB CGS_SLAVE_FILES "slave/*cgs*.cpp")
    list(REMOVE_ITEM SLAVE_C_FILES "${CMAKE_SOURCE_DIR}/slave/slave_crc32_bench.c")
    set(SINGLE_SLAVE_FILES ${SLAVE_FILES})
    list(REMOVE_ITEM SINGLE_SLAVE_FILES ${CGS_SLAVE_FILES})
    list(APPEND SINGLE_SLAVE_FILES "${CMAKE_SOURCE_DIR}/slave/slave_sort_mpi.cpp")
    list(APPEND SINGLE_SLAVE_FILES "${CMAKE_SOURCE_DIR}/slave/slave_markdup_mpi.cpp")
    list(APPEND SINGLE_SLAVE_FILES "${CMAKE_SOURCE_DIR}/slave/slave_fixmate_mpi.cpp")
    list(REMOVE_DUPLICATES SINGLE_SLAVE_FILES)
    # Include slave_zlib source files (compile with -mslave)
    # Only files in slave_zlib directory, not subdirectories
    file(GLOB SLAVE_ZLIB_FILES "slave/slave_zlib/*.c")
    
    # Create object library for slave code with special flags
    add_library(slave_objects OBJECT ${SINGLE_SLAVE_FILES} ${SLAVE_C_FILES} ${SLAVE_ZLIB_FILES})
    target_compile_options(slave_objects PRIVATE ${SLAVE_FLAGS})
    target_include_directories(slave_objects PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/htslib-1.20
        ${CMAKE_SOURCE_DIR}/ext/libdeflate-1.20
        ${CMAKE_SOURCE_DIR}/slave/slave_zlib
    )
    target_compile_definitions(slave_objects PRIVATE PLATFORM_SUNWAY)
    add_library(slave_cgs_objects OBJECT ${SLAVE_FILES} ${SLAVE_C_FILES} ${SLAVE_ZLIB_FILES})
    target_compile_options(slave_cgs_objects PRIVATE ${SLAVE_FLAGS})
    target_include_directories(slave_cgs_objects PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/htslib-1.20
        ${CMAKE_SOURCE_DIR}/ext/libdeflate-1.20
        ${CMAKE_SOURCE_DIR}/slave/slave_zlib
    )
    target_compile_definitions(slave_cgs_objects PRIVATE PLATFORM_SUNWAY ENABLE_CGS_BAM2BAM ENABLE_CGS_BAM2SAM ENABLE_CGS_SAM2BAM)
    
    # Main executable with host code and slave objects
    add_executable(${PROJECT_NAME} ${SINGLE_SRC_FILES} $<TARGET_OBJECTS:slave_objects>)
    # Apply host flags to all source files in src directory
    target_compile_options(${PROJECT_NAME} PRIVATE ${HOST_FLAGS})
    add_executable(RabbitBAM-CGS ${CGS_SRC_FILES} $<TARGET_OBJECTS:slave_cgs_objects>)
    target_compile_options(RabbitBAM-CGS PRIVATE ${HOST_FLAGS})
    target_compile_definitions(RabbitBAM-CGS PRIVATE PLATFORM_SUNWAY ENABLE_CGS_BAM2BAM ENABLE_CGS_BAM2SAM ENABLE_CGS_SAM2BAM)
    add_executable(RabbitBAM-MPI ${MPI_SRC_FILES} $<TARGET_OBJECTS:slave_objects>)
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
