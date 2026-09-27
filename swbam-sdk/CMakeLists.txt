cmake_minimum_required(VERSION 3.13)
include(cmake/SunwayToolchain.cmake)
project(SWBAM VERSION 0.1.0 LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(BUILD_SWBAM_EXAMPLES "Build the two MPI SDK examples" ON)
set(SUNWAY_MPI_ROOT "/usr/sw/mpi/mpi_20220608_SEA" CACHE PATH "Sunway MPI installation root")
set(SWBAM_DEPENDENCY_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../ext" CACHE PATH
    "Directory containing the built Sunway htslib-1.20 and libdeflate-1.20")
set(SWBAM_HTSLIB_ROOT "${SWBAM_DEPENDENCY_ROOT}/htslib-1.20" CACHE PATH "Built HTSlib source directory")
set(SWBAM_LIBDEFLATE_ROOT "${SWBAM_DEPENDENCY_ROOT}/libdeflate-1.20" CACHE PATH "Built libdeflate source directory")
include(cmake/SWBAMDependencies.cmake)
swbam_define_dependencies(
    "${SWBAM_HTSLIB_ROOT}" "${SWBAM_HTSLIB_ROOT}/libhts.a"
    "${SWBAM_LIBDEFLATE_ROOT}" "${SWBAM_LIBDEFLATE_ROOT}/build/libdeflate.a")

# These flags are part of the MPE/CPE parameter-layout and codegen contract.
add_library(swbam_platform INTERFACE)
set_target_properties(swbam_platform PROPERTIES EXPORT_NAME platform)
target_compile_features(swbam_platform INTERFACE cxx_std_11)
target_compile_options(swbam_platform INTERFACE
    -g -O3 -ffast-math -funroll-loops -fstrict-aliasing -faddress_align=64)
target_compile_definitions(swbam_platform INTERFACE PLATFORM_SUNWAY __STDC_LIMIT_MACROS TIMING)
target_include_directories(swbam_platform INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_library(swbam_host INTERFACE)
set_target_properties(swbam_host PROPERTIES EXPORT_NAME host)
target_link_libraries(swbam_host INTERFACE swbam_platform)
target_compile_options(swbam_host INTERFACE -mhost)
add_library(swbam::host ALIAS swbam_host)

add_library(swbam_cpe INTERFACE)
set_target_properties(swbam_cpe PROPERTIES EXPORT_NAME cpe)
target_link_libraries(swbam_cpe INTERFACE swbam_platform swbam::htslib)
target_compile_options(swbam_cpe INTERFACE -mslave -msimd)
add_library(swbam::cpe ALIAS swbam_cpe)

add_library(swbam_io STATIC lib/io/swbam_io.cpp lib/io/swbam_rank_body_sink.cpp)
set_target_properties(swbam_io PROPERTIES EXPORT_NAME io)
target_link_libraries(swbam_io PUBLIC swbam_host swbam::htslib)
add_library(swbam::io ALIAS swbam_io)

add_library(swbam_mpi_runtime STATIC lib/mpi/swbam_mpi_runtime.cpp)
set_target_properties(swbam_mpi_runtime PROPERTIES EXPORT_NAME mpi_runtime)
target_link_libraries(swbam_mpi_runtime PUBLIC swbam_io swbam::mpi_platform)
add_library(swbam::mpi_runtime ALIAS swbam_mpi_runtime)

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
set_target_properties(swbam_cpe_runtime PROPERTIES EXPORT_NAME cpe_runtime)
target_link_libraries(swbam_cpe_runtime PUBLIC swbam_io)
add_library(swbam::cpe_runtime ALIAS swbam_cpe_runtime)

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
foreach(kind core operator support)
    set(target swbam_cpe_${kind}_objects)
    set_target_properties(${target} PROPERTIES EXPORT_NAME cpe_${kind}_objects)
    target_link_libraries(${target} PRIVATE swbam_cpe)
    target_include_directories(${target} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/slave ${CMAKE_CURRENT_SOURCE_DIR}/slave/slave_zlib)
endforeach()

# Export the actual objects: the Sunway hybrid linker needs them directly.
add_library(swbam_cpe_kernels INTERFACE)
set_target_properties(swbam_cpe_kernels PROPERTIES EXPORT_NAME cpe_kernels)
foreach(kind support core operator)
    target_sources(swbam_cpe_kernels INTERFACE
        $<BUILD_INTERFACE:$<TARGET_OBJECTS:swbam_cpe_${kind}_objects>>
        $<INSTALL_INTERFACE:$<TARGET_OBJECTS:swbam::cpe_${kind}_objects>>)
endforeach()
add_library(swbam::cpe_kernels ALIAS swbam_cpe_kernels)

add_library(swbam_sdk INTERFACE)
set_target_properties(swbam_sdk PROPERTIES EXPORT_NAME swbam)
target_link_options(swbam_sdk INTERFACE -g -O3 -ffast-math -funroll-loops
    -fstrict-aliasing -mhybrid -Wl,-zmuldefs)
target_link_libraries(swbam_sdk INTERFACE
    swbam_cpe_runtime swbam_mpi_runtime swbam_io swbam_cpe_kernels
    pthread rt swbam::deflate swbam::htslib z swbam::mpi_platform)
add_library(swbam::swbam ALIAS swbam_sdk)

install(TARGETS swbam_platform swbam_host swbam_cpe swbam_io swbam_mpi_runtime
        swbam_cpe_runtime swbam_cpe_kernels swbam_sdk
        swbam_cpe_support_objects swbam_cpe_core_objects swbam_cpe_operator_objects
    EXPORT SWBAMTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    OBJECTS DESTINATION ${CMAKE_INSTALL_LIBDIR}/swbam/slave)
install(DIRECTORY include/swbam DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Ship the exact adapted host dependencies, not source-tree search paths.
install(DIRECTORY "${SWBAM_HTSLIB_ROOT}/htslib" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h")
install(FILES "${SWBAM_LIBDEFLATE_ROOT}/libdeflate.h" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES "${SWBAM_HTSLIB_ROOT}/libhts.a"
    "${SWBAM_LIBDEFLATE_ROOT}/build/libdeflate.a"
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES "${SWBAM_HTSLIB_ROOT}/LICENSE" DESTINATION ${CMAKE_INSTALL_DATADIR}/SWBAM/licenses/htslib)
install(FILES "${SWBAM_LIBDEFLATE_ROOT}/COPYING" DESTINATION ${CMAKE_INSTALL_DATADIR}/SWBAM/licenses/libdeflate)

set(package_dir ${CMAKE_INSTALL_LIBDIR}/cmake/SWBAM)
configure_package_config_file(cmake/SWBAMConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/SWBAMConfig.cmake
    INSTALL_DESTINATION ${package_dir}
    PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/SWBAMConfigVersion.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(EXPORT SWBAMTargets NAMESPACE swbam:: DESTINATION ${package_dir})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SWBAMConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/SWBAMConfigVersion.cmake
    cmake/SWBAMDependencies.cmake cmake/SunwayToolchain.cmake
    DESTINATION ${package_dir})
if(BUILD_SWBAM_EXAMPLES)
    add_subdirectory(examples)
endif()
