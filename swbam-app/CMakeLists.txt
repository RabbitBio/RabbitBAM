cmake_minimum_required(VERSION 3.13)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR AND NOT CMAKE_TOOLCHAIN_FILE)
    message(FATAL_ERROR "Use -DCMAKE_TOOLCHAIN_FILE=<sdk-prefix>/lib/cmake/SWBAM/SunwayToolchain.cmake")
endif()
project(SWBAMApp LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT TARGET swbam::swbam)
    find_package(SWBAM 0.1 CONFIG REQUIRED)
endif()

add_library(swbam_cpe_algorithm_objects OBJECT
    slave/sort_mpi.cpp slave/collate_mpi.cpp slave/markdup_mpi.cpp slave/fixmate_mpi.cpp)
target_link_libraries(swbam_cpe_algorithm_objects PRIVATE swbam::cpe)
target_include_directories(swbam_cpe_algorithm_objects PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(RabbitBAM-MPI
    main.cpp CmdInfo.cpp
    commands/swbam_bam2bam_mpi.cpp
    commands/swbam_bam2sam_mpi.cpp
    commands/swbam_fixmate_mpi.cpp
    commands/swbam_flagstat_mpi.cpp
    commands/swbam_io_check_mpi.cpp
    commands/swbam_mpi_common.cpp
    commands/swbam_sam2bam_mpi.cpp
    commands/swbam_stats_mpi.cpp
    algorithms/sort/swbam_sort_mpi.cpp
    algorithms/sort/swbam_sort_optimized_mpi.cpp
    algorithms/collate/swbam_collate_mpi.cpp
    algorithms/markdup/swbam_markdup_mpi.cpp
    algorithms/markdup/swbam_markdup_stream_mpi.cpp
    pipelines/swbam_dedup_pipeline_mpi.cpp
    pipelines/swbam_dedup_workflow_mpi.cpp
    $<TARGET_OBJECTS:swbam_cpe_algorithm_objects>)
target_compile_definitions(RabbitBAM-MPI PRIVATE ENABLE_MPI_BAM2BAM)
target_include_directories(RabbitBAM-MPI PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/third_party)
target_link_libraries(RabbitBAM-MPI PRIVATE swbam::swbam)

option(USE_SWLU "Enable SWLU profiling library" OFF)
if(USE_SWLU)
    find_path(SWLU_INCLUDE_DIR swlu.h HINTS ${SWLU_ROOT}/include)
    find_library(SWLU_LIBRARY swlu HINTS ${SWLU_ROOT})
    if(NOT SWLU_INCLUDE_DIR OR NOT SWLU_LIBRARY)
        message(FATAL_ERROR "Set SWLU_ROOT to the SWLU profiling installation")
    endif()
    target_compile_definitions(RabbitBAM-MPI PRIVATE USE_SWLU)
    target_include_directories(RabbitBAM-MPI PRIVATE ${SWLU_INCLUDE_DIR})
    target_link_libraries(RabbitBAM-MPI PRIVATE ${SWLU_LIBRARY})
endif()
install(TARGETS RabbitBAM-MPI RUNTIME DESTINATION bin)
