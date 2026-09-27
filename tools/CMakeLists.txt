# Workspace-only developer probes. CRC probes inspect the bundled codec internals.
add_executable(markdup_verify_fast markdup_verify_fast.cpp)
target_link_libraries(markdup_verify_fast PRIVATE swbam::host)
add_executable(swbam-mpi-memory-probe mpi_memory_capacity_probe.cpp)
add_executable(swbam-mpi-memcpy-probe mpi_memcpy_probe.cpp)
foreach(target swbam-mpi-memory-probe swbam-mpi-memcpy-probe)
    target_link_libraries(${target} PRIVATE swbam::swbam)
endforeach()

add_library(cgs_memory_probe_slave_object OBJECT slave_cgs_memory_probe.cpp)
target_link_libraries(cgs_memory_probe_slave_object PRIVATE swbam::cpe)
target_include_directories(cgs_memory_probe_slave_object PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(swbam-cgs-memory-probe cgs_memory_probe.cpp
    $<TARGET_OBJECTS:cgs_memory_probe_slave_object>)
target_link_libraries(swbam-cgs-memory-probe PRIVATE swbam::host pthread rt m dl)
target_link_options(swbam-cgs-memory-probe PRIVATE -mhybrid)

add_library(crc32_bench_slave_objects OBJECT slave_crc32_bench.c)
target_link_libraries(crc32_bench_slave_objects PRIVATE swbam::cpe)
target_include_directories(crc32_bench_slave_objects PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../swbam-sdk/slave
    ${CMAKE_CURRENT_SOURCE_DIR}/../swbam-sdk/slave/slave_zlib)
add_executable(RabbitBAM-CRC32-Bench bench_crc32.cpp $<TARGET_OBJECTS:crc32_bench_slave_objects>)
target_link_libraries(RabbitBAM-CRC32-Bench PRIVATE swbam::host pthread rt m dl)
target_link_options(RabbitBAM-CRC32-Bench PRIVATE -mhybrid)
