cmake_minimum_required(VERSION 3.18)
project(RabbitBAM LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(HTSlib 1.15 REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_path(LIBDEFLATE_INCLUDE_DIR libdeflate.h
  HINTS "${LIBDEFLATE_ROOT}" ENV LIBDEFLATE_ROOT
  PATH_SUFFIXES include REQUIRED)
find_library(LIBDEFLATE_LIBRARY NAMES deflate
  HINTS "${LIBDEFLATE_ROOT}" ENV LIBDEFLATE_ROOT
  PATH_SUFFIXES lib lib64 REQUIRED)

add_library(rabbitbam_io STATIC
  src/bamsort/rb_sort.cpp
  src/bamindex/rb_index.cpp
  src/rabbitbam/BamTools.cpp
  src/rabbitbam/BamRead.cpp
  src/rabbitbam/BamCompress.cpp
  src/rabbitbam/BamCompleteBlock.cpp
  src/rabbitbam/BamReader.cpp)
target_compile_features(rabbitbam_io PUBLIC cxx_std_11)
target_include_directories(rabbitbam_io PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/bamsort"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/rabbitbam"
  "${HTSlib_INCLUDE_DIR}"
  "${LIBDEFLATE_INCLUDE_DIR}")
target_link_libraries(rabbitbam_io PUBLIC
  "${HTSlib_LIBRARY}" "${LIBDEFLATE_LIBRARY}"
  ZLIB::ZLIB Threads::Threads OpenMP::OpenMP_CXX)

add_executable(rabbitbam src/main.cpp)
target_link_libraries(rabbitbam PRIVATE rabbitbam_io)
install(TARGETS rabbitbam RUNTIME DESTINATION bin)

include(CTest)
if(BUILD_TESTING)
  add_executable(test_bam_io tests/test_bam_io.cpp)
  target_link_libraries(test_bam_io PRIVATE rabbitbam_io)
  add_test(NAME bam_io_roundtrip
    COMMAND test_bam_io "${CMAKE_CURRENT_BINARY_DIR}/bam_io_roundtrip.bam")
  set_tests_properties(bam_io_roundtrip PROPERTIES TIMEOUT 30)
endif()
