cmake_minimum_required(VERSION 3.10)
project(RabbitBAM LANGUAGES C CXX)

# 选平台：默认 x86
set(PLATFORM "x86" CACHE STRING "Choose platform: x86 or sunway")

message(STATUS "Build platform = ${PLATFORM}")

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 公共 include 路径
include_directories(
    ${CMAKE_SOURCE_DIR}/common
    ${CMAKE_SOURCE_DIR}/common/${PLATFORM}
)

# 收集 common 下所有 cpp 文件
file(GLOB COMMON_SRC
    ${CMAKE_SOURCE_DIR}/common/*.cpp
    ${CMAKE_SOURCE_DIR}/common/${PLATFORM}/*.cpp
)

# 平台专属代码
if(PLATFORM STREQUAL "x86")
    message(STATUS "Compiling x86 version...")
    add_subdirectory(x86)
elseif(PLATFORM STREQUAL "sunway")
    message(STATUS "Compiling Sunway version...")
    add_subdirectory(sunway)
else()
    message(FATAL_ERROR "Unknown platform: ${PLATFORM}")
endif()
