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
set(COMMON_FLAGS "-g -O3 -ffast-math -DTIMING")

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

if(PLATFORM_SUNWAY)
    # For Sunway, also compile slave code
    file(GLOB SLAVE_FILES "slave/*.cpp")
    file(GLOB SLAVE_C_FILES   "slave/*.c")
    # Include slave_zlib source files (compile with -mslave)
    # Only files in slave_zlib directory, not subdirectories
    file(GLOB SLAVE_ZLIB_FILES "slave/slave_zlib/*.c")
    
    # Create object library for slave code with special flags
    add_library(slave_objects OBJECT ${SLAVE_FILES}  ${SLAVE_C_FILES} ${SLAVE_ZLIB_FILES})
    target_compile_options(slave_objects PRIVATE ${SLAVE_FLAGS})
    target_include_directories(slave_objects PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/htslib-1.20
        ${CMAKE_SOURCE_DIR}/ext/libdeflate-1.20
        ${CMAKE_SOURCE_DIR}/slave/slave_zlib
    )
    target_compile_definitions(slave_objects PRIVATE PLATFORM_SUNWAY)
    
    # Main executable with host code and slave objects
    add_executable(${PROJECT_NAME} ${SRC_FILES} $<TARGET_OBJECTS:slave_objects>)
    # Apply host flags to all source files in src directory
    target_compile_options(${PROJECT_NAME} PRIVATE ${HOST_FLAGS})
else()
    # For x86, only compile src files
    add_executable(${PROJECT_NAME} ${SRC_FILES})
endif()

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
    pthread
    rt
    deflate
    hts
    z
)

# Link SWLU library if enabled
if(USE_SWLU)
    target_link_libraries(${PROJECT_NAME} PRIVATE swlu)
endif()

# Add static linking option
# target_link_options(${PROJECT_NAME} PRIVATE -static)

# Sunway platform specific link options
if(PLATFORM_SUNWAY)
    target_link_options(${PROJECT_NAME} PRIVATE -mhybrid)
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

