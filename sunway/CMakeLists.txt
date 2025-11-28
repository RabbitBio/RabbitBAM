# # sunway 平台源文件
# file(GLOB SW_SRC
#     ${CMAKE_SOURCE_DIR}/sunway/*.cpp
# )

# add_executable(rabbitbam_sw
#     ${COMMON_SRC}
#     ${SW_SRC}
# )

# target_include_directories(rabbitbam_sw PRIVATE
#     ${CMAKE_SOURCE_DIR}/sunway
# )



# =============================
# Sunway 目录
set(DIR_SRC ${CMAKE_SOURCE_DIR}/sunway)           # MPE 主线程源文件
set(SLAVE_DIR_SRC ${CMAKE_SOURCE_DIR}/sunway/slave)  # CPE 从核源文件

# MPE 主线程 C++ 文件
file(GLOB MPE_CPP_SRC "${DIR_SRC}/*.cpp")
# CPE 从核 C++ 文件
file(GLOB CPE_CPP_SRC "${SLAVE_DIR_SRC}/*.cpp")

# 最终目标可执行文件
add_executable(rabbitbam_sw
    ${COMMON_SRC}      # 公共代码
    ${MPE_CPP_SRC}     # MPE cpp
    ${CPE_CPP_SRC}     # CPE cpp
)

# =============================
# 编译选项
# MPE 主线程源文件
set_source_files_properties(${MPE_CPP_SRC} PROPERTIES COMPILE_FLAGS "-mhost -g -O3 -w")
# CPE cpp 文件
set_source_files_properties(${CPE_CPP_SRC} PROPERTIES COMPILE_FLAGS "-mslave -msimd -g -O3 -w")

# =============================
# include 路径
target_include_directories(rabbitbam_sw PRIVATE
    ${DIR_SRC}
    ${SLAVE_DIR_SRC}
)

# =============================
# 链接库和选项
target_link_libraries(rabbitbam_sw PRIVATE -static -lz -lpthread -lrt)
target_compile_options(rabbitbam_sw PRIVATE -mhybrid)
