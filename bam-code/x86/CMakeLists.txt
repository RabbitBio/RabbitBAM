# x86 平台的源文件
file(GLOB X86_SRC
    ${CMAKE_SOURCE_DIR}/x86/*.cpp
)

# 最终可执行文件（或静态库）
add_executable(rabbitbam_x86
    ${COMMON_SRC}
    ${X86_SRC}
)

target_include_directories(rabbitbam_x86 PRIVATE
    ${CMAKE_SOURCE_DIR}/x86
)



# # ----- HTSLIB / libdeflate 的路径 -----
# set(HTSLIB_INSTALL_PATH "/home/user_home/ylf/someGit/bam_deps/htslib-1.20-install")
# set(LIBDEFLATE_INSTALL_PATH "/home/user_home/ylf/someGit/bam_deps/libdeflate-1.20-install")

# include_directories(
#     ${HTSLIB_INSTALL_PATH}/include
#     ${LIBDEFLATE_INSTALL_PATH}/include
# )

# link_directories(
#     ${HTSLIB_INSTALL_PATH}/lib
#     ${HTSLIB_INSTALL_PATH}/lib64
#     ${LIBDEFLATE_INSTALL_PATH}/lib
#     ${LIBDEFLATE_INSTALL_PATH}/lib64
# )

# set(X86_LIBS hts deflate z pthread)

# # ----- tools 模块 -----
# set(TOOLS_SRC
#     BamTools.cpp BamTools.h
#     Buffer.cpp Buffer.h
#     config.h
# )

# add_library(rabbitbamtools SHARED ${TOOLS_SRC})
# target_link_libraries(rabbitbamtools PUBLIC ${X86_LIBS})

# # ----- read 模块 -----
# set(READ_SRC
#     BamCompleteBlock.cpp BamCompleteBlock.h
#     BamCompress.cpp BamCompress.h
#     BamBlock.cpp BamBlock.h
#     BamRead.cpp BamRead.h
#     BamReader.cpp BamReader.h
# )

# add_library(rabbitbamread SHARED ${READ_SRC})
# target_link_libraries(rabbitbamread PUBLIC rabbitbamtools ${X86_LIBS})

# # ----- write 模块 -----
# set(WRITE_SRC
#     BamWriter.cpp BamWriter.h
#     BamWriteCompress.cpp BamWriteCompress.h
# )

# add_library(rabbitbamwrite SHARED ${WRITE_SRC})
# target_link_libraries(rabbitbamwrite PUBLIC rabbitbamtools ${X86_LIBS})

# # ----- 最终可执行文件 (rabbitbam_x86) -----
# # 这里依赖顶层已经生成的 `COMMON_SRC`
# add_executable(rabbitbam_x86
#     block_mul.cpp
#     ${COMMON_SRC}
# )

# target_link_libraries(rabbitbam_x86
#     rabbitbamtools
#     rabbitbamread
#     rabbitbamwrite
#     ${X86_LIBS}
# )
