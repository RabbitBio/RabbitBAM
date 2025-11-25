# sunway 平台源文件
file(GLOB SW_SRC
    ${CMAKE_SOURCE_DIR}/sunway/*.cpp
)

add_executable(rabbitbam_sw
    ${COMMON_SRC}
    ${SW_SRC}
)

target_include_directories(rabbitbam_sw PRIVATE
    ${CMAKE_SOURCE_DIR}/sunway
)
