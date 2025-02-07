set(SUBPROJ_NAME pa_frame)

if (NOT DEFINED JSONCPP OR ${JSONCPP} STREQUAL "")
    message(FATAL_ERROR "没有在环境变量中设置jsoncpp库目录的路径.")#检查JSONCPP的路径
endif()
include_directories(${JSONCPP}/include) 

file (GLOB FRAME_FILES
    assert.hpp
    command_structures.h
    debug.cpp
    debug.h
    io.cpp
    io.h
    parse.cpp
    parse.h
)

add_library(${SUBPROJ_NAME} STATIC ${FRAME_FILES})

target_link_libraries(pa_frame PRIVATE
    debug ${JSONCPP}/build/lib/Debug/jsoncpp.lib
    optimized ${JSONCPP}/build/lib/Release/jsoncpp.lib
)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
copy_file(${SUBPROJ_NAME} jsoncpp.dll ${JSONCPP}/build/bin/Debug ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG})
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
copy_file(${SUBPROJ_NAME} jsoncpp.dll ${JSONCPP}/build/bin/Release ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE})
endif()