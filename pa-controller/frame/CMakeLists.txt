set(SUBPROJ_NAME pa_frame)

include_directories(${BOOST})

file (GLOB FRAME_FILES
    assert.hpp
    command_structures.h
    debug.cpp
    debug.h
    io.cpp
    io.h
    parse.cpp
    parse.h
    common_program.cpp
    common_program.h
)

add_library(${SUBPROJ_NAME} STATIC ${FRAME_FILES})

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
target_link_libraries(${SUBPROJ_NAME} PRIVATE 
    debug   ${BOOST}/stage/lib/libboost_json-vc143-mt-gd-x64-1_84.lib
    optimized   ${BOOST}/stage/lib/libboost_json-vc143-mt-x64-1_84.lib  
    debug   ${BOOST}/stage/lib/libboost_container-vc143-mt-gd-x64-1_84.lib
    optimized   ${BOOST}/stage/lib/libboost_container-vc143-mt-x64-1_84.lib     
)
elseif()
    message(FATAL_ERROR 尚未编译clang版本的Boost库)
endif()