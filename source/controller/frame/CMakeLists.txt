set(SUBPROJ_NAME frame)

include_directories(${BOOST_PATH}/include/boost-1_87)

file (GLOB FRAME_FILES
    assert.hpp
    action_args.h
    action_args.cpp
    global_args.cpp
    global_args.h
    io.cpp
    io.h
    common_program.cpp
    common_program.h
    process_parallel_executor.cpp
    process_parallel_executor.h
    metric.h
    metric.cpp
)

add_library(${SUBPROJ_NAME} STATIC ${FRAME_FILES})

# 设置要包含的boost库
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set (VC_SUFFIX "vc143")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(LIB_BUILD_TYPE mt-gd)
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(LIB_BUILD_TYPE mt)
    endif()
    set (LIB_SUFFIX "x64-1_87")
    foreach(LIB_NAME json container program_options thread chrono filesystem)
        message(STATUS ${BOOST_PATH}/lib/libboost_${LIB_NAME}-${VC_SUFFIX}-${LIB_BUILD_TYPE}-${LIB_SUFFIX}.lib)
        target_link_libraries(${SUBPROJ_NAME} PRIVATE 
        ${BOOST_PATH}/lib/libboost_${LIB_NAME}-${VC_SUFFIX}-${LIB_BUILD_TYPE}-${LIB_SUFFIX}.lib)
    endforeach()
elseif()
    message(FATAL_ERROR 尚未编译clang版本的Boost库)
endif()