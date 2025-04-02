set(SUBPROJ_NAME evaluator)
include_directories(${CMAKE_SOURCE_DIR}/pa-controller)
include_directories(${BOOST_PATH}/include/boost-1_87)

file (GLOB EVALUATOR_FILES
    evaluator.cpp
    evaluator.h
)

add_executable(${SUBPROJ_NAME} ${EVALUATOR_FILES})

target_link_libraries(${SUBPROJ_NAME} PRIVATE
    pa_frame
)

# 设置要包含的boost库
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set (VC_SUFFIX "vc143")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(LIB_BUILD_TYPE mt-gd)
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(LIB_BUILD_TYPE mt)
    endif()
    set (LIB_SUFFIX "x64-1_87")
    foreach(LIB_NAME json program_options)
        message(STATUS ${BOOST_PATH}/lib/libboost_${LIB_NAME}-${VC_SUFFIX}-${LIB_BUILD_TYPE}-${LIB_SUFFIX}.lib)
        target_link_libraries(${SUBPROJ_NAME} PRIVATE 
        ${BOOST_PATH}/lib/libboost_${LIB_NAME}-${VC_SUFFIX}-${LIB_BUILD_TYPE}-${LIB_SUFFIX}.lib)
    endforeach()
elseif()
    message(FATAL_ERROR 尚未编译clang版本的Boost库)
endif()