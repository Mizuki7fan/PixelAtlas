project(Preprocess)

include_directories(${CMAKE_SOURCE_DIR}/pa-controller)
include_directories(${CMAKE_SOURCE_DIR}/pa-steps/common)
include_directories(SYSTEM ${CGAL_PATH}/include)
include_directories(${BOOST_PATH}/include/boost-1_87)
message(STATUS ${CMAKE_SOURCE_DIR}/pa-controller)

file(GLOB cpp_files
    main.cpp
    process.h
    process.cpp
 )
add_executable(${PROJECT_NAME} ${cpp_files})

target_link_libraries(${PROJECT_NAME} PRIVATE
    pa_frame
    CGALMeshInterface
)