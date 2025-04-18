project(Pixelation)

include_directories(${CMAKE_SOURCE_DIR}/pa-controller)

file(GLOB cpp_files
    main.cpp
    process.h
    process.cpp
 )

 add_executable(${PROJECT_NAME} ${cpp_files})

target_link_libraries(${PROJECT_NAME} PRIVATE
    pa_frame
)