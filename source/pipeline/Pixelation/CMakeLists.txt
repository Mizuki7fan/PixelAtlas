project(Pixelation)

include_directories(${SOURCE_CODE_PATH}/controller)

file(GLOB cpp_files
    main.cpp
    process.h
    process.cpp
    naive_pixelator.h
    naive_pixelator.cpp
 )

 add_executable(${PROJECT_NAME} ${cpp_files})

target_link_libraries(${PROJECT_NAME} PRIVATE
    frame
)