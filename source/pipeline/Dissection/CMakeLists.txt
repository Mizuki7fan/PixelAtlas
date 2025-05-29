project(Dissection)

include_directories(${SOURCE_CODE_PATH}/controller)
include_directories(${SOURCE_CODE_PATH}/algorithm)


file (GLOB cpp_files
    main.cpp
    process.h
    process.cpp
    grid_dissector.h
    grid_dissector.cpp
)

add_executable(${PROJECT_NAME} ${cpp_files})

target_link_libraries(${PROJECT_NAME} PRIVATE
    frame
)