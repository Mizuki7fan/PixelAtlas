project(Pixelation)

include_directories(${SOURCE_CODE_PATH}/controller)
include_directories(${SOURCE_CODE_PATH}/interface)
include_directories(SYSTEM ${CGAL_PATH}/include)
include_directories(${BOOST_PATH}/include/boost-1_87)

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
    CGAL-Interface
)