project(Preprocess)

if (${LINSYSSOLVER} STREQUAL "Eigen")
    include_directories(SYSTEM ${EIGEN_PATH})
    add_definitions(-DUSE_EIGEN)
elseif(${LINSYSSOLVER} STREQUAL "MKL")
    add_definitions(-DUSE_MKL)
endif()

include_directories(${SOURCE_CODE_PATH}/controller)
include_directories(${SOURCE_CODE_PATH}/interface)
include_directories(${SOURCE_CODE_PATH}/algorithm)
include_directories(SYSTEM ${EIGEN_PATH})
include_directories(SYSTEM ${TRIANGLE_PATH})
include_directories(SYSTEM ${CGAL_PATH}/include)
include_directories(${BOOST_PATH}/include/boost-1_87)

file(GLOB cpp_files
    main.cpp
    process.h
    process.cpp
 )
add_executable(${PROJECT_NAME} ${cpp_files})

target_link_libraries(${PROJECT_NAME} PRIVATE
    frame
    EfficientBijectiveParameterization
    Triangle-Interface
    CGAL-Interface
)