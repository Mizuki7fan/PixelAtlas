# 设置项目用到的各种库

# 读取设置系统环境变量
set (MyLib $ENV{MyLibPath})
if(NOT MyLib)
    message(FATAL_ERROR "没有在环境变量中设置本地库目录MyLibPath的值.")
else()
    message("本地库目录为: ${MyLib}")
endif()

set (BOOST_PATH ${MyLib}/boost)
set (CGAL_PATH ${MyLib}/CGAL-6.0.1)
set (EIGEN_PATH ${MyLib}/eigen-3.4.0)
set (TRIANGLE_PATH ${MyLib}/triangle)
set (OpenMesh_PATH ${MyLib}/OpenMesh-11.0)

#${CMAKE_SOURCE_DIR}=项目根目录=G:/PixelAtlas
