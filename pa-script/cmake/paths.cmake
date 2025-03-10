# 设置项目用到的各种库

# 读取设置系统环境变量
set (MyLib $ENV{MyLibPath})
if(NOT MyLib)
    message(FATAL_ERROR "没有在环境变量中设置本地库目录MyLibPath的值.")
else()
    message("本地库目录为: ${MyLib}")
endif()

set (BOOST ${MyLib}/boost)

#${CMAKE_SOURCE_DIR}=项目根目录=G:/PixelAtlas
