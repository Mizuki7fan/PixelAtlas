#include "process.h"
#include "CGALMeshInterface/CGALMeshInterface.h"
#include <frame/global_defs.h>
#include <iostream>
// 函数执行入口
void MainProcess(std::string model_name) {
  cgl::SurfaceMesh3 mesh;
  cgl::MeshOperation::IO::loadOBJ(model_name.c_str(), mesh);
}