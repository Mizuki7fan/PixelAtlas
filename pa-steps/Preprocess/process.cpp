#include "process.h"
#include "EfficientBijectiveParameterization/BiljectivePara.h"
#include <frame/global_defs.h>
#include <iostream>

// 函数执行入口
// void GenereateUVMeshByEfficientBiljectiveParameterization(
//     const cgl::SurfaceMesh3 &tri_mesh, cgl::SurfaceMesh3 &uv_mesh) {
//   std::cout << tri_mesh.num_edges() << uv_mesh.num_edges() << std::endl;
//   std::unique_ptr<BiljectivePara> bi_para =
//       std::make_unique<BiljectivePara>(tri_mesh);

//   bi_para->Load();
//   bi_para->Parameterization();
//   uv_mesh = bi_para->GetResult();
// }

void MainProcess(std::string model_name) {

  Mesh mesh;
  OpenMesh::IO::read_mesh(mesh, model_name);
  BiljectivePara *bil_para = new BiljectivePara(mesh, model_name);

  bil_para->load(); // 读入网格、完成shell的构建
  bil_para->parameterization();
  bil_para->write_obj("uv_mesh.obj");

  delete bil_para;
  bil_para = NULL;

  // system("pause");

  // cgl::SurfaceMesh3 tri_mesh, uv_mesh;
  // cgl::MeshOperation::IO::loadOBJ(model_name.c_str(), tri_mesh);
  // // 网格应该是带天然切缝的, 如果不带则应该生成切缝,
  // // 可以在前序步骤中加入greedycut

  // // 根据tri_mesh生成uv_mesh
  // GenereateUVMeshByEfficientBiljectiveParameterization(tri_mesh, uv_mesh);
  // if (frm::global::DebugLevel() > 1) {
  //   CGAL::IO::write_OBJ("uv_mesh.obj", uv_mesh);
  // }
}