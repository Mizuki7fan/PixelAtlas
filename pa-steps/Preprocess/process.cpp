#include "process.h"
#include "EfficientBijectiveParameterization/BiljectivePara.h"
#include "frame/io.h"
#include <filesystem>
#include <frame/global_defs.h>

namespace fs = std::filesystem;

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

void MainProcess() {
  fs::path model_path = frm::global::CurrFile();
  Mesh mesh;
  OpenMesh::IO::read_mesh(mesh, model_path.string());
  std::unique_ptr<BiljectivePara> bil_para =
      std::make_unique<BiljectivePara>(mesh, model_path.string());

  bil_para->load(); // 读入网格、完成shell的构建
  bil_para->parameterization();

  std::ofstream fout = frm::CreateResultFilestream("uv_mesh.obj");
  bil_para->write_obj(fout);
  fout.close();
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