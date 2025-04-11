#pragma once
#include "ParaFun.h"
#include "ShellData.h"
#include <CGAL-Interface/CGAL-Interface.h>

using namespace Eigen;
using namespace std;

class BiljectivePara {
public:
  BiljectivePara(const cgl::SurfaceMesh3 &mesh);
  ~BiljectivePara();

  void Parameterization();
  void Load();
  void WriteUVMesh(std::ofstream &fout);
  double GetDistortion() { return last_mesh_energy_; }

private:
  void ShellTri(MatrixXi &tri, MatrixXd &pre_pos, MatrixXd &pos, VectorXi &bnd);
  double AdjustWeight(double conv_mesh, double last_energy);

private:
  const cgl::SurfaceMesh3 &mesh_;
  cgl::SurfaceMesh3 uv_mesh_;
  ShellData shell_data_;
  std::shared_ptr<ParaFun> solver_ = nullptr;
  double last_mesh_energy_;
  double convgence_con_rate_ = 1e-5;
  static constexpr int kNumMaxIteration = 5000;
  bool is_initialization_;

  bool weight1_, weight2_;
};