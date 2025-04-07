#pragma once
#include "Parafun.h"
#include "ShellData.h"
#include <CGAL-Interface/CGAL-Interface.h>

using namespace Eigen;
using namespace std;

class BiljectivePara {
public:
  BiljectivePara(const cgl::SurfaceMesh3 &m, string filename);
  ~BiljectivePara();

  void parameterization();
  void load();
  void WriteUVMesh(std::ofstream &fout);
  void shelltri(MatrixXi &tri, MatrixXd &pre_pos, MatrixXd &pos, VectorXi &bnd);
  double adjust_weight(double conv_mesh, double last_energy);
  double GetDistortion() { return last_mesh_energy_; }

protected:
  string modelname;
  const cgl::SurfaceMesh3 &mesh;
  cgl::SurfaceMesh3 uv_mesh_;
  ShellData shell_data;
  std::shared_ptr<Parafun> parafun_solver = nullptr;
  double last_mesh_energy_;
  double convgence_con_rate = 1e-5;
  int MAX_ITER_NUM = 5000;
  bool is_initialization;

  bool weight1;
  bool weight2;
};