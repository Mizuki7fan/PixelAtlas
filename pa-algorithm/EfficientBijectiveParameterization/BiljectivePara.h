#pragma once
#include "Parafun.h"
#include "ShellData.h"
#include <CGALMeshInterface/CGALMeshInterface.h>

class BiljectivePara {
public:
  BiljectivePara(const cgl::SurfaceMesh3 &M);
  ~BiljectivePara() {};

public:
  void Parameterization();
  void Load();
  cgl::SurfaceMesh3 GetResult();

private:
  void ShellTri(Eigen::MatrixXi &tri, Eigen::MatrixXd &pre_pos,
                Eigen::MatrixXd &pos, Eigen::VectorXi &bnd);
  double AdjustWeight(double conv_mesh, double last_energy);

  const cgl::SurfaceMesh3 &tri_mesh_;
  cgl::SurfaceMesh3 uv_mesh_;
  std::unique_ptr<ShellData> shell_data_ = nullptr;
  std::unique_ptr<Parafun> parafun_solver_ = nullptr;
  const double convgence_con_rate_ = 1e-5;
  const int kMaxIterNum = 5000;
  bool is_initialization_;

  bool weight_1_, weight_2_;
};