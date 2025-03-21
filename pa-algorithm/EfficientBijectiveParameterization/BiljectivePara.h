#pragma once
#include "Parafun.h"
#include "ShellData.h"
#include <CGALMeshInterface/CGALMeshInterface.h>

class BiljectivePara {
public:
  BiljectivePara(const cgl::SurfaceMesh3 &M);
  ~BiljectivePara();

  void parameterization();
  void load();
  cgl::SurfaceMesh3 getResult();
  void shelltri(Eigen::MatrixXi &tri, Eigen::MatrixXd &pre_pos,
                Eigen::MatrixXd &pos, Eigen::VectorXi &bnd);
  double adjust_weight(double conv_mesh, double last_energy);

protected:
  std::string modelname;
  const cgl::SurfaceMesh3 &mesh;
  cgl::SurfaceMesh3 uv_mesh;
  ShellData shell_data;
  std::shared_ptr<Parafun> parafun_solver = nullptr;

  double convgence_con_rate = 1e-5;
  int MAX_ITER_NUM = 5000;
  bool is_initialization;

  bool weight1;
  bool weight2;
};