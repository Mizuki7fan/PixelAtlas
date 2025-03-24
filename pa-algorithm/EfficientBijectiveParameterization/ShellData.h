#pragma once

#include <Eigen/Dense>

class ShellData {
public:
  explicit ShellData();
  void AddNewPatch(const Eigen::MatrixXd &V_in, const Eigen::MatrixXi &F_ref,
                   const Eigen::RowVectorXd &center);
  Eigen::MatrixXd w_uv_;     // whole domain uv: mesh + free vertices
  Eigen::MatrixXd w_uv_pre_; // whole domain uv: mesh + free vertices
  double mesh_measure_;      // area or volume
  long v_num_, f_num_;
  long sf_num_;
  double energy_;       // mesh energy
  Eigen::MatrixXi m_T_; // input initial mesh F/T
  Eigen::MatrixXi s_T_; // shell domain tets: shell tets
  Eigen::MatrixXd m_V_; // input initial mesh V
  double shell_factor_ = 10;
  void UpdateShell();
  Eigen::VectorXd m_M_; // mesh area or volume
  Eigen::VectorXd s_M_; // shell area or volume
  Eigen::MatrixXi surface_F_;
  Eigen::VectorXi frame_ids_;
  Eigen::MatrixXd frame_V_;
  std::vector<int> bnd_sizes_;
  Eigen::VectorXi internal_bnd_;

private:
  void MeshImprove();

  int mv_num_, mf_num_, sv_num_;

  //   Eigen::VectorXd w_M_; // area/volume weights for whole

  //   Eigen::VectorXi external_bnd_;

  std::vector<int> component_sizes_; // multi-chart support

  //   const int kDim_ = 2; // dimension for ambient space. Same for mesh/shell
};
