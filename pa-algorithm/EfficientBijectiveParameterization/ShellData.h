#pragma once
#include <Eigen/Dense>

class ShellData {
public:
  ShellData();
  void AddNewPatch(const Eigen::MatrixXd &, const Eigen::MatrixXi &,
                   const Eigen::RowVectorXd &center);
  void UpdateShell();

private:
  void MeshImprove();

public:
  double mesh_measure_;             // area or volume
  Eigen::MatrixXd w_uv_, w_uv_pre_; // whole domain uv: mesh + free vertices
  unsigned int num_shell_faces_, num_vertices_, num_faces_;
  double mesh_energy_;          // mesh energy
  Eigen::MatrixXi shell_faces_; // shell domain tets: shell tets
  Eigen::MatrixXi whole_triangles_;
  Eigen::VectorXi frame_ids;
  Eigen::MatrixXd frame_V;
  Eigen::MatrixXi m_T;      // input initial mesh F/T
  Eigen::VectorXd m_M, s_M; // mesh/shell area or volume
  std::vector<int> bnd_sizes;
  Eigen::VectorXi internal_bnd;
  Eigen::MatrixXd m_V; // input initial mesh V
  double shell_factor = 10;

private:
  int dim = 2; // dimension for ambient space. Same for mesh/shell
  long mv_num, mf_num;
  long sv_num;
  Eigen::MatrixXi w_T;
  Eigen::VectorXd w_M; // area/volume weights for whole
  Eigen::VectorXi external_bnd;
  std::vector<int> component_sizes; // multi-chart support
};
