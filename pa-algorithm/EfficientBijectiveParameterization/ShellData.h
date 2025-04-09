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
  double mesh_measure_; // area or volume
  Eigen::MatrixXd whole_uv_,
      whole_uv_pre_; // whole domain uv: mesh + free vertices
  unsigned int num_shell_faces_, num_vertices_, num_faces_;
  double mesh_energy_;          // mesh energy
  Eigen::MatrixXi shell_faces_; // shell domain tets: shell tets
  Eigen::MatrixXi whole_faces_;
  Eigen::VectorXi frame_ids_;
  Eigen::MatrixXd frame_vertices_;
  Eigen::MatrixXi mesh_faces_;
  Eigen::MatrixXd mesh_vertices_; // input initial mesh F/T
  Eigen::VectorXd mesh_face_area_,
      shell_face_area_; // mesh/shell area or volume
  std::vector<int> bnd_sizes_;
  Eigen::VectorXi internal_bnd_;
  double shell_factor_ = 10;

private:
  // 静态const变量, 该数值不会改变
  static const int kDim = 2; // Dimension for ambient space. Same for mesh/shell
  unsigned int num_mesh_vertices_, num_mesh_faces_;
  //  long num_mesh_vertices_, num_mesh_faces_;
  unsigned int num_shell_vertices_;
  Eigen::MatrixXi w_T;
  Eigen::VectorXd w_M; // area/volume weights for whole
  Eigen::VectorXi external_bnd;
  std::vector<int> component_sizes; // multi-chart support
};
