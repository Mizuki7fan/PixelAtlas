#include "ShellData.h"
#include "Utility.h"
#include <AlgoKit/BoundaryLoop.h>
#include <Triangle-Interface/TriangleInterface.h>
#include <iostream>
#include <numbers>

ShellData::ShellData()
    : mesh_measure_(0),      //
      num_shell_faces_(0),   //
      num_mesh_vertices_(0), //
      num_mesh_faces_(0),    //
      num_shell_vertices_(0) //
{
  whole_uv_.resize(0, 0);
};

void ShellData::UpdateShell() {
  num_mesh_vertices_ = static_cast<unsigned int>(mesh_vertices_.rows());
  num_mesh_faces_ = static_cast<unsigned int>(mesh_faces_.rows());

  num_vertices_ = static_cast<unsigned int>(whole_uv_.rows());
  num_shell_faces_ = static_cast<unsigned int>(shell_faces_.rows());

  num_shell_vertices_ = num_vertices_ - num_mesh_vertices_;
  num_faces_ = num_shell_faces_ + num_mesh_faces_;

  shell_face_area_ = Eigen::VectorXd::Constant(num_shell_faces_, shell_factor_);
}

void ShellData::MeshImprove() {
  Eigen::MatrixXd m_uv = whole_uv_.topRows(num_mesh_vertices_);
  Eigen::MatrixXd V_bnd(internal_bnd_.size(), 2);
  for (int i = 0; i < internal_bnd_.size(); i++)
    V_bnd.row(i) = m_uv.row(internal_bnd_(i));

  if (frame_vertices_.size() == 0) {
    Eigen::VectorXd uv_max = m_uv.colwise().maxCoeff();
    Eigen::VectorXd uv_min = m_uv.colwise().minCoeff();
    Eigen::VectorXd uv_mid = (uv_max + uv_min) / 2.;
    // double raduis1 = (uv_max(0) - uv_min(0)) / 2.;
    double raduis1 = V_bnd.row(0).norm();
    // std::cout << (uv_max(0) - uv_min(0)) / 2 << " " << raduis1 << " "
    //           << uv_mid(0) << " " << uv_mid(1) << std::endl;
    double raduis2 = 1.05 * raduis1; // 1.1 rabbit

    int frame_points = static_cast<int>(internal_bnd_.size() / 5);
    if (frame_points > 200)
      frame_points = 200;
    frame_points = std::max(13, frame_points);
    // std::cout << "shell boundary: " << frame_points << std::endl;

    double delta_angle = 2 * std::numbers::pi / frame_points;
    frame_vertices_.resize(frame_points, 2);
    for (int i = 0; i < frame_points; ++i) {
      frame_vertices_.row(i) << raduis2 * cos(i * delta_angle),
          raduis2 * sin(-i * delta_angle);
    }

    frame_ids_ = Eigen::VectorXi::LinSpaced(
        frame_vertices_.rows(), num_mesh_vertices_,
        num_mesh_vertices_ + static_cast<unsigned int>(frame_vertices_.rows()) -
            1);
  } else {
    for (int i = 0; i < frame_vertices_.rows(); ++i)
      frame_vertices_.row(i) = whole_uv_.row(frame_ids_(i));
  }

  Eigen::MatrixXd V;
  Eigen::MatrixXi E;

  V.resize(V_bnd.rows() + frame_vertices_.rows(), V_bnd.cols());
  V << V_bnd, frame_vertices_;

  E.resize(V.rows(), 2);
  for (int i = 0; i < E.rows(); i++)
    E.row(i) << i, i + 1;
  int acc_bs = 0;
  for (auto bs : bnd_sizes_) {
    E(acc_bs + bs - 1, 1) = acc_bs;
    acc_bs += bs;
  }
  E(V.rows() - 1, 1) = acc_bs;
  if (acc_bs != internal_bnd_.size())
    std::runtime_error("acc_bs != internal_bnd_.size()");

  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(component_sizes_.size(), 2);
  {
    int hole_f = 0;
    int hole_i = 0;
    for (auto cs : component_sizes_) {
      for (int i = 0; i < 3; i++)
        H.row(hole_i) += m_uv.row(mesh_faces_(hole_f, i)); // redoing step 2
      hole_f += cs;
      hole_i++;
    }
  }
  H /= 3.;

  Eigen::MatrixXd uv2;
  GenerateTriangulate(V, E, H, uv2, shell_faces_);
  int bnd_n = static_cast<int>(internal_bnd_.size());

  for (auto i = 0; i < shell_faces_.rows(); i++) {
    for (auto j = 0; j < shell_faces_.cols(); j++) {
      int &x = shell_faces_(i, j);
      if (x < bnd_n)
        x = internal_bnd_(x);
      else
        x += static_cast<int>(m_uv.rows()) - bnd_n;
    }
  }

  whole_faces_.resize(mesh_faces_.rows() + shell_faces_.rows(), 3);
  whole_faces_ << mesh_faces_, shell_faces_;

  whole_uv_.conservativeResize(m_uv.rows() - bnd_n + uv2.rows(), 2);
  whole_uv_.bottomRows(uv2.rows() - bnd_n) =
      uv2.bottomRows(-bnd_n + uv2.rows());

  UpdateShell();
}

void ShellData::AddNewPatch(const Eigen::MatrixXd &V_in,  //
                            const Eigen::MatrixXi &F_ref, //
                            const Eigen::RowVectorXd &center) {
  Eigen::VectorXd M;
  M.resize(F_ref.rows());
  { // 计算每个三角形的面积
    Eigen::Vector3d v0, v1, v2, e01, e02, normal_f;
    double area_f;
    for (int i = 0; i < F_ref.rows(); i++) {
      v0 = V_in.row(F_ref(i, 0));
      v1 = V_in.row(F_ref(i, 1));
      v2 = V_in.row(F_ref(i, 2));
      e01 = v1 - v0;
      e02 = v2 - v0;
      normal_f = e01.cross(e02);
      area_f = normal_f.norm() / 2.0;
      M(i) = area_f;
    }
  }

  mesh_face_area_.conservativeResize(num_mesh_faces_ + M.size());
  mesh_face_area_.bottomRows(M.size()) = M;
  mesh_measure_ += M.sum();

  const Eigen::MatrixXd &V_ref = V_in;
  Eigen::MatrixXd uv_init;
  Eigen::VectorXi bnd;
  Eigen::MatrixXd bnd_uv;
  std::vector<std::vector<int>> all_bnds;
  AlgoKit::GetBoundaryLoop(F_ref, all_bnds);
  int num_holes = static_cast<int>(all_bnds.size()) - 1;

  std::sort(all_bnds.begin(), all_bnds.end(),
            [](auto &a, auto &b) { return a.size() > b.size(); });

  bnd = Eigen::Map<Eigen::VectorXi>(all_bnds[0].data(), all_bnds[0].size());

  MapVerticesToCircle(V_ref, bnd, bnd_uv);
  bnd_uv *= sqrt(M.sum() / std::numbers::pi);
  bnd_uv.rowwise() += center;

  std::cout << "Mesh Measure " << M.sum() << "; number holes " << num_holes
            << std::endl;
  std::cout << "hole: " << num_holes << std::endl;
  if (num_holes == 0) {
    if (bnd.rows() == V_ref.rows()) {
      std::cout << "All vert on boundary" << std::endl;
      uv_init.resize(V_ref.rows(), 2);
      for (int i = 0; i < bnd.rows(); i++)
        uv_init.row(bnd(i)) = bnd_uv.row(i);
    } else {
      Tutte(static_cast<int>(V_ref.rows()), F_ref, bnd, bnd_uv, uv_init);
    }
  } else {
    const Eigen::MatrixXi &F = F_ref;
    const Eigen::MatrixXd &V = V_in;
    Eigen::VectorXi &primary_bnd = bnd;
    // fill holes
    int n_filled_faces = 0;
    int real_F_num = static_cast<int>(F.rows());
    for (int i = 0; i < num_holes; i++)
      n_filled_faces += static_cast<int>(all_bnds[i + 1].size());

    Eigen::MatrixXi F_filled(n_filled_faces + real_F_num, 3);
    F_filled.topRows(real_F_num) = F;

    int new_vert_id = static_cast<int>(V.rows());
    int new_face_id = real_F_num;

    for (int i = 0; i < num_holes; i++) {
      auto it = all_bnds[i + 1].begin();
      auto back = all_bnds[i + 1].end() - 1;
      F_filled.row(new_face_id++) << *it, *back, new_vert_id;
      while (it != back) {
        F_filled.row(new_face_id++) << *(it + 1), *(it), new_vert_id;
        it++;
      }
      new_vert_id++;
    }
    if (new_face_id != F_filled.rows())
      std::runtime_error("new_face_id != F_filled.rows()");
    if (new_vert_id != V.rows() + num_holes)
      std::runtime_error("new_vert_id!=V.rows()+num_holes");

    Tutte(static_cast<int>(V_ref.rows()) + num_holes, F_filled, primary_bnd,
          bnd_uv, uv_init);
    uv_init.conservativeResize(V.rows(), 2);
  }

  // writeObj(uv_init, F_ref, "tutte.obj");

  component_sizes_.push_back(static_cast<int>(F_ref.rows()));

  if (num_mesh_vertices_ == 0) {
    whole_uv_ = uv_init;
  } else {
    Eigen::MatrixXd m_uv = whole_uv_.topRows(num_mesh_vertices_);
    whole_uv_.resize(m_uv.rows() + uv_init.rows(), 2);
    whole_uv_ << m_uv, uv_init;
  }

  for (auto cur_bnd : all_bnds) {
    internal_bnd_.conservativeResize(internal_bnd_.size() + cur_bnd.size());
    internal_bnd_.bottomRows(cur_bnd.size()) =
        Eigen::Map<Eigen::ArrayXi>(cur_bnd.data(), cur_bnd.size()) +
        num_mesh_vertices_;
    bnd_sizes_.push_back(static_cast<int>(cur_bnd.size()));
  }

  mesh_faces_.conservativeResize(num_mesh_faces_ + F_ref.rows(), 3);
  mesh_faces_.bottomRows(F_ref.rows()) = F_ref.array() + num_mesh_vertices_;
  num_mesh_faces_ += static_cast<int>(F_ref.rows());

  mesh_vertices_.conservativeResize(num_mesh_vertices_ + V_ref.rows(), 3);
  mesh_vertices_.bottomRows(V_ref.rows()) = V_ref;
  num_mesh_vertices_ += static_cast<int>(V_ref.rows());

  frame_vertices_.resize(0, 0);

  MeshImprove();
}
