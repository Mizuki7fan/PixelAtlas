#include "ParaFun.h"
#if defined(USE_MKL)
#include <LinSysSolver-Interface/MKLPardisoSolver.h>
#elif defined(USE_EIGEN)
#include <LinSysSolver-Interface/EigenLinSolver.h>
#endif

ParaFun::ParaFun(ShellData &data) : shell_data_(data) {
  solver_ = nullptr;
  is_first_ = true;
  bound_distortion_K_ = 250;
  barrer_coef_ = shell_data_.mesh_measure_ * 1e-8;
  std::cout << barrer_coef_ << std::endl;
}

void ParaFun::AfterMeshImprove() {
  total_num_ = shell_data_.num_vertices_;
  F_N_ = shell_data_.num_faces_;
  V_N_ = shell_data_.num_vertices_;
  F_[0].resize(F_N_);
  F_[1].resize(F_N_);
  F_[2].resize(F_N_);

  position_of_mesh_.resize(2 * total_num_);

  for (size_t i = 0; i < kDim; i++) {
    position_of_mesh_.block(i * total_num_, 0, total_num_, 1) =
        shell_data_.whole_uv_.col(i);
  }

  Init();
}
void ParaFun::Init() {
  for (int i = 0; i < F_N_; ++i) {
    F_[0][i] = shell_data_.whole_faces_(i, 0);
    F_[1][i] = shell_data_.whole_faces_(i, 1);
    F_[2][i] = shell_data_.whole_faces_(i, 2);
  }

  HandleMinTri();
  InitArea();
  SetVirtualTri();
  PreCalculate();

  const double *pos = position_of_mesh_.data();

  x_min = pos[0];
  x_max = pos[0];
  y_min = pos[V_N_];
  y_max = pos[V_N_];
  for (int i = 1; i < V_N_; ++i) {
    if (pos[i] < x_min) {
      x_min = pos[i];
    } else if (pos[i] > x_max) {
      x_max = pos[i];
    }

    if (pos[i + V_N_] < y_min) {
      y_min = pos[i + V_N_];
    } else if (pos[i + V_N_] > y_max) {
      y_max = pos[i + V_N_];
    }
  }

  int frame_num = shell_data_.frame_ids_.size();
  average_length = 0;
  for (int j = 0; j < frame_num - 1; ++j) {
    average_length += (shell_data_.frame_vertices_.row(j) -
                       shell_data_.frame_vertices_.row(j + 1))
                          .norm();
  }
  average_length += (shell_data_.frame_vertices_.row(frame_num - 1) -
                     shell_data_.frame_vertices_.row(0))
                        .norm();
  average_length = average_length / frame_num;
  threhold = average_length / 8.0;
  threhold = 2 * threhold; // dis * 2

  cellx_num = std::ceil((x_max - x_min) / average_length);
  celly_num = std::ceil((y_max - y_min) / average_length);
  cell_points.clear();
  cell_points.resize(cellx_num * celly_num);

  if (solver_ != nullptr)
    solver_ = nullptr;

#if defined(USE_MKL)
  solver_ = std::make_unique<MKLPardisoSolver>();
#elif defined(USE_EIGEN)
  solver_ = std::make_unique<EigenLinSolver>();
#endif
  solver_->ia_ = pardiso_ia_;
  solver_->ja_ = pardiso_ja_;
  solver_->a_.resize(pardiso_ja_.size());
  solver_->nnz_ = pardiso_ja_.size();
  solver_->num_ = 2 * V_N_;
  long time_beg, time_end;
  time_beg = clock();
  solver_->PardisoInit();
  time_end = clock();
  double time_consumption = (time_end - time_beg) / 1000.0;
  time_1_ += time_consumption;

  FunGrid(position_of_mesh_);
  int f0, f1, f2;
  double x0, y0, x1, y1, x2, y2;
  double dis, E_b, energy2 = 0;
  for (int i = 0; i < AV_F_N_; ++i) {
    f0 = V_F0_[AV_ID_[i]];
    f1 = V_F1_[AV_ID_[i]];
    f2 = V_F2_[AV_ID_[i]];

    x0 = pos[f0];
    y0 = pos[f0 + V_N_];

    x1 = pos[f1];
    y1 = pos[f1 + V_N_];

    x2 = pos[f2];
    y2 = pos[f2 + V_N_];

    dis = GetDistance(x0, y0, x1, y1, x2, y2);
    if (dis < 0) {
      std::cout << "distance is zero" << std::endl;
    }
    E_b = (1 - threhold / dis) * (1 - threhold / dis);
    energy2 += E_b;
  }
  energy_barrier_ = energy2;

  if (is_first_)
    is_first_ = false;
}
void ParaFun::InitArea() {
  area_.resize(F_N_);
  int src_t_num = shell_data_.mesh_faces_.rows();
  area_src_.resize(src_t_num);
  for (int i = 0; i < src_t_num; ++i) {
    area_src_[i] = shell_data_.mesh_face_area_(i);
    area_[i] = shell_data_.mesh_face_area_(i);
    if (area_[i] < 1e-15) {
      std::cout << "error" << std::endl;
    }
  }
  for (int i = src_t_num; i < F_N_; ++i) {
    area_[i] = shell_data_.shell_face_area_(i - src_t_num);
    if (area_[i] < 1e-15) {
      std::cout << "error" << std::endl;
    }
  }
}
void ParaFun::SetVirtualTri() {
  Eigen::VectorXi boundary_vertex = shell_data_.frame_ids_;
  BE_N_ = boundary_vertex.size();

  V_F_N_ = BE_N_ * (BE_N_ - 2);
  V_F0_.resize(V_F_N_);
  V_F1_.resize(V_F_N_);
  V_F2_.resize(V_F_N_);
  is_active_.resize(V_F_N_, -1);
  AV_ID_.reserve(V_F_N_);

  boundary_vertexID_.clear();
  boundary_vertexID_.resize(V_N_, -1);

  int id = 0;
  for (int i = 0; i < BE_N_; ++i) {
    boundary_vertexID_[boundary_vertex(i)] = id;
    id++;
  }

  int id_start, id_end, id_mid;
  int k;
  for (int i = 0; i < BE_N_; ++i) {
    id_start = boundary_vertex(i);
    id_end = boundary_vertex((i + 1) % BE_N_);
    k = 0;

    for (int j = 0; j < BE_N_; ++j) {
      if (id_start != boundary_vertex(j) && id_end != boundary_vertex(j)) {
        V_F0_[i * (BE_N_ - 2) + k] = id_start;
        V_F1_[i * (BE_N_ - 2) + k] = id_end;
        V_F2_[i * (BE_N_ - 2) + k] = boundary_vertex(j);
        ++k;
      }
    }
  }
}
void ParaFun::PreCalculate() {
  source_p00_.resize(F_N_);
  source_p01_.resize(F_N_);
  source_p10_.resize(F_N_);
  source_p11_.resize(F_N_);
  if (is_first_) {
    for (int i = 0; i < shell_data_.mesh_faces_.rows(); ++i) {
      double p00, p01, p10, p11;
      LocalCoordinateInverse(i, p00, p01, p10, p11);
      source_p00_[i] = p00;
      source_p01_[i] = p01;
      source_p10_[i] = p10;
      source_p11_[i] = p11;
    }
  }

  for (int i = shell_data_.mesh_faces_.rows(); i < F_N_; ++i) {
    double p00, p01, p10, p11;
    LocalCoordinateInverseScaf(i, p00, p01, p10, p11);
    source_p00_[i] = p00;
    source_p01_[i] = p01;
    source_p10_[i] = p10;
    source_p11_[i] = p11;
  }
  update_p00_ = source_p00_;
  update_p01_ = source_p01_;
  update_p10_ = source_p10_;
  update_p11_ = source_p11_;

  pardiso_ia_.clear();
  pardiso_ia_.reserve(2 * V_N_ + 1);
  pardiso_ja_.clear();
  pardiso_ja_.reserve(8 * V_N_ + BE_N_ * BE_N_ * 2);
  typedef Eigen::Triplet<int> T;
  std::vector<T> tripletlist;
  std::vector<std::set<int>> VV_tmp;
  VV_tmp.resize(V_N_);
  for (size_t i = 0; i < shell_data_.mesh_faces_.rows(); i++) {
    int vid[3];
    for (size_t j = 0; j < shell_data_.mesh_faces_.cols(); j++) {
      vid[j] = shell_data_.mesh_faces_(i, j);
    }
    VV_tmp[vid[0]].insert(vid[1]);
    VV_tmp[vid[0]].insert(vid[2]);
    VV_tmp[vid[1]].insert(vid[0]);
    VV_tmp[vid[1]].insert(vid[2]);
    VV_tmp[vid[2]].insert(vid[0]);
    VV_tmp[vid[2]].insert(vid[1]);
  }

  std::vector<int> s_vid;
  for (size_t i = 0; i < shell_data_.shell_faces_.rows(); i++) {
    s_vid.clear();
    for (size_t j = 0; j < shell_data_.shell_faces_.cols(); j++) {
      int s_id = shell_data_.shell_faces_(i, j);
      s_vid.push_back(s_id);
    }
    VV_tmp[s_vid[0]].insert(s_vid[1]);
    VV_tmp[s_vid[0]].insert(s_vid[2]);
    VV_tmp[s_vid[1]].insert(s_vid[0]);
    VV_tmp[s_vid[1]].insert(s_vid[2]);
    VV_tmp[s_vid[2]].insert(s_vid[0]);
    VV_tmp[s_vid[2]].insert(s_vid[1]);
  }

  std::vector<int> v_vid;
  for (int i = 0; i < shell_data_.frame_ids_.size(); ++i) {
    for (int j = 0; j < shell_data_.frame_ids_.size(); ++j) {
      if (j == i) {
        continue;
      }
      VV_tmp[shell_data_.frame_ids_(i)].insert(shell_data_.frame_ids_(j));
    }
  }

  for (int i = 0; i < V_N_; i++) {
    pardiso_ia_.push_back(pardiso_ja_.size());
    VV_tmp[i].insert(i);
    std::vector<int> row_id;
    for (auto &var : VV_tmp[i]) {
      row_id.push_back(var);
    }
    std::vector<int>::iterator iter =
        std::find(row_id.begin(), row_id.end(), i);
    int dd = 0;
    for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++) {
      pardiso_ja_.push_back(row_id[k]);
      tripletlist.push_back(T(i, row_id[k], dd));
      ++dd;
    }
    for (int k = 0; k < row_id.size(); k++) {
      pardiso_ja_.push_back(row_id[k] + V_N_);
      tripletlist.push_back(T(i, row_id[k] + V_N_, dd));
      ++dd;
    }
  }
  for (int i = V_N_; i < 2 * V_N_; i++) {
    pardiso_ia_.push_back(pardiso_ja_.size());
    std::vector<int> row_id;
    for (auto &var : VV_tmp[i - V_N_]) {
      row_id.push_back(var);
    }
    std::vector<int>::iterator iter =
        std::find(row_id.begin(), row_id.end(), i - V_N_);
    int dd = 0;
    for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++) {
      pardiso_ja_.push_back(row_id[k] + V_N_);
      tripletlist.push_back(T(i, row_id[k] + V_N_, dd));
      ++dd;
    }
  }

  Eigen::SparseMatrix<int> find_id_in_rows;
  find_id_in_rows.resize(2 * V_N_, 2 * V_N_);
  find_id_in_rows.setFromTriplets(tripletlist.begin(), tripletlist.end());
  pardiso_ia_.push_back(pardiso_ja_.size());
  id_h00_.resize(F_N_ + V_F_N_, -1);
  id_h01_.resize(F_N_ + V_F_N_, -1);
  id_h02_.resize(F_N_ + V_F_N_, -1);
  id_h03_.resize(F_N_ + V_F_N_, -1);
  id_h04_.resize(F_N_ + V_F_N_, -1);
  id_h05_.resize(F_N_ + V_F_N_, -1);
  id_h11_.resize(F_N_ + V_F_N_, -1);
  id_h12_.resize(F_N_ + V_F_N_, -1);
  id_h13_.resize(F_N_ + V_F_N_, -1);
  id_h14_.resize(F_N_ + V_F_N_, -1);
  id_h15_.resize(F_N_ + V_F_N_, -1);
  id_h22_.resize(F_N_ + V_F_N_, -1);
  id_h23_.resize(F_N_ + V_F_N_, -1);
  id_h24_.resize(F_N_ + V_F_N_, -1);
  id_h25_.resize(F_N_ + V_F_N_, -1);
  id_h33_.resize(F_N_ + V_F_N_, -1);
  id_h34_.resize(F_N_ + V_F_N_, -1);
  id_h35_.resize(F_N_ + V_F_N_, -1);
  id_h44_.resize(F_N_ + V_F_N_, -1);
  id_h45_.resize(F_N_ + V_F_N_, -1);
  id_h55_.resize(F_N_ + V_F_N_, -1);

  for (int i = 0; i < F_N_; i++) {
    int f0 = F_[0][i];
    int f1 = F_[1][i];
    int f2 = F_[2][i];
    int f3 = f0 + V_N_;
    int f4 = f1 + V_N_;
    int f5 = f2 + V_N_;
    int min01 = std::min(f0, f1);
    int max01 = f0 + f1 - min01;
    int min02 = std::min(f0, f2);
    int max02 = f0 + f2 - min02;
    int min12 = std::min(f1, f2);
    int max12 = f1 + f2 - min12;
    id_h00_[i] = pardiso_ia_[f0];
    id_h01_[i] = pardiso_ia_[min01] + find_id_in_rows.coeff(min01, max01);
    id_h02_[i] = pardiso_ia_[min02] + find_id_in_rows.coeff(min02, max02);
    id_h03_[i] = pardiso_ia_[f0] + find_id_in_rows.coeff(f0, f3);
    id_h04_[i] = pardiso_ia_[f0] + find_id_in_rows.coeff(f0, f4);
    id_h05_[i] = pardiso_ia_[f0] + find_id_in_rows.coeff(f0, f5);
    id_h11_[i] = pardiso_ia_[f1];
    id_h12_[i] = pardiso_ia_[min12] + find_id_in_rows.coeff(min12, max12);
    id_h13_[i] = pardiso_ia_[f1] + find_id_in_rows.coeff(f1, f3);
    id_h14_[i] = pardiso_ia_[f1] + find_id_in_rows.coeff(f1, f4);
    id_h15_[i] = pardiso_ia_[f1] + find_id_in_rows.coeff(f1, f5);
    id_h22_[i] = pardiso_ia_[f2];
    id_h23_[i] = pardiso_ia_[f2] + find_id_in_rows.coeff(f2, f3);
    id_h24_[i] = pardiso_ia_[f2] + find_id_in_rows.coeff(f2, f4);
    id_h25_[i] = pardiso_ia_[f2] + find_id_in_rows.coeff(f2, f5);
    id_h33_[i] = pardiso_ia_[f3];
    id_h34_[i] = pardiso_ia_[min01 + V_N_] +
                 find_id_in_rows.coeff(min01 + V_N_, max01 + V_N_);
    id_h35_[i] = pardiso_ia_[min02 + V_N_] +
                 find_id_in_rows.coeff(min02 + V_N_, max02 + V_N_);
    id_h44_[i] = pardiso_ia_[f4];
    id_h45_[i] = pardiso_ia_[min12 + V_N_] +
                 find_id_in_rows.coeff(min12 + V_N_, max12 + V_N_);
    id_h55_[i] = pardiso_ia_[f5];
  }

  for (int i = F_N_; i < F_N_ + V_F_N_; i++) {
    int f0 = V_F0_[i - F_N_];
    int f1 = V_F1_[i - F_N_];
    int f2 = V_F2_[i - F_N_];
    int f3 = V_F0_[i - F_N_] + V_N_;
    int f4 = V_F1_[i - F_N_] + V_N_;
    int f5 = V_F2_[i - F_N_] + V_N_;
    int min01 = std::min(f0, f1);
    int max01 = f0 + f1 - min01;
    int min02 = std::min(f0, f2);
    int max02 = f0 + f2 - min02;
    int min12 = std::min(f1, f2);
    int max12 = f1 + f2 - min12;
    id_h00_[i] = pardiso_ia_[f0];
    id_h01_[i] = pardiso_ia_[min01] + find_id_in_rows.coeff(min01, max01);
    id_h02_[i] = pardiso_ia_[min02] + find_id_in_rows.coeff(min02, max02);
    id_h03_[i] = pardiso_ia_[f0] + find_id_in_rows.coeff(f0, f3);
    id_h04_[i] = pardiso_ia_[f0] + find_id_in_rows.coeff(f0, f4);
    id_h05_[i] = pardiso_ia_[f0] + find_id_in_rows.coeff(f0, f5);
    id_h11_[i] = pardiso_ia_[f1];
    id_h12_[i] = pardiso_ia_[min12] + find_id_in_rows.coeff(min12, max12);
    id_h13_[i] = pardiso_ia_[f1] + find_id_in_rows.coeff(f1, f3);
    id_h14_[i] = pardiso_ia_[f1] + find_id_in_rows.coeff(f1, f4);
    id_h15_[i] = pardiso_ia_[f1] + find_id_in_rows.coeff(f1, f5);
    id_h22_[i] = pardiso_ia_[f2];
    id_h23_[i] = pardiso_ia_[f2] + find_id_in_rows.coeff(f2, f3);
    id_h24_[i] = pardiso_ia_[f2] + find_id_in_rows.coeff(f2, f4);
    id_h25_[i] = pardiso_ia_[f2] + find_id_in_rows.coeff(f2, f5);
    id_h33_[i] = pardiso_ia_[f3];
    id_h34_[i] = pardiso_ia_[min01 + V_N_] +
                 find_id_in_rows.coeff(min01 + V_N_, max01 + V_N_);
    id_h35_[i] = pardiso_ia_[min02 + V_N_] +
                 find_id_in_rows.coeff(min02 + V_N_, max02 + V_N_);
    id_h44_[i] = pardiso_ia_[f4];
    id_h45_[i] = pardiso_ia_[min12 + V_N_] +
                 find_id_in_rows.coeff(min12 + V_N_, max12 + V_N_);
    id_h55_[i] = pardiso_ia_[f5];
  }
  int nnz = pardiso_ja_.size();
  // std:cout << "density_: " << double(nnz) << std::endl;
  density_ += double(nnz) / 4 / V_N_ / V_N_;
}
void ParaFun::HandleMinTri() {
  double min_bnd_edge_len = std::numeric_limits<double>::infinity();
  int acc_bnd = 0;
  for (int i = 0; i < shell_data_.bnd_sizes_.size(); i++) {
    int current_size = shell_data_.bnd_sizes_[i];

    for (int e = acc_bnd; e < acc_bnd + current_size - 1; e++) {
      min_bnd_edge_len = (std::min)(
          min_bnd_edge_len,
          (shell_data_.whole_uv_.row(shell_data_.internal_bnd_(e)) -
           shell_data_.whole_uv_.row(shell_data_.internal_bnd_(e + 1)))
              .squaredNorm());
    }
    min_bnd_edge_len = (std::min)(
        min_bnd_edge_len,
        (shell_data_.whole_uv_.row(shell_data_.internal_bnd_(acc_bnd)) -
         shell_data_.whole_uv_.row(
             shell_data_.internal_bnd_(acc_bnd + current_size - 1)))
            .squaredNorm());
    acc_bnd += current_size;
  }

  area_threshold_ = min_bnd_edge_len / 4.0;
}

double ParaFun::BPE(bool is_ip_convrate, bool is_slim_convrate) {
  shell_data_.whole_uv_pre_ = shell_data_.whole_uv_;
  if (is_ip_convrate) {
    UpdateSourceSameT();
  }
  bool is_interp = is_ip_convrate && (Intp_T_Min_ < 0.999);
  bool is_slim = is_interp && is_slim_convrate && (change_to_cm_flag_ < 0.99);
  // is_interp = false;
  if (is_slim) {
    SLIM(is_interp);
    // CM(is_interp);
  } else {
    CM(is_interp);
  }
  shell_data_.whole_uv_ =
      Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::ColMajor>>(
          position_of_mesh_.data(), total_num_, kDim);

  return energy_mesh_ + area_.back() * energy_shell_ +
         barrer_coef_ * energy_barrier_;
}
void ParaFun::UpdateSourceSameT() {
  double t_min = 1;
  int geqK = 0;
  int update_fn = shell_data_.mesh_faces_.rows();
  std::vector<double> all_s0;
  all_s0.resize(update_fn);
  std::vector<double> all_s1;
  all_s1.resize(update_fn);
  std::vector<double> all_w00;
  all_w00.resize(update_fn);
  std::vector<double> all_w01;
  all_w01.resize(update_fn);
  std::vector<double> all_w10;
  all_w10.resize(update_fn);
  std::vector<double> all_w11;
  all_w11.resize(update_fn);

  int f0, f1, f2;
  double x0, y0, x1, y1, x2, y2;
  double det;
  double E_d;
  double tt;
  double new_sig0, new_sig1;
  double j00, j01, j10, j11;
  double p00, p01, p10, p11;
  double q00, q01, q10, q11;
  double *position = position_of_mesh_.data();

  for (int i = 0; i < update_fn; ++i) {
    f0 = F_[0][i];
    f1 = F_[1][i];
    f2 = F_[2][i];
    x0 = position[f0];
    y0 = position[f0 + total_num_];
    x1 = position[f1];
    y1 = position[f1 + total_num_];
    x2 = position[f2];
    y2 = position[f2 + total_num_];

    q00 = x1 - x0;
    q01 = x2 - x0;
    q10 = y1 - y0;
    q11 = y2 - y0;
    p00 = source_p00_[i];
    p01 = source_p01_[i];
    p10 = source_p10_[i];
    p11 = source_p11_[i];
    j00 = p00 * q00 + p10 * q01;
    j01 = p01 * q00 + p11 * q01;
    j10 = p00 * q10 + p10 * q11;
    j11 = p01 * q10 + p11 * q11;

    det = j00 * j11 - j01 * j10;
    E_d =
        (1 + 1 / (det * det)) * (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);

    double alpha_0 = j00 + j11;
    double alpha_1 = j10 - j01;
    double beta_0 = j00 - j11;
    double beta_1 = j10 + j01;
    double alpha_norm = 0.5 * sqrt(alpha_0 * alpha_0 + alpha_1 * alpha_1);
    double beta_norm = 0.5 * sqrt(beta_0 * beta_0 + beta_1 * beta_1);

    double sig0 = alpha_norm + beta_norm;
    double sig1 = alpha_norm - beta_norm;
    all_s0[i] = sig0;
    all_s1[i] = sig1;

    if (beta_norm < 1e-15) {
      all_w00[i] = 0.0;
      all_w01[i] = 0.0;
      all_w10[i] = 0.0;
      all_w11[i] = 0.0;
    } else {
      double temp = 1 / (sig1 * sig1 - sig0 * sig0);
      all_w00[i] =
          temp * (j00 * j00 + j10 * j10 - 0.5 * (sig0 * sig0 + sig1 * sig1));
      all_w01[i] = temp * (j00 * j01 + j10 * j11);
      all_w10[i] = temp * (j01 * j00 + j11 * j10);
      all_w11[i] =
          temp * (j01 * j01 + j11 * j11 - 0.5 * (sig0 * sig0 + sig1 * sig1));
    }

    if (E_d <= bound_distortion_K_) {
      geqK++;
    } else {
      tt = NewtonEquation(sig0, sig1, bound_distortion_K_);
      if (tt < t_min) {
        t_min = tt;
      }
    }
  }

  change_to_cm_flag_ = (double)geqK / update_fn;
  update_p00_ = source_p00_;
  update_p01_ = source_p01_;
  update_p10_ = source_p10_;
  update_p11_ = source_p11_;

  for (int i = 0; i < update_fn; ++i) {
    double sig0 = all_s0[i];
    double sig1 = all_s1[i];
    new_sig0 = pow(sig0, t_min - 1);
    new_sig1 = pow(sig1, t_min - 1);

    double delta_new = new_sig1 - new_sig0;
    double plus_new = 0.5 * (new_sig1 + new_sig0);
    double w00 = delta_new * all_w00[i] + plus_new;
    double w01 = delta_new * all_w01[i];
    double w10 = delta_new * all_w10[i];
    double w11 = delta_new * all_w11[i] + plus_new;

    p00 = source_p00_[i];
    p01 = source_p01_[i];
    p10 = source_p10_[i];
    p11 = source_p11_[i];
    update_p00_[i] = p00 * w00 + p01 * w10;
    update_p01_[i] = p00 * w01 + p01 * w11;
    update_p10_[i] = p10 * w00 + p11 * w10;
    update_p11_[i] = p10 * w01 + p11 * w11;
  }

  Intp_T_Min_ = t_min;
}
void ParaFun::SLIM(bool is_interp) {
  double area_now;
  int f0, f1, f2;
  double j00, j01, j10, j11;
  double p00, p01, p10, p11;
  double q00, q01, q10, q11;

  double x0, y0, x1, y1, x2, y2;

  double alpha_norm, beta_norm;

  double alpha_0, alpha_1, beta_0, beta_1;

  double sig0, sig1;

  double det, tr;
  double r0, r1, r2, r3;
  double d00, d01, d02, d10, d11, d12;

  double new_sig0, new_sig1;
  double temp;
  double w00, w01, w10, w11;
  double p1, p2, p3, w1, w2, w3;

  double h00, h01, h02, h03, h04, h05, h11, h12, h13, h14, h15, h22, h23, h24,
      h25, h33, h34, h35, h44, h45, h55;
  double *position = position_of_mesh_.data();

  int nnz = pardiso_ja_.size();
  pardiso_a_.clear();
  pardiso_b_.clear();
  pardiso_a_.resize(nnz, 0.0);
  pardiso_b_.resize(2 * V_N_, 0.0);

  double *tmp_p00;
  double *tmp_p01;
  double *tmp_p10;
  double *tmp_p11;

  if (is_interp) {
    tmp_p00 = update_p00_.data();
    tmp_p01 = update_p01_.data();
    tmp_p10 = update_p10_.data();
    tmp_p11 = update_p11_.data();
  } else {
    tmp_p00 = source_p00_.data();
    tmp_p01 = source_p01_.data();
    tmp_p10 = source_p10_.data();
    tmp_p11 = source_p11_.data();
  }

  int src_t_num = shell_data_.mesh_faces_.rows();
  for (int i = 0; i < src_t_num; ++i) {
    area_now = area_[i];
    f0 = F_[0][i];
    f1 = F_[1][i];
    f2 = F_[2][i];
    x0 = position[f0];
    y0 = position[f0 + total_num_];
    x1 = position[f1];
    y1 = position[f1 + total_num_];
    x2 = position[f2];
    y2 = position[f2 + total_num_];

    q00 = x1 - x0;
    q01 = x2 - x0;
    q10 = y1 - y0;
    q11 = y2 - y0;
    p00 = tmp_p00[i];
    p01 = tmp_p01[i];
    p10 = tmp_p10[i];
    p11 = tmp_p11[i];
    j00 = p00 * q00 + p10 * q01;
    j01 = p01 * q00 + p11 * q01;
    j10 = p00 * q10 + p10 * q11;
    j11 = p01 * q10 + p11 * q11;

    alpha_0 = j00 + j11;
    alpha_1 = j10 - j01;
    beta_0 = j00 - j11;
    beta_1 = j10 + j01;
    alpha_norm = 0.5 * sqrt(alpha_0 * alpha_0 + alpha_1 * alpha_1);
    beta_norm = 0.5 * sqrt(beta_0 * beta_0 + beta_1 * beta_1);

    sig0 = alpha_norm + beta_norm;
    sig1 = alpha_norm - beta_norm;
    new_sig0 =
        sqrt(1 + 1 / sig0 + 1 / (sig0 * sig0) + 1 / (sig0 * sig0 * sig0));
    new_sig1 =
        sqrt(1 + 1 / sig1 + 1 / (sig1 * sig1) + 1 / (sig1 * sig1 * sig1));

    if (abs(sig1 - sig0) < 1e-10) {
      temp = 0;
    } else {
      temp = (new_sig1 - new_sig0) / (sig1 * sig1 - sig0 * sig0);
    }

    w00 = temp * (j00 * j00 + j01 * j01 - 0.5 * (sig0 * sig0 + sig1 * sig1)) +
          0.5 * (new_sig0 + new_sig1);
    w01 = temp * (j00 * j10 + j01 * j11);
    w10 = temp * (j10 * j00 + j11 * j01);
    w11 = temp * (j10 * j10 + j11 * j11 - 0.5 * (sig0 * sig0 + sig1 * sig1)) +
          0.5 * (new_sig0 + new_sig1);
    p1 = p00 * p00 + p01 * p01;
    p2 = p00 * p10 + p01 * p11;
    p3 = p10 * p10 + p11 * p11;
    w1 = w00 * w00 + w10 * w10;
    w2 = w00 * w01 + w10 * w11;
    w3 = w01 * w01 + w11 * w11;

    h00 = area_now * (p1 + p2 + p2 + p3) * w1;
    h01 = -area_now * (p1 + p2) * w1;
    h02 = -area_now * (p2 + p3) * w1;
    h03 = area_now * (p1 + p2 + p2 + p3) * w2;
    h04 = -area_now * (p1 + p2) * w2;
    h05 = -area_now * (p2 + p3) * w2;
    h11 = area_now * p1 * w1;
    h12 = area_now * p2 * w1;
    h13 = -area_now * (p1 + p2) * w2;
    h14 = area_now * p1 * w2;
    h15 = area_now * p2 * w2;
    h22 = area_now * p3 * w1;
    h23 = -area_now * (p2 + p3) * w2;
    h24 = area_now * p2 * w2;
    h25 = area_now * p3 * w2;
    h33 = area_now * (p1 + p2 + p2 + p3) * w3;
    h34 = -area_now * (p1 + p2) * w3;
    h35 = -area_now * (p2 + p3) * w3;
    h44 = area_now * p1 * w3;
    h45 = area_now * p2 * w3;
    h55 = area_now * p3 * w3;

    det = j00 * j11 - j01 * j10;
    if (det < 0) {
      std::cout << i << std::endl;
    }
    tr = (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);
    d00 = -p00 - p10;
    d01 = p00;
    d02 = p10;
    d10 = -p01 - p11;
    d11 = p01;
    d12 = p11;
    r0 =
        area_now * ((1 + 1 / (det * det)) * j00 - tr * j11 / (det * det * det));
    r1 =
        area_now * ((1 + 1 / (det * det)) * j01 + tr * j10 / (det * det * det));
    r2 =
        area_now * ((1 + 1 / (det * det)) * j10 + tr * j01 / (det * det * det));
    r3 =
        area_now * ((1 + 1 / (det * det)) * j11 - tr * j00 / (det * det * det));

    pardiso_b_[f0] -= r0 * d00 + r1 * d10;
    pardiso_b_[f1] -= r0 * d01 + r1 * d11;
    pardiso_b_[f2] -= r0 * d02 + r1 * d12;
    pardiso_b_[f0 + V_N_] -= r2 * d00 + r3 * d10;
    pardiso_b_[f1 + V_N_] -= r2 * d01 + r3 * d11;
    pardiso_b_[f2 + V_N_] -= r2 * d02 + r3 * d12;

    pardiso_a_[id_h00_[i]] += h00;
    pardiso_a_[id_h01_[i]] += h01;
    pardiso_a_[id_h02_[i]] += h02;
    pardiso_a_[id_h03_[i]] += h03;
    pardiso_a_[id_h04_[i]] += h04;
    pardiso_a_[id_h05_[i]] += h05;
    pardiso_a_[id_h11_[i]] += h11;
    pardiso_a_[id_h12_[i]] += h12;
    pardiso_a_[id_h13_[i]] += h13;
    pardiso_a_[id_h14_[i]] += h14;
    pardiso_a_[id_h15_[i]] += h15;
    pardiso_a_[id_h22_[i]] += h22;
    pardiso_a_[id_h23_[i]] += h23;
    pardiso_a_[id_h24_[i]] += h24;
    pardiso_a_[id_h25_[i]] += h25;
    pardiso_a_[id_h33_[i]] += h33;
    pardiso_a_[id_h34_[i]] += h34;
    pardiso_a_[id_h35_[i]] += h35;
    pardiso_a_[id_h44_[i]] += h44;
    pardiso_a_[id_h45_[i]] += h45;
    pardiso_a_[id_h55_[i]] += h55;
  }

  for (int i = src_t_num; i < F_N_; ++i) {
    area_now = area_[i];
    f0 = F_[0][i];
    f1 = F_[1][i];
    f2 = F_[2][i];
    x0 = position[f0];
    y0 = position[f0 + total_num_];
    x1 = position[f1];
    y1 = position[f1 + total_num_];
    x2 = position[f2];
    y2 = position[f2 + total_num_];

    q00 = x1 - x0;
    q01 = x2 - x0;
    q10 = y1 - y0;
    q11 = y2 - y0;
    p00 = tmp_p00[i];
    p01 = tmp_p01[i];
    p10 = tmp_p10[i];
    p11 = tmp_p11[i];
    j00 = p00 * q00 + p10 * q01;
    j01 = p01 * q00 + p11 * q01;
    j10 = p00 * q10 + p10 * q11;
    j11 = p01 * q10 + p11 * q11;
    alpha_0 = j00 + j11;
    alpha_1 = j10 - j01;
    beta_0 = j00 - j11;
    beta_1 = j10 + j01;

    alpha_norm = 0.5 * sqrt(alpha_0 * alpha_0 + alpha_1 * alpha_1);
    beta_norm = 0.5 * sqrt(beta_0 * beta_0 + beta_1 * beta_1);
    sig0 = alpha_norm + beta_norm;
    sig1 = alpha_norm - beta_norm;
    new_sig0 =
        sqrt(1 + 1 / sig0 + 1 / (sig0 * sig0) + 1 / (sig0 * sig0 * sig0));
    new_sig1 =
        sqrt(1 + 1 / sig1 + 1 / (sig1 * sig1) + 1 / (sig1 * sig1 * sig1));

    if (abs(sig1 - sig0) < 1e-10) {
      temp = 0;
    } else {
      temp = (new_sig1 - new_sig0) / (sig1 * sig1 - sig0 * sig0);
    }

    w00 = temp * (j00 * j00 + j01 * j01 - 0.5 * (sig0 * sig0 + sig1 * sig1)) +
          0.5 * (new_sig0 + new_sig1);
    w01 = temp * (j00 * j10 + j01 * j11);
    w10 = temp * (j10 * j00 + j11 * j01);
    w11 = temp * (j10 * j10 + j11 * j11 - 0.5 * (sig0 * sig0 + sig1 * sig1)) +
          0.5 * (new_sig0 + new_sig1);
    p1 = p00 * p00 + p01 * p01;
    p2 = p00 * p10 + p01 * p11;
    p3 = p10 * p10 + p11 * p11;
    w1 = w00 * w00 + w10 * w10;
    w2 = w00 * w01 + w10 * w11;
    w3 = w01 * w01 + w11 * w11;

    h00 = area_now * (p1 + p2 + p2 + p3) * w1;
    h01 = -area_now * (p1 + p2) * w1;
    h02 = -area_now * (p2 + p3) * w1;
    h03 = area_now * (p1 + p2 + p2 + p3) * w2;
    h04 = -area_now * (p1 + p2) * w2;
    h05 = -area_now * (p2 + p3) * w2;
    h11 = area_now * p1 * w1;
    h12 = area_now * p2 * w1;
    h13 = -area_now * (p1 + p2) * w2;
    h14 = area_now * p1 * w2;
    h15 = area_now * p2 * w2;
    h22 = area_now * p3 * w1;
    h23 = -area_now * (p2 + p3) * w2;
    h24 = area_now * p2 * w2;
    h25 = area_now * p3 * w2;
    h33 = area_now * (p1 + p2 + p2 + p3) * w3;
    h34 = -area_now * (p1 + p2) * w3;
    h35 = -area_now * (p2 + p3) * w3;
    h44 = area_now * p1 * w3;
    h45 = area_now * p2 * w3;
    h55 = area_now * p3 * w3;

    det = j00 * j11 - j01 * j10;
    tr = (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);
    d00 = -p00 - p10;
    d01 = p00;
    d02 = p10;
    d10 = -p01 - p11;
    d11 = p01;
    d12 = p11;
    r0 =
        area_now * ((1 + 1 / (det * det)) * j00 - tr * j11 / (det * det * det));
    r1 =
        area_now * ((1 + 1 / (det * det)) * j01 + tr * j10 / (det * det * det));
    r2 =
        area_now * ((1 + 1 / (det * det)) * j10 + tr * j01 / (det * det * det));
    r3 =
        area_now * ((1 + 1 / (det * det)) * j11 - tr * j00 / (det * det * det));

    pardiso_b_[f0] -= r0 * d00 + r1 * d10;
    pardiso_b_[f1] -= r0 * d01 + r1 * d11;
    pardiso_b_[f2] -= r0 * d02 + r1 * d12;
    pardiso_b_[f0 + V_N_] -= r2 * d00 + r3 * d10;
    pardiso_b_[f1 + V_N_] -= r2 * d01 + r3 * d11;
    pardiso_b_[f2 + V_N_] -= r2 * d02 + r3 * d12;

    pardiso_a_[id_h00_[i]] += h00;
    pardiso_a_[id_h01_[i]] += h01;
    pardiso_a_[id_h02_[i]] += h02;
    pardiso_a_[id_h03_[i]] += h03;
    pardiso_a_[id_h04_[i]] += h04;
    pardiso_a_[id_h05_[i]] += h05;
    pardiso_a_[id_h11_[i]] += h11;
    pardiso_a_[id_h12_[i]] += h12;
    pardiso_a_[id_h13_[i]] += h13;
    pardiso_a_[id_h14_[i]] += h14;
    pardiso_a_[id_h15_[i]] += h15;
    pardiso_a_[id_h22_[i]] += h22;
    pardiso_a_[id_h23_[i]] += h23;
    pardiso_a_[id_h24_[i]] += h24;
    pardiso_a_[id_h25_[i]] += h25;
    pardiso_a_[id_h33_[i]] += h33;
    pardiso_a_[id_h34_[i]] += h34;
    pardiso_a_[id_h35_[i]] += h35;
    pardiso_a_[id_h44_[i]] += h44;
    pardiso_a_[id_h45_[i]] += h45;
    pardiso_a_[id_h55_[i]] += h55;
  }

  double len01, len0i, len1i, coef1, coef2, dis;
  int i;
  double s0, s1, e0, e1, p0;
  double h_u, uu, aa, bb, cc;
  double a1x0, a1x1, a1x2, a1x3, a1x4, a1x5, a2x0, a2x1, a2x2, a2x3, a2x4, a2x5;

  AV_F_N_H_ = AV_F_N_;
  V_F0_H_ = V_F0_;
  V_F1_H_ = V_F1_;
  V_F2_H_ = V_F2_;
  for (int ii = 0; ii < AV_F_N_; ++ii) {
    f0 = V_F0_[AV_ID_[ii]];
    f1 = V_F1_[AV_ID_[ii]];
    f2 = V_F2_[AV_ID_[ii]];
    x0 = position[f0];
    y0 = position[f0 + V_N_];
    x1 = position[f1];
    y1 = position[f1 + V_N_];
    x2 = position[f2];
    y2 = position[f2 + V_N_];

    s0 = x0 + x1;
    s1 = y0 + y1;
    e0 = x1 + x0;
    e1 = y1 + y0;
    p0 = x2;
    p1 = y2;
    len01 = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    len0i = sqrt((p0 - s0) * (p0 - s0) + (p1 - s1) * (p1 - s1));
    len1i = sqrt((p0 - e0) * (p0 - e0) + (p1 - e1) * (p1 - e1));
    dis = len0i + len1i - len01;

    h_u = 2 * barrer_coef_ * threhold * (dis - threhold) / dis / dis / dis;
    uu = 2 * barrer_coef_ * threhold * (3 * threhold - 2 * dis) / dis / dis /
         dis / dis;

    a1x0 = (x0 - x2) / len0i - (x0 - x1) / len01;
    a1x1 = (x1 - x2) / len1i - (x1 - x0) / len01;
    a1x2 = (x2 - x0) / len0i + (x2 - x1) / len1i;
    a1x3 = (y0 - y2) / len0i - (y0 - y1) / len01;
    a1x4 = (y1 - y2) / len1i - (y1 - y0) / len01;
    a1x5 = (y2 - y0) / len0i + (y2 - y1) / len1i;

    coef1 = len01 * len01;
    coef2 = h_u / coef1 / len01;

    aa = (x0 - x1) * (x0 - x1);
    bb = (x1 - x0) * (y1 - y0);
    cc = (y0 - y1) * (y0 - y1);
    h00 = coef2 * (aa - coef1);
    h01 = coef2 * (coef1 - aa);
    h02 = 0.0;
    h03 = coef2 * bb;
    h04 = -coef2 * bb;
    h05 = 0.0;
    h11 = coef2 * (aa - coef1);
    h12 = 0.0;
    h13 = -coef2 * bb;
    h14 = coef2 * bb;
    h15 = 0.0;
    h22 = 0.0;
    h23 = 0.0;
    h24 = 0.0;
    h25 = 0.0;
    h33 = coef2 * (cc - coef1);
    h34 = coef2 * (coef1 - cc);
    h35 = 0.0;
    h44 = coef2 * (cc - coef1);
    h45 = 0.0;
    h55 = 0.0;

    h00 += uu * a1x0 * a1x0;
    h01 += uu * a1x0 * a1x1;
    h02 += uu * a1x0 * a1x2;
    h03 += uu * a1x0 * a1x3;
    h04 += uu * a1x0 * a1x4;
    h05 += uu * a1x0 * a1x5;
    h11 += uu * a1x1 * a1x1;
    h12 += uu * a1x1 * a1x2;
    h13 += uu * a1x1 * a1x3;
    h14 += uu * a1x1 * a1x4;
    h15 += uu * a1x1 * a1x5;
    h22 += uu * a1x2 * a1x2;
    h23 += uu * a1x2 * a1x3;
    h24 += uu * a1x2 * a1x4;
    h25 += uu * a1x2 * a1x5;
    h33 += uu * a1x3 * a1x3;
    h34 += uu * a1x3 * a1x4;
    h35 += uu * a1x3 * a1x5;
    h44 += uu * a1x4 * a1x4;
    h45 += uu * a1x4 * a1x5;
    h55 += uu * a1x5 * a1x5;

    pardiso_b_[f0] -= h_u * (a1x0);
    pardiso_b_[f1] -= h_u * (a1x1);
    pardiso_b_[f2] -= h_u * (a1x2);
    pardiso_b_[f0 + V_N_] -= h_u * (a1x3);
    pardiso_b_[f1 + V_N_] -= h_u * (a1x4);
    pardiso_b_[f2 + V_N_] -= h_u * (a1x5);

    i = F_N_ + AV_ID_[ii];
    pardiso_a_[id_h00_[i]] += h00;
    pardiso_a_[id_h01_[i]] += h01;
    pardiso_a_[id_h02_[i]] += h02;
    pardiso_a_[id_h03_[i]] += h03;
    pardiso_a_[id_h04_[i]] += h04;
    pardiso_a_[id_h05_[i]] += h05;
    pardiso_a_[id_h11_[i]] += h11;
    pardiso_a_[id_h12_[i]] += h12;
    pardiso_a_[id_h13_[i]] += h13;
    pardiso_a_[id_h14_[i]] += h14;
    pardiso_a_[id_h15_[i]] += h15;
    pardiso_a_[id_h22_[i]] += h22;
    pardiso_a_[id_h23_[i]] += h23;
    pardiso_a_[id_h24_[i]] += h24;
    pardiso_a_[id_h25_[i]] += h25;
    pardiso_a_[id_h33_[i]] += h33;
    pardiso_a_[id_h34_[i]] += h34;
    pardiso_a_[id_h35_[i]] += h35;
    pardiso_a_[id_h44_[i]] += h44;
    pardiso_a_[id_h45_[i]] += h45;
    pardiso_a_[id_h55_[i]] += h55;
  }

  pardiso_a_[0] += 1.0;
  pardiso_a_[pardiso_ia_[V_N_]] += 1.0;
  solver_->a_ = pardiso_a_;
  solver_->rhs_ = pardiso_b_;
  long time_beg, time_end;
  time_beg = clock();
  solver_->Factorize();
  time_end = clock();
  double time_consumption = (time_end - time_beg) / 1000.0;
  time_2_ += time_consumption;
  time_beg = clock();
  solver_->PardisoSolver();
  time_end = clock();
  time_consumption = (time_end - time_beg) / 1000.0;
  time_3_ += time_consumption;
  // std::cout << "numerical factorize" << time_consumption << std::endl;

  std::vector<double> result_d = solver_->result_;
  Eigen::VectorXd negative_grad(2 * total_num_), d(2 * total_num_);
  for (int i = 0; i < V_N_; i++) {
    negative_grad(i) = pardiso_b_[i];
    negative_grad(i + total_num_) = pardiso_b_[i + V_N_];
    d(i) = result_d[i];
    d(i + total_num_) = result_d[i + V_N_];
  }
  solver_->FreeNumericalFactorizationMemory();

  double temp_t;
  MaxStep(position_of_mesh_, d, temp_t);
  double alpha = std::min(1.0, 0.8 * temp_t);
  BacktrackingLineSearch(position_of_mesh_, d, negative_grad, alpha, is_interp);
  // std::cout << "slim step: " << alpha << std::endl;
  position_of_mesh_ += alpha * d;
  EnergySource();
}
void ParaFun::CM(bool is_interp) {
  double area_now;
  int f0, f1, f2;
  double j00, j01, j10, j11;
  double p00, p01, p10, p11;
  double q00, q01, q10, q11;

  double x0, y0, x1, y1, x2, y2;
  double hi_0, hi_1;
  double alpha_0, alpha_1, beta_0, beta_1;
  double ss1, ss2, sig0, sig1;
  double alpha_norm, beta_norm;
  double h_u, h_v, walpha, wbeta;
  double a1x0, a1x1, a1x2, a1x3, a1x4, a1x5, a2x0, a2x1, a2x2, a2x3, a2x4, a2x5;

  double aa, bb, cc;
  double uu, vv, uv;
  double u, v;
  double h00, h01, h02, h03, h04, h05, h11, h12, h13, h14, h15, h22, h23, h24,
      h25, h33, h34, h35, h44, h45, h55;
  double *position = position_of_mesh_.data();
  int nnz = pardiso_ja_.size();
  pardiso_a_.clear();
  pardiso_b_.clear();
  pardiso_a_.resize(nnz, 0.0);
  pardiso_b_.resize(2 * V_N_, 0.0);

  double *tmp_p00;
  double *tmp_p01;
  double *tmp_p10;
  double *tmp_p11;
  if (is_interp) {
    tmp_p00 = update_p00_.data();
    tmp_p01 = update_p01_.data();
    tmp_p10 = update_p10_.data();
    tmp_p11 = update_p11_.data();
  } else {
    tmp_p00 = source_p00_.data();
    tmp_p01 = source_p01_.data();
    tmp_p10 = source_p10_.data();
    tmp_p11 = source_p11_.data();
  }

  for (int i = 0; i < shell_data_.mesh_faces_.rows(); ++i) {
    area_now = area_[i];
    f0 = F_[0][i];
    f1 = F_[1][i];
    f2 = F_[2][i];
    x0 = position[f0];
    y0 = position[f0 + V_N_];
    x1 = position[f1];
    y1 = position[f1 + V_N_];
    x2 = position[f2];
    y2 = position[f2 + V_N_];

    q00 = x1 - x0;
    q01 = x2 - x0;
    q10 = y1 - y0;
    q11 = y2 - y0;
    p00 = tmp_p00[i];
    p01 = tmp_p01[i];
    p10 = tmp_p10[i];
    p11 = tmp_p11[i];
    j00 = p00 * q00 + p10 * q01;
    j01 = p01 * q00 + p11 * q01;
    j10 = p00 * q10 + p10 * q11;
    j11 = p01 * q10 + p11 * q11;

    alpha_0 = j00 + j11;
    alpha_1 = j10 - j01;
    beta_0 = j00 - j11;
    beta_1 = j10 + j01;
    alpha_norm = 0.5 * sqrt(alpha_0 * alpha_0 + alpha_1 * alpha_1);
    beta_norm = 0.5 * sqrt(beta_0 * beta_0 + beta_1 * beta_1);
    /*if (beta_norm < 1e-30)
    {
            beta_norm = 1e-10;
    }*/
    ss1 = (p00) * (p00 + p10) + (p01) * (p01 + p11);
    ss2 = (p10) * (p00 + p10) + (p11) * (p01 + p11);

    double h1 = p00 * p00 + p01 * p01;
    double h2 = p00 * p10 + p01 * p11;
    double h3 = p10 * p10 + p11 * p11;
    double h4 = p00 * p11 - p01 * p10;

    a1x0 = alpha_0 * (-p00 - p10) + alpha_1 * (p01 + p11);
    a1x1 = alpha_0 * p00 - alpha_1 * p01;
    a1x2 = alpha_0 * p10 - alpha_1 * p11;
    a1x3 = alpha_0 * (-p01 - p11) + alpha_1 * (-p00 - p10);
    a1x4 = alpha_0 * p01 + alpha_1 * p00;
    a1x5 = alpha_0 * p11 + alpha_1 * p10;
    a2x0 = beta_0 * (-p00 - p10) + beta_1 * (-p01 - p11);
    a2x1 = beta_0 * p00 + beta_1 * p01;
    a2x2 = beta_0 * p10 + beta_1 * p11;
    a2x3 = beta_0 * (p01 + p11) + beta_1 * (-p00 - p10);
    a2x4 = -beta_0 * p01 + beta_1 * p00;
    a2x5 = -beta_0 * p11 + beta_1 * p10;
    sig0 = alpha_norm + beta_norm;
    sig1 = alpha_norm - beta_norm;

    hi_0 = 2 + 6 * 1 / (sig0 * sig0 * sig0 * sig0);
    hi_1 = 2 + 6 * 1 / (sig1 * sig1 * sig1 * sig1);
    aa = 0.25 / alpha_norm;
    bb = 0.25 / beta_norm;
    uu = aa * aa * (area_now * hi_0 + area_now * hi_1);
    vv = bb * bb * (area_now * hi_0 + area_now * hi_1);
    uv = aa * bb * (area_now * hi_0 - area_now * hi_1);
    h_u = area_now * (2 * sig0 - 2 * 1 / (sig0 * sig0 * sig0));
    h_v = area_now * (2 * sig1 - 2 * 1 / (sig1 * sig1 * sig1));

    walpha = h_u + h_v;
    wbeta = h_u - h_v;
    double hwa1 = (walpha * 0.25 / alpha_norm);
    double hwa2 =
        -(walpha * 0.25 * 0.25 / (alpha_norm * alpha_norm * alpha_norm));
    double hwb1 = (wbeta * 0.25 / beta_norm);
    double hwb2 = -(wbeta * 0.25 * 0.25 / (beta_norm * beta_norm * beta_norm));

    h00 = uu * a1x0 * a1x0 + vv * a2x0 * a2x0 + uv * a1x0 * a2x0 +
          uv * a2x0 * a1x0;
    h01 = uu * a1x0 * a1x1 + vv * a2x0 * a2x1 + uv * a1x0 * a2x1 +
          uv * a2x0 * a1x1;
    h02 = uu * a1x0 * a1x2 + vv * a2x0 * a2x2 + uv * a1x0 * a2x2 +
          uv * a2x0 * a1x2;
    h03 = uu * a1x0 * a1x3 + vv * a2x0 * a2x3 + uv * a1x0 * a2x3 +
          uv * a2x0 * a1x3;
    h04 = uu * a1x0 * a1x4 + vv * a2x0 * a2x4 + uv * a1x0 * a2x4 +
          uv * a2x0 * a1x4;
    h05 = uu * a1x0 * a1x5 + vv * a2x0 * a2x5 + uv * a1x0 * a2x5 +
          uv * a2x0 * a1x5;
    h11 = uu * a1x1 * a1x1 + vv * a2x1 * a2x1 + uv * a1x1 * a2x1 +
          uv * a2x1 * a1x1;
    h12 = uu * a1x1 * a1x2 + vv * a2x1 * a2x2 + uv * a1x1 * a2x2 +
          uv * a2x1 * a1x2;
    h13 = uu * a1x1 * a1x3 + vv * a2x1 * a2x3 + uv * a1x1 * a2x3 +
          uv * a2x1 * a1x3;
    h14 = uu * a1x1 * a1x4 + vv * a2x1 * a2x4 + uv * a1x1 * a2x4 +
          uv * a2x1 * a1x4;
    h15 = uu * a1x1 * a1x5 + vv * a2x1 * a2x5 + uv * a1x1 * a2x5 +
          uv * a2x1 * a1x5;
    h22 = uu * a1x2 * a1x2 + vv * a2x2 * a2x2 + uv * a1x2 * a2x2 +
          uv * a2x2 * a1x2;
    h23 = uu * a1x2 * a1x3 + vv * a2x2 * a2x3 + uv * a1x2 * a2x3 +
          uv * a2x2 * a1x3;
    h24 = uu * a1x2 * a1x4 + vv * a2x2 * a2x4 + uv * a1x2 * a2x4 +
          uv * a2x2 * a1x4;
    h25 = uu * a1x2 * a1x5 + vv * a2x2 * a2x5 + uv * a1x2 * a2x5 +
          uv * a2x2 * a1x5;
    h33 = uu * a1x3 * a1x3 + vv * a2x3 * a2x3 + uv * a1x3 * a2x3 +
          uv * a2x3 * a1x3;
    h34 = uu * a1x3 * a1x4 + vv * a2x3 * a2x4 + uv * a1x3 * a2x4 +
          uv * a2x3 * a1x4;
    h35 = uu * a1x3 * a1x5 + vv * a2x3 * a2x5 + uv * a1x3 * a2x5 +
          uv * a2x3 * a1x5;
    h44 = uu * a1x4 * a1x4 + vv * a2x4 * a2x4 + uv * a1x4 * a2x4 +
          uv * a2x4 * a1x4;
    h45 = uu * a1x4 * a1x5 + vv * a2x4 * a2x5 + uv * a1x4 * a2x5 +
          uv * a2x4 * a1x5;
    h55 = uu * a1x5 * a1x5 + vv * a2x5 * a2x5 + uv * a1x5 * a2x5 +
          uv * a2x5 * a1x5;

    if (walpha >= 0) {
      h00 += hwa1 * (ss1 + ss2) + hwa2 * a1x0 * a1x0;
      h01 += hwa1 * (-ss1) + hwa2 * a1x0 * a1x1;
      h02 += hwa1 * (-ss2) + hwa2 * a1x0 * a1x2;
      h03 += hwa2 * a1x0 * a1x3;
      h04 += hwa1 * (h4) + hwa2 * a1x0 * a1x4;
      h05 += hwa1 * (-h4) + hwa2 * a1x0 * a1x5;
      h11 += hwa1 * (h1) + hwa2 * a1x1 * a1x1;
      h12 += hwa1 * (h2) + hwa2 * a1x1 * a1x2;
      h13 += hwa1 * (-h4) + hwa2 * a1x1 * a1x3;
      h14 += hwa2 * a1x1 * a1x4;
      h15 += hwa1 * (h4) + hwa2 * a1x1 * a1x5;
      h22 += hwa1 * (h3) + hwa2 * a1x2 * a1x2;
      h23 += hwa1 * (h4) + hwa2 * a1x2 * a1x3;
      h24 += hwa1 * (-h4) + hwa2 * a1x2 * a1x4;
      h25 += hwa2 * a1x2 * a1x5;
      h33 += hwa1 * (ss1 + ss2) + hwa2 * a1x3 * a1x3;
      h34 += hwa1 * (-ss1) + hwa2 * a1x3 * a1x4;
      h35 += hwa1 * (-ss2) + hwa2 * a1x3 * a1x5;
      h44 += hwa1 * (h1) + hwa2 * a1x4 * a1x4;
      h45 += hwa1 * (h2) + hwa2 * a1x4 * a1x5;
      h55 += hwa1 * (h3) + hwa2 * a1x5 * a1x5;
    }
    h00 += hwb1 * (ss1 + ss2) + hwb2 * a2x0 * a2x0;
    h01 += hwb1 * (-ss1) + hwb2 * a2x0 * a2x1;
    h02 += hwb1 * (-ss2) + hwb2 * a2x0 * a2x2;
    h03 += hwb2 * a2x0 * a2x3;
    h04 += hwb1 * (-h4) + hwb2 * a2x0 * a2x4;
    h05 += hwb1 * (h4) + hwb2 * a2x0 * a2x5;
    h11 += hwb1 * (h1) + hwb2 * a2x1 * a2x1;
    h12 += hwb1 * (h2) + hwb2 * a2x1 * a2x2;
    h13 += hwb1 * (h4) + hwb2 * a2x1 * a2x3;
    h14 += hwb2 * a2x1 * a2x4;
    h15 += hwb1 * (-h4) + hwb2 * a2x1 * a2x5;
    h22 += hwb1 * (h3) + hwb2 * a2x2 * a2x2;
    h23 += hwb1 * (-h4) + hwb2 * a2x2 * a2x3;
    h24 += hwb1 * (h4) + hwb2 * a2x2 * a2x4;
    h25 += hwb2 * a2x2 * a2x5;
    h33 += hwb1 * (ss1 + ss2) + hwb2 * a2x3 * a2x3;
    h34 += hwb1 * (-ss1) + hwb2 * a2x3 * a2x4;
    h35 += hwb1 * (-ss2) + hwb2 * a2x3 * a2x5;
    h44 += hwb1 * (h1) + hwb2 * a2x4 * a2x4;
    h45 += hwb1 * (h2) + hwb2 * a2x4 * a2x5;
    h55 += hwb1 * (h3) + hwb2 * a2x5 * a2x5;

    u = aa * walpha;
    v = bb * wbeta;
    pardiso_b_[f0] -= (u * a1x0 + v * a2x0);
    pardiso_b_[f1] -= (u * a1x1 + v * a2x1);
    pardiso_b_[f2] -= (u * a1x2 + v * a2x2);
    pardiso_b_[f0 + V_N_] -= (u * a1x3 + v * a2x3);
    pardiso_b_[f1 + V_N_] -= (u * a1x4 + v * a2x4);
    pardiso_b_[f2 + V_N_] -= (u * a1x5 + v * a2x5);

    pardiso_a_[id_h00_[i]] += h00;
    pardiso_a_[id_h01_[i]] += h01;
    pardiso_a_[id_h02_[i]] += h02;
    pardiso_a_[id_h03_[i]] += h03;
    pardiso_a_[id_h04_[i]] += h04;
    pardiso_a_[id_h05_[i]] += h05;
    pardiso_a_[id_h11_[i]] += h11;
    pardiso_a_[id_h12_[i]] += h12;
    pardiso_a_[id_h13_[i]] += h13;
    pardiso_a_[id_h14_[i]] += h14;
    pardiso_a_[id_h15_[i]] += h15;
    pardiso_a_[id_h22_[i]] += h22;
    pardiso_a_[id_h23_[i]] += h23;
    pardiso_a_[id_h24_[i]] += h24;
    pardiso_a_[id_h25_[i]] += h25;
    pardiso_a_[id_h33_[i]] += h33;
    pardiso_a_[id_h34_[i]] += h34;
    pardiso_a_[id_h35_[i]] += h35;
    pardiso_a_[id_h44_[i]] += h44;
    pardiso_a_[id_h45_[i]] += h45;
    pardiso_a_[id_h55_[i]] += h55;
  }

  for (int i = shell_data_.mesh_faces_.rows(); i < F_N_; i++) {
    area_now = area_[i];
    f0 = F_[0][i];
    f1 = F_[1][i];
    f2 = F_[2][i];
    x0 = position[f0];
    y0 = position[f0 + V_N_];
    x1 = position[f1];
    y1 = position[f1 + V_N_];
    x2 = position[f2];
    y2 = position[f2 + V_N_];

    q00 = x1 - x0;
    q01 = x2 - x0;
    q10 = y1 - y0;
    q11 = y2 - y0;
    p00 = tmp_p00[i];
    p01 = tmp_p01[i];
    p10 = tmp_p10[i];
    p11 = tmp_p11[i];
    j00 = p00 * q00 + p10 * q01;
    j01 = p01 * q00 + p11 * q01;
    j10 = p00 * q10 + p10 * q11;
    j11 = p01 * q10 + p11 * q11;
    alpha_0 = j00 + j11;
    alpha_1 = j10 - j01;
    beta_0 = j00 - j11;
    beta_1 = j10 + j01;

    alpha_norm = 0.5 * sqrt(alpha_0 * alpha_0 + alpha_1 * alpha_1);
    beta_norm = 0.5 * sqrt(beta_0 * beta_0 + beta_1 * beta_1);
    if (beta_norm < 1e-30) {
      beta_norm = 1e-10;
    }

    ss1 = (p00) * (p00 + p10) + (p01) * (p01 + p11);
    ss2 = (p10) * (p00 + p10) + (p11) * (p01 + p11);
    double h1 = p00 * p00 + p01 * p01;
    double h2 = p00 * p10 + p01 * p11;
    double h3 = p10 * p10 + p11 * p11;
    double h4 = p00 * p11 - p01 * p10;

    a1x0 = alpha_0 * (-p00 - p10) + alpha_1 * (p01 + p11);
    a1x1 = alpha_0 * p00 - alpha_1 * p01;
    a1x2 = alpha_0 * p10 - alpha_1 * p11;
    a1x3 = alpha_0 * (-p01 - p11) + alpha_1 * (-p00 - p10);
    a1x4 = alpha_0 * p01 + alpha_1 * p00;
    a1x5 = alpha_0 * p11 + alpha_1 * p10;
    a2x0 = beta_0 * (-p00 - p10) + beta_1 * (-p01 - p11);
    a2x1 = beta_0 * p00 + beta_1 * p01;
    a2x2 = beta_0 * p10 + beta_1 * p11;
    a2x3 = beta_0 * (p01 + p11) + beta_1 * (-p00 - p10);
    a2x4 = -beta_0 * p01 + beta_1 * p00;
    a2x5 = -beta_0 * p11 + beta_1 * p10;
    sig0 = alpha_norm + beta_norm;
    sig1 = alpha_norm - beta_norm;

    hi_0 = 2 + 6 * 1 / (sig0 * sig0 * sig0 * sig0);
    hi_1 = 2 + 6 * 1 / (sig1 * sig1 * sig1 * sig1);
    aa = 0.25 / alpha_norm;
    bb = 0.25 / beta_norm;
    uu = aa * aa * (area_now * hi_0 + area_now * hi_1);
    vv = bb * bb * (area_now * hi_0 + area_now * hi_1);
    uv = aa * bb * (area_now * hi_0 - area_now * hi_1);
    h_u = area_now * (2 * sig0 - 2 * 1 / (sig0 * sig0 * sig0));
    h_v = area_now * (2 * sig1 - 2 * 1 / (sig1 * sig1 * sig1));
    walpha = h_u + h_v;
    wbeta = h_u - h_v;

    double hwa1 = (walpha * 0.25 / alpha_norm);
    double hwa2 =
        -(walpha * 0.25 * 0.25 / (alpha_norm * alpha_norm * alpha_norm));
    double hwb1 = (wbeta * 0.25 / beta_norm);
    double hwb2 = -(wbeta * 0.25 * 0.25 / (beta_norm * beta_norm * beta_norm));
    h00 = uu * a1x0 * a1x0 + vv * a2x0 * a2x0 + uv * a1x0 * a2x0 +
          uv * a2x0 * a1x0;
    h01 = uu * a1x0 * a1x1 + vv * a2x0 * a2x1 + uv * a1x0 * a2x1 +
          uv * a2x0 * a1x1;
    h02 = uu * a1x0 * a1x2 + vv * a2x0 * a2x2 + uv * a1x0 * a2x2 +
          uv * a2x0 * a1x2;
    h03 = uu * a1x0 * a1x3 + vv * a2x0 * a2x3 + uv * a1x0 * a2x3 +
          uv * a2x0 * a1x3;
    h04 = uu * a1x0 * a1x4 + vv * a2x0 * a2x4 + uv * a1x0 * a2x4 +
          uv * a2x0 * a1x4;
    h05 = uu * a1x0 * a1x5 + vv * a2x0 * a2x5 + uv * a1x0 * a2x5 +
          uv * a2x0 * a1x5;
    h11 = uu * a1x1 * a1x1 + vv * a2x1 * a2x1 + uv * a1x1 * a2x1 +
          uv * a2x1 * a1x1;
    h12 = uu * a1x1 * a1x2 + vv * a2x1 * a2x2 + uv * a1x1 * a2x2 +
          uv * a2x1 * a1x2;
    h13 = uu * a1x1 * a1x3 + vv * a2x1 * a2x3 + uv * a1x1 * a2x3 +
          uv * a2x1 * a1x3;
    h14 = uu * a1x1 * a1x4 + vv * a2x1 * a2x4 + uv * a1x1 * a2x4 +
          uv * a2x1 * a1x4;
    h15 = uu * a1x1 * a1x5 + vv * a2x1 * a2x5 + uv * a1x1 * a2x5 +
          uv * a2x1 * a1x5;
    h22 = uu * a1x2 * a1x2 + vv * a2x2 * a2x2 + uv * a1x2 * a2x2 +
          uv * a2x2 * a1x2;
    h23 = uu * a1x2 * a1x3 + vv * a2x2 * a2x3 + uv * a1x2 * a2x3 +
          uv * a2x2 * a1x3;
    h24 = uu * a1x2 * a1x4 + vv * a2x2 * a2x4 + uv * a1x2 * a2x4 +
          uv * a2x2 * a1x4;
    h25 = uu * a1x2 * a1x5 + vv * a2x2 * a2x5 + uv * a1x2 * a2x5 +
          uv * a2x2 * a1x5;
    h33 = uu * a1x3 * a1x3 + vv * a2x3 * a2x3 + uv * a1x3 * a2x3 +
          uv * a2x3 * a1x3;
    h34 = uu * a1x3 * a1x4 + vv * a2x3 * a2x4 + uv * a1x3 * a2x4 +
          uv * a2x3 * a1x4;
    h35 = uu * a1x3 * a1x5 + vv * a2x3 * a2x5 + uv * a1x3 * a2x5 +
          uv * a2x3 * a1x5;
    h44 = uu * a1x4 * a1x4 + vv * a2x4 * a2x4 + uv * a1x4 * a2x4 +
          uv * a2x4 * a1x4;
    h45 = uu * a1x4 * a1x5 + vv * a2x4 * a2x5 + uv * a1x4 * a2x5 +
          uv * a2x4 * a1x5;
    h55 = uu * a1x5 * a1x5 + vv * a2x5 * a2x5 + uv * a1x5 * a2x5 +
          uv * a2x5 * a1x5;

    if (walpha >= 0) {
      h00 += hwa1 * (ss1 + ss2) + hwa2 * a1x0 * a1x0;
      h01 += hwa1 * (-ss1) + hwa2 * a1x0 * a1x1;
      h02 += hwa1 * (-ss2) + hwa2 * a1x0 * a1x2;
      h03 += hwa2 * a1x0 * a1x3;
      h04 += hwa1 * (h4) + hwa2 * a1x0 * a1x4;
      h05 += hwa1 * (-h4) + hwa2 * a1x0 * a1x5;
      h11 += hwa1 * (h1) + hwa2 * a1x1 * a1x1;
      h12 += hwa1 * (h2) + hwa2 * a1x1 * a1x2;
      h13 += hwa1 * (-h4) + hwa2 * a1x1 * a1x3;
      h14 += hwa2 * a1x1 * a1x4;
      h15 += hwa1 * (h4) + hwa2 * a1x1 * a1x5;
      h22 += hwa1 * (h3) + hwa2 * a1x2 * a1x2;
      h23 += hwa1 * (h4) + hwa2 * a1x2 * a1x3;
      h24 += hwa1 * (-h4) + hwa2 * a1x2 * a1x4;
      h25 += hwa2 * a1x2 * a1x5;
      h33 += hwa1 * (ss1 + ss2) + hwa2 * a1x3 * a1x3;
      h34 += hwa1 * (-ss1) + hwa2 * a1x3 * a1x4;
      h35 += hwa1 * (-ss2) + hwa2 * a1x3 * a1x5;
      h44 += hwa1 * (h1) + hwa2 * a1x4 * a1x4;
      h45 += hwa1 * (h2) + hwa2 * a1x4 * a1x5;
      h55 += hwa1 * (h3) + hwa2 * a1x5 * a1x5;
    }
    h00 += hwb1 * (ss1 + ss2) + hwb2 * a2x0 * a2x0;
    h01 += hwb1 * (-ss1) + hwb2 * a2x0 * a2x1;
    h02 += hwb1 * (-ss2) + hwb2 * a2x0 * a2x2;
    h03 += hwb2 * a2x0 * a2x3;
    h04 += hwb1 * (-h4) + hwb2 * a2x0 * a2x4;
    h05 += hwb1 * (h4) + hwb2 * a2x0 * a2x5;
    h11 += hwb1 * (h1) + hwb2 * a2x1 * a2x1;
    h12 += hwb1 * (h2) + hwb2 * a2x1 * a2x2;
    h13 += hwb1 * (h4) + hwb2 * a2x1 * a2x3;
    h14 += hwb2 * a2x1 * a2x4;
    h15 += hwb1 * (-h4) + hwb2 * a2x1 * a2x5;
    h22 += hwb1 * (h3) + hwb2 * a2x2 * a2x2;
    h23 += hwb1 * (-h4) + hwb2 * a2x2 * a2x3;
    h24 += hwb1 * (h4) + hwb2 * a2x2 * a2x4;
    h25 += hwb2 * a2x2 * a2x5;
    h33 += hwb1 * (ss1 + ss2) + hwb2 * a2x3 * a2x3;
    h34 += hwb1 * (-ss1) + hwb2 * a2x3 * a2x4;
    h35 += hwb1 * (-ss2) + hwb2 * a2x3 * a2x5;
    h44 += hwb1 * (h1) + hwb2 * a2x4 * a2x4;
    h45 += hwb1 * (h2) + hwb2 * a2x4 * a2x5;
    h55 += hwb1 * (h3) + hwb2 * a2x5 * a2x5;

    u = aa * walpha;
    v = bb * wbeta;
    pardiso_b_[f0] -= (u * a1x0 + v * a2x0);
    pardiso_b_[f1] -= (u * a1x1 + v * a2x1);
    pardiso_b_[f2] -= (u * a1x2 + v * a2x2);
    pardiso_b_[f0 + V_N_] -= (u * a1x3 + v * a2x3);
    pardiso_b_[f1 + V_N_] -= (u * a1x4 + v * a2x4);
    pardiso_b_[f2 + V_N_] -= (u * a1x5 + v * a2x5);

    pardiso_a_[id_h00_[i]] += h00;
    pardiso_a_[id_h01_[i]] += h01;
    pardiso_a_[id_h02_[i]] += h02;
    pardiso_a_[id_h03_[i]] += h03;
    pardiso_a_[id_h04_[i]] += h04;
    pardiso_a_[id_h05_[i]] += h05;
    pardiso_a_[id_h11_[i]] += h11;
    pardiso_a_[id_h12_[i]] += h12;
    pardiso_a_[id_h13_[i]] += h13;
    pardiso_a_[id_h14_[i]] += h14;
    pardiso_a_[id_h15_[i]] += h15;
    pardiso_a_[id_h22_[i]] += h22;
    pardiso_a_[id_h23_[i]] += h23;
    pardiso_a_[id_h24_[i]] += h24;
    pardiso_a_[id_h25_[i]] += h25;
    pardiso_a_[id_h33_[i]] += h33;
    pardiso_a_[id_h34_[i]] += h34;
    pardiso_a_[id_h35_[i]] += h35;
    pardiso_a_[id_h44_[i]] += h44;
    pardiso_a_[id_h45_[i]] += h45;
    pardiso_a_[id_h55_[i]] += h55;
  }

  double len01, len0i, len1i, coef1, coef2, dis;
  int i;
  double s0, s1, e0, e1, p0, p1;
  AV_F_N_H_ = AV_F_N_;
  V_F0_H_ = V_F0_;
  V_F1_H_ = V_F1_;
  V_F2_H_ = V_F2_;
  for (int ii = 0; ii < AV_F_N_; ++ii) {
    f0 = V_F0_[AV_ID_[ii]];
    f1 = V_F1_[AV_ID_[ii]];
    f2 = V_F2_[AV_ID_[ii]];
    x0 = position[f0];
    y0 = position[f0 + V_N_];
    x1 = position[f1];
    y1 = position[f1 + V_N_];
    x2 = position[f2];
    y2 = position[f2 + V_N_];

    s0 = x0;
    s1 = y0;
    e0 = x1;
    e1 = y1;
    p0 = x2;
    p1 = y2;
    len01 = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    len0i = sqrt((p0 - s0) * (p0 - s0) + (p1 - s1) * (p1 - s1));
    len1i = sqrt((p0 - e0) * (p0 - e0) + (p1 - e1) * (p1 - e1));
    dis = len0i + len1i - len01;

    h_u = 2 * barrer_coef_ * threhold * (dis - threhold) / dis / dis / dis;
    uu = 2 * barrer_coef_ * threhold * (3 * threhold - 2 * dis) / dis / dis /
         dis / dis;
    a1x0 = (x0 - x2) / len0i - (x0 - x1) / len01;
    a1x1 = (x1 - x2) / len1i - (x1 - x0) / len01;
    a1x2 = (x2 - x0) / len0i + (x2 - x1) / len1i;
    a1x3 = (y0 - y2) / len0i - (y0 - y1) / len01;
    a1x4 = (y1 - y2) / len1i - (y1 - y0) / len01;
    a1x5 = (y2 - y0) / len0i + (y2 - y1) / len1i;

    coef1 = len01 * len01;
    coef2 = h_u / coef1 / len01;
    aa = (x0 - x1) * (x0 - x1);
    bb = (x1 - x0) * (y1 - y0);
    cc = (y0 - y1) * (y0 - y1);
    h00 = coef2 * (aa - coef1);
    h01 = coef2 * (coef1 - aa);
    h02 = 0.0;
    h03 = coef2 * bb;
    h04 = -coef2 * bb;
    h05 = 0.0;
    h11 = coef2 * (aa - coef1);
    h12 = 0.0;
    h13 = -coef2 * bb;
    h14 = coef2 * bb;
    h15 = 0.0;
    h22 = 0.0;
    h23 = 0.0;
    h24 = 0.0;
    h25 = 0.0;
    h33 = coef2 * (cc - coef1);
    h34 = coef2 * (coef1 - cc);
    h35 = 0.0;
    h44 = coef2 * (cc - coef1);
    h45 = 0.0;
    h55 = 0.0;

    h00 += uu * a1x0 * a1x0;
    h01 += uu * a1x0 * a1x1;
    h02 += uu * a1x0 * a1x2;
    h03 += uu * a1x0 * a1x3;
    h04 += uu * a1x0 * a1x4;
    h05 += uu * a1x0 * a1x5;
    h11 += uu * a1x1 * a1x1;
    h12 += uu * a1x1 * a1x2;
    h13 += uu * a1x1 * a1x3;
    h14 += uu * a1x1 * a1x4;
    h15 += uu * a1x1 * a1x5;
    h22 += uu * a1x2 * a1x2;
    h23 += uu * a1x2 * a1x3;
    h24 += uu * a1x2 * a1x4;
    h25 += uu * a1x2 * a1x5;
    h33 += uu * a1x3 * a1x3;
    h34 += uu * a1x3 * a1x4;
    h35 += uu * a1x3 * a1x5;
    h44 += uu * a1x4 * a1x4;
    h45 += uu * a1x4 * a1x5;
    h55 += uu * a1x5 * a1x5;

    pardiso_b_[f0] -= h_u * (a1x0);
    pardiso_b_[f1] -= h_u * (a1x1);
    pardiso_b_[f2] -= h_u * (a1x2);
    pardiso_b_[f0 + V_N_] -= h_u * (a1x3);
    pardiso_b_[f1 + V_N_] -= h_u * (a1x4);
    pardiso_b_[f2 + V_N_] -= h_u * (a1x5);

    i = F_N_ + AV_ID_[ii];
    pardiso_a_[id_h00_[i]] += h00;
    pardiso_a_[id_h01_[i]] += h01;
    pardiso_a_[id_h02_[i]] += h02;
    pardiso_a_[id_h03_[i]] += h03;
    pardiso_a_[id_h04_[i]] += h04;
    pardiso_a_[id_h05_[i]] += h05;
    pardiso_a_[id_h11_[i]] += h11;
    pardiso_a_[id_h12_[i]] += h12;
    pardiso_a_[id_h13_[i]] += h13;
    pardiso_a_[id_h14_[i]] += h14;
    pardiso_a_[id_h15_[i]] += h15;
    pardiso_a_[id_h22_[i]] += h22;
    pardiso_a_[id_h23_[i]] += h23;
    pardiso_a_[id_h24_[i]] += h24;
    pardiso_a_[id_h25_[i]] += h25;
    pardiso_a_[id_h33_[i]] += h33;
    pardiso_a_[id_h34_[i]] += h34;
    pardiso_a_[id_h35_[i]] += h35;
    pardiso_a_[id_h44_[i]] += h44;
    pardiso_a_[id_h45_[i]] += h45;
    pardiso_a_[id_h55_[i]] += h55;
  }

  pardiso_a_[0] += 1.0;
  pardiso_a_[pardiso_ia_[V_N_]] += 1.0;
  solver_->a_ = pardiso_a_;
  solver_->rhs_ = pardiso_b_;
  long time_beg, time_end;
  time_beg = clock();
  solver_->Factorize();
  time_end = clock();
  double time_consumption = (time_end - time_beg) / 1000.0;
  time_2_ += time_consumption;
  time_beg = clock();
  solver_->PardisoSolver();
  time_end = clock();
  time_consumption = (time_end - time_beg) / 1000.0;
  time_3_ += time_consumption;
  // std::cout << "numerical factorize" << time_consumption << std::endl;

  std::vector<double> result_d = solver_->result_;
  Eigen::VectorXd negative_grad(2 * V_N_), d(2 * V_N_);
  for (int i = 0; i < 2 * V_N_; i++) {
    negative_grad(i) = pardiso_b_[i];
    d(i) = result_d[i];
  }
  solver_->FreeNumericalFactorizationMemory();

  double temp_t;
  MaxStep(position_of_mesh_, d, temp_t);
  double alpha = 0.95 * temp_t;
  // double alpha = 0.8*temp_t;
  BacktrackingLineSearch(position_of_mesh_, d, negative_grad, alpha);
  position_of_mesh_ += alpha * d;
  // std::cout <<is_interp<< " cm step: " << alpha << std::endl;
  EnergySource();
}

void ParaFun::MaxStep(const Eigen::VectorXd &xx, const Eigen::VectorXd &qq,
                      double &step) {
  double temp_t = std::numeric_limits<double>::infinity();
  int f0, f1, f2;
  double a, b, c, tt, tt1, tt2;
  double x0, x1, x2, y0, y1, y2, u0, u1, u2, v0, v1, v2;
  const double *x = xx.data();

  for (int i = 0; i < F_N_; ++i) {
    f0 = F_[0][i];
    f1 = F_[1][i];
    f2 = F_[2][i];

    x0 = x[f0];
    x1 = x[f1];
    x2 = x[f2];
    y0 = x[f0 + V_N_];
    y1 = x[f1 + V_N_];
    y2 = x[f2 + V_N_];
    u0 = qq[f0];
    u1 = qq[f1];
    u2 = qq[f2];
    v0 = qq[f0 + V_N_];
    v1 = qq[f1 + V_N_];
    v2 = qq[f2 + V_N_];

    a = (u1 - u0) * (v2 - v0) - (v1 - v0) * (u2 - u0);
    b = (u1 - u0) * (y2 - y0) + (x1 - x0) * (v2 - v0) - (y1 - y0) * (u2 - u0) -
        (x2 - x0) * (v1 - v0);
    c = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
    tt = GetSmallestPosQuadZero(a, b, c);

    if (temp_t > tt) {
      temp_t = tt;
    }
  }

  temp_t = 0.95 * temp_t;
  const double *dir = qq.data();
  DetectTmax(xx, qq, temp_t);
  double a2, b2, c2, a3, b3, c3;
  double root3[3];
  double d, q, r, r2, q3, A, B, tt3, tt4;
  double q1, q2, p1, p2, D, sqD, y, m2, m, n, alpha, beta, gamma, delta;
  double x0_up, x1_up, x2_up, y0_up, y1_up, y2_up;

  for (int i = 0; i < AV_F_N_; ++i) {
    is_active_[AV_ID_[i]] = -1;

    f0 = V_F0_[AV_ID_[i]];
    f1 = V_F1_[AV_ID_[i]];
    f2 = V_F2_[AV_ID_[i]];
    x0 = x[f0];
    x1 = x[f1];
    x2 = x[f2];
    y0 = x[f0 + V_N_];
    y1 = x[f1 + V_N_];
    y2 = x[f2 + V_N_];
    u0 = dir[f0];
    u1 = dir[f1];
    u2 = dir[f2];
    v0 = dir[f0 + V_N_];
    v1 = dir[f1 + V_N_];
    v2 = dir[f2 + V_N_];

    // no triangle flip
    a = (u1 - u0) * (v2 - v0) - (v1 - v0) * (u2 - u0);
    b2 = (y1 - y0) * (u2 - u0) + (x2 - x0) * (v1 - v0);
    b = (u1 - u0) * (y2 - y0) + (x1 - x0) * (v2 - v0) - b2;
    c = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
    tt = std::numeric_limits<double>::infinity();
    if (b * b - 4 * a * c >= 0) {
      tt1 = 1 / (2 * a) * (-b + sqrt(b * b - 4 * a * c));
      tt2 = 1 / (2 * a) * (-b - sqrt(b * b - 4 * a * c));
      if (tt1 > 0 && tt2 > 0) {
        tt = std::min(tt1, tt2);
      }
      if (tt1 > 0 && tt2 < 0) {
        tt = tt1;
      }
      if (tt1 < 0 && tt2 > 0) {
        tt = tt2;
      }
    }

    if (temp_t > tt) {
      x0_up = x0 + tt * u0;
      x1_up = x1 + tt * u1;
      x2_up = x2 + tt * u2;
      y0_up = y0 + tt * v0;
      y1_up = y1 + tt * v1;
      y2_up = y2 + tt * v2;

      double dot1 =
          (x1_up - x0_up) * (x2_up - x0_up) + (y1_up - y0_up) * (y2_up - y0_up);
      double dot2 =
          (x0_up - x1_up) * (x2_up - x1_up) + (y0_up - y1_up) * (y2_up - y1_up);

      if (dot1 > 0 && dot2 > 0) {
        temp_t = tt;
      }
    }
  }

  // std::cout << step << std::endl;

  bool check_isintersection = false;
  do {
    check_isintersection = false;
    check_isintersection = CheckIntersection(xx + temp_t * qq);
    if (check_isintersection) {
      temp_t = 0.8 * temp_t;
    }
  } while (check_isintersection);

  step = temp_t;
}
void ParaFun::DetectTmax(const Eigen::VectorXd &x, const Eigen::VectorXd &d,
                         double &tmax) {
  const double *pos = x.data();
  std::vector<double> x_update(2 * V_N_);
  for (int i = 0; i < 2 * V_N_; ++i) {
    x_update[i] = x[i] + tmax * d[i];
  }
  x_min = pos[0];
  x_max = pos[0];
  y_min = pos[V_N_];
  y_max = pos[V_N_];

  double x_min_up = x_update[0];
  double x_max_up = x_update[0];
  double y_min_up = x_update[V_N_];
  double y_max_up = x_update[V_N_];

  for (int i = 1; i < V_N_; ++i) {
    if (pos[i] < x_min) {
      x_min = pos[i];
    } else if (pos[i] > x_max) {
      x_max = pos[i];
    }

    if (pos[i + V_N_] < y_min) {
      y_min = pos[i + V_N_];
    } else if (pos[i + V_N_] > y_max) {
      y_max = pos[i + V_N_];
    }

    if (x_update[i] < x_min_up) {
      x_min_up = x_update[i];
    } else if (x_update[i] > x_max_up) {
      x_max_up = x_update[i];
    }

    if (x_update[i + V_N_] < y_min_up) {
      y_min_up = x_update[i + V_N_];
    } else if (x_update[i + V_N_] > y_max_up) {
      y_max_up = x_update[i + V_N_];
    }
  }
  x_min = std::min(x_min, x_min_up);
  x_max = std::max(x_max, x_max_up);
  y_min = std::min(y_min, y_min_up);
  y_max = std::max(y_max, y_max_up);

  lengthgrid_x = (x_max - x_min) / (cellx_num - 1);
  x_max = x_min + cellx_num * lengthgrid_x;
  lengthgrid_y = (y_max - y_min) / (celly_num - 1);
  y_max = y_min + celly_num * lengthgrid_y;

  for (int j = 0; j < cell_points.size(); ++j) {
    cell_points[j].clear();
  }

  Eigen::VectorXi boundary_vertex = shell_data_.frame_ids_;
  int id, id_x_min, id_x_max, id_y_min, id_y_max;
  double s0, s1, e0, e1, p0, p1, s0_, s1_, e0_, e1_;
  double l_x_min, l_x_max, l_y_min, l_y_max;

  AV_ID_.clear();
  AV_F_N_ = 0;
  int id_start, id_end;
  Eigen::Vector4d localx, localy;
  double a, len, b;

  for (int i = 0; i < BE_N_; ++i) {
    id = boundary_vertex(i);
    s0 = pos[id];
    s1 = pos[id + V_N_];
    e0 = x_update[id];
    e1 = x_update[id + V_N_];

    l_x_min = std::min(s0, e0);
    l_x_max = std::max(s0, e0);
    l_y_min = std::min(s1, e1);
    l_y_max = std::max(s1, e1);

    id_x_min = std::floor((l_x_min - x_min) / lengthgrid_x);
    id_x_max = std::floor((l_x_max - x_min) / lengthgrid_x);
    id_y_min = std::floor((l_y_min - y_min) / lengthgrid_y);
    id_y_max = std::floor((l_y_max - y_min) / lengthgrid_y);
    if (id_y_max > celly_num - 1) {
      id_y_max = celly_num - 1;
    }
    if (id_x_max > cellx_num - 1) {
      id_x_max = cellx_num - 1;
    }

    for (int j_y = id_y_min; j_y < id_y_max + 1; ++j_y) {
      for (int j_x = id_x_min; j_x < id_x_max + 1; ++j_x) {
        int id_grid = j_x + j_y * cellx_num;
        cell_points[id_grid].push_back(id);
      }
    }
  }

  for (int i = 0; i < BE_N_; ++i) {
    id_start = boundary_vertex(i);
    id_end = boundary_vertex((i + 1) % BE_N_);
    s0 = pos[id_start];
    s1 = pos[id_start + V_N_];
    e0 = pos[id_end];
    e1 = pos[id_end + V_N_];
    s0_ = x_update[id_start];
    s1_ = x_update[id_start + V_N_];
    e0_ = x_update[id_end];
    e1_ = x_update[id_end + V_N_];

    localx[0] = s0;
    localx[1] = e0;
    localx[2] = s0_;
    localx[3] = e0_;
    localy[0] = s1;
    localy[1] = e1;
    localy[2] = s1_;
    localy[3] = e1_;
    l_x_min = localx.minCoeff();
    l_x_max = localx.maxCoeff();
    l_y_min = localy.minCoeff();
    l_y_max = localy.maxCoeff();

    len = sqrt((l_x_max - l_x_min) * (l_x_max - l_x_min) +
               (l_y_max - l_y_min) * (l_y_max - l_y_min));

    id_x_min = std::floor((l_x_min - x_min) / lengthgrid_x);
    id_x_max = std::floor((l_x_max - x_min) / lengthgrid_x);
    id_y_min = std::floor((l_y_min - y_min) / lengthgrid_y);
    id_y_max = std::floor((l_y_max - y_min) / lengthgrid_y);

    if (id_y_max > celly_num - 1) {
      id_y_max = celly_num - 1;
    }
    if (id_x_max > cellx_num - 1) {
      id_x_max = cellx_num - 1;
    }
    if (id_x_min < 0) {
      id_x_min = 0;
    }
    if (id_y_min < 0) {
      id_y_min = 0;
    }

    for (int j_y = id_y_min; j_y < id_y_max + 1; ++j_y) {
      for (int j_x = id_x_min; j_x < id_x_max + 1; ++j_x) {
        int id_grid = j_x + j_y * cellx_num;
        int size = cell_points[id_grid].size();
        for (int k = 0; k < size; ++k) {
          int id_mid = cell_points[id_grid][k];
          if (id_start == id_mid || id_end == id_mid) {
            continue;
          }

          if (BE_N_ - 1 == i) {
            if (is_active_[i * (BE_N_ - 2) + boundary_vertexID_[id_mid] - 1] <
                0) {
              is_active_[i * (BE_N_ - 2) + boundary_vertexID_[id_mid] - 1] = 1;
              AV_ID_.push_back(i * (BE_N_ - 2) + boundary_vertexID_[id_mid] -
                               1);
            }
          } else if (boundary_vertexID_[id_mid] <
                     boundary_vertexID_[id_start]) {
            if (is_active_[i * (BE_N_ - 2) + boundary_vertexID_[id_mid]] < 0) {
              is_active_[i * (BE_N_ - 2) + boundary_vertexID_[id_mid]] = 1;
              AV_ID_.push_back(i * (BE_N_ - 2) + boundary_vertexID_[id_mid]);
            }

          } else {
            if (is_active_[i * (BE_N_ - 2) + boundary_vertexID_[id_mid] - 2] <
                0) {
              is_active_[i * (BE_N_ - 2) + boundary_vertexID_[id_mid] - 2] = 1;
              AV_ID_.push_back(i * (BE_N_ - 2) + boundary_vertexID_[id_mid] -
                               2);
            }
          }
        }
      }
    }
  }

  AV_F_N_ = AV_ID_.size();
  //	std::cout << AV_F_N_ << std::endl;
}
double ParaFun::GetSmallestPosQuadZero(double a, double b, double c) {
  using namespace std;
  double t1, t2;
  if (std::abs(a) < 1.0e-10) {
    a *= 1e6;
    b *= 1e6;
    c *= 1e6;
  }
  if (std::abs(a) > 1.0e-10) {
    double delta_in = pow(b, 2) - 4 * a * c;
    if (delta_in <= 0) {
      return INFINITY;
    }

    double delta = sqrt(delta_in); // delta >= 0
    if (b >= 0)                    // avoid subtracting two similar numbers
    {
      double bd = -b - delta;
      t1 = 2 * c / bd;
      t2 = bd / (2 * a);
    } else {
      double bd = -b + delta;
      t1 = bd / (2 * a);
      t2 = (2 * c) / bd;
    }

    assert(std::isfinite(t1));
    assert(std::isfinite(t2));

    if (a < 0)
      std::swap(t1, t2); // make t1 > t2
                         // return the smaller positive root if it exists,
                         // otherwise return infinity
    if (t1 > 0) {
      return t2 > 0 ? t2 : t1;
    } else {
      return INFINITY;
    }
  } else {
    if (b == 0)
      return INFINITY; // just to avoid divide-by-zero
    t1 = -c / b;
    return t1 > 0 ? t1 : INFINITY;
  }
}
bool ParaFun::CheckIntersection(const Eigen::VectorXd &pos) {
  Eigen::VectorXi boundary_vertex = shell_data_.frame_ids_;
  BE_N_ = boundary_vertex.size();

  int id_start, id_end, id_before, id_mid1, id_mid2;
  for (int i = 0; i < BE_N_; ++i) {
    id_start = boundary_vertex(i);
    id_end = boundary_vertex((i + 1) % BE_N_);
    id_before = boundary_vertex((i - 1 + BE_N_) % BE_N_);

    for (int j = 0; j < BE_N_; ++j) {
      if (id_before == boundary_vertex(j))
        continue;
      if (id_start == boundary_vertex(j))
        continue;
      if (id_end == boundary_vertex(j))
        continue;

      id_mid1 = boundary_vertex(j);
      id_mid2 = boundary_vertex((j + 1) % BE_N_);
      double x0 = pos(id_start);
      double y0 = pos(id_start + V_N_);
      double x1 = pos(id_end);
      double y1 = pos(id_end + V_N_);

      double x2 = pos(id_mid1);
      double y2 = pos(id_mid1 + V_N_);
      double x3 = pos(id_mid2);
      double y3 = pos(id_mid2 + V_N_);

      double cross1 = (y2 - y1) * (x2 - x0) - (y2 - y0) * (x2 - x1);
      double cross2 = (y3 - y1) * (x3 - x0) - (y3 - y0) * (x3 - x1);
      double cross3 = (y0 - y3) * (x0 - x2) - (y0 - y2) * (x0 - x3);
      double cross4 = (y1 - y3) * (x1 - x2) - (y1 - y2) * (x1 - x3);
      if (cross1 * cross2 < 0 && cross3 * cross4 < 0) {
        return true;
      }
    }
  }
  return false;
}

void ParaFun::BacktrackingLineSearch(const Eigen::VectorXd &x,
                                     const Eigen::VectorXd &d,
                                     const Eigen::VectorXd &negetive_grad,
                                     double &alpha, bool is_interp) {
  double h = 0.5;
  double tt = -(negetive_grad.transpose() * d)(0, 0);
  double c = 0.2;
  double ex;
  Energy(x, ex, is_interp, false);
  ex = ex + barrer_coef_ * energy_barrier_;
  double e;
  Eigen::VectorXd x_new = x + alpha * d;
  Energy(x_new, e, is_interp);
  while (e > ex + alpha * c * tt) {
    alpha = h * alpha;
    x_new = x + alpha * d;
    Energy(x_new, e, is_interp);
  }
}
void ParaFun::Energy(const Eigen::VectorXd &position, double &energyupdate,
                     bool is_interp, bool is_whole) {
  double energy = 0;

  int f0, f1, f2;
  double x0, y0, x1, y1, x2, y2;
  double det, E_d, max_E_d = 0;
  double j00, j01, j10, j11;
  double p00, p01, p10, p11;
  double q00, q01, q10, q11;
  const double *pos = position.data();

  double *tmp_p00;
  double *tmp_p01;
  double *tmp_p10;
  double *tmp_p11;

  if (is_interp) {
    tmp_p00 = update_p00_.data();
    tmp_p01 = update_p01_.data();
    tmp_p10 = update_p10_.data();
    tmp_p11 = update_p11_.data();
  } else {
    tmp_p00 = source_p00_.data();
    tmp_p01 = source_p01_.data();
    tmp_p10 = source_p10_.data();
    tmp_p11 = source_p11_.data();
  }

  for (int i = 0; i < F_N_; ++i) {
    f0 = F_[0][i];
    f1 = F_[1][i];
    f2 = F_[2][i];

    x0 = pos[f0];
    y0 = pos[f0 + total_num_];

    x1 = pos[f1];
    y1 = pos[f1 + total_num_];

    x2 = pos[f2];
    y2 = pos[f2 + total_num_];

    q00 = x1 - x0;
    q01 = x2 - x0;
    q10 = y1 - y0;
    q11 = y2 - y0;

    p00 = tmp_p00[i];
    p01 = tmp_p01[i];
    p10 = tmp_p10[i];
    p11 = tmp_p11[i];

    j00 = p00 * q00 + p10 * q01;
    j01 = p01 * q00 + p11 * q01;
    j10 = p00 * q10 + p10 * q11;
    j11 = p01 * q10 + p11 * q11;

    det = j00 * j11 - j01 * j10;
    E_d =
        (1 + 1 / (det * det)) * (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);

    energy += area_[i] * E_d;
  }

  energyupdate = energy;

  if (is_whole) {
    FunGrid(position);
    double dis, E_b, energy2 = 0;
    for (int i = 0; i < AV_F_N_; ++i) {
      f0 = V_F0_[AV_ID_[i]];
      f1 = V_F1_[AV_ID_[i]];
      f2 = V_F2_[AV_ID_[i]];

      x0 = pos[f0];
      y0 = pos[f0 + V_N_];

      x1 = pos[f1];
      y1 = pos[f1 + V_N_];

      x2 = pos[f2];
      y2 = pos[f2 + V_N_];

      dis = GetDistance(x0, y0, x1, y1, x2, y2);
      if (dis < 0) {
        std::cout << "distance is zero" << std::endl;
      }
      E_b = (1 - threhold / dis) * (1 - threhold / dis);
      energy2 += E_b;
    }
    energyupdate = energy + barrer_coef_ * energy2;
    energy_barrier_ = energy2;
    // std::cout << AV_F_N_ << std::endl;
  }
}

void ParaFun::EnergySource() {
  double end_e_area = 0;

  int f0, f1, f2;
  double x0, y0, x1, y1, x2, y2;
  double det, E_1, E_2, max_E_d = 0;

  double j00, j01, j10, j11;
  double p00, p01, p10, p11;
  double q00, q01, q10, q11;

  const double *pos = position_of_mesh_.data();
  for (int i = 0; i < shell_data_.mesh_faces_.rows(); ++i) {
    f0 = F_[0][i];
    f1 = F_[1][i];
    f2 = F_[2][i];

    x0 = pos[f0];
    y0 = pos[f0 + total_num_];

    x1 = pos[f1];
    y1 = pos[f1 + total_num_];

    x2 = pos[f2];
    y2 = pos[f2 + total_num_];

    q00 = x1 - x0;
    q01 = x2 - x0;
    q10 = y1 - y0;
    q11 = y2 - y0;

    p00 = source_p00_[i];
    p01 = source_p01_[i];
    p10 = source_p10_[i];
    p11 = source_p11_[i];

    j00 = p00 * q00 + p10 * q01;
    j01 = p01 * q00 + p11 * q01;
    j10 = p00 * q10 + p10 * q11;
    j11 = p01 * q10 + p11 * q11;

    det = j00 * j11 - j01 * j10;

    E_1 = (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);
    E_2 = 1.0 / (det * det) * E_1;

    /*if ((E_1 + E_2) > max_E_d)
    {
            max_E_d = E_1 + E_2;
    }*/

    end_e_area += ((E_1 + E_2) * area_[i]);
  }
  // std::cout << max_E_d << std::endl;
  energy_mesh_ = end_e_area;

  end_e_area = 0;
  for (int i = shell_data_.mesh_faces_.rows(); i < F_N_; ++i) {
    f0 = F_[0][i];
    f1 = F_[1][i];
    f2 = F_[2][i];

    x0 = pos[f0];
    y0 = pos[f0 + total_num_];

    x1 = pos[f1];
    y1 = pos[f1 + total_num_];

    x2 = pos[f2];
    y2 = pos[f2 + total_num_];

    q00 = x1 - x0;
    q01 = x2 - x0;
    q10 = y1 - y0;
    q11 = y2 - y0;

    p00 = source_p00_[i];
    p01 = source_p01_[i];
    p10 = source_p10_[i];
    p11 = source_p11_[i];

    j00 = p00 * q00 + p10 * q01;
    j01 = p01 * q00 + p11 * q01;
    j10 = p00 * q10 + p10 * q11;
    j11 = p01 * q10 + p11 * q11;

    det = j00 * j11 - j01 * j10;

    E_1 = (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);
    E_2 = 1.0 / (det * det) * E_1;

    end_e_area += (E_1 + E_2);
  }

  energy_shell_ = end_e_area;
}
double ParaFun::ComputeEnergy(const Eigen::MatrixXd &x, bool whole) {
  double end_e_one_temp = 0, end_e_area = 0;

  int f0, f1, f2;
  double x0, y0, x1, y1, x2, y2;
  double det, E_1, E_2;

  double j00, j01, j10, j11;
  double p00, p01, p10, p11;
  double q00, q01, q10, q11;

  const double *pos = x.data();
  int src_t_num = shell_data_.mesh_faces_.rows();

  for (int i = 0; i < src_t_num; ++i) {
    f0 = F_[0][i];
    f1 = F_[1][i];
    f2 = F_[2][i];

    x0 = pos[f0];
    y0 = pos[f0 + total_num_];

    x1 = pos[f1];
    y1 = pos[f1 + total_num_];

    x2 = pos[f2];
    y2 = pos[f2 + total_num_];

    q00 = x1 - x0;
    q01 = x2 - x0;
    q10 = y1 - y0;
    q11 = y2 - y0;

    p00 = source_p00_[i];
    p01 = source_p01_[i];
    p10 = source_p10_[i];
    p11 = source_p11_[i];

    j00 = p00 * q00 + p10 * q01;
    j01 = p01 * q00 + p11 * q01;
    j10 = p00 * q10 + p10 * q11;
    j11 = p01 * q10 + p11 * q11;

    det = j00 * j11 - j01 * j10;

    E_1 = (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);
    E_2 = 1.0 / (det * det) * E_1;

    end_e_one_temp += E_1;
    end_e_one_temp += E_2;
    end_e_area += ((E_1 + E_2) * area_src_[i]);
  }

  if (whole) {
    for (int i = src_t_num; i < F_N_; ++i) {
      f0 = F_[0][i];
      f1 = F_[1][i];
      f2 = F_[2][i];

      x0 = pos[f0];
      y0 = pos[f0 + total_num_];

      x1 = pos[f1];
      y1 = pos[f1 + total_num_];

      x2 = pos[f2];
      y2 = pos[f2 + total_num_];

      q00 = x1 - x0;
      q01 = x2 - x0;
      q10 = y1 - y0;
      q11 = y2 - y0;

      p00 = source_p00_[i];
      p01 = source_p01_[i];
      p10 = source_p10_[i];
      p11 = source_p11_[i];

      j00 = p00 * q00 + p10 * q01;
      j01 = p01 * q00 + p11 * q01;
      j10 = p00 * q10 + p10 * q11;
      j11 = p01 * q10 + p11 * q11;

      det = j00 * j11 - j01 * j10;

      E_1 = (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);
      E_2 = 1.0 / (det * det) * E_1;

      end_e_one_temp += E_1;
      end_e_one_temp += E_2;
      end_e_area += ((E_1 + E_2) * area_[i]);
    }

    end_e_area += barrer_coef_ * energy_barrier_;
  }
  // cout << "compute energy with scaf "<<whole<<" "<< end_e_area << endl;
  // cout << "compute energy with scaf/area " << whole << " " <<
  // end_e_area/shell_data_.mesh_measure_ << endl;

  return end_e_area;
}

void ParaFun::LocalCoordinateInverse(int i, double &p00, double &p01,
                                     double &p10, double &p11) {
  /*int f0 = F_[0][i];
  int f1 = F_[1][i];
  int f2 = F_[2][i];

  Vector3d x_(shell_data_.mesh_vertices_(f1, 0) - shell_data_.mesh_vertices_(f0,
  0), shell_data_.mesh_vertices_(f1, 1) - shell_data_.mesh_vertices_(f0, 1),
  shell_data_.mesh_vertices_(f1, 2)
  - shell_data_.mesh_vertices_(f0, 2)); double x1_0 = x_.norm(); x_ /= x1_0;
  Vector3d l_(shell_data_.mesh_vertices_(f2, 0) - shell_data_.mesh_vertices_(f0,
  0), shell_data_.mesh_vertices_(f2, 1) - shell_data_.mesh_vertices_(f0, 1),
  shell_data_.mesh_vertices_(f2, 2)
  - shell_data_.mesh_vertices_(f0, 2));

  Vector3d n_ = x_.cross(l_);
  n_.normalize();
  Vector3d y_ = n_.cross(x_);
  double x2_0 = l_.dot(x_);
  double y2_0 = l_.dot(y_);

  p00 = 1 / x1_0;
  p01 = -x2_0 / (x1_0*y2_0);
  p10 = 0;
  p11 = 1 / y2_0;*/

  int f0 = F_[0][i];
  int f1 = F_[1][i];
  int f2 = F_[2][i];
  double area_tri = area_[i];
  double area_min = 1e-15;
  if (area_tri > area_min) {
    Eigen::Vector3d t_(shell_data_.mesh_vertices_(f0, 0),
                       shell_data_.mesh_vertices_(f0, 1),
                       shell_data_.mesh_vertices_(f0, 2));
    Eigen::Vector3d u_(shell_data_.mesh_vertices_(f1, 0),
                       shell_data_.mesh_vertices_(f1, 1),
                       shell_data_.mesh_vertices_(f1, 2));
    Eigen::Vector3d v_(shell_data_.mesh_vertices_(f2, 0),
                       shell_data_.mesh_vertices_(f2, 1),
                       shell_data_.mesh_vertices_(f2, 2));
    Eigen::Vector3d x_(
        shell_data_.mesh_vertices_(f1, 0) - shell_data_.mesh_vertices_(f0, 0),
        shell_data_.mesh_vertices_(f1, 1) - shell_data_.mesh_vertices_(f0, 1),
        shell_data_.mesh_vertices_(f1, 2) - shell_data_.mesh_vertices_(f0, 2));
    double x1_0 = x_.norm();
    x_ /= x1_0;
    Eigen::Vector3d l_(
        shell_data_.mesh_vertices_(f2, 0) - shell_data_.mesh_vertices_(f0, 0),
        shell_data_.mesh_vertices_(f2, 1) - shell_data_.mesh_vertices_(f0, 1),
        shell_data_.mesh_vertices_(f2, 2) - shell_data_.mesh_vertices_(f0, 2));
    Eigen::Vector3d n_ = x_.cross(l_);
    n_.normalize();
    Eigen::Vector3d y_ = n_.cross(x_);
    double x2_0 = l_.dot(x_);
    double y2_0 = l_.dot(y_);
    p00 = 1 / x1_0;
    p01 = -x2_0 / (x1_0 * y2_0);
    p10 = 0;
    p11 = 1 / y2_0;
  } else {
    // cout << "area too small!!!!!!!!!!!!! " << endl;
    double h = sqrt((2 * area_min) / sqrt(3.0));
    double x1_0 = h;
    double x2_0 = h / 2.0;
    double y2_0 = sqrt(3.0) * h / 2.0;
    p00 = 1 / x1_0;
    p01 = -x2_0 / (x1_0 * y2_0);
    p10 = 0;
    p11 = 1 / y2_0;
  }
}
void ParaFun::LocalCoordinateInverseScaf(int i, double &p00, double &p01,
                                         double &p10, double &p11) {
  int f0 = F_[0][i];
  int f1 = F_[1][i];
  int f2 = F_[2][i];

  Eigen::Vector2d x_(
      shell_data_.whole_uv_(f1, 0) - shell_data_.whole_uv_(f0, 0),
      shell_data_.whole_uv_(f1, 1) - shell_data_.whole_uv_(f0, 1));
  Eigen::Vector2d l_(
      shell_data_.whole_uv_(f2, 0) - shell_data_.whole_uv_(f0, 0),
      shell_data_.whole_uv_(f2, 1) - shell_data_.whole_uv_(f0, 1));

  double area_tri = std::abs(x_(0) * l_(1) - x_(1) * l_(0));
  double x1_0, x2_0, y2_0;
  if (area_tri > area_threshold_) {
    x1_0 = x_.norm();
    x_ /= x1_0;
    Eigen::Vector2d y_(-x_(1), x_(0));
    x2_0 = l_.dot(x_);
    y2_0 = l_.dot(y_);
  } else {
    // cout << "area too small!!!!!!!!!!!!! " << endl;
    double h = sqrt((2 * area_threshold_) / sqrt(3.0));
    x1_0 = h;
    x2_0 = h / 2.0;
    y2_0 = sqrt(3.0) * h / 2.0;
  }
  p00 = 1 / x1_0;
  p01 = -x2_0 / (x1_0 * y2_0);
  p10 = 0;
  p11 = 1 / y2_0;
}

double ParaFun::NewtonEquation(const double &a, const double &b,
                               const double &K) {
  double tt = 1;
  double E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) +
               pow(1 / b, 2 * tt) - K;
  while (abs(E_d) > 1e-5) {
    tt = tt - 1 /
                  (2 * log(a) * pow(a, 2 * tt) + 2 * log(b) * pow(b, 2 * tt) +
                   2 * log(1 / a) * pow(1 / a, 2 * tt) +
                   2 * log(1 / b) * pow(1 / b, 2 * tt)) *
                  (pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) +
                   pow(1 / b, 2 * tt) - K);
    E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) +
          pow(1 / b, 2 * tt) - K;
  }
  return tt;
}
void ParaFun::AdjustShellWeight(double new_weight) {
  shell_data_.shell_factor_ = new_weight;
  shell_data_.UpdateShell();
  InitArea();
}

double ParaFun::GetDistance(double s0, double s1, double e0, double e1,
                            double p0, double p1) {
  double s0_ = s0, s1_ = s1;
  double e0_ = e0, e1_ = e1;

  double len01 = std::sqrt((s0 - e0) * (s0 - e0) + (s1 - e1) * (s1 - e1));
  double len0i = std::sqrt((p0 - s0_) * (p0 - s0_) + (p1 - s1_) * (p1 - s1_));
  double len1i = std::sqrt((p0 - e0_) * (p0 - e0_) + (p1 - e1_) * (p1 - e1_));

  return len0i + len1i - len01;
}

void ParaFun::FunGrid(const Eigen::VectorXd &x) {
  const double *pos = x.data();

  x_min = pos[0];
  x_max = pos[0];
  y_min = pos[V_N_];
  y_max = pos[V_N_];
  for (int i = 1; i < V_N_; ++i) {
    if (pos[i] < x_min) {
      x_min = pos[i];
    } else if (pos[i] > x_max) {
      x_max = pos[i];
    }

    if (pos[i + V_N_] < y_min) {
      y_min = pos[i + V_N_];
    } else if (pos[i + V_N_] > y_max) {
      y_max = pos[i + V_N_];
    }
  }

  lengthgrid_x = (x_max - x_min) / (cellx_num - 1);
  x_max = x_min + cellx_num * lengthgrid_x;
  lengthgrid_y = (y_max - y_min) / (celly_num - 1);
  y_max = y_min + celly_num * lengthgrid_y;

  for (int j = 0; j < cell_points.size(); ++j) {
    cell_points[j].clear();
  }
  Eigen::VectorXi boundary_vertex = shell_data_.frame_ids_;
  int bound_num = BE_N_;
  double l_x_min, l_x_max, l_y_min, l_y_max, b, len;
  int id_x_min, id_x_max, id_y_min, id_y_max;
  int id_start, id_end;
  double s0, s1, e0, e1, p0, p1, s0_, s1_, e0_, e1_;

  AV_ID_.clear();
  AV_F_N_ = 0;
  double dis;

  for (int i = 0; i < bound_num; ++i) {
    int id = boundary_vertex(i);
    int x_i = std::floor((pos[id] - x_min) / lengthgrid_x);
    int y_i = std::floor((pos[id + V_N_] - y_min) / lengthgrid_y);
    if (y_i > celly_num - 1) {
      y_i = celly_num - 1;
    }
    if (x_i > cellx_num - 1) {
      x_i = cellx_num - 1;
    }
    cell_points[y_i * cellx_num + x_i].push_back(id);
  }

  for (int i = 0; i < bound_num; ++i) {
    id_start = boundary_vertex(i);
    id_end = boundary_vertex((i + 1) % BE_N_);
    s0 = pos[id_start];
    s1 = pos[id_start + V_N_];
    e0 = pos[id_end];
    e1 = pos[id_end + V_N_];
    len = sqrt((s0 - e0) * (s0 - e0) + (s1 - e1) * (s1 - e1));
    b = sqrt(threhold * len / 2 + threhold * threhold / 4);
    s0_ = s0;
    s1_ = s1;
    e0_ = e0;
    e1_ = e1;
    l_x_min = std::min(s0_, e0_);
    l_x_max = std::max(s0_, e0_);
    l_y_min = std::min(s1_, e1_);
    l_y_max = std::max(s1_, e1_);

    // local box
    l_x_min = std::max(l_x_min - b, x_min);
    l_x_max = std::min(l_x_max + b, x_max);
    l_y_min = std::max(l_y_min - b, y_min);
    l_y_max = std::min(l_y_max + b, y_max);

    id_x_min = std::floor((l_x_min - x_min) / lengthgrid_x);
    id_x_max = std::floor((l_x_max - x_min) / lengthgrid_x);
    id_y_min = std::floor((l_y_min - y_min) / lengthgrid_y);
    id_y_max = std::floor((l_y_max - y_min) / lengthgrid_y);

    if (id_y_max > celly_num - 1) {
      id_y_max = celly_num - 1;
    }
    if (id_x_max > cellx_num - 1) {
      id_x_max = cellx_num - 1;
    }

    for (int j_y = id_y_min; j_y < id_y_max + 1; ++j_y) {
      for (int j_x = id_x_min; j_x < id_x_max + 1; ++j_x) {
        int id_cell = j_x + j_y * cellx_num;
        int size = cell_points[id_cell].size();
        for (int k = 0; k < size; ++k) {
          int id_p = cell_points[id_cell][k];

          if (id_start == id_p || id_end == id_p) {
            continue;
          }

          p0 = pos[id_p];
          p1 = pos[id_p + V_N_];
          dis = sqrt((p0 - s0_) * (p0 - s0_) + (p1 - s1_) * (p1 - s1_)) +
                sqrt((p0 - e0_) * (p0 - e0_) + (p1 - e1_) * (p1 - e1_)) - len;
          if (dis > threhold) {
            continue;
          }
          if (dis < 0) {
            std::cout << "distance error" << std::endl;
          }

          if (BE_N_ - 1 == i) {
            AV_ID_.push_back(i * (BE_N_ - 2) + boundary_vertexID_[id_p] - 1);
          } else if (id_p < id_start) {
            AV_ID_.push_back(i * (BE_N_ - 2) + boundary_vertexID_[id_p]);
          } else {
            AV_ID_.push_back(i * (BE_N_ - 2) + boundary_vertexID_[id_p] - 2);
          }
        }
      }
    }
  }

  AV_F_N_ = AV_ID_.size();
  // std::cout << AV_F_N_ << std::endl;
}

ParaFun::~ParaFun() {}
