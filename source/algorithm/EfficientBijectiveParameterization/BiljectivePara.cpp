#include "BiljectivePara.h"

BiljectivePara::BiljectivePara(const cgl::SurfaceMesh3 &m) : mesh_(m) {
  is_initialization_ = false;
  weight1_ = true;
  weight2_ = true;
  last_mesh_energy_ = std::numeric_limits<double>::max();
}

BiljectivePara::~BiljectivePara() {}

void BiljectivePara::Parameterization() {
  if (!is_initialization_)
    Load();

  int iteration_count = 0;
  double conv_rate_mesh = 1.0;
  double conv_rate_all = 1.0;
  long time_beg, time_end;

  solver_->AfterMeshImprove();
  last_mesh_energy_ = solver_->ComputeEnergy(shell_data_.whole_uv_, false) /
                      shell_data_.mesh_measure_;
  solver_->AdjustShellWeight((last_mesh_energy_)*shell_data_.mesh_measure_ /
                             (shell_data_.num_shell_faces_) / 1000.0);
  // solver_->AdjustShellWeight((last_mesh_energy)*shell_data_.mesh_measure_
  // / (shell_data_.sf_num) / 10.0);
  // solver_->AdjustShellWeight((last_mesh_energy)*shell_data_.mesh_measure_
  // / (shell_data_.sf_num) / 100000.0);
  double last_all_energy = solver_->ComputeEnergy(shell_data_.whole_uv_, true);

  // std::cout << "last_mesh_energy:" << last_mesh_energy << endl;
  // std::cout << "last_all_energy:" << last_all_energy << endl;

  time_beg = clock();
  time_end = clock();
  double time = (time_end - time_beg) / 1000.0;
  std::cout << time << " " << last_mesh_energy_ << endl;
  for (int i = 0; i < kNumMaxIteration; i++) {
    iteration_count++;

    bool is_ip_convrate = true;
    if (conv_rate_mesh < 0.01)
      is_ip_convrate = false;
    bool is_slim_convrate = false;
    if (conv_rate_mesh > 0.1)
      is_slim_convrate = true;

    shell_data_.mesh_energy_ = solver_->BPE(is_ip_convrate, is_slim_convrate);
    double current_mesh_energy =
        solver_->energy_mesh_ / shell_data_.mesh_measure_;
    double current_all_energy = shell_data_.mesh_energy_;
    double mesh_energy_decrease = last_mesh_energy_ - current_mesh_energy;
    double all_energy_decrease = last_all_energy - current_all_energy;
    time_end = clock();
    double time = (time_end - time_beg) / 1000.0;
    std::cout << time << " " << current_mesh_energy << "	"
              << solver_->AV_F_N_H_ << endl;
    // cout << "Mesh Energy: " << current_mesh_energy << "\tEnergy Decrease: "
    // << mesh_energy_decrease << endl; cout << "All Energy:" <<
    // current_all_energy << "\tEnergy Decrease" << all_energy_decrease << endl;
    conv_rate_mesh = abs(mesh_energy_decrease) / last_mesh_energy_;
    conv_rate_all = abs(all_energy_decrease) / last_all_energy;
    last_mesh_energy_ = current_mesh_energy;
    last_all_energy = AdjustWeight(conv_rate_mesh, last_mesh_energy_);

    if (conv_rate_all < convgence_con_rate_) {
      break;
    }
  }

  time_end = clock();
  double time_consumption = (time_end - time_beg) / 1000.0;
  std::cout << "iter: " << iteration_count << std::endl;
  std::cout << "time: " << time_consumption << std::endl;
  std::cout << "energy:" << last_mesh_energy_ << std::endl;
  std::cout << "per time:" << time_consumption / iteration_count << std::endl;
  cout << solver_->time_1_ << " " << solver_->time_2_ << " " << solver_->time_3_
       << endl;
  cout << solver_->time_1_ << " " << solver_->time_2_ / iteration_count << " "
       << solver_->time_3_ / iteration_count << endl;
  cout << solver_->density_ << endl;
}

double BiljectivePara::AdjustWeight(double conv_mesh, double last_mesh_energy) {
  if (conv_mesh > 1e-3 && weight1_) {
    // solver_->AdjustShellWeight((last_mesh_energy)*shell_data_.mesh_measure_
    // / (shell_data_.sf_num) / 10.0);
    solver_->AdjustShellWeight((last_mesh_energy)*shell_data_.mesh_measure_ /
                               (shell_data_.num_shell_faces_) / 1000.0);
    // solver_->AdjustShellWeight(4.0*shell_data_.mesh_measure_ /
    // (shell_data_.sf_num) / 100000.0);
  } else if (conv_mesh > 1e-5 && weight2_) {
    // solver_->AdjustShellWeight((last_mesh_energy)*shell_data_.mesh_measure_
    // / (shell_data_.sf_num) / 10.0);
    solver_->AdjustShellWeight(4.0 * shell_data_.mesh_measure_ /
                               (shell_data_.num_shell_faces_) / 100000.0);
    // solver_->AdjustShellWeight(4.0*shell_data_.mesh_measure_ /
    // (shell_data_.sf_num) / 1000.0);
    weight1_ = false;
  } else {
    // solver_->AdjustShellWeight((last_mesh_energy)*shell_data_.mesh_measure_
    // / (shell_data_.sf_num) / 10.0);
    solver_->AdjustShellWeight(4.0 * shell_data_.mesh_measure_ /
                               shell_data_.num_shell_faces_ / 100000000.0);
    // solver_->AdjustShellWeight(4.0*shell_data_.mesh_measure_ /
    // (shell_data_.sf_num) / 1000.0);
    weight2_ = false;
  }

  double last_all_energy = solver_->energy_mesh_ +
                           solver_->area_.back() * solver_->energy_shell_ +
                           solver_->barrer_coef_ * solver_->energy_barrier_;

  // std::cout << conv_mesh <<"   "<<solver_->area.back() <<"    "<<
  // last_all_energy << std::endl;

  return last_all_energy;
}

void BiljectivePara::Load() {
  is_initialization_ = true;

  int v_num = mesh_.num_vertices();
  int f_num = mesh_.num_faces();
  Eigen::MatrixXd V(v_num, 3);
  Eigen::MatrixXi F(f_num, 3);

  for (CGAL::SM_Vertex_index vertex : mesh_.vertices()) {
    cgl::Point3II pos = mesh_.point(vertex);
    for (int j = 0; j < 3; ++j) {
      V(vertex.idx(), j) = pos[j];
    }
  }

  for (const CGAL::SM_Face_index &face : mesh_.faces()) {
    CGAL::SM_Vertex_index face_vertex[3];
    CGAL::SM_Halfedge_index halfedge = mesh_.halfedge(face);
    face_vertex[0] = mesh_.target(halfedge);
    face_vertex[1] = mesh_.target(mesh_.next(halfedge));
    face_vertex[2] = mesh_.target(mesh_.next(mesh_.next(halfedge)));
    F(face.idx(), 0) = face_vertex[0].idx();
    F(face.idx(), 1) = face_vertex[1].idx();
    F(face.idx(), 2) = face_vertex[2].idx();
  }

  shell_data_ = ShellData();
  shell_data_.AddNewPatch(V, F, Eigen::RowVector2d(0, 0));
  solver_.reset(new ParaFun(shell_data_));

  uv_mesh_ = mesh_;
  for (CGAL::SM_Vertex_index vertex : uv_mesh_.vertices()) {
    auto p = shell_data_.whole_uv_.row(vertex.idx());
    cgl::Point3II pos(p(0), p(1), 0.0);
    uv_mesh_.point(vertex) = pos;
  }
}

void BiljectivePara::WriteUVMesh(std::ofstream &of_obj) {
  for (CGAL::SM_Vertex_index vertex : uv_mesh_.vertices()) {
    auto pos = shell_data_.whole_uv_.row(vertex.idx());
    uv_mesh_.point(vertex) = cgl::Point3II(pos[0], pos[1], 0.0);
  }

  // 将uv_mesh放缩到(0,1)范围, 且中心为(0.5,0.5)
  cgl::MeshOperation::UV::NormalizeUVToUnitSquare(uv_mesh_);
  CGAL::IO::write_OBJ(of_obj, uv_mesh_);
}

void BiljectivePara::ShellTri(MatrixXi &tri, MatrixXd &pre_pos, MatrixXd &pos,
                              VectorXi &bnd) {
  tri = shell_data_.shell_faces_;
  pre_pos = shell_data_.whole_uv_pre_;
  pos = shell_data_.whole_uv_;
  bnd = shell_data_.frame_ids_;
}
