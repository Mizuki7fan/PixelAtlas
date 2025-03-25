#include "BiljectivePara.h"

BiljectivePara::BiljectivePara(const cgl::SurfaceMesh3 &M)
    : tri_mesh_(M),              //
      is_initialization_(false), //
      weight_1_(true),           //
      weight_2_(true)            //
{}

void BiljectivePara::Load() {
  is_initialization_ = true;
  int v_num = tri_mesh_.num_vertices();
  int f_num = tri_mesh_.num_faces();

  Eigen::MatrixXd V(v_num, 3);
  Eigen::MatrixXi F(f_num, 3);

  for (CGAL::SM_Vertex_index vertex : tri_mesh_.vertices()) {
    cgl::Point3II point = tri_mesh_.point(vertex);
    for (int j = 0; j < 3; ++j)
      V(vertex.idx(), j) = point[j];
  }

  for (CGAL::SM_Face_index face : tri_mesh_.faces()) {
    CGAL::SM_Vertex_index face_vertex[3];
    CGAL::SM_Halfedge_index halfedge = tri_mesh_.halfedge(face);
    face_vertex[0] = tri_mesh_.target(halfedge);
    face_vertex[1] = tri_mesh_.target(tri_mesh_.next(halfedge));
    face_vertex[2] = tri_mesh_.target(tri_mesh_.next(tri_mesh_.next(halfedge)));
    F(face.idx(), 0) = face_vertex[0].idx();
    F(face.idx(), 1) = face_vertex[1].idx();
    F(face.idx(), 2) = face_vertex[2].idx();
  }

  shell_data_ = std::make_unique<ShellData>();
  shell_data_->AddNewPatch(V, F, Eigen::RowVector2d(0, 0));
  parafun_solver_ = std::make_unique<Parafun>(shell_data_);

  uv_mesh_ = tri_mesh_;
  for (CGAL::SM_Vertex_index vertex : tri_mesh_.vertices()) {
    auto p = shell_data_->w_uv_.row(vertex.idx());
    uv_mesh_.point(vertex) = cgl::Point3II(p[0], p[1], 0.0);
  }
}

void BiljectivePara::Parameterization() {
  if (!is_initialization_)
    Load();

  int iteration_count = 0;
  double conv_rate_mesh = 1.0;
  double conv_rate_all = 1.0;
  long time_beg, time_end;

  parafun_solver_->AfterMeshImprove();
  double last_mesh_energy =
      parafun_solver_->ComputeEnergy(shell_data_->w_uv_, false) /
      shell_data_->mesh_measure_;
  parafun_solver_->AdjustShellWeight(
      (last_mesh_energy)*shell_data_->mesh_measure_ / (shell_data_->sf_num_) /
      1000.0);
  double last_all_energy =
      parafun_solver_->ComputeEnergy(shell_data_->w_uv_, true);

  time_beg = clock();
  time_end = clock();
  double time = (time_end - time_beg) / 1000.0;
  std::cout << time << " " << last_mesh_energy << std::endl;
  for (int i = 0; i < kMaxIterNum; i++) {
    iteration_count++;

    bool is_ip_convrate = true;
    if (conv_rate_mesh < 0.01)
      is_ip_convrate = false;
    bool is_slim_convrate = false;
    if (conv_rate_mesh > 0.1)
      is_slim_convrate = true;

    shell_data_->energy_ =
        parafun_solver_->BPE(is_ip_convrate, is_slim_convrate);
    double current_mesh_energy =
        parafun_solver_->energy_mesh_ / shell_data_->mesh_measure_;
    double current_all_energy = shell_data_->energy_;
    double mesh_energy_decrease = last_mesh_energy - current_mesh_energy;
    double all_energy_decrease = last_all_energy - current_all_energy;
    time_end = clock();
    double time_2 = (time_end - time_beg) / 1000.0;
    std::cout << time_2 << " " << current_mesh_energy << "	"
              << parafun_solver_->AV_F_N_H_ << std::endl;
    // cout << "Mesh Energy: " << current_mesh_energy << "\tEnergy Decrease: "
    // << mesh_energy_decrease << endl; cout << "All Energy:" <<
    // current_all_energy << "\tEnergy Decrease" << all_energy_decrease << endl;
    conv_rate_mesh = std::abs(mesh_energy_decrease) / last_mesh_energy;
    conv_rate_all = std::abs(all_energy_decrease) / last_all_energy;
    last_mesh_energy = current_mesh_energy;
    last_all_energy = AdjustWeight(conv_rate_mesh, last_mesh_energy);

    if (conv_rate_all < convgence_con_rate_)
      break;
  }

  time_end = clock();
  double time_consumption = (time_end - time_beg) / 1000.0;
  std::cout << "iter: " << iteration_count << std::endl;
  std::cout << "time: " << time_consumption << std::endl;
  std::cout << "energy:" << last_mesh_energy << std::endl;
  std::cout << "per time:" << time_consumption / iteration_count << std::endl;
  std::cout << parafun_solver_->time_1_ << " " << parafun_solver_->time_2_
            << " " << parafun_solver_->time_3_ << std::endl;
  std::cout << parafun_solver_->time_1_ << " "
            << parafun_solver_->time_2_ / iteration_count << " "
            << parafun_solver_->time_3_ / iteration_count << std::endl;
}

cgl::SurfaceMesh3 BiljectivePara::GetResult() {
  for (CGAL::SM_Vertex_index vertex : uv_mesh_.vertices()) {
    auto p = shell_data_->w_uv_.row(vertex.idx());
    uv_mesh_.point(vertex) = cgl::Point3II(p[0], p[1], 0.0);
  }
  return uv_mesh_;
}

double BiljectivePara::AdjustWeight(double conv_mesh, double last_mesh_energy) {
  if (conv_mesh > 1e-3 && weight_1_) {
    // parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure
    // / (shell_data.sf_num) / 10.0);
    parafun_solver_->AdjustShellWeight(
        (last_mesh_energy)*shell_data_->mesh_measure_ / (shell_data_->sf_num_) /
        1000.0);
    // parafun_solver->adjust_shell_weight(4.0*shell_data.mesh_measure /
    // (shell_data.sf_num) / 100000.0);
  } else if (conv_mesh > 1e-5 && weight_2_) {
    // parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure
    // / (shell_data.sf_num) / 10.0);
    parafun_solver_->AdjustShellWeight(4.0 * shell_data_->mesh_measure_ /
                                       (shell_data_->sf_num_) / 100000.0);
    // parafun_solver->adjust_shell_weight(4.0*shell_data.mesh_measure /
    // (shell_data.sf_num) / 1000.0);
    weight_1_ = false;
  } else {
    // parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure
    // / (shell_data.sf_num) / 10.0);
    parafun_solver_->AdjustShellWeight(4.0 * shell_data_->mesh_measure_ /
                                       shell_data_->sf_num_ / 100000000.0);
    // parafun_solver->adjust_shell_weight(4.0*shell_data.mesh_measure /
    // (shell_data.sf_num) / 1000.0);
    weight_2_ = false;
  }

  double last_all_energy =
      parafun_solver_->energy_mesh_ +
      parafun_solver_->area_.back() * parafun_solver_->energy_shell_ +
      parafun_solver_->barrier_coef_ * parafun_solver_->energy_barrier_;

  // std::cout << conv_mesh <<"   "<<parafun_solver->area.back() <<"    "<<
  // last_all_energy << std::endl;

  return last_all_energy;
};