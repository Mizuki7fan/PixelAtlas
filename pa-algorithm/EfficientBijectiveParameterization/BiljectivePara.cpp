#include "BiljectivePara.h"

BiljectivePara::BiljectivePara(const cgl::SurfaceMesh3 &m, string filename)
    : mesh(m) {
  string file_str = filename;
  modelname = file_str.substr(file_str.find_last_of('/') + 1);
  modelname.replace(modelname.end() - 4, modelname.end(), "");
  is_initialization = false;
  weight1 = true;
  weight2 = true;
  last_mesh_energy_ = std::numeric_limits<double>::max();
}

BiljectivePara::~BiljectivePara() {}

void BiljectivePara::parameterization() {
  if (!is_initialization) {
    load();
  }

  int iteration_count = 0;
  double conv_rate_mesh = 1.0;
  double conv_rate_all = 1.0;
  long time_beg, time_end;
  bool is_first = true;
  bool is_second = true;

  parafun_solver->AfterMeshImprove();
  last_mesh_energy_ =
      parafun_solver->ComputeEnergy(shell_data.whole_uv_, false) /
      shell_data.mesh_measure_;
  parafun_solver->adjust_shell_weight(
      (last_mesh_energy_)*shell_data.mesh_measure_ /
      (shell_data.num_shell_faces_) / 1000.0);
  // parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure_
  // / (shell_data.sf_num) / 10.0);
  // parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure_
  // / (shell_data.sf_num) / 100000.0);
  double last_all_energy =
      parafun_solver->ComputeEnergy(shell_data.whole_uv_, true);

  // std::cout << "last_mesh_energy:" << last_mesh_energy << endl;
  // std::cout << "last_all_energy:" << last_all_energy << endl;

  time_beg = clock();
  time_end = clock();
  double time = (time_end - time_beg) / 1000.0;
  std::cout << time << " " << last_mesh_energy_ << endl;
  for (int i = 0; i < MAX_ITER_NUM; i++) {
    iteration_count++;

    bool is_ip_convrate = true;
    if (conv_rate_mesh < 0.01)
      is_ip_convrate = false;
    bool is_slim_convrate = false;
    if (conv_rate_mesh > 0.1)
      is_slim_convrate = true;

    shell_data.mesh_energy_ =
        parafun_solver->BPE(is_ip_convrate, is_slim_convrate);
    double current_mesh_energy =
        parafun_solver->energy_mesh / shell_data.mesh_measure_;
    double current_all_energy = shell_data.mesh_energy_;
    double mesh_energy_decrease = last_mesh_energy_ - current_mesh_energy;
    double all_energy_decrease = last_all_energy - current_all_energy;
    time_end = clock();
    double time = (time_end - time_beg) / 1000.0;
    std::cout << time << " " << current_mesh_energy << "	"
              << parafun_solver->AV_F_N_H << endl;
    // cout << "Mesh Energy: " << current_mesh_energy << "\tEnergy Decrease: "
    // << mesh_energy_decrease << endl; cout << "All Energy:" <<
    // current_all_energy << "\tEnergy Decrease" << all_energy_decrease << endl;
    conv_rate_mesh = abs(mesh_energy_decrease) / last_mesh_energy_;
    conv_rate_all = abs(all_energy_decrease) / last_all_energy;
    last_mesh_energy_ = current_mesh_energy;
    last_all_energy = adjust_weight(conv_rate_mesh, last_mesh_energy_);

    if (conv_rate_all < convgence_con_rate) {
      break;
    }
  }

  time_end = clock();
  double time_consumption = (time_end - time_beg) / 1000.0;
  std::cout << "iter: " << iteration_count << std::endl;
  std::cout << "time: " << time_consumption << std::endl;
  std::cout << "energy:" << last_mesh_energy_ << std::endl;
  std::cout << "per time:" << time_consumption / iteration_count << std::endl;
  cout << parafun_solver->time1 << " " << parafun_solver->time2 << " "
       << parafun_solver->time3 << endl;
  cout << parafun_solver->time1 << " "
       << parafun_solver->time2 / iteration_count << " "
       << parafun_solver->time3 / iteration_count << endl;
  cout << parafun_solver->density << endl;
}

double BiljectivePara::adjust_weight(double conv_mesh,
                                     double last_mesh_energy) {
  if (conv_mesh > 1e-3 && weight1) {
    // parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure_
    // / (shell_data.sf_num) / 10.0);
    parafun_solver->adjust_shell_weight(
        (last_mesh_energy)*shell_data.mesh_measure_ /
        (shell_data.num_shell_faces_) / 1000.0);
    // parafun_solver->adjust_shell_weight(4.0*shell_data.mesh_measure_ /
    // (shell_data.sf_num) / 100000.0);
  } else if (conv_mesh > 1e-5 && weight2) {
    // parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure_
    // / (shell_data.sf_num) / 10.0);
    parafun_solver->adjust_shell_weight(4.0 * shell_data.mesh_measure_ /
                                        (shell_data.num_shell_faces_) /
                                        100000.0);
    // parafun_solver->adjust_shell_weight(4.0*shell_data.mesh_measure_ /
    // (shell_data.sf_num) / 1000.0);
    weight1 = false;
  } else {
    // parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure_
    // / (shell_data.sf_num) / 10.0);
    parafun_solver->adjust_shell_weight(4.0 * shell_data.mesh_measure_ /
                                        shell_data.num_shell_faces_ /
                                        100000000.0);
    // parafun_solver->adjust_shell_weight(4.0*shell_data.mesh_measure_ /
    // (shell_data.sf_num) / 1000.0);
    weight2 = false;
  }

  double last_all_energy =
      parafun_solver->energy_mesh +
      parafun_solver->area.back() * parafun_solver->energy_shell +
      parafun_solver->barrer_coef * parafun_solver->energy_barrier;

  // std::cout << conv_mesh <<"   "<<parafun_solver->area.back() <<"    "<<
  // last_all_energy << std::endl;

  return last_all_energy;
}

void BiljectivePara::load() {
  is_initialization = true;

  int v_num = mesh.num_vertices();
  int f_num = mesh.num_faces();
  Eigen::MatrixXd V(v_num, 3);
  Eigen::MatrixXi F(f_num, 3);

  for (CGAL::SM_Vertex_index vertex : mesh.vertices()) {
    cgl::Point3II pos = mesh.point(vertex);
    for (int j = 0; j < 3; ++j) {
      V(vertex.idx(), j) = pos[j];
    }
  }

  for (const CGAL::SM_Face_index &face : mesh.faces()) {
    CGAL::SM_Vertex_index face_vertex[3];
    CGAL::SM_Halfedge_index halfedge = mesh.halfedge(face);
    face_vertex[0] = mesh.target(halfedge);
    face_vertex[1] = mesh.target(mesh.next(halfedge));
    face_vertex[2] = mesh.target(mesh.next(mesh.next(halfedge)));
    F(face.idx(), 0) = face_vertex[0].idx();
    F(face.idx(), 1) = face_vertex[1].idx();
    F(face.idx(), 2) = face_vertex[2].idx();
  }

  shell_data = ShellData();
  shell_data.AddNewPatch(V, F, Eigen::RowVector2d(0, 0));
  parafun_solver.reset(new ParaFun(shell_data));

  uv_mesh_ = mesh;
  for (CGAL::SM_Vertex_index vertex : uv_mesh_.vertices()) {
    auto p = shell_data.whole_uv_.row(vertex.idx());
    cgl::Point3II pos(p(0), p(1), 0.0);
    uv_mesh_.point(vertex) = pos;
  }
}

void BiljectivePara::WriteUVMesh(std::ofstream &of_obj) {
  for (CGAL::SM_Vertex_index vertex : uv_mesh_.vertices()) {
    auto pos = shell_data.whole_uv_.row(vertex.idx());
    uv_mesh_.point(vertex) = cgl::Point3II(pos[0], pos[1], 0.0);
  }

  CGAL::IO::write_OBJ(of_obj, uv_mesh_);
}

void BiljectivePara::shelltri(MatrixXi &tri, MatrixXd &pre_pos, MatrixXd &pos,
                              VectorXi &bnd) {
  tri = shell_data.shell_faces_;
  pre_pos = shell_data.whole_uv_pre_;
  pos = shell_data.whole_uv_;
  bnd = shell_data.frame_ids_;
}
