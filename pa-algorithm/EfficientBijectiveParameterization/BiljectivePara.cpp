#include "BiljectivePara.h"

BiljectivePara::BiljectivePara(Mesh &m, string filename) : mesh(m) {
  string file_str = filename;
  modelname = file_str.substr(file_str.find_last_of('/') + 1);
  modelname.replace(modelname.end() - 4, modelname.end(), "");
  is_initialization = false;
  weight1 = true;
  weight2 = true;
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

  parafun_solver->after_mesh_improve();
  double last_mesh_energy =
      parafun_solver->compute_energy(shell_data.w_uv, false) /
      shell_data.mesh_measure;
  parafun_solver->adjust_shell_weight(
      (last_mesh_energy)*shell_data.mesh_measure / (shell_data.sf_num) /
      1000.0);
  // parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure
  // / (shell_data.sf_num) / 10.0);
  // parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure
  // / (shell_data.sf_num) / 100000.0);
  double last_all_energy =
      parafun_solver->compute_energy(shell_data.w_uv, true);

  // std::cout << "last_mesh_energy:" << last_mesh_energy << endl;
  // std::cout << "last_all_energy:" << last_all_energy << endl;

  time_beg = clock();
  time_end = clock();
  double time = (time_end - time_beg) / 1000.0;
  std::cout << time << " " << last_mesh_energy << endl;
  for (int i = 0; i < MAX_ITER_NUM; i++) {
    iteration_count++;

    bool is_ip_convrate = true;
    if (conv_rate_mesh < 0.01)
      is_ip_convrate = false;
    bool is_slim_convrate = false;
    if (conv_rate_mesh > 0.1)
      is_slim_convrate = true;

    shell_data.energy = parafun_solver->BPE(is_ip_convrate, is_slim_convrate);
    double current_mesh_energy =
        parafun_solver->energy_mesh / shell_data.mesh_measure;
    double current_all_energy = shell_data.energy;
    double mesh_energy_decrease = last_mesh_energy - current_mesh_energy;
    double all_energy_decrease = last_all_energy - current_all_energy;
    time_end = clock();
    double time = (time_end - time_beg) / 1000.0;
    std::cout << time << " " << current_mesh_energy << "	"
              << parafun_solver->AV_F_N_H << endl;
    // cout << "Mesh Energy: " << current_mesh_energy << "\tEnergy Decrease: "
    // << mesh_energy_decrease << endl; cout << "All Energy:" <<
    // current_all_energy << "\tEnergy Decrease" << all_energy_decrease << endl;
    conv_rate_mesh = abs(mesh_energy_decrease) / last_mesh_energy;
    conv_rate_all = abs(all_energy_decrease) / last_all_energy;
    last_mesh_energy = current_mesh_energy;
    last_all_energy = adjust_weight(conv_rate_mesh, last_mesh_energy);

    if (conv_rate_all < convgence_con_rate)
    //	if (current_mesh_energy < 4.0544)
    {
      break;
    }
  }

  time_end = clock();
  double time_consumption = (time_end - time_beg) / 1000.0;
  std::cout << "iter: " << iteration_count << std::endl;
  std::cout << "time: " << time_consumption << std::endl;
  std::cout << "energy:" << last_mesh_energy << std::endl;
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
    // parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure
    // / (shell_data.sf_num) / 10.0);
    parafun_solver->adjust_shell_weight(
        (last_mesh_energy)*shell_data.mesh_measure / (shell_data.sf_num) /
        1000.0);
    // parafun_solver->adjust_shell_weight(4.0*shell_data.mesh_measure /
    // (shell_data.sf_num) / 100000.0);
  } else if (conv_mesh > 1e-5 && weight2) {
    // parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure
    // / (shell_data.sf_num) / 10.0);
    parafun_solver->adjust_shell_weight(4.0 * shell_data.mesh_measure /
                                        (shell_data.sf_num) / 100000.0);
    // parafun_solver->adjust_shell_weight(4.0*shell_data.mesh_measure /
    // (shell_data.sf_num) / 1000.0);
    weight1 = false;
  } else {
    // parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure
    // / (shell_data.sf_num) / 10.0);
    parafun_solver->adjust_shell_weight(4.0 * shell_data.mesh_measure /
                                        shell_data.sf_num / 100000000.0);
    // parafun_solver->adjust_shell_weight(4.0*shell_data.mesh_measure /
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

  int v_num = mesh.n_vertices();
  int f_num = mesh.n_faces();
  Eigen::MatrixXd V(v_num, 3);
  Eigen::MatrixXi F(f_num, 3);

  for (int i = 0; i < v_num; ++i) {
    auto vertex = mesh.vertex_handle(i);
    Mesh::Point pos = mesh.point(vertex);
    for (int j = 0; j < 3; ++j) {
      V(i, j) = pos[j];
    }
  }

  for (int i = 0; i < f_num; ++i) {
    auto face = mesh.face_handle(i);
    int dd = 0;
    for (Mesh::FaceVertexIter it2 = mesh.fv_begin(face);
         it2 != mesh.fv_end(face); ++it2) {
      auto vertex_ = it2.handle();
      switch (dd) {
      case 0:
        F(i, 0) = vertex_.idx();
        break;
      case 1:
        F(i, 1) = vertex_.idx();
        break;
      case 2:
        F(i, 2) = vertex_.idx();
        break;
      default:
        break;
      }
      dd++;
    }
  }

  shell_data = ShellData();
  shell_data.add_new_patch(V, F, Eigen::RowVector2d(0, 0));
  parafun_solver.reset(new Parafun(shell_data));

  for (int i = 0; i < v_num; ++i) {
    auto vertex = mesh.vertex_handle(i);
    auto p = shell_data.w_uv.row(i);
    Mesh::Point pos(p(0), p(1), 0.0);
    mesh.set_point(vertex, pos);
  }
  //   Mesh::VertexIter v_it = mesh.vertices_begin();
  //   OpenMesh::Vec3d nor;
  //   for (v_it; v_it != mesh.vertices_end(); ++v_it) {
  //     mesh.set_normal(v_it, OpenMesh::Vec3d(0, 0, 0.95));
  //   }
  //   Mesh::ConstFaceIter fIt(mesh.faces_begin()), fEnd(mesh.faces_end());
  //   for (; fIt != fEnd; ++fIt) {
  //     mesh.set_normal(fIt, OpenMesh::Vec3d(0, 0, 0.95));
  //   }
}

void BiljectivePara::write_obj(std::string outstr) {
  Eigen::MatrixXd &V_in = shell_data.w_uv;
  int v_num = mesh.n_vertices();
  for (int i = 0; i < v_num; ++i) {
    auto vertex = mesh.vertex_handle(i);
    Mesh::Point pos(V_in(i, 0), V_in(i, 1), 0.0);
    mesh.set_point(vertex, pos);
  }

  Eigen::MatrixXi &F_ref = shell_data.m_T;

  ofstream of_obj(outstr, ios::trunc);

  if (V_in.cols() == 3) {
    for (size_t vid = 0; vid < V_in.rows(); vid++) {
      of_obj << "v " << V_in(vid, 0) << " " << V_in(vid, 1) << " "
             << V_in(vid, 2) << endl;
    }
  } else if (V_in.cols() == 2) {
    for (size_t vid = 0; vid < V_in.rows(); vid++) {
      of_obj << "v " << V_in(vid, 0) << " " << V_in(vid, 1) << " " << 0.0
             << endl;
    }
  }

  for (size_t fi = 0; fi < F_ref.rows(); fi++) {
    of_obj << "f " << F_ref(fi, 0) + 1 << " " << F_ref(fi, 1) + 1 << " "
           << F_ref(fi, 2) + 1 << endl;
  }
  of_obj.close();
}

void BiljectivePara::shelltri(MatrixXi &tri, MatrixXd &pre_pos, MatrixXd &pos,
                              VectorXi &bnd) {
  tri = shell_data.s_T;
  pre_pos = shell_data.w_uv_pre;
  pos = shell_data.w_uv;
  bnd = shell_data.frame_ids;
}
