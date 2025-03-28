#include "Common.h"
#if defined(USE_MKL)
#include <LinSysSolver-Interface/MKLPardisoSolver.h>
#elif defined(USE_EIGEN)
#include <LinSysSolver-Interface/EigenLinSolver.h>
#endif
#include <numbers>
#include <set>

/*加入判断边界线是否闭合,return true means open curve; return false means closed
curve. In closed case, the first and last item of order_boundary are same.*/
void Tutte(const int &num_vertices,              //
           const Eigen::MatrixXi &face_vertices, //
           const Eigen::VectorXi &bnd,           //
           const Eigen::MatrixXd &bnd_uv,        //
           Eigen::MatrixXd &uv_init) {
  int num_faces = static_cast<int>(face_vertices.rows());
  uv_init.resize(num_vertices, 2);

  std::set<int> bound_idxes;
  for (std::size_t i = 0; i < static_cast<std::size_t>(bnd.size()); i++) {
    bound_idxes.insert(bnd(i));
    uv_init.row(bnd(i)) << bnd_uv(i, 0), bnd_uv(i, 1);
  }

  std::vector<std::set<int>> VV_tmp(num_vertices);
  for (std::size_t i = 0; i < num_faces; i++) {
    VV_tmp[face_vertices[0]].insert(face_vertices[1]);
    VV_tmp[face_vertices[0]].insert(face_vertices[2]);

    VV_tmp[face_vertices[1]].insert(face_vertices[0]);
    VV_tmp[face_vertices[1]].insert(face_vertices[2]);

    VV_tmp[face_vertices[2]].insert(face_vertices[0]);
    VV_tmp[face_vertices[2]].insert(face_vertices[1]);
  }

  std::unique_ptr<Solver> solver = nullptr;

#if defined(USE_MKL)
  solver = std::unique_ptr<MKLPardisoSolver>();
#elif defined(USE_EIGEN)
  solver = std::make_unique<EigenLinSolver>();
#endif

  std::vector<double> solver_tu;
  std::vector<double> solver_tv;

  solver->ia_.reserve(num_vertices + 1);
  solver->ja_.reserve(8 * num_vertices);
  solver->a_.reserve(8 * num_vertices);
  solver_tu.resize(num_vertices, 0.0);
  solver_tv.resize(num_vertices, 0.0);

  for (size_t i = 0; i < num_vertices; i++) {
    solver->ia_.push_back(solver->ja_.size());

    if (bound_idxes.count(i) > 0) {
      solver->ja_.push_back(i);
      solver->a_.push_back(1.0);

      solver_tu[i] = uv_init(i, 0);
      solver_tv[i] = uv_init(i, 1);

    } else {
      solver->ja_.push_back(i);
      solver->a_.push_back(VV_tmp[i].size());
      std::vector<int> row_id;
      row_id.reserve(VV_tmp[i].size());
      double bu = 0.0;
      double bv = 0.0;
      for (auto &vv_id : VV_tmp[i]) {
        if (bound_idxes.count(vv_id) > 0) {
          bu += uv_init(vv_id, 0);
          bv += uv_init(vv_id, 1);
        } else {
          if (vv_id > i) {
            row_id.push_back(vv_id);
          }
        }
      }
      std::sort(row_id.begin(), row_id.end(), std::less<int>());
      for (size_t j = 0; j < row_id.size(); j++) {
        solver->ja_.push_back(row_id[j]);
        solver->a_.push_back(-1.0);
      }
      solver_tu[i] = bu;
      solver_tv[i] = bv;
    }
  }
  solver->ia_.push_back(solver->ja_.size());

  solver->nnz_ = solver->ja_.size();
  solver->num_ = num_vertices;

  solver->pardiso_init();
  solver->rhs_ = solver_tu;

  solver->factorize();
  solver->pardiso_solver();

  for (std::size_t i = 0; i < num_vertices; i++)
    uv_init(i, 0) = solver->result_[i];

  solver->rhs_ = solver_tv;
  solver->pardiso_solver();
  for (std::size_t i = 0; i < num_vertices; i++)
    uv_init(i, 1) = solver->result_[i];
}

void preCalc_pardiso(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                     Solver &pardiso) {

  int V_N = V.rows();
  int F_N = F.rows();
  pardiso.ia_.clear();
  pardiso.ia_.reserve(2 * V_N + 1);
  pardiso.ja_.clear();
  pardiso.ja_.reserve(8 * V_N);

  std::vector<std::set<int>> VV_tmp;
  VV_tmp.resize(V_N);
  for (size_t i = 0; i < F_N; i++) {
    int vid[3];

    for (size_t j = 0; j < F.cols(); j++) {
      vid[j] = F(i, j);
    }
    VV_tmp[vid[0]].insert(vid[1]);
    VV_tmp[vid[0]].insert(vid[2]);

    VV_tmp[vid[1]].insert(vid[0]);
    VV_tmp[vid[1]].insert(vid[2]);

    VV_tmp[vid[2]].insert(vid[0]);
    VV_tmp[vid[2]].insert(vid[1]);
  }

  for (int i = 0; i < V_N; i++) {
    pardiso.ia_.push_back(pardiso.ja_.size());
    VV_tmp[i].insert(i);
    vector<int> row_id;
    for (auto &var : VV_tmp[i]) {
      row_id.push_back(var);
    }

    vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);

    int dd = 0;
    for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++) {
      pardiso.ja_.push_back(row_id[k]);
      ++dd;
    }
    for (int k = 0; k < row_id.size(); k++) {
      pardiso.ja_.push_back(row_id[k] + V_N);
      ++dd;
    }
  }
  for (int i = V_N; i < 2 * V_N; i++) {
    pardiso.ia_.push_back(pardiso.ja_.size());
    vector<int> row_id;
    for (auto &var : VV_tmp[i - V_N]) {
      row_id.push_back(var);
    }
    vector<int>::iterator iter =
        std::find(row_id.begin(), row_id.end(), i - V_N);

    int dd = 0;
    for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++) {
      pardiso.ja_.push_back(row_id[k] + V_N);
      ++dd;
    }
  }
  pardiso.ia_.push_back(pardiso.ja_.size());
}

void map_vertices_to_circle(const Eigen::MatrixXd &V,
                            const Eigen::VectorXi &bnd, Eigen::MatrixXd &UV) {
  // Get sorted list of boundary vertices
  std::vector<int> interior, map_ij;
  map_ij.resize(V.rows());

  std::vector<bool> isOnBnd(V.rows(), false);
  for (int i = 0; i < bnd.size(); i++) {
    isOnBnd[bnd[i]] = true;
    map_ij[bnd[i]] = i;
  }

  for (int i = 0; i < (int)isOnBnd.size(); i++) {
    if (!isOnBnd[i]) {
      map_ij[i] = interior.size();
      interior.push_back(i);
    }
  }

  // Map boundary to unit circle
  std::vector<double> len(bnd.size());
  len[0] = 0.;

  for (int i = 1; i < bnd.size(); i++) {
    len[i] = len[i - 1] + (V.row(bnd[i - 1]) - V.row(bnd[i])).norm();
  }
  double total_len =
      len[len.size() - 1] + (V.row(bnd[0]) - V.row(bnd[bnd.size() - 1])).norm();

  UV.resize(bnd.size(), 2);

  for (int i = 0; i < bnd.size(); i++) {
    double frac = len[i] * 2. * std::numbers::pi / total_len;
    // double frac = i * 2. * M_PI / (bnd.size());
    UV.row(map_ij[bnd[i]]) << cos(frac), sin(frac);
  }
}

void boundary_loop(const Eigen::MatrixXi &F_ref,
                   std::vector<std::vector<int>> &boundaryloop) {
  std::vector<std::vector<int>> boundaryEdges;
  std::vector<std::vector<int>> edges;
  int n_fvs = F_ref.cols();

  for (int it = 0; it < F_ref.rows(); it++) {
    for (int i = 0; i < n_fvs; i++) {
      int var = F_ref(it, i);
      int var_n = F_ref(it, (i + 1) % n_fvs);
      if (var > var_n)
        std::swap(var, var_n);
      std::vector<int> edge(4);
      edge[0] = var;
      edge[1] = var_n;
      edge[2] = it;
      edge[3] = i;
      edges.emplace_back(edge);
    }
  }
  std::sort(edges.begin(), edges.end());
  int i = 1;
  for (; i < edges.size();) {
    auto &r1 = edges[i - 1];
    auto &r2 = edges[i];
    if ((r1[0] == r2[0]) && (r1[1] == r2[1])) {
      i += 2;
    } else {
      boundaryEdges.emplace_back(edges[i - 1]);
      i++;
    }
  }
  if (i == edges.size())
    boundaryEdges.emplace_back(edges.back());

  for (auto &var : boundaryEdges) {
    var[0] = F_ref(var[2], var[3]);
    var[1] = F_ref(var[2], (var[3] + 1) % n_fvs);
  }
  int ev0 = boundaryEdges.front()[0];
  int ev1 = boundaryEdges.front()[1];

  vector<int> visited;
  visited.resize(boundaryEdges.size(), 0);
  visited[0] = 1;
  vector<int> loop0;
  loop0.push_back(ev1);
  while (ev1 != ev0) {
    for (int i = 1; i < boundaryEdges.size(); i++) {
      if (visited[i] == 1)
        continue;
      if (boundaryEdges[i][0] == ev1) {
        visited[i] = 1;
        ev1 = boundaryEdges[i][1];
        loop0.push_back(ev1);
        break;
      }
    }
  }
  boundaryloop.emplace_back(loop0);
}
