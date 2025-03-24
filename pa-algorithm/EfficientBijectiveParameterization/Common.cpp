#include "Common.h"
#if defined(USE_MKL)
#include "LinSysSolverInterface/MKLPardisoSolver.h"
#elif defined(USE_EIGEN)
#include "LinSysSolverInterface/EigenLinSolver.h"
#endif
#include "frame/assert.hpp"
#include <iostream>
#include <numbers>

// 使得cut成为有序点列
void OrderTheBoundary(std::vector<std::list<int>> &order_boundary,
                      const std::vector<int> &nextlocation) {
  int n_es = static_cast<int>(nextlocation.size()) / 2;
  std::set<int> isused;
  for (int i = 0; i < n_es; i++) {
    isused.insert(i);
  }

  int s_1 = 0;
  int s_2 = 1;

  // bool flag = true;
  while (!isused.empty()) {
    s_1 = *(isused.begin()) * 2;
    s_2 = s_1 + 1;
    std::list<int> list_1, list_2;
    if (GrowFromP(s_1, isused, list_1, nextlocation)) {
      GrowFromP(s_2, isused, list_2, nextlocation);
      for (auto it = list_2.begin(); it != list_2.end(); it++) {
        list_1.push_front(*it);
      }
    }
    order_boundary.emplace_back(list_1);
  }
}

/*加入判断边界线是否闭合,return true means open curve; return false means closed
curve. In closed case, the first and last item of order_boundary are same.*/
bool GrowFromP(int p, std::set<int> &isused, std::list<int> &order_boundary,
               const std::vector<int> &nextlocation) {
  int loc = nextlocation[p];
  isused.erase(p / 2);
  if (loc == -1) {
    order_boundary.push_back(p);
    return true;
  } else {
    order_boundary.push_back(loc);
  }

  while (true) {
    isused.erase(loc / 2);
    if (loc % 2 == 1) {
      if ((loc - 1) == p) {
        order_boundary.push_back(nextlocation[p]);
        return false;
      }

      if (nextlocation[loc - 1] != -1) {
        order_boundary.push_back(nextlocation[loc - 1]);
        loc = nextlocation[loc - 1];
      } else {
        order_boundary.push_back(loc - 1);
        break;
      }
    } else {
      if ((loc + 1) == p) {
        order_boundary.push_back(nextlocation[p]);
        return false;
      }
      if (nextlocation[loc + 1] != -1) {
        order_boundary.push_back(nextlocation[loc + 1]);
        loc = nextlocation[loc + 1];
      } else {
        order_boundary.push_back(loc + 1);
        break;
      }
    }
  }

  return true;
}

void Tutte(int V_N, const Eigen::MatrixXi &F, const Eigen::VectorXi &bnd,
           const Eigen::MatrixXd &bnd_uv, Eigen::MatrixXd &uv_init) {
  int F_N = static_cast<int>(F.rows());

  uv_init.resize(V_N, 2);
  std::set<int> bound_ids;
  for (int i = 0; i < bnd.size(); i++) {
    bound_ids.insert(bnd(i));
    uv_init.row(bnd(i)) << bnd_uv(i, 0), bnd_uv(i, 1);
  }

  std::vector<std::set<int>> VV_tmp;
  VV_tmp.resize(V_N);
  for (int i = 0; i < F_N; i++) {
    int vid[3];

    for (int j = 0; j < F.cols(); j++) {
      vid[j] = F(i, j);
    }
    VV_tmp[vid[0]].insert(vid[1]);
    VV_tmp[vid[0]].insert(vid[2]);

    VV_tmp[vid[1]].insert(vid[0]);
    VV_tmp[vid[1]].insert(vid[2]);

    VV_tmp[vid[2]].insert(vid[0]);
    VV_tmp[vid[2]].insert(vid[1]);
  }

  std::unique_ptr<Solver> solver = nullptr;
#if defined(USE_MKL)
  solver = std::make_unique<MKLPardisoSolver>();
#elif defined(USE_EIGEN)
  solver = std::make_unique<EigenLinSolver>();
#endif
  PA_ASSERT(solver != nullptr);

  std::vector<double> pardiso_tu;
  std::vector<double> pardiso_tv;

  solver->ia_.reserve(V_N + 1);
  solver->ja_.reserve(8 * V_N);
  solver->a_.reserve(8 * V_N);
  pardiso_tu.resize(V_N, 0.0);
  pardiso_tv.resize(V_N, 0.0);

  for (int i = 0; i < V_N; i++) {
    solver->ia_.push_back(static_cast<int>(solver->ja_.size()));

    if (bound_ids.count(i) > 0) {
      solver->ja_.push_back(i);
      solver->a_.push_back(1.0);

      pardiso_tu[i] = uv_init(i, 0);
      pardiso_tv[i] = uv_init(i, 1);

    } else {
      solver->ja_.push_back(i);
      solver->a_.push_back(static_cast<int>(VV_tmp[i].size()));
      std::vector<int> row_id;
      row_id.reserve(VV_tmp[i].size());
      double bu = 0.0;
      double bv = 0.0;
      for (auto &vv_id : VV_tmp[i]) {
        if (bound_ids.count(vv_id) > 0) {
          bu += uv_init(vv_id, 0);
          bv += uv_init(vv_id, 1);
        } else {
          if (vv_id > i) {
            row_id.push_back(vv_id);
          }
        }
      }
      sort(row_id.begin(), row_id.end(), std::less<int>());
      for (std::size_t j = 0; j < row_id.size(); j++) {
        solver->ja_.push_back(row_id[j]);
        solver->a_.push_back(-1.0);
      }
      pardiso_tu[i] = bu;
      pardiso_tv[i] = bv;
    }
  }
  solver->ia_.push_back(static_cast<int>(solver->ja_.size()));

  solver->nnz_ = static_cast<int>(solver->ja_.size());
  solver->num_ = V_N;

  solver->pardiso_init();
  solver->rhs_ = pardiso_tu;

  solver->factorize();
  solver->pardiso_solver();

  for (int i = 0; i < V_N; i++) {
    uv_init(i, 0) = solver->result_[i];
  }

  solver->rhs_ = pardiso_tv;
  solver->pardiso_solver();
  for (int i = 0; i < V_N; i++) {
    uv_init(i, 1) = solver->result_[i];
  }
}

void PreCalcPardiso(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                    std::unique_ptr<Solver> &solver) {
  int V_N = static_cast<int>(V.rows());
  int F_N = static_cast<int>(F.rows());
  solver->ia_.clear();
  solver->ia_.reserve(2 * V_N + 1);
  solver->ja_.clear();
  solver->ja_.reserve(8 * V_N);

  std::vector<std::set<int>> VV_tmp;
  VV_tmp.resize(V_N);
  for (int i = 0; i < F_N; i++) {
    int vid[3];

    for (int j = 0; j < F.cols(); j++) {
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
    solver->ia_.push_back(static_cast<int>(solver->ja_.size()));
    VV_tmp[i].insert(i);
    std::vector<int> row_id;
    for (auto &var : VV_tmp[i]) {
      row_id.push_back(var);
    }

    std::vector<int>::iterator iter =
        std::find(row_id.begin(), row_id.end(), i);

    for (std::size_t k = std::distance(row_id.begin(), iter); k < row_id.size();
         k++) {
      solver->ja_.push_back(row_id[k]);
    }
    for (std::size_t k = 0; k < row_id.size(); k++) {
      solver->ja_.push_back(row_id[k] + V_N);
    }
  }
  for (int i = V_N; i < 2 * V_N; i++) {
    solver->ia_.push_back(static_cast<int>(solver->ja_.size()));
    std::vector<int> row_id;
    for (auto &var : VV_tmp[i - V_N]) {
      row_id.push_back(var);
    }
    std::vector<int>::iterator iter =
        std::find(row_id.begin(), row_id.end(), i - V_N);

    for (std::size_t k = std::distance(row_id.begin(), iter); k < row_id.size();
         k++) {
      solver->ja_.push_back(row_id[k] + V_N);
    }
  }
  solver->ia_.push_back(static_cast<int>(solver->ja_.size()));
}

void BoundaryLoop(const Eigen::MatrixXi &F_ref,
                  std::vector<std::vector<int>> &boundaryloop) {
  std::vector<std::vector<int>> boundaryEdges;
  std::vector<std::vector<int>> edges;
  int n_fvs = static_cast<int>(F_ref.cols());

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
  std::size_t i = 1;
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

  std::vector<int> visited;
  visited.resize(boundaryEdges.size(), 0);
  visited[0] = 1;
  std::vector<int> loop0;
  loop0.push_back(ev1);
  while (ev1 != ev0) {
    for (std::size_t j = 1; j < boundaryEdges.size(); j++) {
      if (visited[j] == 1)
        continue;
      if (boundaryEdges[j][0] == ev1) {
        visited[j] = 1;
        ev1 = boundaryEdges[j][1];
        loop0.push_back(ev1);
        break;
      }
    }
  }
  boundaryloop.emplace_back(loop0);
}

void MapVerticesToCircle(const Eigen::MatrixXd &V,   //
                         const Eigen::VectorXi &bnd, //
                         Eigen::MatrixXd &UV) {
  // Get sorted list of boundary vertices
  std::vector<int> interior, map_ij;
  map_ij.resize(V.rows());

  std::vector<bool> isOnBnd(V.rows(), false);
  for (int i = 0; i < bnd.size(); i++) {
    isOnBnd[bnd[i]] = true;
    map_ij[bnd[i]] = i;
  }

  for (std::size_t i = 0; i < isOnBnd.size(); i++) {
    if (!isOnBnd[i]) {
      map_ij[i] = static_cast<int>(interior.size());
      interior.push_back(static_cast<int>(i));
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
