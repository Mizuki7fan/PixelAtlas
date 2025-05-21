#include "Utility.h"
#if defined(USE_MKL)
#include <LinSysSolver-Interface/MKLPardisoSolver.h>
#elif defined(USE_EIGEN)
#include <LinSysSolver-Interface/EigenLinSolver.h>
#endif
#include <numbers>
#include <unordered_set>

/*加入判断边界线是否闭合,return true means open curve; return false means closed
curve. In closed case, the first and last item of order_boundary are same.*/
void Tutte(const int &num_vertices,              //
           const Eigen::MatrixXi &face_vertices, // 三角面片顶点索引
           const Eigen::VectorXi &bnd,           // 边界顶点索引
           const Eigen::MatrixXd &bnd_uv,        // 边界顶点预设UV坐标
           Eigen::MatrixXd &uv_init) {           // 输出：所有顶点UV坐标
  std::size_t num_faces = static_cast<std::size_t>(face_vertices.rows());
  uv_init.resize(num_vertices, 2);

  std::unordered_set<int> bound_idxes;
  for (int i = 0; i < bnd.size(); i++) {
    bound_idxes.insert(bnd(i));
    uv_init.row(bnd(i)) << bnd_uv(i, 0), bnd_uv(i, 1);
  }

  std::vector<std::unordered_set<int>> VV_tmp(num_vertices);
  for (std::size_t i = 0; i < num_faces; i++) {
    VV_tmp[face_vertices(i, 0)].insert(face_vertices(i, 1));
    VV_tmp[face_vertices(i, 0)].insert(face_vertices(i, 2));

    VV_tmp[face_vertices(i, 1)].insert(face_vertices(i, 0));
    VV_tmp[face_vertices(i, 1)].insert(face_vertices(i, 2));

    VV_tmp[face_vertices(i, 2)].insert(face_vertices(i, 0));
    VV_tmp[face_vertices(i, 2)].insert(face_vertices(i, 1));
  }

  std::unique_ptr<Solver> solver = nullptr;

#if defined(USE_MKL)
  solver = std::unique_ptr<MKLPardisoSolver>();
#elif defined(USE_EIGEN)
  solver = std::make_unique<EigenLinSolver>();
#endif

  if (solver == nullptr) {
    throw std::runtime_error("未设置USE_MKL或者USE_EIGEN");
  }

  std::vector<double> solver_tu;
  std::vector<double> solver_tv;

  solver->ia_.reserve(num_vertices + 1);
  solver->ja_.reserve(8 * num_vertices);
  solver->a_.reserve(8 * num_vertices);
  solver_tu.resize(num_vertices, 0.0);
  solver_tv.resize(num_vertices, 0.0);

  for (int i = 0; i < num_vertices; i++) {
    solver->ia_.push_back(static_cast<int>(solver->ja_.size()));

    if (bound_idxes.count(i) > 0) {
      solver->ja_.push_back(i);
      solver->a_.push_back(1.0);

      solver_tu[i] = uv_init(i, 0);
      solver_tv[i] = uv_init(i, 1);

    } else {
      solver->ja_.push_back(i);
      solver->a_.push_back(static_cast<double>(VV_tmp[i].size()));
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
  solver->ia_.push_back(static_cast<int>(solver->ja_.size()));

  solver->nnz_ = static_cast<int>(solver->ja_.size());
  solver->num_ = num_vertices;

  solver->PardisoInit();
  solver->rhs_ = solver_tu;

  solver->Factorize();
  solver->PardisoSolver();

  for (int i = 0; i < num_vertices; i++)
    uv_init(i, 0) = solver->result_[i];

  solver->rhs_ = solver_tv;
  solver->PardisoSolver();
  for (int i = 0; i < num_vertices; i++)
    uv_init(i, 1) = solver->result_[i];
}

void MapVerticesToCircle(const Eigen::MatrixXd &V,   //
                         const Eigen::VectorXi &bnd, //
                         Eigen::MatrixXd &UV) {
  // Get sorted list of boundary vertices
  std::vector<int> interior, map_ij(V.rows());

  std::vector<bool> is_on_bnd(V.rows(), false);
  for (int i = 0; i < bnd.size(); i++) {
    is_on_bnd[bnd[i]] = true;
    map_ij[bnd[i]] = i;
  }

  for (std::size_t i = 0; i < is_on_bnd.size(); i++) {
    if (!is_on_bnd[i]) {
      map_ij[i] = static_cast<int>(interior.size());
      interior.push_back(static_cast<int>(i));
    }
  }

  // Map boundary to unit circle
  std::vector<double> len(bnd.size());
  len[0] = 0.;

  for (int i = 1; i < bnd.size(); i++)
    len[i] = len[i - 1] + (V.row(bnd[i - 1]) - V.row(bnd[i])).norm();

  double total_len =
      len[len.size() - 1] + (V.row(bnd[0]) - V.row(bnd[bnd.size() - 1])).norm();

  UV.resize(bnd.size(), 2);

  for (int i = 0; i < bnd.size(); i++) {
    double frac = len[i] * 2. * std::numbers::pi / total_len;
    UV.row(map_ij[bnd[i]]) << cos(frac), sin(frac);
  }
}