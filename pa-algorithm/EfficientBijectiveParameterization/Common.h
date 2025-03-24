#pragma once
#include "LinSysSolverInterface/Solver.h"
#include <Eigen/Core>
#include <set>

void OrderTheBoundary(std::vector<std::list<int>> &order_boundary,
                      const std::vector<int> &nextlocation);
bool GrowFromP(int p, std::set<int> &isused, std::list<int> &order_boundary,
               const std::vector<int> &nextlocation);
void Tutte(int V_N, const Eigen::MatrixXi &F, const Eigen::VectorXi &bnd,
           const Eigen::MatrixXd &bnd_uv, Eigen::MatrixXd &uv_init);
void PreCalcPardiso(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                    std::unique_ptr<Solver> &solver);
void BoundaryLoop(const Eigen::MatrixXi &F_ref,
                  std::vector<std::vector<int>> &boundaryEdges);
void MapVerticesToCircle(const Eigen::MatrixXd &V, const Eigen::VectorXi &bnd,
                         Eigen::MatrixXd &UV);