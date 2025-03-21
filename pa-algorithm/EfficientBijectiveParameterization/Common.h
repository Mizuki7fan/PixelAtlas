#pragma once
#include "LinSysSolverInterface/Solver.h"
#include <Eigen/Core>
#include <list>
#include <math.h>
#include <set>
#include <vector>

void orderTheBoundary(std::vector<std::list<int>> &order_boundary,
                      const std::vector<int> &nextlocation);
bool growFromP(int p, std::set<int> &isused, std::list<int> &order_boundary,
               const std::vector<int> &nextlocation);

void map_vertices_to_circle(const Eigen::MatrixXd &V,
                            const Eigen::VectorXi &bnd, Eigen::MatrixXd &UV);

void Tutte(std::size_t V_N, const Eigen::MatrixXi &F,
           const Eigen::VectorXi &bnd, const Eigen::MatrixXd &bnd_uv,
           Eigen::MatrixXd &uv_init);
void preCalc_pardiso(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                     std::unique_ptr<Solver> &solver);
void writeObj(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
              const std::string &outfile);
void boundary_loop(const Eigen::MatrixXi &F_ref,
                   std::vector<std::vector<int>> &boundaryEdges);