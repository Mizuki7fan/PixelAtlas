#pragma once

#include <Eigen/Core>
#include <LinSysSolver-Interface/Solver.h>

using namespace std;

void Tutte(const int &V_N,                //
           const Eigen::MatrixXi &F,      //
           const Eigen::VectorXi &bnd,    //
           const Eigen::MatrixXd &bnd_uv, //
           Eigen::MatrixXd &uv_init);

void MapVerticesToCircle(const Eigen::MatrixXd &V,   //
                         const Eigen::VectorXi &bnd, //
                         Eigen::MatrixXd &UV);

void preCalc_pardiso(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                     Solver &pardiso);
void boundary_loop(const Eigen::MatrixXi &F_ref,
                   std::vector<std::vector<int>> &boundaryEdges);