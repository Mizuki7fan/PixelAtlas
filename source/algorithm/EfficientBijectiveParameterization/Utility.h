#pragma once
#include <Eigen/Core>

using namespace std;

void Tutte(const int &V_N,                //
           const Eigen::MatrixXi &F,      //
           const Eigen::VectorXi &bnd,    //
           const Eigen::MatrixXd &bnd_uv, //
           Eigen::MatrixXd &uv_init);

void MapVerticesToCircle(const Eigen::MatrixXd &V,   //
                         const Eigen::VectorXi &bnd, //
                         Eigen::MatrixXd &UV);