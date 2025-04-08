#pragma once
#include <Eigen/Core>
#include <triangle.h>
// triangle库已经有了名为triangle()的函数, 这里避免重名,
// 改为GenerateTriangulate()
void GenerateTriangulate(const Eigen::MatrixXd &bnd_pts,
                         const double &area_threshold, Eigen::MatrixXd &pts,
                         Eigen::MatrixXi &FV);
void GenerateTriangulate(const Eigen::MatrixXd &bnd_pts,
                         const Eigen::MatrixXi &E, const Eigen::MatrixXd &hole,
                         Eigen::MatrixXd &pts, Eigen::MatrixXi &FV);