#pragma once
#include <array>
#include <math.h>

// 基础几何
// 点线关系/点面关系/计算重心坐标/四元数等等
namespace AlgoKit {
int Orient2D(const double &px, const double &py,   //
             const double &p0x, const double &p0y, //
             const double &p1x, const double &p1y);

bool CheckPointInTriangle2D(const double &px, const double &py,   //
                            const double &p0x, const double &p0y, //
                            const double &p1x, const double &p1y, //
                            const double &p2x, const double &p2y);

// 检测box和三角形的相交
bool CheckBoxTriangleIntersection(const std::array<double, 2> &bbmin,
                                  const std::array<double, 2> &bbmax,
                                  const std::array<double, 2> &v0, //
                                  const std::array<double, 2> &v1, //
                                  const std::array<double, 2> &v2);
}; // namespace AlgoKit