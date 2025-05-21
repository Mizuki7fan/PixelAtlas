#pragma once
#include <Eigen/Core>

namespace AlgoKit {
void GetBoundaryLoop(const Eigen::MatrixXi &F_ref, //
                     std::vector<std::vector<int>> &boundaryEdges);
}
