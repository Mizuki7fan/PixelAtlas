#include "BoundaryLoop.h"

namespace AlgoKit {
void GetBoundaryLoop(const Eigen::MatrixXi &face_vertices, //
                     std::vector<std::vector<int>> &boundary_loop) {
  std::vector<std::vector<int>> boundary_edges;
  std::vector<std::vector<int>> all_edges;
  // 每个面含有几个顶点
  constexpr int kNumFaceVertices = 3;

  for (int face_idx = 0; face_idx < face_vertices.rows(); face_idx++) {
    for (int i = 0; i < kNumFaceVertices; i++) {
      int vertex_idx = face_vertices(face_idx, i);
      int next_vertex_idx = face_vertices(face_idx, (i + 1) % kNumFaceVertices);
      if (vertex_idx > next_vertex_idx)
        std::swap(vertex_idx, next_vertex_idx);
      std::vector<int> edge(4);
      edge[0] = vertex_idx;
      edge[1] = next_vertex_idx;
      edge[2] = face_idx;
      edge[3] = i;
      all_edges.emplace_back(edge);
    }
  }
  std::sort(all_edges.begin(), all_edges.end());
  int edge_idx = 1;
  for (; edge_idx < static_cast<int>(all_edges.size());) {
    auto &r1 = all_edges[edge_idx - 1];
    auto &r2 = all_edges[edge_idx];
    if ((r1[0] == r2[0]) && (r1[1] == r2[1])) {
      edge_idx += 2;
    } else {
      boundary_edges.emplace_back(all_edges[edge_idx - 1]);
      edge_idx++;
    }
  }
  if (edge_idx == static_cast<int>(all_edges.size()))
    boundary_edges.emplace_back(all_edges.back());

  for (auto &var : boundary_edges) {
    var[0] = face_vertices(var[2], var[3]);
    var[1] = face_vertices(var[2], (var[3] + 1) % kNumFaceVertices);
  }
  int ev0 = boundary_edges.front()[0];
  int ev1 = boundary_edges.front()[1];

  std::vector<int> visited(boundary_edges.size(), 0);
  visited[0] = 1;
  std::vector<int> loop0;
  loop0.push_back(ev1);
  while (ev1 != ev0) {
    for (std::size_t i = 1; i < boundary_edges.size(); i++) {
      if (visited[i] == 1)
        continue;
      if (boundary_edges[i][0] == ev1) {
        visited[i] = 1;
        ev1 = boundary_edges[i][1];
        loop0.push_back(ev1);
        break;
      }
    }
  }
  boundary_loop.emplace_back(loop0);
}
} // namespace AlgoKit