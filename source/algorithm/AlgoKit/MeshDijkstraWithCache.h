#pragma once
#include <CGAL-Interface/CGAL-Interface.h>
#include <queue>

struct DijkstraNode {
  int idx;
  double distance;
  DijkstraNode(int idx, double distance) : idx(idx), distance(distance) {}
  bool operator<(const DijkstraNode &other) const {
    return distance > other.distance;
  }
};

class MeshDijkstraWithCache {
public:
  explicit MeshDijkstraWithCache(const cgl::SurfaceMesh3 &mesh);

  // 算这两个vertex之间的距离
  std::pair<double, bool>
  GetVertexVertexDistnace(const CGAL::SM_Vertex_index vertex_0,
                          const CGAL::SM_Vertex_index vertex_1);
  std::pair<double, std::vector<CGAL::SM_Halfedge_index>>
  GetVertexPath(const CGAL::SM_Vertex_index &vertex_0,
                const CGAL::SM_Vertex_index &vertex_1);

private:
  const cgl::SurfaceMesh3 &mesh_;
  // 记录缓存
  std::vector<std::vector<double>> vertex_distance_;
  std::vector<std::vector<int>> vertex_is_visited_;
  std::vector<std::priority_queue<DijkstraNode>> vertex_dijkstra_que_;

private:
  std::vector<std::vector<int>> VV;
  std::vector<std::vector<int>> VE;
  std::vector<std::vector<int>> VVp;
  std::vector<double> EL;
};
