#include "MeshDijkstraWithCache.h"

MeshDijkstraWithCache::MeshDijkstraWithCache(const cgl::SurfaceMesh3 &mesh)
    : mesh_(mesh) {
  vertex_distance_.resize(mesh_.num_vertices());
  vertex_is_visited_.resize(mesh.num_vertices());
  vertex_dijkstra_que_.resize(mesh.num_vertices());

  VV.resize(mesh.num_vertices());
  VE.resize(mesh.num_vertices());
  VVp.resize(mesh.num_vertices());

  for (auto vertex : mesh.vertices()) {
    VV[vertex.idx()].reserve(6);
    VE[vertex.idx()].reserve(6);
    for (auto halfedge : mesh.halfedges_around_target(mesh.halfedge(vertex))) {
      VV[vertex.idx()].push_back(mesh.source(halfedge).idx());
      VE[vertex.idx()].push_back(mesh.edge(halfedge).idx());
    }
  }

  EL.resize(mesh.num_edges());
  for (auto edge : mesh.edges()) {
    auto halfedge = mesh.halfedge(edge);
    auto vertex_0 = mesh.source(halfedge);
    auto vertex_1 = mesh.target(halfedge);
    EL[edge.idx()] = CGAL::sqrt(
        (mesh.point(vertex_0) - mesh.point(vertex_1)).squared_length());
  }
}

std::pair<double, bool> MeshDijkstraWithCache::GetVertexVertexDistnace(
    const CGAL::SM_Vertex_index vertex_0,
    const CGAL::SM_Vertex_index vertex_1) {
  int idx_0, idx_1;
  if (vertex_0.idx() < vertex_1.idx()) {
    idx_0 = vertex_0.idx();
    idx_1 = vertex_1.idx();
  } else {
    idx_0 = vertex_1.idx();
    idx_1 = vertex_0.idx();
  }

  // 如果是已经计算获得的值, 则直接返回
  if (vertex_is_visited_[idx_0].size() != 0)
    if (vertex_is_visited_[idx_0][idx_1] == 1) {
      //  std::cout << std::format("访问缓存命中: {} -> {}\n", idx_0, idx_1);
      return {vertex_distance_[idx_0][idx_1], true};
    }

  std::vector<int> &is_visited = vertex_is_visited_[idx_0];
  std::vector<double> &distance = vertex_distance_[idx_0];
  std::vector<int> &vertex_prev = VVp[idx_0];
  std::priority_queue<DijkstraNode> &que = vertex_dijkstra_que_[idx_0];
  if (is_visited.size() == 0) {
    // 首次对idx_0展开搜索
    is_visited.resize(mesh_.num_vertices(), 0);
    distance.resize(mesh_.num_vertices(), DBL_MAX);
    distance[idx_0] = 0.0;
    que.push(DijkstraNode(idx_0, 0.0));
    vertex_prev.resize(mesh_.num_vertices(), -1);
    vertex_prev[idx_0] = idx_0;
  }

  while (!que.empty()) {
    DijkstraNode node = std::move(que.top());
    que.pop();
    // 这里缓存机制有问题, 如果输入的idx_0和idx_1是相同的, 则que就会为0
    if (is_visited[node.idx] == 1)
      continue;
    is_visited[node.idx] = 1;
    for (std::size_t u = 0; u < VV[node.idx].size(); ++u) {
      int &vertex_idx = VV[node.idx][u];
      int &edge_idx = VE[node.idx][u];
      if (is_visited[vertex_idx])
        continue;
      if (distance[node.idx] + EL[edge_idx] < distance[vertex_idx]) {
        distance[vertex_idx] = distance[node.idx] + EL[edge_idx];
        que.push(DijkstraNode(vertex_idx, distance[vertex_idx]));
        vertex_prev[vertex_idx] = node.idx;
      }
    }
    if (node.idx == idx_1)
      return {distance[node.idx], false};
  }
  return {-1, false};
}

std::pair<double, std::vector<CGAL::SM_Halfedge_index>>
MeshDijkstraWithCache::GetVertexPath(const CGAL::SM_Vertex_index &vertex_0,
                                     const CGAL::SM_Vertex_index &vertex_1) {
  GetVertexVertexDistnace(vertex_0, vertex_1);
  int idx_0, idx_1;
  if (vertex_0.idx() < vertex_1.idx()) {
    idx_0 = vertex_0.idx();
    idx_1 = vertex_1.idx();
  } else {
    idx_0 = vertex_1.idx();
    idx_1 = vertex_0.idx();
  }

  std::vector<CGAL::SM_Halfedge_index> path;
  path.reserve(mesh_.num_edges());

  std::vector<int> &v_p = VVp[idx_0];
  int curr_idx = idx_1;

  while (v_p[curr_idx] != curr_idx) {
    int prev_idx = v_p[curr_idx];
    for (auto halfedge : mesh_.halfedges_around_target(
             mesh_.halfedge(CGAL::SM_Vertex_index(curr_idx)))) {
      if (mesh_.source(halfedge).idx() == static_cast<std::size_t>(prev_idx)) {
        path.push_back(halfedge);
        break;
      }
    }
    curr_idx = prev_idx;
  }
  return {vertex_distance_[idx_0][idx_1], path};
}
