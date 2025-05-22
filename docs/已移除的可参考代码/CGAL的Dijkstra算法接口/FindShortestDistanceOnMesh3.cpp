#include "FindShortestDistanceOnMesh3.h"
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/numeric.hpp>

DijkstraAllTargetReachGuardingVisitor::DijkstraAllTargetReachGuardingVisitor(
    const AnySinglePassRange<CGAL::SM_Vertex_index> &target_vtxs) {
  boost::range::insert(unreached_target_vtxs_, target_vtxs);
}

void DijkstraAllTargetReachGuardingVisitor::examine_vertex(
    CGAL::SM_Vertex_index vertex, const cgl::SurfaceMesh3 &) {
  if (unreached_target_vtxs_.empty())
    return;
  unreached_target_vtxs_.erase(vertex);
  if (unreached_target_vtxs_.empty())
    throw(AllTargetsReachedException());
}

void FindShortestPathsOnMesh3(
    const cgl::SurfaceMesh3 &mesh, CGAL::SM_Vertex_index source,
    const AnySinglePassRange<CGAL::SM_Vertex_index> &targets,
    std::vector<CGAL::SM_Vertex_index> &map_vtx_to_prev_vtx,
    std::vector<double> &map_vtx_to_distance) {

  DijkstraAllTargetReachGuardingVisitor targets_visitor(targets);

  map_vtx_to_prev_vtx.resize(mesh.num_vertices());
  map_vtx_to_distance.resize(mesh.num_vertices());

  try {
    boost::dijkstra_shortest_paths(
        mesh, source,
        boost::predecessor_map(PredecessorMap(map_vtx_to_prev_vtx.begin()))
            .distance_map(DistanceMap(map_vtx_to_distance.begin()))
            .visitor(targets_visitor));
  } catch (
      const DijkstraAllTargetReachGuardingVisitor::AllTargetsReachedException
          &) {
  }
}

std::pair<double, std::vector<CGAL::SM_Halfedge_index>>
FindShortestHalfedgePathOnMesh3(const cgl::SurfaceMesh3 &mesh,
                                cgl::SurfaceMesh3::vertex_index source,
                                cgl::SurfaceMesh3::vertex_index target) {

  std::vector<CGAL::SM_Vertex_index> map_vtx_to_prev_vtx;
  std::vector<double> map_vtx_to_distance;
  FindShortestPathsOnMesh3(mesh, source,
                           std::vector<cgl::SurfaceMesh3::vertex_index>{target},
                           map_vtx_to_prev_vtx, map_vtx_to_distance);

  if (map_vtx_to_prev_vtx[target.idx()] == target)
    return {}; // No path found
  auto curr_vtx = target;
  std::vector<cgl::SurfaceMesh3::halfedge_index> halfedges_of_shortest_path;
  while (curr_vtx != source) {
    auto prev_vtx = map_vtx_to_prev_vtx[curr_vtx.idx()];

    for (auto halfedge :
         mesh.halfedges_around_target(mesh.halfedge(curr_vtx))) {
      if (mesh.source(halfedge) == prev_vtx) {
        halfedges_of_shortest_path.push_back(halfedge);
        break;
      }
    }

    curr_vtx = prev_vtx;
  }

  return {map_vtx_to_distance[target.idx()], halfedges_of_shortest_path};
}