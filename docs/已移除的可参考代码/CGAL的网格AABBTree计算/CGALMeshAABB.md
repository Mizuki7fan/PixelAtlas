#include <CGAL/Polygon_mesh_processing/locate.h>
using AABBFaceGraphPrimitive =
    CGAL::AABB_face_graph_triangle_primitive<cgl::SurfaceMesh3>;
using AABBFaceGraphTraits =
    CGAL::AABB_traits_3<cgl::KernelII, AABBFaceGraphPrimitive>;
using AABBTriTree = CGAL::AABB_tree<AABBFaceGraphTraits>;
  AABBTriTree aabb_uv_mesh_;

    CGAL::Polygon_mesh_processing::build_AABB_tree(uv_mesh_, aabb_uv_mesh_);
  aabb_uv_mesh_.accelerate_distance_queries();


  for (std::size_t i = 0; i < grid_->V.size(); ++i) {
    std::array<double, 2> &coord = grid_->V[i].coord;
    cgl::Point3II query(coord[0], coord[1], 0);
    auto location = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(
        query, aabb_uv_mesh_, uv_mesh_);
    CGAL::SM_Face_index face = location.first;