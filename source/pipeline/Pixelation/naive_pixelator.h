#pragma once
#include <CGAL-Interface/CGAL-Interface.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <HierarchicalPixelGrid/HierarchicalPixelGrid.h>

using AABBFaceGraphPrimitive =
    CGAL::AABB_face_graph_triangle_primitive<cgl::SurfaceMesh3>;
using AABBFaceGraphTraits =
    CGAL::AABB_traits_3<cgl::KernelII, AABBFaceGraphPrimitive>;
using AABBTriTree = CGAL::AABB_tree<AABBFaceGraphTraits>;

class NaivePixelator {
public:
  NaivePixelator(const cgl::SurfaceMesh3 &uv_mesh, int grid_size);
  void run();

private:
  std::unique_ptr<HierarchicalPixelGrid> grid_ = nullptr;
  const cgl::SurfaceMesh3 &uv_mesh_;
  AABBTriTree aabb_uv_mesh_;

  const int grid_size_;
};