#include "naive_pixelator.h"
#include "frame/global_args.h"
#include "frame/io.h"
#include <AlgoKit/BasicGeometry.h>

NaivePixelator::NaivePixelator(const cgl::SurfaceMesh3 &uv_mesh, int grid_size)
    : uv_mesh_(uv_mesh), grid_size_(grid_size) {
  // 构建uv_mesh的AABB
  CGAL::Polygon_mesh_processing::build_AABB_tree(uv_mesh_, aabb_uv_mesh_);
  aabb_uv_mesh_.accelerate_distance_queries();
  // 可能不需要使用AABB数据结构, 而是算每个三角面片所影响的grid格点,
  // 这样效率更高, 以后可以替换
};

void NaivePixelator::run() {
  // create_grid();
  grid_ = std::make_unique<HierarchicalPixelGrid>(grid_size_);
  std::vector<uint8_t> grid_vertex_state(grid_->V.size(), 0);

  if (global::DebugLevel() > 0) {
    std::ofstream filestream = frm::CreateDebugFilestream("uv_mesh.obj");
    CGAL::IO::write_OBJ(filestream, uv_mesh_);
    filestream.close();
  }
  for (std::size_t i = 0; i < grid_->V.size(); ++i) {
    std::array<double, 2> &coord = grid_->V[i].coord;
    cgl::Point3II query(coord[0], coord[1], 0);
    auto location = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(
        query, aabb_uv_mesh_, uv_mesh_);
    CGAL::SM_Face_index face = location.first;
    CGAL::SM_Vertex_index face_vertex[3];
    CGAL::SM_Halfedge_index halfedge = uv_mesh_.halfedge(face);
    face_vertex[0] = uv_mesh_.target(halfedge);
    face_vertex[1] = uv_mesh_.target(uv_mesh_.next(halfedge));
    face_vertex[2] = uv_mesh_.target(uv_mesh_.next(uv_mesh_.next(halfedge)));
    // 判定点是否在面内
    bool is_inside = AlgoKit::CheckPointInTriangle2D(
        coord[0], coord[1], //
        uv_mesh_.point(face_vertex[0]).x(),
        uv_mesh_.point(face_vertex[0]).y(), //
        uv_mesh_.point(face_vertex[1]).x(),
        uv_mesh_.point(face_vertex[1]).y(), //
        uv_mesh_.point(face_vertex[2]).x(), uv_mesh_.point(face_vertex[2]).y());
    if (is_inside)
      grid_->V[i].state = GridVertex::STATE::in_mesh;
    else
      grid_->V[i].state = GridVertex::STATE::outof_mesh;
  }

  for (std::size_t i = 0; i < grid_->F.size(); ++i) {
    GridFace &face = grid_->F[i];
    if (face.ldVertex->state == GridVertex::STATE::in_mesh ||
        face.luVertex->state == GridVertex::STATE::in_mesh ||
        face.rdVertex->state == GridVertex::STATE::in_mesh ||
        face.ruVertex->state == GridVertex::STATE::in_mesh) {
      face.state = GridFace::STATE::in_mesh;
    }
  }

  // 还有一种情况:对于部分非常狭长的三角形, 可能刺过某个grid格点,
  if (global::DebugLevel() > 0) {
    std::ofstream filestream = frm::CreateDebugFilestream("grid.findvef");
    grid_->PrintFindVEF(filestream);
    filestream.close();

    filestream = frm::CreateDebugFilestream("grid.obj");
    grid_->PrintQuadMeshOBJ(filestream);
    filestream.close();
  }
}
