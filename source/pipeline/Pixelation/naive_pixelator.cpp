#include "naive_pixelator.h"
#include "frame/global_args.h"
#include "frame/io.h"
#include <AlgoKit/BasicGeometry.h>

NaivePixelator::NaivePixelator(const cgl::SurfaceMesh3 &uv_mesh, int grid_size)
    : uv_mesh_(uv_mesh), grid_size_(grid_size) {};

void NaivePixelator::run() {
  // create_grid();
  grid_ = std::make_unique<HierarchicalPixelGrid>(grid_size_);
  std::vector<uint8_t> grid_vertex_state(grid_->V.size(), 0);

  if (global::DebugLevel() > 0) {
    std::ofstream filestream = frm::CreateDebugFilestream("uv_mesh.obj");
    CGAL::IO::write_OBJ(filestream, uv_mesh_);
    filestream.close();
  }

  std::vector<GridElement> map_vtx_to_grid_element(uv_mesh_.num_vertices(),
                                                   std::monostate());

  for (CGAL::SM_Vertex_index vertex : uv_mesh_.vertices()) {
    std::array<double, 2> coord_ = {uv_mesh_.point(vertex).x(),
                                    uv_mesh_.point(vertex).y()};
    map_vtx_to_grid_element[vertex.idx()] = grid_->LocateGridElement(coord_);
  }

  for (CGAL::SM_Edge_index eh : uv_mesh_.edges()) {
    std::array<CGAL::SM_Vertex_index, 2> edge_vertex;
    edge_vertex[0] = uv_mesh_.source(uv_mesh_.halfedge(eh, 0));
    edge_vertex[1] = uv_mesh_.source(uv_mesh_.halfedge(eh, 1));

    // 根据grid_element0确定边所划过的面的上下界
    std::array<int, 2> bbMin = {grid_size_, grid_size_}, bbMax = {0, 0};
    auto UpdateEdgeGridBBox = [&](const GridFace *face) {
      if (face == NULL)
        return;
      bbMin[0] = (bbMin[0] > face->X_) ? face->X_ : bbMin[0];
      bbMin[1] = (bbMin[1] > face->Y_) ? face->Y_ : bbMin[1];
      bbMax[0] = (bbMax[0] < face->X_) ? face->X_ : bbMax[0];
      bbMax[1] = (bbMax[1] < face->Y_) ? face->Y_ : bbMax[1];
    };

    for (int i = 0; i < 2; ++i) {
      GridElement &element = map_vtx_to_grid_element[edge_vertex[i].idx()];
      if (std::holds_alternative<GridFace>(element)) {
        GridFace &face = std::get<GridFace>(element);
        UpdateEdgeGridBBox(&face);
      } else if (std::holds_alternative<GridEdge>(element)) {
        GridHalfEdge *lhalfedge = std::get<GridEdge>(element).HalfEdge0;
        UpdateEdgeGridBBox(lhalfedge->Face);
        GridHalfEdge *rhalfedge = std::get<GridEdge>(element).HalfEdge1;
        UpdateEdgeGridBBox(rhalfedge->Face);
      } else if (std::holds_alternative<GridVertex>(element)) {
        GridVertex &vertex = std::get<GridVertex>(element);
        UpdateEdgeGridBBox(vertex.ldFace);
        UpdateEdgeGridBBox(vertex.rdFace);
        UpdateEdgeGridBBox(vertex.luFace);
        UpdateEdgeGridBBox(vertex.ruFace);
      }
    }

    if (global::DebugLevel() > 0) {
      std::cout << std::format("bbox: [{} {}]x[{} {}]", bbMin[0], bbMin[1],
                               bbMax[0], bbMax[1])
                << std::endl;
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
