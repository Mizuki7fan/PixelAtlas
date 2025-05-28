#include "naive_pixelator.h"
#include "frame/global_args.h"
#include "frame/io.h"
#include "frame/pa_assert.hpp"
#include <AlgoKit/IntersectionDetection2D.h>

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

  if (global::DebugLevel() > 0) {
    dbg_file = frm::CreateDebugFilestream(std::format("eh.findvef"));
  }

  for (CGAL::SM_Edge_index eh : uv_mesh_.edges()) {
    std::array<CGAL::SM_Vertex_index, 2> edge_vertex;
    edge_vertex[0] = uv_mesh_.source(uv_mesh_.halfedge(eh, 0));
    edge_vertex[1] = uv_mesh_.source(uv_mesh_.halfedge(eh, 1));

    // 根据grid_element0确定边所划过的点的坐标
    std::array<int, 2> bbMin = {grid_size_, grid_size_}, bbMax = {0, 0};
    for (int i = 0; i < 2; ++i) {
      GridElement &element = map_vtx_to_grid_element[edge_vertex[i].idx()];
      if (std::holds_alternative<GridFace>(element)) {
        GridFace &face = std::get<GridFace>(element);
        bbMin[0] = (bbMin[0] > face.ldVertex->X) ? face.ldVertex->X : bbMin[0];
        bbMin[1] = (bbMin[1] > face.ldVertex->Y) ? face.ldVertex->Y : bbMin[1];
        bbMax[0] = (bbMax[0] < face.ruVertex->X) ? face.ruVertex->X : bbMax[0];
        bbMax[1] = (bbMax[1] < face.ruVertex->Y) ? face.ruVertex->Y : bbMax[1];
      } else if (std::holds_alternative<GridEdge>(element)) {
        GridEdge &edge = std::get<GridEdge>(element);
        GridVertex vertex[2];
        if (edge.HalfEdge0->BeginVertex != NULL) {
          vertex[0] = *edge.HalfEdge0->BeginVertex;
          vertex[1] = *edge.HalfEdge0->EndVertex;
        } else {
          vertex[0] = *edge.HalfEdge1->BeginVertex;
          vertex[1] = *edge.HalfEdge1->EndVertex;
        }
        for (int i = 0; i < 2; ++i) {
          bbMin[0] = (bbMin[0] > vertex[i].X) ? vertex[i].X : bbMin[0];
          bbMin[1] = (bbMin[1] > vertex[i].Y) ? vertex[i].Y : bbMin[1];
          bbMax[0] = (bbMax[0] < vertex[i].X) ? vertex[i].X : bbMax[0];
          bbMax[1] = (bbMax[1] < vertex[i].Y) ? vertex[i].Y : bbMax[1];
        }
      } else if (std::holds_alternative<GridVertex>(element)) {
        GridVertex &vertex = std::get<GridVertex>(element);
        bbMin[0] = (bbMin[0] > vertex.X) ? vertex.X : bbMin[0];
        bbMin[1] = (bbMin[1] > vertex.Y) ? vertex.Y : bbMin[1];
        bbMax[0] = (bbMax[0] < vertex.X) ? vertex.X : bbMax[0];
        bbMax[1] = (bbMax[1] < vertex.Y) ? vertex.Y : bbMax[1];
      }
    }

    if (global::DebugLevel() > 0) {
      int area = (bbMax[0] - bbMin[0]) * (bbMax[1] - bbMin[1]);
      if (area > 4) {
        dbg_file << "PE" << std::endl;
        dbg_file << std::format("{} {} 0 {} {} 0",
                                uv_mesh_.point(edge_vertex[0]).x(),
                                uv_mesh_.point(edge_vertex[0]).y(),
                                uv_mesh_.point(edge_vertex[1]).x(),
                                uv_mesh_.point(edge_vertex[1]).y())
                 << std::endl;
        GridVertex vertex[4];
        vertex[0] = grid_->Vertex(bbMin[0], bbMin[1]);
        vertex[1] = grid_->Vertex(bbMax[0], bbMin[1]);
        vertex[2] = grid_->Vertex(bbMax[0], bbMax[1]);
        vertex[3] = grid_->Vertex(bbMin[0], bbMax[1]);
        for (int i = 0; i < 4; ++i) {
          dbg_file << std::format(
                          "{} {} 0 {} {} 0", //
                          vertex[i].coord[0] / static_cast<double>(grid_size_),
                          vertex[i].coord[1] / static_cast<double>(grid_size_),
                          vertex[(i + 1) % 4].coord[0] /
                              static_cast<double>(grid_size_),
                          vertex[(i + 1) % 4].coord[1] /
                              static_cast<double>(grid_size_))
                   << std::endl;
        }
      }
    }

    std::array<double, 2> coord_0 = {uv_mesh_.point(edge_vertex[0]).x(),
                                     uv_mesh_.point(edge_vertex[0]).y()};
    std::array<double, 2> coord_1 = {uv_mesh_.point(edge_vertex[1]).x(),
                                     uv_mesh_.point(edge_vertex[1]).y()};

    // 根据bbMin和bbMax确定grid的面, 然后比较该面与当前edge是否相交
    std::array<int, 2> range = {bbMax[0] - bbMin[0], bbMax[1] - bbMin[1]};
    for (int i = 0; i < range[0]; ++i)
      for (int j = 0; j < range[1]; ++j) {
        const GridVertex &vertex = grid_->Vertex(bbMin[0] + i, bbMin[0] + j);
        GridFace *face = vertex.ruFace;
        PA_ASSERT_WITH_MSG(face != NULL, "GridFace is NULL");
        if (face->state == GridFace::STATE::in_mesh)
          continue;
        if (AlgoKit::CheckBoxLineSegmentIntersection(
                {face->ldVertex->coord[0] / static_cast<double>(grid_size_),
                 face->ldVertex->coord[1] /
                     static_cast<double>(grid_size_)}, // bbmin
                {face->ruVertex->coord[0] / static_cast<double>(grid_size_),
                 face->ruVertex->coord[1] /
                     static_cast<double>(grid_size_)}, // bbmax
                coord_0, coord_1)) {
          face->state = GridFace::STATE::in_mesh;
        }
      };
  }

  if (global::DebugLevel() > 0) {
    dbg_file.close();
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
