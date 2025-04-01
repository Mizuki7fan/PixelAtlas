#include "CGAL-Interface.h"

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/centroid.h>

#include <frame/assert.hpp>
namespace PMP = CGAL::Polygon_mesh_processing;

namespace cgl {
std::array<CGAL::SM_Vertex_index, 3>
MeshOperation::GetFaceVertices(const cgl::SurfaceMesh3 &mesh, //
                               const CGAL::SM_Face_index &face) {
  std::array<CGAL::SM_Vertex_index, 3> face_vertex;
  CGAL::SM_Halfedge_index halfedge = mesh.halfedge(face);
  face_vertex[0] = mesh.target(halfedge);
  face_vertex[1] = mesh.target(mesh.next(halfedge));
  face_vertex[2] = mesh.target(mesh.next(mesh.next(halfedge)));
  return face_vertex;
}

std::array<CGAL::SM_Face_index, 3>
MeshOperation::GetNeighbourFaces(const cgl::SurfaceMesh3 &mesh, //
                                 const CGAL::SM_Face_index &face) {
  std::array<CGAL::SM_Face_index, 3> neighbour_faces;
  CGAL::SM_Halfedge_index halfedge = mesh.halfedge(face);
  CGAL::SM_Halfedge_index next_halfedge = mesh.next(halfedge);
  CGAL::SM_Halfedge_index next_next_halfedge = mesh.next(next_halfedge);

  neighbour_faces[0] = mesh.face(mesh.opposite(halfedge));
  neighbour_faces[1] = mesh.face(mesh.opposite(next_halfedge));
  neighbour_faces[2] = mesh.face(mesh.opposite(next_next_halfedge));

  return neighbour_faces;
}

CGAL::SM_Edge_index
MeshOperation::GetEdgeByTwoVertex(const cgl::SurfaceMesh3 &mesh,        //
                                  const CGAL::SM_Vertex_index vertex_0, //
                                  const CGAL::SM_Vertex_index vertex_1) {
  // 遍历与 vertex_0 相邻的所有半边
  for (auto halfedge :
       CGAL::halfedges_around_target(mesh.halfedge(vertex_0), mesh)) {
    // 获取半边的源顶点
    auto source_vertex = mesh.source(halfedge);
    // 检查源顶点是否是 vertex_1
    if (source_vertex == vertex_1) {
      // 返回对应的边
      return mesh.edge(halfedge);
    }
  }
  // 如果没有找到，则返回 null_edge()
  return mesh.null_edge();
}

CGAL::SM_Face_index
MeshOperation::GetFaceByThreeVertex(const cgl::SurfaceMesh3 &mesh,        //
                                    const CGAL::SM_Vertex_index vertex_0, //
                                    const CGAL::SM_Vertex_index vertex_1, //
                                    const CGAL::SM_Vertex_index vertex_2) {
  for (auto halfedge : mesh.halfedges_around_target(mesh.halfedge(vertex_0))) {
    CGAL::SM_Halfedge_index next_next_halfedge = mesh.next(mesh.next(halfedge));
    if (mesh.source(next_next_halfedge) == vertex_1 &&
        mesh.target(next_next_halfedge) == vertex_2)
      return mesh.face(halfedge);
    if (mesh.source(next_next_halfedge) == vertex_2 &&
        mesh.target(next_next_halfedge) == vertex_1)
      return mesh.face(halfedge);
  }
  return mesh.null_face();
}

bool MeshOperation::IO::loadOFF(const char *filename, SurfaceMesh3 &Mesh) {
  CGAL::IO::read_OFF(filename, Mesh);
  return true;
}

bool MeshOperation::IO::loadOBJ(const char *filename, SurfaceMesh3 &Mesh) {
  bool load_OK = CGAL::IO::read_OBJ(filename, Mesh);
  PA_ASSERT(load_OK);
  return true;
}

bool MeshOperation::IO::loadVTK(const char *filename, SurfaceMesh3 &Mesh) {

  Mesh.clear();
  std::ifstream ifs(filename);
  if (!ifs.is_open())
    return false;
  std::string line;
  std::istringstream iss;
  std::string flag;

  std::vector<cgl::Point3II> points;              // 顶点位置
  std::vector<std::array<unsigned int, 3>> faces; // 面

  while (std::getline(ifs, line)) {
    iss.clear();
    iss.str(line);
    iss >> flag;

    if (flag == "POINTS") {
      std::size_t num_points;
      iss >> num_points;
      points.resize(num_points);

      for (std::size_t i = 0; i < num_points; ++i) {
        std::getline(ifs, line);
        iss.clear();
        iss.str(line);
        double x, y, z;
        iss >> x >> y >> z;
        points[i] = cgl::Point3II(x, y, z);
      }
    } else if (flag == "CELLS") {
      std::size_t num_faces;
      iss >> num_faces;
      faces.resize(num_faces);
      for (std::size_t i = 0; i < num_faces; ++i) {
        std::getline(ifs, line);
        iss.clear();
        iss.str(line);
        unsigned int v0, v1, v2;
        int num_face_vertices = -1;
        iss >> num_face_vertices >> v0 >> v1 >> v2;
        faces[i] = {{v0, v1, v2}};
      }
    }
  }

  for (const auto &pt : points)
    Mesh.add_vertex(pt);
  for (auto face_vertex_idx : faces)
    Mesh.add_face(CGAL::SM_Vertex_index(face_vertex_idx[0]),
                  CGAL::SM_Vertex_index(face_vertex_idx[1]),
                  CGAL::SM_Vertex_index(face_vertex_idx[2]));

  auto mesh_full_path =
      Mesh.add_property_map<CGAL::SM_Vertex_index, std::string>("m:filepath")
          .first;
  mesh_full_path[CGAL::SM_Vertex_index(0)] = filename;
  return true;
}
} // namespace cgl
