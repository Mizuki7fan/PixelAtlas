#include "CGALMeshInterface.h"

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/centroid.h>

#include <frame/assert.hpp>
namespace PMP = CGAL::Polygon_mesh_processing;

namespace cgl {
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
