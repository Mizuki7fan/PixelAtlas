#pragma once

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
namespace cgl {

using KernelII = CGAL::Simple_cartesian<double>;
typedef KernelII::Point_3 Point3II;
typedef KernelII::Segment_3 Seg3II;
typedef KernelII::Vector_3 Vec3II;
typedef KernelII::Plane_3 Plane3II;
typedef KernelII::Line_3 Line3II;
typedef KernelII::Triangle_3 Triangle3II;

typedef CGAL::Surface_mesh<Point3II> SurfaceMesh3;

class MeshOperation {
public:
  class IO {
  public:
    static bool loadOFF(const char *filename, SurfaceMesh3 &Mesh);
    static bool loadOBJ(const char *filename, SurfaceMesh3 &Mesh);
    static bool loadVTK(const char *filename, SurfaceMesh3 &Mesh);
  };
};
}; // namespace cgl
