#pragma once
#include <CGAL-Interface/CGAL-Interface.h>
class NaivePixelator {
public:
  NaivePixelator(const cgl::SurfaceMesh3 &uv_mesh, int grid_size);
  void run();

private:
  const cgl::SurfaceMesh3 &uv_mesh;
  const int grid_size;
};