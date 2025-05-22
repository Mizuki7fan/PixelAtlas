#pragma once
#include <CGAL-Interface/CGAL-Interface.h>
#include <HierarchicalPixelGrid/HierarchicalPixelGrid.h>
class NaivePixelator {
public:
  NaivePixelator(const cgl::SurfaceMesh3 &uv_mesh, int grid_size);
  void run();

private:
  std::unique_ptr<HierarchicalPixelGrid> grid = nullptr;
  const cgl::SurfaceMesh3 &uv_mesh;
  const int grid_size;
};