#include "naive_pixelator.h"

NaivePixelator::NaivePixelator(const cgl::SurfaceMesh3 &uv_mesh, int grid_size)
    : uv_mesh(uv_mesh), grid_size(grid_size) {};

void NaivePixelator::run() {
  // create_grid();
  grid = std::make_unique<HierarchicalPixelGrid>(grid_size);
}
