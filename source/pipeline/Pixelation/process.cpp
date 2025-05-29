#include "process.h"
#include "naive_pixelator.h"
#include <CGAL-Interface/CGAL-Interface.h>
#include <frame/global_args.h>
#include <frame/io.h>
#include <frame/metric.h>

namespace fs = std::filesystem;
void MainProcess() {
  fs::path instance_path = global::InstanceFullPath();
  // 读取超参数grid_size
  int grid_size =
      std::get<int>(global::ActionHyperParameters().at("grid_size"));
  cgl::SurfaceMesh3 uv_mesh;
  std::ifstream uv_mesh_file = frm::GetInputFilestream("uv_mesh.obj");
  CGAL::IO::read_OBJ(uv_mesh_file, uv_mesh);

  NaivePixelator pixelator(uv_mesh, grid_size);
  pixelator.run();

  std::ofstream grid_file = frm::CreateOutputFilestream("grid.grid");
  pixelator.write_grid(grid_file);
  grid_file.close();
}