#include "process.h"
#include "naive_pixelator.h"
#include <AlgoKit/MeshDijkstraWithCache.h>
#include <CGAL-Interface/CGAL-Interface.h>
#include <frame/global_args.h>
#include <frame/io.h>
#include <frame/metric.h>

namespace fs = std::filesystem;
using GA = frm::GlobalArguments;
void MainProcess() {
  fs::path instance_path = global::InstanceFullPath();
  // 读取超参数grid_size
  int grid_size =
      std::get<int>(global::ActionHyperParameters().at("grid_size"));
  cgl::SurfaceMesh3 uv_mesh;
  std::ifstream uv_mesh_file = frm::GetInputFilestream("uv_mesh.obj");
  CGAL::IO::read_OBJ(uv_mesh_file, uv_mesh);

  MeshDijkstraWithCache mesh_cache(uv_mesh);

  std::size_t num_hit = 0, num_total = 0;
  for (std::size_t i = 0; i < uv_mesh.num_vertices() - 1; ++i) {
    std::cout << "i: " << i << std::endl;
    for (std::size_t j = i + 1; j < uv_mesh.num_vertices(); ++j) {
      std::pair<double, bool> hit = mesh_cache.GetVertexVertexDistnace(
          CGAL::SM_Vertex_index(i), CGAL::SM_Vertex_index(j));
      if (hit.second)
        num_hit++;
      num_total++;
    }
  }

  std::cout << std::format("{}/{}命中率: {:.2f}%\n", num_hit, num_total,
                           100.0 * num_hit / num_total);

  NaivePixelator pixelator(uv_mesh, grid_size);
  pixelator.run();
}