#include "process.h"
#include <CGAL-Interface/CGAL-Interface.h>
#include <EfficientBijectiveParameterization/BiljectivePara.h>
#include <filesystem>
#include <frame/global_args.h>
#include <frame/io.h>
#include <frame/metric.h>


namespace fs = std::filesystem;
using GA = frm::GlobalArguments;

void MainProcess() {
  fs::path instance_path = global::InstanceFullPath();
  cgl::SurfaceMesh3 tri_mesh;
  CGAL::IO::read_OBJ(instance_path.string(), tri_mesh);
  std::unique_ptr<BiljectivePara> bil_para =
      std::make_unique<BiljectivePara>(tri_mesh);

  bil_para->Load(); // 读入网格、完成shell的构建
  bil_para->Parameterization();

  std::ofstream fout = frm::CreateOutputFilestream("uv_mesh.obj");
  bil_para->WriteUVMesh(fout);
  fout.close();

  // 写入metrics
  std::unordered_map<std::string, std::string> metrics;
  metrics["distortion"] = "double";
  metrics["length"] = "double";

  std::unordered_map<std::string, frm::ValueType> metric_values;
  metric_values["distortion"] = bil_para->GetDistortion();
  metric_values["length"] = 4.00;

  fout = frm::CreateMetricsFilestream();
  frm::WriteMetricJsonFile(fout, metrics, metric_values);
  fout.close();
}