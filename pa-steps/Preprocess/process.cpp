#include "process.h"
#include "EfficientBijectiveParameterization/BiljectivePara.h"
#include <CGAL-Interface/CGAL-Interface.h>
#include <filesystem>
#include <frame/global_defs.h>
#include <frame/io.h>
#include <frame/metric.h>


namespace fs = std::filesystem;

void MainProcess() {
  fs::path model_path = frm::global::CurrFile();
  cgl::SurfaceMesh3 tri_mesh;
  CGAL::IO::read_OBJ(model_path.string(), tri_mesh);
  std::unique_ptr<BiljectivePara> bil_para =
      std::make_unique<BiljectivePara>(tri_mesh, model_path.string());

  bil_para->load(); // 读入网格、完成shell的构建
  bil_para->parameterization();

  std::ofstream fout = frm::CreateResultFilestream("uv_mesh.obj");
  bil_para->WriteUVMesh(fout);
  fout.close();

  // 写入metrics
  std::unordered_map<std::string, std::string> metrics;
  metrics["distortion"] = "double";
  std::unordered_map<std::string, frm::MetricValue> metric_values;
  metric_values["distortion"] = bil_para->GetDistortion();

  fout = frm::CreateMetricsFilestream();
  frm::WriteMetricJsonFile(fout, metrics, metric_values);
  fout.close();
}