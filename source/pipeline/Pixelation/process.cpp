#include "process.h"
#include "naive_pixelator.h"
#include <frame/global_defs.h>
#include <frame/io.h>
#include <frame/metric.h>

namespace fs = std::filesystem;
void MainProcess() {
  std::string work_name = frm::global::GlobalArgs().work_name;
  std::cout << work_name << std::endl;
  fs::path instance_path = frm::global::InstancePath();
}