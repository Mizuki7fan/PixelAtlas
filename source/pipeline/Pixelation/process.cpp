#include "process.h"
#include "naive_pixelator.h"
#include <frame/global_args.h>
#include <frame/io.h>
#include <frame/metric.h>

namespace fs = std::filesystem;
using GA = frm::GlobalArguments;
void MainProcess() {
  fs::path instance_path = global::InstanceFullPath();
  std::cout << "pixelation main process" << std::endl;
}