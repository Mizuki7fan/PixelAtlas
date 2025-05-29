#include "process.h"
#include <frame/global_args.h>
#include <frame/io.h>
#include <iostream>

namespace fs = std::filesystem;

void MainProcess() {
  fs::path instance_path = global::InstanceFullPath();
  std::cout << instance_path << std::endl;
}