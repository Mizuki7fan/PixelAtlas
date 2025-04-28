#include "io.h"
#include "global_args.h"
#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

namespace frm {
using GA = GlobalArguments;

std::ofstream CreateOutputFilestream(const std::string &path) {
  // 生成输出的路径
  const fs::path &curr_result_dir = GA::I().ActionResultDir();
  std::string path_str;
  if (GA::I().UseIndividualInstanceDir()) {
    path_str = std::format("{}/{}", curr_result_dir.string(), path);
  } else {
    path_str =
        std::format("{}/{}-{}", curr_result_dir.string(),
                    GA::I().InstanceFullPath().filename().string(), path);
  }
  if (GA::I().DebugLevel() > 0) {
    std::cout << "CreateOutputFilestream: " << path_str << std::endl;
  }
  return std::ofstream(path_str);
}

std::ofstream CreateDebugFilestream(const std::string &path) {
  // 生成输出的路径
  return std::ofstream(GA::I().ActionDebugDir() / path);
}

std::ofstream CreateMetricsFilestream() {
  // 生成输出的路径
  std::string path_str;
  if (GA::I().UseIndividualInstanceDir()) {
    path_str = std::format("{}/{}", GA::I().ActionResultDir().string(),
                           "metrics.json");
  } else {
    path_str =
        std::format("{}/{}-metrics.json", GA::I().ActionResultDir().string(),
                    GA::I().InstanceFullPath().filename().string());
  }
  if (GA::I().DebugLevel() > 0) {
    std::cout << "CreateMetricsFilestream: " << path_str << std::endl;
  }

  std::ofstream fout(path_str);
  return fout;
}
} // namespace frm