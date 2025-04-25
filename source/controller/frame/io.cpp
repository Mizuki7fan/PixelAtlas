#include "io.h"
#include "global_defs.h"
#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

namespace frm {

std::ofstream CreateOutputFilestream(const std::string &path) {
  // 生成输出的路径
  fs::path curr_result_dir = global::ActionResultDir();
  std::string path_str;
  if (global::UseIndividualInstanceDir()) {
    path_str = std::format("{}/{}", curr_result_dir.string(), path);
  } else {
    path_str = std::format("{}/{}-{}", curr_result_dir.string(),
                           global::InstancePath().filename().string(), path);
  }
  if (global::DebugLevel() > 0) {
    std::cout << "CreateOutputFilestream: " << path_str << std::endl;
  }
  return std::ofstream(path_str);
}

std::ofstream CreateDebugFilestream(const std::string &path) {
  // 生成输出的路径
  fs::path curr_debug_dir = global::ActionDebugDir();
  return std::ofstream(global::ActionDebugDir() / path);
}

std::ofstream CreateMetricsFilestream() {
  // 生成输出的路径
  std::string path_str;
  if (global::UseIndividualInstanceDir()) {
    path_str = std::format("{}/{}", global::ActionResultDir().string(),
                           "metrics.json");
  } else {
    path_str =
        std::format("{}/{}-metrics.json", global::ActionResultDir().string(),
                    global::InstancePath().filename().string());
  }
  if (global::DebugLevel() > 0) {
    std::cout << "CreateMetricsFilestream: " << path_str << std::endl;
  }

  std::ofstream fout(path_str);
  return fout;
}
} // namespace frm