#include "io.h"
#include <iostream>

namespace frm {

std::ofstream CreateResultFilestream(const std::string &path) {
  // 生成输出的路径
  fs::path curr_result_dir = global::CurrResultDir();
  std::string path_str;
  if (global::UseIndividualModelDir()) {
    path_str = std::format("{}/{}", curr_result_dir.string(), path);
  } else {
    path_str = std::format("{}/{}-{}", curr_result_dir.string(),
                           global::CurrFile().filename().string(), path);
  }
  if (global::DebugLevel() > 0) {
    std::cout << "CreateResultFilestream: " << path_str << std::endl;
  }
  return std::ofstream(path_str);
}

std::ofstream CreateDebugFilestream(const std::string &path) {
  // 生成输出的路径
  fs::path curr_debug_dir = global::CurrDebugDir();
  return std::ofstream(global::CurrDebugDir() / path);
}
} // namespace frm