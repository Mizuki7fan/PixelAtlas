#include "io.h"
#include "global_defs.h"
#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

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

std::ofstream CreateMetricsFilestreamBegin() {
  // 生成输出的路径
  std::string path_str;
  if (global::UseIndividualModelDir()) {
    path_str =
        std::format("{}/{}", global::CurrResultDir().string(), "metrics.json");
  } else {
    path_str =
        std::format("{}/{}-metrics.json", global::CurrResultDir().string(),
                    global::CurrFile().filename().string());
  }
  if (global::DebugLevel() > 0) {
    std::cout << "CreateMetricsFilestream: " << path_str << std::endl;
  }

  std::ofstream fout(path_str);
  fout << "{\n";

  return fout;
}

void WriteMetrics(std::ofstream &fout, const std::string &name, double value) {
  fout << "  \"" << name << "\": " << value << "\n";
}

void WriteMetricsFilestreamEnd(std::ofstream &fout) {
  fout << "}\n";
  fout.close();
}
} // namespace frm