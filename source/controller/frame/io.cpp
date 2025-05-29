#include "io.h"
#include "global_args.h"
#include "pa_assert.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

namespace frm {
using GA = GlobalArguments;

std::ofstream CreateOutputFilestream(const std::string &path) {
  if (global::DebugLevel() > 1) {
    const GA &ga = GlobalArguments::I();
    std::cout << "path: " << path << std::endl;
    for (auto str : global::ActionOutputs())
      std::cout << str << std::endl;
    std::cout << global::ActionOutputs().count(path) << std::endl;
  }
  PA_ASSERT_WITH_MSG(global::ActionOutputs().count(path) != 0,
                     std::format("{}不符合本工具的输出", path));

  // 生成输出的路径
  const fs::path &curr_result_dir = global::ActionResultDir();
  std::string path_str;
  if (global::UseIndividualInstanceDir()) {
    path_str = std::format("{}/{}", curr_result_dir.string(), path);
  } else {
    path_str =
        std::format("{}/{}-{}", curr_result_dir.string(),
                    global::InstanceFullPath().filename().string(), path);
  }
  if (global::DebugLevel() > 0) {
    std::cout << "CreateOutputFilestream: " << path_str << std::endl;
  }
  return std::ofstream(path_str);
}

std::ofstream CreateDebugFilestream(const std::string &path) {
  // 生成输出的路径
  const fs::path &curr_debug_dir = global::ActionDebugDir();
  std::string path_str;
  if (global::UseIndividualInstanceDir()) {
    path_str = std::format("{}/{}", curr_debug_dir.string(), path);
  } else {
    path_str =
        std::format("{}/{}-{}", curr_debug_dir.string(),
                    global::InstanceFullPath().filename().string(), path);
  }
  return std::ofstream(path_str);
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
                    global::InstanceFullPath().filename().string());
  }
  if (global::DebugLevel() > 0) {
    std::cout << "CreateMetricsFilestream: " << path_str << std::endl;
  }

  std::ofstream fout(path_str);
  return fout;
}

ValueType GetHyperParameter(std::string name) {
  // 读取当前工具的超参数
  PA_ASSERT_WITH_MSG(global::ActionHyperParameters().count(name) != 0,
                     "超参数不存在");
  return global::ActionHyperParameters().at(name);
}

std::ifstream GetInputFilestream(const std::string &input_name) {
  PA_ASSERT_WITH_MSG(global::BatchInstanceRegex().empty(), "非单例模式");
  PA_ASSERT_WITH_MSG(global::ActionInputs().count(input_name) != 0,
                     "input_name不是设计的输入文件");
  std::size_t action_idx = global::ActionInputs().at(input_name);
  fs::path input_full_path =
      global::ActionResultDir(action_idx) /
      std::format("{}-{}", global::SingleInstanceName(), input_name);
  PA_ASSERT_WITH_MSG(fs::exists(input_full_path), "前置输入文件不存在");

  return std::ifstream(input_full_path);
}
} // namespace frm