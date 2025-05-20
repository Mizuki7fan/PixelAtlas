#include "common_program.h"
#include "pa_assert.hpp"
#include "process_parallel_executor.h"
#include <iostream>
#include <regex>

namespace frm {

using GA = GlobalArguments;
GlobalArguments &ga = GlobalArguments::I();

bool CommonProgram::PrepareWorkingDirectory() {
  // 检查work文件夹是否存在
  std::cout << global::WorkDir() << std::endl;
  if (!fs::exists(global::WorkDir()))
    fs::create_directories(global::WorkDir());

  // 检查并刷新本步骤的工作目录
  if (global::CleanActionCache())
    if (fs::exists(global::ActionDir()))
      fs::remove_all(global::ActionDir());

  if (!fs::exists(global::ActionDir()))
    fs::create_directories(global::ActionDir());

  if (!fs::exists(global::ActionDebugDir())) {
    fs::create_directories(global::ActionDebugDir());
  }

  if (!fs::exists(global::ActionLogDir()))
    fs::create_directories(global::ActionLogDir());

  if (!fs::exists(global::ActionResultDir()))
    fs::create_directories(global::ActionResultDir());

  return true;
}

// 确定本次运行的run_gargets
bool CommonProgram::SelectRunTargets() {
  run_targets_.clear();
  // 从dataset中选择本次运行的文件
  if (!global::SingleInstanceName()
           .empty()) { // 如果字符串g_single_instance不为空
    // 验证该文件的前置输入是合法的
    if (CheckRunTargetInputValidity(global::SingleInstanceName()))
      run_targets_.push_back(global::SingleInstanceName());
  } else if (!global::BatchInstanceRegex().empty()) {
    std::smatch match;
    const std::string &file_regex = global::BatchInstanceRegex();
    std::regex pattern(file_regex);
    for (const auto &entry :
         std::filesystem::directory_iterator(global::DatasetDir())) {
      if (entry.is_directory())
        continue;
      std::string filename = entry.path().filename().string();
      if (std::regex_match(filename, match, pattern))
        if (CheckRunTargetInputValidity(filename))
          run_targets_.push_back(entry.path().filename().string());
    }
  }
  return true;
}

bool CommonProgram::CheckRunTargetInputValidity(
    const std::string &instance_name) const {
  fs::path instance_full_path = global::DatasetDir() / instance_name;
  if (!fs::exists(instance_full_path))
    return false;

  if (global::CurrActionIdx() != 0) {
    // 需要检查当前工具的输入依赖文件是否都存在
    for (const auto &[input_name, action_idx] : global::ActionInputs()) {
      fs::path input_full_path =
          global::ActionResultDir(action_idx) /
          std::format("{}-{}", instance_name, input_name);
      if (!fs::exists(input_full_path))
        return false;
    }
  }
  return true;
}

int CommonProgram::Run(const std::function<void()> &func) const {
  if (run_targets_.size() == 0)
    return 0;
  // 如果只处理单个文件, 则直接调用func
  if (!global::SingleInstanceName().empty()) {
    func();
    return 0;
  } else {
    ProcessParallelExecutor process_parallel_executor(
        func, run_targets_, global::NumParallelCnt(),
        global::CurrActionPath().string());
    process_parallel_executor.Exec();
    return 0;
  }
  return 1;
}

CommonProgram::CommonProgram(int argc, char *argv[]) {
  try {
    GA::I().Initialize(argc, argv, {});
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    exit(1);
  }
  // 准备当前工具运行需要的文件夹
  PrepareWorkingDirectory();
  // 选择本次运行需要处理的模型
  SelectRunTargets();
};

} // namespace frm
