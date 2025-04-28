#include "common_program.h"
#include "assert.hpp"
#include "process_parallel_executor.h"
#include <iostream>
#include <regex>

namespace frm {

using GA = GlobalArguments;

bool CommonProgram::PrepareWorkingDirectory() {
  // 检查work文件夹是否存在

  if (!fs::exists(GA::I().WorkDir()))
    fs::create_directories(GA::I().WorkDir());

  // 检查并刷新本步骤的工作目录
  if (GA::I().CleanActionCache())
    if (fs::exists(GA::I().ActionDir()))
      fs::remove_all(GA::I().ActionDir());

  if (!fs::exists(GA::I().ActionDir()))
    fs::create_directories(GA::I().ActionDir());

  if (!fs::exists(GA::I().ActionDebugDir())) {
    fs::create_directories(GA::I().ActionDebugDir());
  }

  if (!fs::exists(GA::I().ActionLogDir()))
    fs::create_directories(GA::I().ActionLogDir());

  if (!fs::exists(GA::I().ActionResultDir()))
    fs::create_directories(GA::I().ActionResultDir());

  return true;
}

// 确定本次运行的run_gargets
bool CommonProgram::SelectRunTargets() {
  run_targets_.clear();
  GlobalArguments &ga = GA::I();
  // 从dataset中选择本次运行的文件
  if (!GA::I()
           .SingleInstanceName()
           .empty()) { // 如果字符串g_single_instance不为空
    // 验证该文件的前置输入是合法的
    if (CheckRunTargetInputValidity(GA::I().InstanceFullPath()))
      run_targets_.push_back(GA::I().InstanceFullPath());
  } else if (!GA::I().BatchInstanceRegex().empty()) {
    std::smatch match;
    const std::string &file_regex = GA::I().BatchInstanceRegex();
    std::regex pattern(file_regex);
    for (const auto &entry :
         std::filesystem::directory_iterator(GA::I().DatasetDir())) {
      if (entry.is_directory())
        continue;
      std::string filename = entry.path().filename().string();
      if (std::regex_match(filename, match, pattern))
        if (CheckRunTargetInputValidity(entry))
          run_targets_.push_back(entry.path().string());
    }
  }
  return true;
}

bool CommonProgram::CheckRunTargetInputValidity(const fs::path &instance_path) {
  if (GA::I().CleanActionCache() == 0)
    return true;

  // std::string instance_name = instance_path.filename().string();

  // for (const std::string &input_name :
  //      all_action_list_[curr_action_idx_].inputs) {
  //   std::size_t input_file_action_idx =
  //       map_input_file_to_action_idx[input_name];
  //   fs::path input_file_path =
  //       g_work_dir /
  //       std::format("{}_{}", input_file_action_idx,
  //                   all_action_list_[input_file_action_idx].name) /
  //       "result" / std::format("{}-{}", instance_name, input_name);
  //   if (!fs::exists(input_file_path))
  //     return false;
  // }

  return true;
  // 计算当前步骤输入模型的
  // std::unordered_map<std::string, std::size_t> map_input_file_to_step_idx;
  // for (std::size_t action_idx = 0; action_idx < curr_action_idx_;
  //      ++action_idx) {
  // }
}

int CommonProgram::Run(const std::function<void()> &func) const {
  if (run_targets_.size() == 0)
    return 0;
  // 如果只处理单个文件, 则直接调用func
  if (!GA::I().SingleInstanceName().empty()) {
    func();
    return 0;
  } else {
    ProcessParallelExecutor process_parallel_executor(
        func, run_targets_, GA::I().NumParallelCnt(),
        GA::I().CurrActionPath().string());
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
