#pragma once
#include "command_structures.h"
#include "parse.h"
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace frm {
// 全局变量

class CommonProgram {
public:
  CommonProgram(int argc, char *argv[]);
  int Run(const std::function<void(std::string model_name)> &func) const;

private:                                // functions
  std::size_t GetCurrentProgramIndex(); // 取当前的cmd的idx
  // std::unordered_set<std::size_t>
  // GetCurrentProgramDependencies(); // 取当前program所依赖的前置program
  bool PrepareWorkingDirectory(); // 准备work文件夹

  bool
  SelectRunTargets(); // 根据依赖的前置工具和文件通配符, 确定本次运行所有的目标

private:
  std::vector<StepArguments> all_step_list;
  std::unordered_map<std::string, std::size_t> map_step_name_to_step_idx;
  std::vector<std::string> all_args;
  std::vector<std::string> all_models;
  std::size_t curr_cmd_idx;
};

fs::path GetGlobalWorkDir();
fs::path GetGlobalCurrCmdDir();
int GetDebugLevel();
} // namespace frm
