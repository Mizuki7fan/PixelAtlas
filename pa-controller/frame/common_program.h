#pragma once
#include "command_structures.h"
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

private:                          // functions
  bool PrepareWorkingDirectory(); // 准备work文件夹

  bool
  SelectRunTargets(); // 根据依赖的前置工具和文件通配符, 确定本次运行所有的目标

private:
  std::vector<StepArguments> all_step_list;
  std::unordered_map<std::string, std::size_t> map_step_name_to_step_idx;
  std::size_t curr_cmd_idx;
  std::string curr_cmd_name;
  fs::path curr_cmd_path;

  // 通过正则表达式筛选出来的例子
  std::vector<std::string> run_targets;
};

fs::path GetGlobalWorkDir();
fs::path GetGlobalCurrCmdDir();
int GetDebugLevel();
std::string GetDataset();
int GetNumParallelCnt();
bool UseIndividualModelDir();
int GetMaxTimeElapsed();
std::string GetParallelLevel();
} // namespace frm
