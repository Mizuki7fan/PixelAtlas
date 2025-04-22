#pragma once
#include "action_arguments.h"
#include "global_defs.h"
#include <filesystem>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem; // 确保命名空间别名在全局或正确作用域
namespace frm {
// 全局变量

int GetDebugLevel();
fs::path GetActionDebugDir();
fs::path GetActionResultDir();
bool GetUseIndividualInstanceDir();
int GetMaxTimeElapsed();
std::string GetDatasetName();
std::string GetWorkName();
fs::path GetInstancePath();

class CommonProgram {
public:
  CommonProgram(int argc, char *argv[]);
  int Run(const std::function<void()> &func) const;

private:                          // functions
  bool PrepareWorkingDirectory(); // 准备work文件夹
  bool PrepareWorkingDirectoryForIndividualRunning();
  bool
  SelectRunTargets(); // 根据依赖的前置工具和文件通配符, 确定本次运行所有的目标
  // 检测当前run_target的前置输入是否完整
  bool CheckRunTargetInputValidity(const fs::path &run_target);

private:
  std::vector<ActionArguments> all_action_list_;
  std::unordered_map<std::string, std::size_t> map_action_name_to_action_idx;
  // 当前步骤的输入文件所对应的action_idx
  std::unordered_map<std::string, std::size_t> map_input_file_to_action_idx;
  std::size_t curr_action_idx_;
  fs::path curr_action_path_;

  // 通过正则表达式筛选出来的例子
  std::vector<fs::path> run_targets_;
};
} // namespace frm
