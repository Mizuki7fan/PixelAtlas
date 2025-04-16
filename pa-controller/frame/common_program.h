#pragma once
#include "command_structures.h"
#include <filesystem>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem; // 确保命名空间别名在全局或正确作用域
namespace frm {
// 全局变量

int GetDebugLevel();
fs::path GetCurrDebugDir();
fs::path GetCurrResultDir();
fs::path GetCurrFile();
std::string GetDatasetStr(); // 取dataset路径
int GetMaxTimeElapsed();
bool GetUseIndividualModelDir();
std::string GetRunName();

class CommonProgram {
public:
  CommonProgram(int argc, char *argv[]);
  int Run(const std::function<void()> &func) const;

private:                          // functions
  bool PrepareWorkingDirectory(); // 准备work文件夹
  bool PrepareWorkingDirectoryForIndividualRunning();
  bool
  SelectRunTargets(); // 根据依赖的前置工具和文件通配符, 确定本次运行所有的目标
  bool SelectRunTargetsOfFirstCmd();
  bool SelectRunTargetsOfFollowingCmd();

private:
  std::vector<StepArguments> all_step_list;
  std::unordered_map<std::string, std::size_t> map_step_name_to_step_idx;
  std::size_t curr_cmd_idx;
  std::string curr_cmd_name;
  fs::path curr_cmd_path;

  // 通过正则表达式筛选出来的例子
  std::vector<fs::path> run_targets;
};
} // namespace frm
