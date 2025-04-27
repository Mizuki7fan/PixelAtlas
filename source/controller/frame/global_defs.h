#pragma once
#include "action_args.h"
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

namespace frm {
struct GlobalArguments {
  // 暂时列举所有的, 之后部分参数可以用函数现场计算
  std::string work_name = "";
  fs::path work_dir;
  fs::path action_dir;
  fs::path action_debug_dir;
  fs::path action_log_dir;
  fs::path action_result_dir;
  int debug_level = 0;
  fs::path instance_path;
  std::string batch_instance_regex;
  std::string single_instance;
  std::string dataset_name;
  ActionArguments action_args;
  int num_parallel_cnt = 1;
  bool use_individual_instance_dir = false;
  int max_time_elapsed = 0;
  bool clean_action_dir = false;
};

} // namespace frm
namespace frm::global {

int DebugLevel();
fs::path ActionDebugDir();
fs::path ActionResultDir();
std::string DatasetName();
fs::path InstancePath();
int MaxTimeElapsed();
bool UseIndividualInstanceDir();
std::string WorkName();
const ActionArguments &ActionArgs();
const GlobalArguments &GlobalArgs();

} // namespace frm::global