#include "common_program.h"
#include "assert.hpp"
#include "parse.h"
#include "process_parallel_executor.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <regex>

namespace po = boost::program_options;

namespace frm {
// 本批次运行的名字, 例如work_baseline;
static std::string g_work_name = "";
std::string GetWorkName() { return g_work_name; }

// work目录 本批次运行的目录
static fs::path g_work_dir;

// 工具的目录, 例如0_Preprocess
static fs::path g_action_dir;

// 工具的调试目录
static fs::path g_action_debug_dir;
fs::path GetActionDebugDir() { return g_action_debug_dir; }

// 当前日志输出目录
static fs::path g_action_log_dir;

// 当前结果输出目录
static fs::path g_action_result_dir;
fs::path GetActionResultDir() { return g_action_result_dir; }

// 调试等级
static int g_debug_level = 0;
int GetDebugLevel() { return g_debug_level; };

// 当前处理的文件名
static fs::path g_instance_path = "";
fs::path GetInstancePath() { return g_instance_path; }

// 批量执行时的文件名通配符
static std::string g_batch_instance_regex = ".*.*";

// 单一文件执行时的文件名
static std::string g_single_instance = "alien.obj";

// 数据集名
static std::string g_dataset_name = ".";
std::string GetDatasetName() { return g_dataset_name; };

// 并行数
static int g_num_parallel_cnt = 1;

// 是否每个文件独立新建文件夹
static bool g_use_individual_instance_dir = false;
bool GetUseIndividualInstanceDir() { return g_use_individual_instance_dir; };

// 程序执行最大用时
static int g_max_time_elapsed = 1800;
int GetMaxTimeElapsed() { return g_max_time_elapsed; };

static bool g_clean_action_cache = false;

bool CommonProgram::PrepareWorkingDirectory() {
  // 检查work文件夹是否存在
  g_work_dir = fs::current_path() / std::format("../work_{}", g_work_name);
  g_action_dir =
      g_work_dir / std::format("{}_{}", curr_action_idx_,
                               all_action_list_[curr_action_idx_].name);

  if (!fs::exists(g_work_dir))
    fs::create_directories(g_work_dir);

  // 检查并刷新本步骤的工作目录
  if (g_clean_action_cache)
    if (fs::exists(g_action_dir))
      fs::remove_all(g_action_dir);

  if (!fs::exists(g_action_dir))
    fs::create_directories(g_action_dir);

  g_action_debug_dir = g_action_dir / "debug";
  if (!fs::exists(g_action_debug_dir)) {
    fs::create_directories(g_action_debug_dir);
  }

  g_action_log_dir = g_action_dir / "log";
  if (!fs::exists(g_action_log_dir))
    fs::create_directories(g_action_log_dir);

  g_action_result_dir = g_action_dir / "result";
  if (!fs::exists(g_action_result_dir))
    fs::create_directories(g_action_result_dir);

  return true;
}

bool CommonProgram::PrepareWorkingDirectoryForIndividualRunning() {
  if (!g_use_individual_instance_dir || run_targets_.size() != 1)
    return true;

  // debug
  g_action_debug_dir = g_action_debug_dir / run_targets_[0].filename();
  if (!fs::exists(g_action_debug_dir))
    fs::create_directories(g_action_debug_dir);

  // log
  g_action_log_dir = g_action_log_dir / run_targets_[0].filename();
  if (!fs::exists(g_action_log_dir))
    fs::create_directories(g_action_log_dir);

  // result
  g_action_result_dir = g_action_result_dir / run_targets_[0].filename();
  if (!fs::exists(g_action_result_dir))
    fs::create_directories(g_action_result_dir);

  return true;
}

// 确定本次运行的run_gargets
bool CommonProgram::SelectRunTargets() {
  // 匹配规则
  fs::path project_root_dir = fs::current_path() / ".." / "..";
  fs::path project_asset_dir = project_root_dir / "asset";
  // 指定数据集
  fs::path run_dataset_dir = project_asset_dir / g_dataset_name;
  run_targets_.clear();

  // 从dataset中选择本次运行的文件
  if (!g_single_instance.empty()) { // 如果字符串g_single_instance不为空
    fs::path g_instance_full_path = run_dataset_dir / g_single_instance;
    PA_ASSERT_WITH_MSG(fs::exists(g_instance_full_path),
                       std::format("{}不存在", g_instance_full_path.string()));
    run_targets_.push_back(g_instance_full_path);
  } else if (!g_batch_instance_regex.empty()) {
    std::smatch match;
    std::string &file_regex = g_batch_instance_regex;
    std::regex pattern(file_regex);
    for (const auto &entry :
         std::filesystem::directory_iterator(run_dataset_dir)) {
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

bool CommonProgram::CheckRunTargetInputValidity(const fs::path &path) {
  if (curr_action_idx_ == 0)
    return true;
  return true;
  // 计算当前步骤输入模型的
  // std::unordered_map<std::string, std::size_t> map_input_file_to_step_idx;
  // for (std::size_t action_idx = 0; action_idx < curr_action_idx_;
  //      ++action_idx) {
  // }
}

// bool CommonProgram::SelectRunTargetsOfFirstCmd() {
//   // 匹配规则
//   fs::path project_root_dir = fs::current_path() / ".." / "..";
//   fs::path project_asset_dir = project_root_dir / "asset";
//   // 指定数据集
//   fs::path run_dataset_dir = project_asset_dir / g_dataset_name;
//   run_targets_.clear();

//   if (!g_single_filename.empty()) {
//     g_curr_file = run_dataset_dir / g_single_filename;
//     PA_ASSERT(fs::exists(g_curr_file));
//     run_targets_.push_back(g_curr_file);
//   } else if (!g_batch_filename_regex_.empty()) {
//     // 按照正则表达式匹配
//     std::smatch match;
//     std::string &file_regex = g_batch_filename_regex_;
//     std::regex pattern(file_regex);

//     // 如果当前工具为第一个工具
//     if (curr_cmd_idx == 0) {
//       for (const auto &entry :
//            std::filesystem::directory_iterator(run_dataset_dir)) {
//         if (entry.is_directory())
//           continue;
//         //
//         std::string filename = entry.path().filename().string();
//         if (std::regex_match(filename, match, pattern)) {
//           run_targets_.push_back(entry.path().string());
//         }
//       }
//       if (g_debug_level_arg > 1) {
//         std::cout << std::format("共成功匹配{}个例子:", run_targets_.size())
//                   << std::endl;
//         for (auto run_file : run_targets)
//           std::cout << run_file << std::endl;
//       }
//     }
//   }
//   return true;
// }

// bool CommonProgram::SelectRunTargetsOfFollowingCmd() {
//   // 算当前工具的输入文件所依赖的前序步骤
//   run_targets_.clear();

//   std::unordered_map<std::string, std::size_t> map_input_file_to_step_idx;
//   for (std::size_t step_idx = 0; step_idx < curr_cmd_idx; ++step_idx) {
//     std::cout << all_step_list[step_idx].step_name << std::endl;
//     for (const std::string &input_file_name :
//          all_step_list[step_idx].output_files) {
//       std::cout << input_file_name << std::endl;
//       map_input_file_to_step_idx[input_file_name] = step_idx;
//     }
//   }

//   fs::path project_root_dir = fs::current_path() / ".." / "..";

//   if (!g_single_filename.empty()) {
//     for (auto input_file_name : all_step_list[curr_cmd_idx].input_files) {
//       PA_ASSERT_WITH_MSG(
//           map_input_file_to_step_idx.count(input_file_name) != 0,
//           std::format("当前输入文件{}的前序步骤不存在", input_file_name));
//       std::size_t curr_input_file_output_cmd_idx =
//           map_input_file_to_step_idx.at(input_file_name);
//       const fs::path &curr_working_dir = g_curr_working_dir / "..";
//       fs::path curr_input_file_dir =
//           curr_working_dir /
//           std::format("{}_{}", curr_input_file_output_cmd_idx,
//                       all_step_list[curr_input_file_output_cmd_idx].step_name);
//       fs::path curr_input_file =
//           curr_input_file_dir / "result" /
//           std::format("{}-{}", g_single_filename, input_file_name);
//       if (!fs::exists(curr_input_file))
//         return false;
//       else {
//         g_input_file_full_path[input_file_name] = curr_input_file;
//       }
//     }
//     run_targets_.push_back(g_single_filename);
//   } else if (!g_batch_filename_regex_.empty()) {
//     // 填充run_targets
//     std::unordered_map<std::string, bool> map_model_name_is_valid;
//     for (auto input_file_name : all_step_list[curr_cmd_idx].input_files) {
//       // 遍历该命令的所有输入文件
//       PA_ASSERT_WITH_MSG(
//           map_input_file_to_step_idx.count(input_file_name) != 0,
//           std::format("当前输入文件{}的前序步骤不存在", input_file_name));
//       std::size_t curr_input_file_output_cmd_idx =
//           map_input_file_to_step_idx.at(input_file_name);
//       const fs::path &curr_working_dir = g_curr_working_dir / "..";
//       fs::path curr_input_file_dir =
//           curr_working_dir /
//           std::format("{}_{}", curr_input_file_output_cmd_idx,
//                       all_step_list[curr_input_file_output_cmd_idx].step_name)
//                       /
//           "result";
//       if (map_model_name_is_valid.empty()) {
//         for (const auto &entry : fs::directory_iterator(curr_input_file_dir))
//         {
//           if (entry.is_directory())
//             continue;
//           std::string filename = entry.path().filename().string();
//           if (filename.size() > input_file_name.size()) {
//             std::size_t suffix_pos = filename.size() -
//             input_file_name.size(); if (filename.substr(suffix_pos) ==
//             input_file_name) {
//               std::string model_name = filename.substr(0, suffix_pos - 1);
//               if (map_model_name_is_valid.count(model_name) == 0) {
//                 // 需要检测model_name是否满足通配条件

//                 map_model_name_is_valid[model_name] = true;
//                 run_targets_.push_back(model_name);
//               }
//             }
//           }
//         }
//       } else {
//         // 搜索该目录下所有包含input_file_name的文件, 并识别模型名
//         for (auto &[model_name, is_valid] : map_model_name_is_valid) {
//           if (!is_valid)
//             continue;
//           fs::path curr_input_file =
//               curr_input_file_dir /
//               std::format("{}-{}", model_name, input_file_name);
//           if (fs::exists(curr_input_file)) {
//             map_model_name_is_valid[model_name] = true;
//             run_targets_.push_back(model_name);
//           } else {
//             map_model_name_is_valid[model_name] = false;
//           }
//         }
//       }
//     }
//   }
//   return true;
// }

int CommonProgram::Run(const std::function<void()> &func) const {
  if (run_targets_.size() == 0)
    return 0;
  // 如果只处理单个文件, 则直接调用func
  if (!g_single_instance.empty()) {
    std::cout << "run target: " << run_targets_[0] << std::endl;
    g_instance_path = run_targets_[0];
    func();
    return 0;
  } else {
    ProcessParallelExecutor process_parallel_executor(
        func, run_targets_, g_num_parallel_cnt, curr_action_path_.string());
    process_parallel_executor.Exec();
    return 0;
  }
  return 1;
}

CommonProgram::CommonProgram(int argc, char *argv[]) {
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "帮助信息") //
      ("debug,d", po::value<int>(&g_debug_level)->default_value(0),
       "调试等级") //
      ("batch",
       po::value<std::string>(&g_batch_instance_regex)->default_value(""),
       "批量执行") //
      ("single", po::value<std::string>(&g_single_instance)->default_value(""),
       "单一执行") //
      ("dataset", po::value<std::string>(&g_dataset_name)->default_value(""),
       "数据集") //
      ("parallel,p", po::value<int>(&g_num_parallel_cnt)->default_value(1),
       "并行数") //
      ("use_individual_instance_dir",
       po::bool_switch(&g_use_individual_instance_dir)->default_value(false),
       "是否每个instance独立建立文件夹") //
      ("max_time_elapsed,t",
       po::value<int>(&g_max_time_elapsed)->default_value(1800),
       "最大耗时") //
      ("clean", po::bool_switch(&g_clean_action_cache)->default_value(false),
       "是否清理单步缓存") //
      ("work_name", po::value<std::string>(&g_work_name)->default_value(""),
       "设置运行批次名");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
      std::cout << desc << std::endl;
      exit(0);
    }
    po::notify(vm);
  } catch (const po::error &e) {
    std::cerr << "error: " << e.what() << std::endl;
    std::cerr << desc << std::endl;
    PA_ASSERT_WITH_MSG(0, "输入异常");
  }

  PA_FALSE_ASSERT_WITH_MSG(g_batch_instance_regex.empty() &&
                               g_single_instance.empty(),
                           "single和batch都未指定需要处理的文件");
  PA_FALSE_ASSERT_WITH_MSG(!g_batch_instance_regex.empty() &&
                               !g_single_instance.empty(),
                           "single和batch都指定了需要处理的文件");

  // 加载all-step.json文件
  all_action_list_ = LoadAllStepList();
  // 记录每个action_name的idx
  for (std::size_t action_idx = 0; action_idx < all_action_list_.size();
       ++action_idx)
    map_action_name_to_action_idx[all_action_list_[action_idx].name] =
        action_idx;
  curr_action_path_ = argv[0];
  std::string curr_action_name = curr_action_path_.stem().string();
  curr_action_idx_ = std::numeric_limits<std::size_t>::max();
  if (map_action_name_to_action_idx.contains(curr_action_name))
    curr_action_idx_ = map_action_name_to_action_idx.at(curr_action_name);

  PA_ASSERT_WITH_MSG(curr_action_idx_ !=
                         std::numeric_limits<std::size_t>::max(),
                     std::format("无效的指令: {}", argv[0]));

  // 统计当前步骤inputs的前置步骤的idx
  for (std::string input_file_name :
       all_action_list_[curr_action_idx_].inputs) {
    for (std::size_t i = 0; i < curr_action_idx_; ++i)
      for (auto output_file_name : all_action_list_[i].outputs)
        if (output_file_name == input_file_name)
          map_input_file_to_action_idx[input_file_name] = i;
  }

  for (std::string input_file_name : all_action_list_[curr_action_idx_].inputs)
    PA_ASSERT_WITH_MSG(map_input_file_to_action_idx.contains(input_file_name),
                       std::format("{}的前置输入文件不存在", input_file_name));

  // 准备当前工具运行需要的文件夹
  PrepareWorkingDirectory();
  // 选择本次运行需要处理的模型
  SelectRunTargets();
  // 为单独文件单独建立文件夹
  PrepareWorkingDirectoryForIndividualRunning();
};

} // namespace frm
