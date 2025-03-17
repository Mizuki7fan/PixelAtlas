#include "common_program.h"
#include "assert.hpp"
#include "global_defs.h"
#include "parse.h"
#include "process_parallel_executor.h"
#include "thread_parallel_executor.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <regex>

namespace po = boost::program_options;

namespace frm {
// work目录
static fs::path g_working_dir;
fs::path GetGlobalWorkDir() { return g_working_dir; };

// 当前工具目录
static fs::path g_curr_working_dir;
fs::path GetCurrWorkingDir() { return g_curr_working_dir; };

// 当前调试输出目录
static fs::path g_curr_debug_dir;
fs::path GetCurrDebugDir() { return g_curr_debug_dir; };

// 当前日志输出目录
static fs::path g_curr_log_dir;
fs::path GetCurrLogDir() { return g_curr_log_dir; };

// 当前结果输出目录
static fs::path g_curr_result_dir;
fs::path GetResultDir() { return g_curr_result_dir; };

// 调试等级
static int g_debug_level_arg = 0;
int GetDebugLevel() { return g_debug_level_arg; };

// 批量执行时的文件名通配符
static std::string g_batch_filename_regex_ = ".*.*";
// 单一文件执行时的文件名
static std::string g_single_filename = "alien.obj";

// 数据集名
static std::string g_dataset_str = ".";
std::string GetDataset() { return g_dataset_str; };

// 并行数
static int g_num_parallel_cnt = 1;
int GetNumParallelCnt() { return g_num_parallel_cnt; };

// 是否每个文件独立新建文件夹
static bool g_use_individual_model_dir = false;
bool UseIndividualModelDir() { return g_use_individual_model_dir; };

// 程序执行最大用时
static int g_max_time_elapsed = 1800;
int GetMaxTimeElapsed() { return g_max_time_elapsed; };

// 并行级别: 默认为进程级的并行, 也可以设置为线程级并行(thread), 不并行
static std::string g_parallel_level = "process";
std::string GetParallelLevel() { return g_parallel_level; };

bool CommonProgram::PrepareWorkingDirectory() {
  // 检查work文件夹是否存在
  g_working_dir = fs::current_path() / "../work";
  // 如果没有work文件夹则新建
  if (!fs::exists(g_working_dir))
    fs::create_directories(g_working_dir);

  // 检查并刷新本步骤的工作目录
  g_curr_working_dir =
      g_working_dir /
      std::format("{}_{}", curr_cmd_idx, all_step_list[curr_cmd_idx].step_name);
  if (fs::exists(g_curr_working_dir))
    fs::remove_all(g_curr_working_dir);
  fs::create_directories(g_curr_working_dir);

  g_curr_debug_dir = g_curr_working_dir / "debug";
  fs::create_directories(g_curr_debug_dir);

  g_curr_log_dir = g_curr_working_dir / "log";
  fs::create_directories(g_curr_log_dir);

  g_curr_result_dir = g_curr_working_dir / "result";
  fs::create_directories(g_curr_result_dir);

  return true;
}

bool CommonProgram::SelectRunTargets() {
  // 匹配规则
  std::smatch match;
  std::string &file_regex = g_batch_filename_regex_;
  std::regex pattern(file_regex);

  // 如果当前工具为第一个工具
  if (curr_cmd_idx == 0) {
    fs::path project_root_dir = fs::current_path() / ".." / "..";
    fs::path project_asset_dir = project_root_dir / "asset";
    // 指定数据集
    fs::path run_dataset_dir = project_asset_dir / g_dataset_str;

    run_targets.clear();
    for (const auto &entry :
         std::filesystem::directory_iterator(run_dataset_dir)) {
      if (entry.is_directory())
        continue;
      //
      std::string filename = entry.path().filename().string();
      if (std::regex_match(filename, match, pattern)) {
        run_targets.push_back(entry.path().string());
      }
    }
    if (g_debug_level_arg > 1) {
      std::cout << std::format("共成功匹配{}个例子:", run_targets.size())
                << std::endl;
      for (auto run_file : run_targets)
        std::cout << run_file << std::endl;
    }
  }

  return true;
}

int CommonProgram::Run(const std::function<void(std::string)> &func) const {
  if (run_targets.size() == 0)
    return 0;
  // 如果只处理单个文件, 则直接调用func
  if (!g_single_filename.empty()) {
    func(g_single_filename);
    return 0;
  } else if (g_num_parallel_cnt == 1) { // 如果并行数为1, 则依次执行
    for (auto run_file : run_targets)
      func(run_file);
  }
  if (g_parallel_level == "thread") {
    ThreadParallelExecutor thread_parallel_executor(func, run_targets,
                                                    g_num_parallel_cnt);
    thread_parallel_executor.Exec();
    return 0;
  } else if (g_parallel_level == "process") {
    // 执行进程级并行
    ProcessParallelExecutor process_parallel_executor(
        func, run_targets, g_num_parallel_cnt, curr_cmd_path.string());
    process_parallel_executor.Exec();
    return 0;
  }

  return 1;
}

CommonProgram::CommonProgram(int argc, char *argv[]) {
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "帮助信息") //
      ("debug,d", po::value<int>(&g_debug_level_arg)->default_value(0),
       "调试等级") //
      ("batch",
       po::value<std::string>(&g_batch_filename_regex_)->default_value(""),
       "文件名正则") //
      ("single", po::value<std::string>(&g_single_filename)->default_value(""),
       "单一执行文件") //
      ("dataset", po::value<std::string>(&g_dataset_str)->default_value(""),
       "数据集名称") //
      ("parallel,p", po::value<int>(&g_num_parallel_cnt)->default_value(1),
       "并行执行数")     //
      ("parallel_level", // 并行级别, 取值必须为thread或者process
       po::value<std::string>(&g_parallel_level)
           ->default_value("process")              // 默认进程级并行
           ->notifier([](const std::string &val) { // 参数校验lambda
             if (val != "process" && val != "thread") {
               throw po::validation_error(
                   po::validation_error::invalid_option_value, "parallel_level",
                   val);
             }
           }),
       "并行级别 (process|thread)") //
      ("use_individual_model_dir",
       po::bool_switch(&g_use_individual_model_dir)->default_value(false),
       "是否每个文件独立建立文件夹") //
      ("max_time_elapsed,t",
       po::value<int>(&g_max_time_elapsed)->default_value(1800), "最大耗时");

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

  if (g_batch_filename_regex_.empty() && g_single_filename.empty()) {
    std::cerr << "error: 未指定需要处理的文件" << std::endl;
    std::cerr << desc << std::endl;
    exit(0);
  } else if (!g_batch_filename_regex_.empty() && !g_single_filename.empty()) {
    std::cerr << "不能既指定batch运行, 又指定单个文件" << std::endl;
    std::cerr << std::format("batch: {}, single: {}", g_batch_filename_regex_,
                             g_single_filename);
    exit(0);
  }

  // 加载all-step.json文件
  all_step_list = LoadAllStepList();
  for (std::size_t cmd_idx = 0; cmd_idx < all_step_list.size(); ++cmd_idx)
    map_step_name_to_step_idx[all_step_list[cmd_idx].step_name] = cmd_idx;

  // 获取当前工具的索引
  curr_cmd_path = argv[0];                       // 取exe路径
  curr_cmd_name = curr_cmd_path.stem().string(); // 取exe名
  if (map_step_name_to_step_idx.contains(curr_cmd_name))
    curr_cmd_idx = map_step_name_to_step_idx.at(curr_cmd_name);
  if (curr_cmd_idx == std::numeric_limits<std::size_t>::max())
    std::cerr << "invalid command: " << argv[0] << std::endl;

  // 准备当前工具运行需要的文件夹
  PrepareWorkingDirectory();
  // 选择本次运行需要处理的模型
  SelectRunTargets();
};

} // namespace frm
