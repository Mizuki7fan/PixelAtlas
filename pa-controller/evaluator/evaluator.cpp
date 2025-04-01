#include "evaluator.h"
#include <boost/program_options.hpp>
#include <filesystem>
#include <frame/assert.hpp>
#include <frame/command_structures.h>
#include <frame/parse.h>
#include <iostream>

namespace po = boost::program_options;
namespace fs = std::filesystem;

// 结果评估模块
// 功能: 输入两次work的名称A,B和工具名称, 在run/evaluate_A_B文件夹下生成评估结果
int main(int argc, char *argv[]) {
  int debug_level_arg = 0;
  std::string work_name[2];
  int cmd_idx = -1;

  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "帮助信息")(
      "debug,d", po::value<int>(&debug_level_arg)->default_value(0),
      "调试等级") //
      ("cmd_idx", po::value<int>(&cmd_idx)->default_value(-1),
       "工具名") //
      ("work1", po::value<std::string>(&work_name[0])->default_value("base"),
       "work_1_name") //
      ("work2", po::value<std::string>(&work_name[1])->default_value("curr"),
       "work_2_name"); //
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

  std::vector<StepArguments> all_step_list = frm::LoadAllStepList();
  PA_ASSERT_WITH_MSG(cmd_idx < all_step_list.size(), "输入cmd_idx无效");

  fs::path project_root_dir = fs::current_path() / ".." / "..";
  fs::path project_run_dir = project_root_dir / "run";

  for (std::size_t work_idx = 0; work_idx < 2; ++work_idx) {
    fs::path curr_work_dir =
        project_run_dir / std::format("work_{}", work_name[work_idx]);
    fs::path curr_cmd_result_dir =
        curr_work_dir /
        std::format("{}_{}", cmd_idx, all_step_list[cmd_idx].step_name) /
        "result";
    // 需要区分是否use_individual_model_dir
    // 要么全是文件夹, 要么不包含文件夹
    std::size_t num_directories = 0, num_files = 0;
    for (const auto &entry : fs::directory_iterator(curr_cmd_result_dir)) {
      if (entry.is_directory()) {
        ++num_directories;
      } else if (entry.is_regular_file()) {
        ++num_files;
      }
    }
    PA_ASSERT_WITH_MSG(num_directories * num_files == 0,
                       "文件夹和文件个数不能同时不为0");
  }

  return 0;
}