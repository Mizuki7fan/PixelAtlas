#include "evaluator.h"
#include <boost/program_options.hpp>
#include <filesystem>
namespace fs = std::filesystem;
namespace po = boost::program_options;

// 结果评估模块
// 功能: 输入两次work的名称A,B和工具名称,
// 在run/evaluate_A_B文件夹下生成评估结果
int main(int argc, char *argv[]) {
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

  // 存放两个work的结果信息
  SampleMetric sample_metric;

  for (std::size_t work_idx = 0; work_idx < 2; ++work_idx) {
    fs::path curr_work_dir =
        project_run_dir / std::format("work_{}", work_name[work_idx]);
    fs::path curr_cmd_result_dir =
        curr_work_dir /
        std::format("{}_{}", cmd_idx, all_step_list[cmd_idx].step_name) /
        "result";
    // 需要区分是否use_individual_model_dir
    // 要么全是文件夹, 要么不包含文件夹
    std::vector<fs::path> all_directories;
    std::vector<fs::path> all_files;

    for (const auto &entry : fs::directory_iterator(curr_cmd_result_dir)) {
      if (entry.is_directory()) {
        all_directories.push_back(entry.path());
      } else if (entry.is_regular_file()) {
        all_files.push_back(entry.path());
      }
    }
    PA_ASSERT_WITH_MSG(all_directories.size() * all_files.size() == 0,
                       "文件夹和文件个数不能同时不为0");
    if (all_files.size() != 0) {
      // 不为每个例子单独建一个文件夹
      for (const fs::path &file_name : all_files) {
        std::cout << "file_name: " << file_name << std::endl;
        if (file_name.filename().string().find("-metrics.json") ==
            std::string::npos)
          continue;
        sample_metric.FromJson(file_name.string(),
                               all_step_list[cmd_idx].metrics, work_idx);
      }
    } else if (all_directories.size() != 0) {
    }
  }

  return 0;
}