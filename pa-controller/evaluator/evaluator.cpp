#include "evaluator.h"
#include <boost/json.hpp>
#include <boost/program_options.hpp>
#include <frame/assert.hpp>
#include <frame/parse.h>
#include <fstream>
#include <iostream>

namespace po = boost::program_options;
Evaluator::Evaluator(int argc, char *argv[]) {
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "帮助信息")(
      "debug,d", po::value<int>(&debug_level_arg)->default_value(0),
      "调试等级") //
      ("cmd_idx", po::value<int>(&cmd_idx_)->default_value(-1),
       "工具名") //
      ("work1", po::value<std::string>(&work_name_[0])->default_value("base"),
       "work_1_name") //
      ("work2", po::value<std::string>(&work_name_[1])->default_value("curr"),
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
}

void Evaluator::LoadData() {
  all_step_list_ = frm::LoadAllStepList();
  PA_ASSERT_WITH_MSG(cmd_idx_ != -1, "输入cmd_idx无效");

  fs::path project_root_dir = fs::current_path() / ".." / "..";
  fs::path project_run_dir = project_root_dir / "run";

  for (std::size_t work_idx = 0; work_idx < 2; ++work_idx) {
    fs::path curr_work_dir =
        project_run_dir / std::format("work_{}", work_name_[work_idx]);
    fs::path curr_cmd_result_dir =
        curr_work_dir /
        std::format("{}_{}", cmd_idx_, all_step_list_[cmd_idx_].step_name) /
        "result";

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
    if (all_files.size() != 0)
      LoadDataFromNoIndividualModelDir(all_files, static_cast<int>(work_idx));
    else if (all_directories.size() != 0)
      LoadDataFromIndividualModelDir();
  }
}

void Evaluator::LoadDataFromNoIndividualModelDir(
    const std::vector<fs::path> &all_files, int curr_work_idx) {
  // 读入信息放入metric_data_中
  std::unordered_map<std::string, frm::Metric> &metric_data =
      sample_metric_data_[curr_work_idx];

  for (const fs::path &file_name : all_files) {
    std::cout << "file_name: " << file_name << std::endl;
    std::string model_name = file_name.filename().string();
    std::size_t pos = model_name.find("-metrics.json");
    if (pos == std::string::npos)
      continue;
    model_name = model_name.substr(0, pos);
    std::cout << "model_name: " << model_name << std::endl;
    std::ifstream ifs(file_name);
    frm::LoadMetricJsonFile(ifs, all_step_list_[cmd_idx_].metrics,
                            metric_data[model_name]);
    ifs.close();
  }
}

void Evaluator::LoadDataFromIndividualModelDir() {}

void Evaluator::PrintData() {
  std::cout << "model_name\t";
  for (auto metric_name : all_step_list_[cmd_idx_].metrics)
    std::cout << std::format("work_0_{}\twork_1_{}\t", metric_name.first,
                             metric_name.first);
  std::cout << std::endl;

  std::set<std::string> all_samples;
  for (std::size_t i = 0; i < 2; ++i)
    for (auto data : sample_metric_data_[i])
      all_samples.insert(data.first);

  for (auto sample : all_samples) {
    std::cout << sample << "\t";
    const frm::Metric &metric_data_0 = sample_metric_data_[0][sample];
    const frm::Metric &metric_data_1 = sample_metric_data_[1][sample];
    for (auto [name, value] : metric_data_0) {
      std::visit(frm::PrintMetricValue{}, value);
      std::cout << "\t";
    }
    for (auto [name, value] : metric_data_1) {
      std::visit(frm::PrintMetricValue{}, value);
      std::cout << "\t";
    }
    std::cout << std::endl;
  }
}