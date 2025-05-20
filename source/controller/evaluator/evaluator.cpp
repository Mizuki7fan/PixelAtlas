#include "evaluator.h"
#include <boost/json.hpp>
#include <boost/program_options.hpp>
#include <frame/pa_assert.hpp>
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
  all_action_list_ = frm::LoadAllActionList();
  PA_ASSERT_WITH_MSG(cmd_idx_ != -1, "输入cmd_idx无效");

  fs::path project_root_dir = fs::current_path() / ".." / "..";
  fs::path project_run_dir = project_root_dir / RUN_DIR_NAME;

  for (std::size_t work_idx = 0; work_idx < 2; ++work_idx) {
    fs::path curr_work_dir =
        project_run_dir / std::format("work_{}", work_name_[work_idx]);
    fs::path curr_cmd_result_dir =
        curr_work_dir /
        std::format("{}_{}", cmd_idx_, all_action_list_[cmd_idx_].name) /
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
    frm::LoadMetricJsonFile(ifs, all_action_list_[cmd_idx_].metrics,
                            metric_data[model_name]);
    ifs.close();
  }
}

void Evaluator::LoadDataFromIndividualModelDir() {}

void Evaluator::PrintData() {
  // 表头格式化
  constexpr int kColWidth = 20; // 增加列宽到20字符
  constexpr int kPairSpace = 4; // 组间额外间距

  std::cout << std::format("{:<{}}", "model_name", kColWidth);

  for (auto metric_name : all_action_list_[cmd_idx_].metrics) {
    // 使用复合格式字符串实现双左对齐
    std::cout << std::format(
        "{:<{}}",
        std::format("\"{}\"-{}  \"{}\"-{}", // 添加两个空格分隔
                    work_name_[0], metric_name.first, work_name_[1],
                    metric_name.first),
        kColWidth * 2 + kPairSpace); // 双倍列宽+间距
  }
  std::cout << "\n";

  std::set<std::string> all_samples;
  for (std::size_t i = 0; i < 2; ++i)
    for (auto data : sample_metric_data_[i])
      all_samples.insert(data.first);

  // 数据行格式化
  for (auto sample : all_samples) {
    std::cout << std::format("{:<{}}", sample, kColWidth);

    const frm::Metric &metric_data_0 = sample_metric_data_[0].count(sample) != 0
                                           ? sample_metric_data_[0][sample]
                                           : frm::Metric();
    const frm::Metric &metric_data_1 = sample_metric_data_[1].count(sample) != 0
                                           ? sample_metric_data_[1][sample]
                                           : frm::Metric();

    // 使用visitor模式统一数值格式
    auto print_value = [](const auto &value) {
      return std::visit(
          [](auto &&arg) -> std::string {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_floating_point_v<T>) {
              return std::format("{:.10f}", arg); // 统一保留5位小数
            } else if constexpr (std::is_integral_v<T>) {
              return std::format("{}", arg);
            }
            return "";
          },
          value);
    };

    for (auto [name, value] : metric_data_0) {
      std::cout << std::format("{:>{}}", print_value(value), kColWidth);
    }
    for (auto [name, value] : metric_data_1) {
      std::cout << std::format("{:>{}}", print_value(value), kColWidth);
    }
    std::cout << "\n";
  }
}

void Evaluator::AnalyseDataDifference() {
  const std::unordered_map<std::string, std::string> &metrics =
      all_action_list_[cmd_idx_].metrics;
  std::set<std::string> all_samples;
  for (std::size_t i = 0; i < 2; ++i)
    for (auto data : sample_metric_data_[i])
      all_samples.insert(data.first);

  map_property_to_num_valid_data_.clear();
  map_property_to_squared_difference_.clear();
  map_property_to_max_difference_.clear();

  for (auto sample : all_samples) {
    const frm::Metric &metric_data_0 = sample_metric_data_[0].count(sample) != 0
                                           ? sample_metric_data_[0][sample]
                                           : frm::Metric();
    const frm::Metric &metric_data_1 = sample_metric_data_[1].count(sample) != 0
                                           ? sample_metric_data_[1][sample]
                                           : frm::Metric();
    // 使用try_emplace优化 - 只在键不存在时插入
    for (auto [name, value] : metric_data_0) {
      map_property_to_num_valid_data_.try_emplace(name, 0, 0);
      map_property_to_squared_difference_.try_emplace(name, 0.0);
      map_property_to_max_difference_.try_emplace(name, "", 0.0);
    }
    for (auto [name, value] : metric_data_1) {
      map_property_to_num_valid_data_.try_emplace(name, 0, 0);
      map_property_to_squared_difference_.try_emplace(name, 0.0);
      map_property_to_max_difference_.try_emplace(name, "", 0.0);
    }
  }

  for (auto sample : all_samples) {
    const frm::Metric &metric_data_0 = sample_metric_data_[0].count(sample) != 0
                                           ? sample_metric_data_[0][sample]
                                           : frm::Metric();
    const frm::Metric &metric_data_1 = sample_metric_data_[1].count(sample) != 0
                                           ? sample_metric_data_[1][sample]
                                           : frm::Metric();
    for (auto property_info : map_property_to_squared_difference_) {
      // 只计算double或者int类型的数值差异
      if (metrics.at(property_info.first) != "DOUBLE" &&
          metrics.at(property_info.first) != "INT")
        continue;
      // 如果该属性不齐全, 则记录为无效数值
      if (metric_data_0.count(property_info.first) == 0 ||
          metric_data_1.count(property_info.first) == 0) {
        map_property_to_num_valid_data_.at(property_info.first).second++;
        continue;
      }
      // 否则记录该属性的差异
      map_property_to_num_valid_data_.at(property_info.first).first++;
      const frm::ValueType &value_0 = metric_data_0.at(property_info.first);
      const frm::ValueType &value_1 = metric_data_1.at(property_info.first);
      double difference = 0;
      if (std::holds_alternative<double>(value_0) &&
          std::holds_alternative<double>(value_1))
        difference = std::get<double>(value_0) - std::get<double>(value_1);
      else if (std::holds_alternative<int>(value_0) &&
               std::holds_alternative<int>(value_1))
        difference = std::get<int>(value_0) - std::get<int>(value_1);
      // 更新squared_difference
      map_property_to_squared_difference_.at(property_info.first) +=
          std::pow(difference, 2);
      // 更新max_difference
      if (std::abs(difference) >
          map_property_to_max_difference_.at(property_info.first).second)
        map_property_to_max_difference_.at(property_info.first) =
            std::pair<std::string, double>(sample, std::abs(difference));
    }
  }
}

void Evaluator::PrintDataDifference() {
  for (auto property_info : map_property_to_squared_difference_) {
    std::size_t num_valid_data =
        map_property_to_num_valid_data_.at(property_info.first).first;
    std::size_t num_invalid_data =
        map_property_to_num_valid_data_.at(property_info.first).second;
    double squared_difference = std::sqrt(
        map_property_to_squared_difference_.at(property_info.first) /
        map_property_to_num_valid_data_.at(property_info.first).first);
    std::pair<std::string, double> max_difference =
        map_property_to_max_difference_.at(property_info.first);
    std::cout << std::format(
                     "metric name: \"{}\",\nvalid data: {}/{},\navg squared "
                     "difference: {},\n"
                     "max difference: ({}, {})",
                     property_info.first, //
                     num_valid_data,
                     num_valid_data + num_invalid_data, //
                     squared_difference,                //
                     max_difference.first, max_difference.second)
              << std::endl;
  }
}
