#pragma once
#include <array>
#include <filesystem>
#include <frame/action_args.h>
#include <frame/metric.h>
#include <string>

namespace fs = std::filesystem;
static int debug_level_arg = 0;

class Evaluator {
public:
  Evaluator(int argc, char *argv[]);
  void LoadData();
  void PrintData();
  void AnalyseDataDifference();
  void PrintDataDifference();

private:
  void
  LoadDataFromNoIndividualModelDir(const std::vector<fs::path> &all_files, //
                                   int curr_work_idx);
  void LoadDataFromIndividualModelDir();
  std::vector<frm::ActionArguments> all_action_list_;
  std::array<std::unordered_map<std::string, frm::Metric>, 2>
      sample_metric_data_;
  std::string work_name_[2];
  int cmd_idx_ = -1;

  std::unordered_map<std::string, std::pair<std::size_t, std::size_t>>
      map_property_to_num_valid_data_;
  std::unordered_map<std::string, double> map_property_to_squared_difference_;
  std::unordered_map<std::string, std::pair<std::string, double>>
      map_property_to_max_difference_;
};
