#pragma once
#include <array>
#include <filesystem>
#include <frame/command_structures.h>
#include <frame/metric.h>
#include <string>

namespace fs = std::filesystem;
static int debug_level_arg = 0;

class Evaluator {
public:
  Evaluator(int argc, char *argv[]);
  void LoadData();
  void PrintData();

private:
  void
  LoadDataFromNoIndividualModelDir(const std::vector<fs::path> &all_files, //
                                   int curr_work_idx);
  void LoadDataFromIndividualModelDir();
  std::vector<StepArguments> all_step_list_;
  std::array<std::unordered_map<std::string, frm::Metric>, 2>
      sample_metric_data_;
  std::string work_name_[2];
  int cmd_idx_ = -1;
};
