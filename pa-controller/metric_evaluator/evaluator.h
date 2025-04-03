#pragma once
#include "metric_io.h"
#include <filesystem>
#include <frame/command_structures.h>
#include <string>
#include <unordered_map>

namespace fs = std::filesystem;
static int debug_level_arg = 0;

struct SampleMetric {
public:
  void FromJson(const std::string &json_str, //
                const std::unordered_map<std::string, std::string> &metrics,
                int work_idx);
  std::unordered_map<std::string, MetricValue> sample_metrics_work_0,
      sample_metrics_work_1;

private:
  std::string sample_name;
};

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
  std::unordered_map<std::string, SampleMetric> metric_data_;

  std::string work_name_[2];
  int cmd_idx_ = -1;
};
