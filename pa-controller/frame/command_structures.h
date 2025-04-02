#pragma once
#include <string>
#include <unordered_map>
#include <unordered_set>

// all-steps.json文件的内容
struct StepArguments { // 表示步骤的各种参数
  int step_idx = -1;   //
  std::string step_name = "";
  std::unordered_set<std::string> input_files;  //
  std::unordered_set<std::string> output_files; //
  std::unordered_map<std::string, std::string> metrics;
};