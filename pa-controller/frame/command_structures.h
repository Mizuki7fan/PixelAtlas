#pragma once
#include <iostream>
#include <string.h>
#include <unordered_map>
#include <unordered_set>

struct InputArguments {                    // 表示控制程序的输入指令
  int num_parallel_cnt = 1;                // 并行数
  int num_debug_level = 0;                 // 调试等级
  std::string filename_regex_str = ".*.*"; // 文件名通配符
  std::string dataset_str = ".";           // 数据集名
  int max_time_elapsed = 1800;             // 用时
};

// all-steps.json文件的内容
struct StepArguments { // 表示步骤的各种参数
  int step_idx = -1;   //
  std::string step_name = "";
  std::unordered_set<std::string> input_files;  //
  std::unordered_set<std::string> output_files; //
};