#pragma once
#include "io.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <variant>

namespace frm {
// all-steps.json文件的内容
struct ActionArguments { // 表示步骤的各种参数
  int idx = -1;          //
  std::string name = "";
  std::unordered_set<std::string> inputs;               // 当前action的输入文件
  std::unordered_set<std::string> outputs;              // 当前action的输出文件
  std::unordered_map<std::string, std::string> metrics; // 当前action的度量
  std::unordered_map<std::string, ValueType>
      hpyer_parameters; // 当前action的输入参数
};

std::vector<ActionArguments> LoadAllActionList();
} // namespace frm
