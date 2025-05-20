#pragma once

#include "io.h"
#include <iostream>
#include <string>
#include <unordered_map>

namespace frm {
struct PrintValueType {
  void operator()(double d) const { std::cout << d; }
  void operator()(int i) const { std::cout << i; }
  void operator()(const std::string &s) const { std::cout << s; }
};
// 度量, 包括度量名和值
using Metric = std::unordered_map<std::string, ValueType>;

// 读入记录metric的json文件
void LoadMetricJsonFile(
    std::ifstream &json_str,                                     //
    const std::unordered_map<std::string, std::string> &metrics, //
    std::unordered_map<std::string, ValueType> &metric_values);
// 写记录metric的json文件
void WriteMetricJsonFile(
    std::ofstream &json_str,                                     //
    const std::unordered_map<std::string, std::string> &metrics, //
    const std::unordered_map<std::string, ValueType> &metric_values);
}; // namespace frm
