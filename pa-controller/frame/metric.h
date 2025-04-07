#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <variant>

namespace frm {
struct PrintMetricValue {
  void operator()(double d) const { std::cout << d; }
  void operator()(int i) const { std::cout << i; }
  void operator()(const std::string &s) const { std::cout << s; }
};

// 度量的值, 允许double/int/std::string类型
using MetricValue = std::variant<double, int, std::string>;
// 度量, 包括度量名和值
using Metric = std::unordered_map<std::string, MetricValue>;

// 读入记录metric的json文件
void LoadMetricJsonFile(
    std::ifstream &json_str,                                     //
    const std::unordered_map<std::string, std::string> &metrics, //
    std::unordered_map<std::string, MetricValue> &metric_values);
// 写记录metric的json文件
void WriteMetricJsonFile(
    std::ofstream &json_str,                                     //
    const std::unordered_map<std::string, std::string> &metrics, //
    const std::unordered_map<std::string, MetricValue> &metric_values);
}; // namespace frm
