#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <variant>

struct PrintMetricValue {
  void operator()(double d) const { std::cout << d; }
  void operator()(int i) const { std::cout << i; }
  void operator()(const std::string &s) const { std::cout << s; }
};

using MetricValue = std::variant<double, int, std::string>;

// 读入记录metric的json文件
void LoadMetricJsonFile(
    std::ifstream &json_str,                                     //
    const std::unordered_map<std::string, std::string> &metrics, //
    std::unordered_map<std::string, MetricValue> &metric_values);

void WriteMetricJsonFile(
    std::ofstream &json_str,                                     //
    const std::unordered_map<std::string, std::string> &metrics, //
    const std::unordered_map<std::string, MetricValue> &metric_values);
