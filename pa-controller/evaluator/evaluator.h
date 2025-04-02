#pragma once
#include <string>
#include <unordered_map>
#include <variant>

using MetricValue = std::variant<double, int, std::string>;
using MetricName = std::string;

class SampleMetric {
public:
  void FromJson(const std::string &json_str, //
                const std::unordered_map<std::string, std::string> &metrics,
                int work_idx);

private:
  std::string sample_name;
  std::unordered_map<MetricName, MetricValue> sample_metrics_work_0,
      sample_metrics_work_1;
};
// 记录模型所有生成结果的数据结构
