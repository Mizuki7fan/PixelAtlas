#include "metric.h"
#include <boost/json.hpp>
#include <fstream>
namespace frm {
void LoadMetricJsonFile(
    std::ifstream &json_str,                                     //
    const std::unordered_map<std::string, std::string> &metrics, //
    std::unordered_map<std::string, MetricValue> &metric_values) {
  // metrics标记每个度量的名字和类型
  // metric_values存储每个度量的值

  std::string json_file_str((std::istreambuf_iterator<char>(json_str)),
                            (std::istreambuf_iterator<char>()));
  auto jv = boost::json::parse(json_file_str);
  const boost::json::object &json_root = jv.as_object();
  for (const auto &[key, val] : json_root) {
    if (metrics.count(key) == 0)
      continue;
    if (metrics.at(key) == "double" || metrics.at(key) == "DOUBLE")
      metric_values[key] = val.as_double();
    else if (metrics.at(key) == "int" || metrics.at(key) == "INT")
      metric_values[key] = static_cast<int>(val.as_int64());
    else if (metrics.at(key) == "string" || metrics.at(key) == "STRING")
      metric_values[key] = static_cast<std::string>(val.as_string());
  }
}

void WriteMetricJsonFile(
    std::ofstream &json_file,                                    //
    const std::unordered_map<std::string, std::string> &metrics, //
    const std::unordered_map<std::string, MetricValue> &metric_values) {
  boost::json::object obj;
  for (auto [metric_name, metric_value] : metric_values) {
    std::cout << metric_name << std::endl;
    std::cout << metrics.at(metric_name) << std::endl;
    if (metrics.count(metric_name) == 0)
      continue;
    if (metrics.at(metric_name) == "double" ||
        metrics.at(metric_name) == "DOUBLE") {
      if (std::holds_alternative<double>(metric_value))
        obj[metric_name] = std::get<double>(metric_value);
    } else if (metrics.at(metric_name) == "int" ||
               metrics.at(metric_name) == "INT") {
      if (std::holds_alternative<int>(metric_value))
        obj[metric_name] = std::get<int>(metric_value);
    } else if (metrics.at(metric_name) == "string" ||
               metrics.at(metric_name) == "STRING") {
      if (std::holds_alternative<std::string>(metric_value))
        obj[metric_name] = std::get<std::string>(metric_value);
    }
  }

  std::string json_str = boost::json::serialize(obj);
  json_file << json_str;
}
}; // namespace frm
