#include "evaluator.h"
#include <boost/json.hpp>

#include <filesystem>
#include <frame/assert.hpp>
#include <frame/command_structures.h>
#include <frame/global_defs.h>
#include <frame/parse.h>
#include <fstream>
#include <iostream>

static int debug_level_arg = 0;
int DebugLevel() { return debug_level_arg; }

void SampleMetric::FromJson(
    const std::string &json_file,                                //
    const std::unordered_map<std::string, std::string> &metrics, //
    int work_idx) {

  std::unordered_map<MetricName, MetricValue> &sample_metrics =
      work_idx == 0 ? sample_metrics_work_0 : sample_metrics_work_1;
  std::ifstream file(json_file);
  if (!file.is_open()) {
    throw std::runtime_error("");
  }

  std::string json_file_str((std::istreambuf_iterator<char>(file)),
                            (std::istreambuf_iterator<char>()));
  if (debug_level_arg > 0) {
    std::cout << json_file << std::endl;
    std::cout << json_file_str << std::endl;
  }
  // boost:json::parse的输入是一个字符串, 该字符串是文件内容;
  auto jv = boost::json::parse(json_file_str);
  if (!jv.is_object()) {
    throw std::runtime_error("Parsed value is not a JSON object.");
  }
  const boost::json::object &json_root = jv.as_object();
  for (const auto &[key, val] : json_root) {
    std::cout << key << std::endl;
    // 识别读入的key, 从而确定val的类型
    if (debug_level_arg > 0) {
      std::cout << metrics.count(key) << " " << metrics.at(key) << std::endl;
    }
    if (metrics.count(key) == 0)
      continue;
    if (metrics.at(key) == "double" || metrics.at(key) == "DOUBLE")
      sample_metrics[key] = val.as_double();
    else if (metrics.at(key) == "int" || metrics.at(key) == "INT")
      sample_metrics[key] = static_cast<int>(val.as_int64());
    else if (metrics.at(key) == "string" || metrics.at(key) == "STRING")
      sample_metrics[key] = static_cast<std::string>(val.as_string());
  }
}