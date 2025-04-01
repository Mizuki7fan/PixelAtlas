#include "parse.h"
#include "global_defs.h"
#include <boost/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

namespace frm {
// 该函数应该放到parse.h中
std::vector<std::string> SplitString(const std::string &full_string,
                                     const std::string &delimiter) {
  std::vector<std::string> parts;
  size_t start = 0;
  size_t end = full_string.find(delimiter);
  while (end != std::string::npos) {
    parts.push_back(full_string.substr(start, end - start));
    start = end + delimiter.length();
    end = full_string.find(delimiter, start);
  }
  parts.push_back(full_string.substr(start, end));
  return parts;
}

std::vector<StepArguments> LoadAllStepList() {
  // 使用boost/json解析json文件
  std::vector<StepArguments> all_step_list;

  fs::path steplist_file_path =
      fs::current_path() / ".." / ".." / STEPLIST_FILE;

  if (global::DebugLevel() > 0) {
    std::cout << "steplist_file: " << steplist_file_path.string() << std::endl;
  }

  std::ifstream file(steplist_file_path);

  PA_ASSERT_WITH_MSG(file.is_open(), "没有找到all_step_list.json文件");
  // 使用std::istreambuf_iterator将整个文件流转换为字符串
  std::string str((std::istreambuf_iterator<char>(file)),
                  (std::istreambuf_iterator<char>()));
  auto jv = boost::json::parse(str);

  // 验证jv是否是对象类型（即{}包裹的数据）
  if (!jv.is_object()) {
    throw std::runtime_error("Parsed value is not a JSON object.");
  }

  const boost::json::object &json_root = jv.as_object();
  all_step_list.resize(json_root.size());
  for (const auto &[key, val] : json_root) {
    PA_ASSERT(val.is_array());
    std::string step_name = key.data(); // 取步骤的名字
    if (val.is_array()) {
      const boost::json::array &step_array = val.get_array();
      for (auto step_value : step_array) {
        const boost::json::object &step_obj = step_value.get_object();
        int step_index = static_cast<int>(step_obj.at("step_idx").as_int64());
        StepArguments &curr_step = all_step_list[step_index];
        curr_step.step_idx = step_index;
        curr_step.step_name = step_name;

        // 读输入文件名
        const boost::json::array &step_input_files =
            step_obj.at("input_files").get_array();
        for (auto value : step_input_files)
          curr_step.input_files.insert(value.as_string().c_str());

        // 读输出文件名
        const boost::json::array &step_output_files =
            step_obj.at("output_files").get_array();
        for (auto value : step_output_files)
          curr_step.output_files.insert(value.as_string().c_str());
        const boost::json::array &step_metrics =
            step_obj.at("metrics").get_array();
        for (auto value : step_metrics)
          curr_step.metrics.push_back(value.as_string().c_str());
      }
    }
  }
  return all_step_list;
}

} // namespace frm