#include "parse.h"
#include "assert.hpp"
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

std::vector<ActionArguments> LoadAllStepList() {
  // 使用boost/json解析json文件
  std::vector<ActionArguments> all_action_list;

  fs::path all_action_list_file_path =
      fs::current_path() / ".." / ".." / ACTIONLIST_FILE;

  std::cout << "action_list_file: " << all_action_list_file_path.string()
            << std::endl;

  std::ifstream file(all_action_list_file_path);

  PA_ASSERT_WITH_MSG(
      file.is_open(),
      std::format("未找到{}文件", all_action_list_file_path.string()));
  // 使用std::istreambuf_iterator将整个文件流转换为字符串
  std::string str((std::istreambuf_iterator<char>(file)),
                  (std::istreambuf_iterator<char>()));
  auto jv = boost::json::parse(str);

  // 验证jv是否是对象类型（即{}包裹的数据）
  if (!jv.is_object()) {
    throw std::runtime_error("Parsed value is not a JSON object.");
  }

  const boost::json::object &json_root = jv.as_object();
  all_action_list.resize(json_root.size());
  for (const auto &[key, val] : json_root) {
    PA_ASSERT(val.is_array());
    std::string action_name = key.data(); // 取步骤的名字
    if (val.is_array()) {
      const boost::json::array &action_array = val.get_array();
      for (auto action_value : action_array) {
        const boost::json::object &action_obj = action_value.get_object();
        int action_index = static_cast<int>(action_obj.at("idx").as_int64());
        ActionArguments &curr_action = all_action_list[action_index];
        curr_action.idx = action_index;
        curr_action.name = action_name;

        // 读输入文件名
        if (action_obj.contains("input")) {
          const boost::json::array &action_input_files =
              action_obj.at("input").get_array();
          for (auto value : action_input_files)
            curr_action.inputs.insert(value.as_string().c_str());
        }

        // 读输出文件名
        if (action_obj.contains("output")) {
          const boost::json::array &action_output_files =
              action_obj.at("output").get_array();
          for (auto value : action_output_files)
            curr_action.outputs.insert(value.as_string().c_str());
        }

        if (action_obj.contains("metrics")) {
          const boost::json::array &action_metrics =
              action_obj.at("metrics").get_array();
          for (auto value : action_metrics) {
            const boost::json::object &metric_obj = value.get_object();
            std::string metric_type = metric_obj.at("type").as_string().c_str();
            std::string metric_name = metric_obj.at("name").as_string().c_str();
            curr_action.metrics[metric_name] = metric_type;
          }
        }
      }
    }
  }
  return all_action_list;
}

} // namespace frm