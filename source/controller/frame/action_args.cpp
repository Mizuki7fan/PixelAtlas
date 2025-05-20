#include "action_args.h"
#include "pa_assert.hpp"
#include <boost/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

namespace frm {
std::vector<ActionArguments> LoadAllActionList() {
  // 使用boost/json解析json文件
  std::vector<ActionArguments> all_action_list;

  fs::path all_action_list_file_path =
      fs::current_path() / ".." / ".." / ACTIONLIST_FILE_NAME;

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
        // 读入工具的参数
        if (action_obj.contains("hyper-parameter")) {
          const boost::json::array &action_parameter =
              action_obj.at("hyper-parameter").get_array();
          for (auto value : action_parameter) {
            const boost::json::object &parameter_obj = value.get_object();
            std::string parameter_type =
                parameter_obj.at("type").as_string().c_str();
            std::string parameter_name =
                parameter_obj.at("name").as_string().c_str();
            const auto &parameter_value = parameter_obj.at("value");
            if (parameter_type == "INT") {
              PA_ASSERT_WITH_MSG(
                  parameter_value.is_number(),
                  std::format("参数{}数据类型异常", parameter_name));
              // as_int64返回的是int64_t, 即long long, 需要强制类型转换
              curr_action.hpyer_parameters.emplace(
                  parameter_name, static_cast<int>(parameter_value.as_int64()));
            } else if (parameter_type == "DOUBLE") {
              PA_ASSERT_WITH_MSG(
                  parameter_value.is_double(),
                  std::format("参数{}数据类型异常", parameter_name));
              curr_action.hpyer_parameters.emplace(
                  parameter_name, parameter_obj.at("value").as_double());
            } else if (parameter_type == "STRING") {
              PA_ASSERT_WITH_MSG(
                  parameter_value.is_string(),
                  std::format("参数{}数据类型异常", parameter_name));
              curr_action.hpyer_parameters.emplace(
                  parameter_name,
                  parameter_obj.at("value").as_string().c_str());
            }
          }
        }
      }
    }
  }
  return all_action_list;
}
} // namespace frm
