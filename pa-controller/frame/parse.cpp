#include "parse.h"
#include <boost/json.hpp>
#include <fstream>
#include <sstream>
#include <string>

namespace frm {
// 该函数应该放到parse.h中
void ParseMultiStepArgs(std::stringstream &all_input_args,
                        std::size_t &begin_cmd_idx, std::size_t &end_cmd_idx,
                        InputArguments &input_args) {
  std::string tag, value;
  std::string exe_name;
  all_input_args >> exe_name;
  while (all_input_args >> tag >> value) {
    if (tag == "-b") { //
      begin_cmd_idx = std::stoi(value);
    } else if (tag == "-e") {
      end_cmd_idx = std::stoi(value);
    } else if (tag == "-p") {
      input_args.num_parallel_cnt = std::stoi(value);
    } else if (tag == "-d") {
      input_args.num_debug_level = std::stoi(value);
    } else if (tag == "-f") {
      input_args.filename_regex_str = value;
    } else if (tag == "-ds") {
      input_args.dataset_str = value;
    } else if (tag == "-t") {
      input_args.max_time_elapsed = std::stoi(value);
    } else {
      PA_ASSERT_WITH_MSG(0, std::format("错误的tag: {}", tag));
    }
  }
  PA_ASSERT(begin_cmd_idx <= end_cmd_idx);
  std::cout
      << std::format(
             "{} start_cmd_idx: {}, end_cmd_idx: {}, "
             "parallel_cnt: {}, debug_level: {}, file_regex: {}, dataset: {}",
             exe_name, begin_cmd_idx, end_cmd_idx, input_args.num_parallel_cnt,
             input_args.num_debug_level, input_args.filename_regex_str,
             input_args.dataset_str)
      << std::endl;
}

InputArguments ParseSingleStepArgs(std::stringstream &all_input_args) {
  InputArguments all_cmd_args;
  std::string tag, value;
  std::string exe_name;
  all_input_args >> exe_name;
  while (all_input_args >> tag >> value) {
    if (tag == "-p") {
      all_cmd_args.num_parallel_cnt = std::stoi(value);
    } else if (tag == "-d") {
      all_cmd_args.num_debug_level = std::stoi(value);
    } else if (tag == "-f") {
      all_cmd_args.filename_regex_str = value;
    } else if (tag == "-ds") {
      all_cmd_args.dataset_str = value;
    } else if (tag == "-t") {
      all_cmd_args.max_time_elapsed = std::stoi(value);
    } else {
      PA_ASSERT_WITH_MSG(0, std::format("错误的tag: {}", tag));
    }
  }
  std::cout
      << std::format(
             "{} , "
             "parallel_cnt: {}, debug_level: {}, file_regex: {}, dataset: {}",
             exe_name, all_cmd_args.num_parallel_cnt,
             all_cmd_args.num_debug_level, all_cmd_args.filename_regex_str,
             all_cmd_args.dataset_str)
      << std::endl;

  return all_cmd_args;
}

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
  std::ifstream file(STEPLIST_FILE);

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
      }
    }
  }
  return all_step_list;
}

} // namespace frm