#include "parse.h"
#include <fstream>
#include <json/json.h> //jsoncpp
#include <sstream>
#include <string>

namespace frm {
// 该函数应该放到parse.h中
InputCommandArguments ParseInputArgs(std::stringstream &all_input_args) {
  InputCommandArguments all_cmd_args;
  std::string tag, value;
  std::string exe_name;
  all_input_args >> exe_name;
  while (all_input_args >> tag >> value) {
    if (tag == "-b") { //
      all_cmd_args.begin_step_idx = std::stoi(value);
    } else if (tag == "-e") {
      all_cmd_args.end_step_idx = std::stoi(value);
    } else if (tag == "-p") {
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
  PA_ASSERT(all_cmd_args.begin_step_idx <= all_cmd_args.end_step_idx);
  std::cout
      << std::format(
             "{} start_cmd_idx: {}, end_cmd_idx: {}, "
             "parallel_cnt: {}, debug_level: {}, file_regex: {}, dataset: {}",
             exe_name, all_cmd_args.begin_step_idx, all_cmd_args.end_step_idx,
             all_cmd_args.num_parallel_cnt, all_cmd_args.num_debug_level,
             all_cmd_args.filename_regex_str, all_cmd_args.dataset_str)
      << std::endl;

  return all_cmd_args;
}

std::vector<StepArguments> ParseAllStepList() {
  std::vector<StepArguments> all_step_list;
  std::ifstream in_stream(STEPLIST_FILE);
  Json::CharReaderBuilder reader_builder;
  Json::Value json_root;
  std::string errs;
  bool parsing_successful_flag =
      Json::parseFromStream(reader_builder, in_stream, &json_root, &errs);
  std::cout << std::format("parsing_successful: {}", parsing_successful_flag)
            << std::endl;
  PA_ASSERT(parsing_successful_flag);
  all_step_list.resize(json_root.size());

  for (auto &json_cmd : json_root.getMemberNames()) {

    const Json::Value &cmd_info = json_root[json_cmd][0];
    int curr_step_idx = cmd_info["step_idx"].asInt();
    StepArguments &curr_step = all_step_list[curr_step_idx];

    curr_step.step_name = json_cmd;     // 命令的名字
    curr_step.step_idx = curr_step_idx; // 命令的序号
    for (const auto &input_file : cmd_info["input_files"]) {
      curr_step.input_files.insert(input_file.asString()); // 命令的输入文件
    }
    for (const auto &output_file : cmd_info["output_files"]) {
      curr_step.output_files.insert(output_file.asString()); // 命令的输出文件
    }
  }

  return all_step_list;
}

// 根据指令筛选本次运行需要执行的所有命令
void ParseRunCommandListsOfSingleStep(std::stringstream &all_input_args_stream,
                                      std::vector<std::string> &run_cmd_list) {
  InputCommandArguments all_input_cmd_args =
      ParseInputArgs(all_input_args_stream);
  std::vector<StepArguments> all_step_list = ParseAllStepList();
  // 已知: 工具名, 其他参数, 和step_list, 返回需要执行的每一条命令
  std::cout << run_cmd_list.size() << std::endl;
}
} // namespace frm