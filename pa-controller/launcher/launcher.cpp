#include <boost/process.hpp>
#include <filesystem>
#include <format>
#include <frame/assert.hpp>
#include <frame/command_structures.h>
#include <frame/parse.h>
#include <iostream>
#include <regex>
#include <sstream>

namespace fs = std::filesystem;
namespace bp = boost::process;

int main(int argc, char *argv[]) {
  // 查看项目的工作目录, 该值通过launch.json的cwd字段指定
  std::cout << "工作目录: " << fs::current_path() << std::endl;
  PA_ASSERT_WITH_MSG(fs::exists(STEPLIST_FILE), "宏STEPLIST_FILE不存在");
  PA_ASSERT_WITH_MSG(fs::exists("jsoncpp.dll"), "jsoncpp.dll不存在");

  // 1. 读取并解析输入命令行参数

  std::cout << "argc: " << argc << std::endl;
  std::cout << "输入参数: " << std::endl;
  for (int i = 0; i < argc; ++i)
    std::cout << argv[i] << " ";
  std::cout << std::endl;
  std::stringstream all_input_cmd_args_stream;
  for (int i = 0; i < argc; ++i)
    all_input_cmd_args_stream << argv[i] << " ";
  // 解析命令行参数
  InputCommandArguments all_input_cmd_args =
      frm::ParseInputArgs(all_input_cmd_args_stream);

  // 2. 读取并解析所有工具参数
  std::vector<StepArguments> all_step_list = frm::LoadAllStepList();

  // 3. 读取并处理work文件夹
  fs::path work_dir = fs::current_path() / ".." / "work"; // 拼接
  std::cout << work_dir << std::endl;
  if (fs::exists(work_dir)) {
    for (fs::directory_iterator iter(work_dir);
         iter != fs::directory_iterator(); ++iter) {
      std::string dir_name = iter->path().filename().string();
      std::cout << dir_name << std::endl;

      std::smatch match; // string-match
      std::regex pattern(R"((\d+)-(.+))");
      //\d表示任意一个数字字符, 等效于[0-9]
      //\d+表示该模式可重复多次
      //()表示捕获分组, 用()括起来的内容会被作为一个整体, 通过match[xx]即可获取
      //.表示任意单个字符
      //+表示该模式可重复多次

      if (std::regex_match(dir_name, match, pattern)) {
        // 匹配成功
        std::cout << std::format("{}匹配成功, 拆成{}和{}", dir_name,
                                 match[1].str(), match[2].str())
                  << std::endl;
        int cmd_idx = std::stoi(match[1].str());
        if (cmd_idx >= all_input_cmd_args.begin_step_idx &&
            cmd_idx <= all_input_cmd_args.end_step_idx) {
          std::cout << std::format("删除路径: {}", iter->path().string())
                    << std::endl;
          fs::remove_all(iter->path()); // 删除文件夹
        }
      } else {
        std::cout << std::format("{}匹配失败", dir_name) << std::endl;
      }
    }
  } else {
    fs::create_directories(work_dir);
  }

  // 4. 进程管理
  for (int cmd_idx = all_input_cmd_args.begin_step_idx;
       cmd_idx <= all_input_cmd_args.end_step_idx; ++cmd_idx) {
    // 依次执行每个步骤
    std::cout << cmd_idx << " " << all_step_list.size() << std::endl;
    StepArguments &curr_step = all_step_list[cmd_idx];
    std::string cmd_filename = curr_step.step_name + ".exe";
    std::cout << "cmd_name: " << cmd_filename << std::endl;
    std::vector<std::string> args;
    args.push_back(std::format("-j{}", all_input_cmd_args.num_parallel_cnt));
    args.push_back(std::format("-d{}", all_input_cmd_args.num_debug_level));
    args.push_back(std::format("-f{}", all_input_cmd_args.filename_regex_str));
    args.push_back(std::format("-t{}", all_input_cmd_args.max_time_elapsed));
    args.push_back(std::format("-ds{}", all_input_cmd_args.dataset_str));

    bp::child process;
    auto env = boost::this_process::environment();
    try {
      process = bp::child(cmd_filename, bp::args(args), env);
      process.wait();
    } catch (const std::exception &e) {
      std::cout << std::format("进程{}异常: {}", cmd_filename, e.what())
                << std::endl;
      exit(0);
    } catch (...) {
      std::cout << std::format("进程{}异常", cmd_filename) << std::endl;
      exit(0);
    }
  }
}