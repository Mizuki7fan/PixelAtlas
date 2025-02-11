#include "common_program.h"
#include <filesystem>
#include <iostream>
// #include <queue>
#include <regex>

namespace fs = std::filesystem;

namespace frm {

std::size_t CommonProgram::GetCurrentProgramIndex() {

  fs::path program_path = all_args[0];                     // 取exe路径
  std::string program_name = program_path.stem().string(); // 取exe名

  // 取当前program的名字
  if (map_step_name_to_step_idx.contains(program_name))
    curr_cmd_idx = map_step_name_to_step_idx.at(program_name);

  return curr_cmd_idx;
}

// 取当前工具的前置依赖工具
std::unordered_set<std::size_t> CommonProgram::GetCurrentProgramDependencies() {
  std::unordered_set<std::size_t> denpendcies;
  // std::queue<std::size_t> queue_denpendcies;
  // queue_denpendcies.push(curr_cmd_idx);

  // while (!queue_denpendcies.empty()) {
  //   std::size_t cmd_idx = queue_denpendcies.front();
  //   denpendcies.insert(cmd_idx);
  //   // 取当前cmd_idx的依赖idx
  //   for (auto depend_cmd_name : all_step_list[cmd_idx].dependencies) {
  //     std::size_t depend_cmd_idx =
  //     map_step_name_to_step_idx[depend_cmd_name]; PA_ASSERT(depend_cmd_idx <
  //     cmd_idx); queue_denpendcies.push(depend_cmd_idx);
  //   }
  // };
  return denpendcies;
}

bool CommonProgram::PrepareWorkingDirectory() {
  // 准备当前项目运行需要依赖的文件夹
  fs::path root_work_dir = fs::current_path() / "..work";
  std::unordered_map<std::size_t, std::string> existing_subdirs;

  if (!fs::exists(root_work_dir)) {
    fs::create_directories(root_work_dir); // 创建work文件夹
  } else {
    for (fs::directory_iterator iter(root_work_dir);
         iter != fs::directory_iterator(); ++iter) {
      std::string subdir_name = iter->path().filename().string();
      std::vector<std::string> parts =
          SplitString(subdir_name, "_"); // 根据"_"的前后拆分字符串
      if (parts.size() != 2)
        continue;
      std::size_t cmd_idx = static_cast<std::size_t>(std::stoi(parts[0]));
      existing_subdirs[cmd_idx] = parts[1];
    }
  }
  // 扫描work文件夹下的所有子文件夹

  std::unordered_set<std::size_t> program_dependencies =
      GetCurrentProgramDependencies();
  // 确定当前工具依赖的前置工具
  // 检测相应的文件夹是否存在

  // 当前工具运行时所缺乏的前置工具
  std::unordered_set<std::size_t> lacking_programs;

  for (auto depend_cmd_idx : program_dependencies) {
    std::cout << depend_cmd_idx << std::endl;
    if (!existing_subdirs.contains(depend_cmd_idx))
      lacking_programs.insert(depend_cmd_idx);
    else if (existing_subdirs[depend_cmd_idx] !=
             all_step_list[depend_cmd_idx].step_name)
      lacking_programs.insert(depend_cmd_idx);
  }

  if (lacking_programs.size() != 0) {
    std::cout << "前置程序的结果缺失:" << std::endl;
    for (auto lacking : lacking_programs)
      std::cout << std::format("{}_{}", lacking,
                               all_step_list[lacking].step_name)
                << std::endl;
    return false;
  }

  // 否则: 新建当前工具的结果文件夹, 并统计本次运行一共需要跑的例子数量
  fs::path curr_working_dir =
      root_work_dir /
      std::format("{}_{}", curr_cmd_idx, all_step_list[curr_cmd_idx].step_name);
  if (fs::exists(curr_working_dir)) {
    fs::remove_all(curr_working_dir);
  }
  fs::create_directories(curr_working_dir);

  // 取当前项目需要运行的所有文件

  return true;
}

bool CommonProgram::SelectRunTargets() {
  //
  std::smatch match;
  std::string &file_regex = input_args.filename_regex_str;
  std::regex pattern(file_regex);

  if (curr_cmd_idx == 0) {
    fs::path project_root_dir = fs::current_path() / ".." / "..";
    fs::path project_asset_dir = project_root_dir / "asset";
    fs::path run_data_dir = project_asset_dir / input_args.dataset_str;

    std::vector<std::string> run_files;
    for (const auto &entry :
         std::filesystem::directory_iterator(run_data_dir)) {
      if (entry.is_directory())
        continue;
      //
      std::string filename = entry.path().filename().string();
      if (std::regex_match(filename, match, pattern)) {
        run_files.push_back(filename);
      }
    }
    std::cout << std::format("共成功匹配{}个例子:", run_files.size())
              << std::endl;
    for (auto run_file : run_files)
      std::cout << run_file << std::endl;
  }

  return true;
}

CommonProgram::CommonProgram(int argc, char *argv[]) {
  for (int i = 0; i < argc; ++i) {
    all_args.push_back(argv[i]);
    std::cout << all_args.back() << " ";
  }
  std::cout << std::endl;

  // 加载all-step.json文件
  all_step_list = LoadAllStepList();
  for (std::size_t cmd_idx = 0; cmd_idx < all_step_list.size(); ++cmd_idx)
    map_step_name_to_step_idx[all_step_list[cmd_idx].step_name] = cmd_idx;

  curr_cmd_idx = std::numeric_limits<std::size_t>::max();

  curr_cmd_idx = GetCurrentProgramIndex();
  if (curr_cmd_idx == std::numeric_limits<std::size_t>::max())
    std::cerr << "invalid command: " << argv[0] << std::endl;

  std::stringstream all_input_cmd_args_stream;
  for (int i = 0; i < argc; ++i)
    all_input_cmd_args_stream << argv[i] << " ";

  input_args = ParseSingleStepArgs(all_input_cmd_args_stream);

  PrepareWorkingDirectory();
  SelectRunTargets();
};

int CommonProgram::Run(const std::function<void(std::string)> &func) const {
  // 清空当前步骤的结果

  // 每个工具将要执行的函数传到这里, 该函数负责进行输入控制以及并行控制
  std::string model_name = "alien.obj";
  func(model_name);
  return 1;
}
} // namespace frm
