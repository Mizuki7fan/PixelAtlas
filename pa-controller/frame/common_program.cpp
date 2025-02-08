#include "common_program.h"
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

namespace frm {

bool CommonProgram::CheckProgramNameValidity() const {
  fs::path program_path = all_args[0];                     // 取exe路径
  std::string program_name = program_path.stem().string(); // 取exe名
  for (std::size_t step_idx = 0; step_idx < all_step_list.size(); ++step_idx)
    if (all_step_list[step_idx].step_name == program_name)
      return true;
  return false;
}

CommonProgram::CommonProgram(int argc, char *argv[]) {
  for (int i = 0; i < argc; ++i) {
    all_args.push_back(argv[i]);
    std::cout << all_args.back() << " ";
  }
  std::cout << std::endl;

  // 加载all-step.json文件
  all_step_list = LoadAllStepList();
  PA_ASSERT(CheckProgramNameValidity());
};

int CommonProgram::Run(const std::function<void(std::string)> &func) const {

  // 每个工具将要执行的函数传到这里, 该函数负责进行输入控制以及并行控制
  std::string model_name = "alien.obj";
  func(model_name);
  return 1;
}
} // namespace frm
