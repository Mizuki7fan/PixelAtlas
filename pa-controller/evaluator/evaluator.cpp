#include "evaluator.h"
#include <boost/program_options.hpp>
#include <frame/assert.hpp>
#include <iostream>

namespace po = boost::program_options;
// 结果评估模块
// 功能: 输入两次work的名称A,B和工具名称, 在run/evaluate_A_B文件夹下生成评估结果
int main(int argc, char *argv[]) {
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "帮助信息"); //

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
      std::cout << desc << std::endl;
      exit(0);
    }
    po::notify(vm);
  } catch (const po::error &e) {
    std::cerr << "error: " << e.what() << std::endl;
    std::cerr << desc << std::endl;
    PA_ASSERT_WITH_MSG(0, "输入异常");
  }
  return 0;
}