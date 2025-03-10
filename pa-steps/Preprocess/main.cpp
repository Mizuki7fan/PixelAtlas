#include <frame/common_program.h>
#include <frame/parse.h>
#include <iostream>

// 函数执行入口
void MainProcess(std::string model_name) {
  std::cout << model_name << std::endl;
}

int main(int argc, char *argv[]) {

  frm::CommonProgram common_program(argc, argv);

  common_program.RunThreadParallel(MainProcess);

  std::cout << frm::GetGlobalWorkDir() << std::endl;
  return 0;
}