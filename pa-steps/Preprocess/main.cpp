#include <frame/common_program.h>
#include <frame/parse.h>
#include <iostream>

void MainProcess(std::string model_name) {
  std::cout << model_name << std::endl;
}

int main(int argc, char *argv[]) {

  frm::CommonProgram common_program(argc, argv);

  common_program.Run(MainProcess);
  return 0;
}