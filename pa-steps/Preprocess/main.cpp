#include <frame/common_program.h>
#include <frame/parse.h>
#include <iostream>

static void MainProcess() { std::cout << "MainProcess" << std::endl; }

int main(int argc, char *argv[]) {

  frm::CommonProgram common_program(argc, argv);

  common_program.Run(MainProcess);
  return 0;
}