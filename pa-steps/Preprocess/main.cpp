#include "process.h"
#include <frame/common_program.h>
#include <iostream>

int main(int argc, char *argv[]) {
  std::cout << argv[0] << " " << argc << std::endl;
  for (int i = 0; i < argc; ++i)
    std::cout << argv[i] << " ";
  std::cout << std::endl;
  frm::CommonProgram common_program(argc, argv);

  common_program.Run(MainProcess);

  return 0;
}