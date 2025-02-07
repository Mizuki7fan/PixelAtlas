#include "common_program.h"
#include <iostream>

namespace frm {
CommonProgram::CommonProgram(int argc, char *argv[]) {
  for (int i = 0; i < argc; ++i) {
    all_args.push_back(argv[i]);
    std::cout << all_args.back() << " ";
  }
  std::cout << std::endl;
};

int CommonProgram::Run(const std::function<void()> &func) const {
  func();
  return 1;
}
} // namespace frm
