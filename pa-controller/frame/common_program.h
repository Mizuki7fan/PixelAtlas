#pragma once
#include "command_structures.h"
#include "parse.h"
#include <functional>
#include <string>
#include <vector>

namespace frm {
class CommonProgram {
public:
  CommonProgram(int argc, char *argv[]);
  int Run(const std::function<void(std::string model_name)> &func) const;

private: // functions
  bool CheckProgramNameValidity() const;

private:
  std::vector<StepArguments> all_step_list;
  std::vector<std::string> all_args;
  std::vector<std::string> all_models;
};
} // namespace frm
