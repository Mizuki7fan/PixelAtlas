#pragma once
#include <functional>
#include <string>
#include <vector>

namespace frm {
class CommonProgram {
public:
  CommonProgram(int argc, char *argv[]);
  int Run(const std::function<void()> &func) const;

private:
  std::vector<std::string> all_args;
};
} // namespace frm
