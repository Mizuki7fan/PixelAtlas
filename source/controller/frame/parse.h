#pragma once
#include "action_arguments.h"

namespace frm {
std::vector<std::string> SplitString(const std::string &full_string,
                                     const std::string &delimiter);

std::vector<ActionArguments> LoadAllStepList();
} // namespace frm
