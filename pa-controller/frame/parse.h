#pragma once

#include "assert.hpp"
#include "command_structures.h"

namespace frm {
std::vector<std::string> SplitString(const std::string &full_string,
                                     const std::string &delimiter);

std::vector<StepArguments> LoadAllStepList();
} // namespace frm
