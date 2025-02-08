#pragma once

#include "assert.hpp"
#include "command_structures.h"

namespace frm {
InputCommandArguments ParseInputArgs(std::stringstream &all_input_args);

std::vector<StepArguments> LoadAllStepList();
} // namespace frm
