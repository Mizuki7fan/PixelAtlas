#pragma once

#include "assert.hpp"
#include "command_structures.h"

namespace frm {
InputCommandArguments ParseInputArgs(std::stringstream &all_input_args);

std::vector<StepArguments> ParseAllStepList();

// 读取每个step需要批量执行的所有指令
void ParseRunCommandListsOfSingleStep(std::stringstream &all_input_args_stream,
                                      std::vector<std::string> &run_cmd_list);

} // namespace frm
