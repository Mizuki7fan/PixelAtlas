#pragma once

#include "assert.hpp"
#include "command_structures.h"

namespace frm {
void ParseMultiStepArgs(std::stringstream &all_input_args,
                        std::size_t &begin_cmd_idx, std::size_t &end_cmd_idx,
                        InputArguments &input_args);
// 解析单步的输入参数
InputArguments ParseSingleStepArgs(std::stringstream &all_input_args);

std::vector<std::string> SplitString(const std::string &full_string,
                                     const std::string &delimiter);

std::vector<StepArguments> LoadAllStepList();
} // namespace frm
