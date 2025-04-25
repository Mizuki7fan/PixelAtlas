#pragma once
#include "action_args.h"
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

namespace frm::global {

int DebugLevel();
fs::path ActionDebugDir();
fs::path ActionResultDir();
std::string DatasetName();
fs::path InstancePath();
int MaxTimeElapsed();
bool UseIndividualInstanceDir();
std::string WorkName();
const ActionArguments &ActionArgs();
} // namespace frm::global