#pragma once
#include "global_defs.h"
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;
namespace frm {
std::ofstream CreateResultFilestream(const std::string &path);
std::ofstream CreateDebugFilestream(const std::string &path);
} // namespace frm
