#pragma once
#include <fstream>

namespace frm {
std::ofstream CreateOutputFilestream(const std::string &path);
std::ofstream CreateDebugFilestream(const std::string &path);
std::ofstream CreateMetricsFilestream();
} // namespace frm
