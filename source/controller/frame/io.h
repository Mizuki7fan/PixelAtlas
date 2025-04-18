#pragma once
#include <fstream>

namespace frm {
std::ofstream CreateResultFilestream(const std::string &path);
std::ofstream CreateDebugFilestream(const std::string &path);
std::ofstream CreateMetricsFilestream();
} // namespace frm
