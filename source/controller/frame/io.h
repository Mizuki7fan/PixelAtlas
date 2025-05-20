#pragma once
#include <fstream>
#include <variant>
namespace frm {
using ValueType = std::variant<double, int, std::string>;

std::ofstream CreateOutputFilestream(const std::string &path);
std::ofstream CreateDebugFilestream(const std::string &path);
std::ofstream CreateMetricsFilestream();
// 读超参数
ValueType GetHyperParameter(std::string name);
std::ifstream GetInputFilestream(const std::string &path);

} // namespace frm
