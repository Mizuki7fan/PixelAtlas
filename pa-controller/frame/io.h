#pragma once
#include <fstream>

namespace frm {
std::ofstream CreateResultFilestream(const std::string &path);
std::ofstream CreateDebugFilestream(const std::string &path);

std::ofstream CreateMetricsFilestreamBegin();
void WriteMetricsFilestreamEnd(std::ofstream &fout);

// 写double类型的度量
void WriteMetrics(std::ofstream &fout, const std::string &name, double value);

} // namespace frm
