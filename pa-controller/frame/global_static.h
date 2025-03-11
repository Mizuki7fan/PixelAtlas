#pragma once
#include <string>

namespace frm::global {
int DebugLevel();
std::string DataSet();
int NumParallelCount();
bool UseIndividualModelDir();
int MaxTimeElapsed();
std::string ParallelLevel();
} // namespace frm::global