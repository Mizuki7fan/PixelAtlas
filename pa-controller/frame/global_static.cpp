#include "global_static.h"
#include "common_program.h"

namespace frm::global {
int DebugLevel() { return frm::GetDebugLevel(); }
std::string DataSet() { return frm::GetDataset(); }
int NumParallelCount() { return frm::GetNumParallelCnt(); }
bool UseIndividualModelDir() { return frm::UseIndividualModelDir(); }
int MaxTimeElapsed() { return frm::GetMaxTimeElapsed(); }
std::string ParallelLevel() { return frm::GetParallelLevel(); }
} // namespace frm::global