#include "global_defs.h"
#include "common_program.h"

namespace frm::global {

int DebugLevel() { return frm::GetDebugLevel(); }

fs::path ActionDebugDir() { return frm::GetActionDebugDir(); }
fs::path ActionResultDir() { return frm::GetActionResultDir(); }
fs::path InstancePath() { return frm::GetInstancePath(); }
std::string DatasetName() { return frm::GetDatasetName(); }
int MaxTimeElapsed() { return frm::GetMaxTimeElapsed(); }
bool UseIndividualInstanceDir() { return frm::GetUseIndividualInstanceDir(); }
std::string WorkName() { return frm::GetWorkName(); }
} // namespace frm::global