#include "global_defs.h"
#include "common_program.h"

namespace frm::global {
int DebugLevel() { return frm::GetDebugLevel(); }

fs::path CurrDebugDir() { return frm::GetCurrDebugDir(); }
fs::path CurrResultDir() { return frm::GetCurrResultDir(); }
fs::path CurrFile() { return frm::GetCurrFile(); }
std::string DatasetStr() { return frm::GetDatasetStr(); }
int MaxTimeElapsed() { return frm::GetMaxTimeElapsed(); }
bool UseIndividualModelDir() { return frm::GetUseIndividualModelDir(); }
std::string RunName() { return frm::GetRunName(); }
} // namespace frm::global