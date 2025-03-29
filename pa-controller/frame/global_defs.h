#pragma once
#include <filesystem>
namespace fs = std::filesystem;

namespace frm::global {

int DebugLevel();
fs::path CurrDebugDir();
fs::path CurrResultDir();
std::string DatasetStr();
fs::path CurrFile();
int MaxTimeElapsed();
bool UseIndividualModelDir();
} // namespace frm::global