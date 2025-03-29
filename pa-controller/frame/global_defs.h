#pragma once
#include <filesystem>
#include <string>
namespace fs = std::filesystem;

namespace frm::global {

int DebugLevel();
fs::path CurrDebugDir();
fs::path CurrResultDir();
std::string DatasetStr();
fs::path CurrFile();
int MaxTimeElapsed();
bool UseIndividualModelDir();
std::string RunName();
} // namespace frm::global