#pragma once
#include <filesystem>
#include <functional>
#include <string>

namespace fs = std::filesystem;
namespace frm {
class ProcessParallelExecutor {
public:
  explicit ProcessParallelExecutor(
      const std::function<void()> &func,        // 需要并行处理的函数
      const std::vector<fs::path> &run_targets, // 需要并行处理的目标
      const std::size_t num_parallel_cnt,       // 并行数
      const std::string exe_path);              // exe路径
  bool Exec();

private:
  const std::function<void()> m_func;
  const std::vector<fs::path> m_run_targets;
  const std::size_t m_num_parallel_cnt;
  const std::string m_exe_path;
};
} // namespace frm
