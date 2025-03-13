#pragma once
#include <functional>
#include <string>

namespace frm {
class ProcessParallelExecutor {
public:
  explicit ProcessParallelExecutor(
      const std::function<void(std::string)> &func, // 需要并行处理的函数
      const std::vector<std::string> &run_targets,  // 需要并行处理的目标
      const std::size_t num_parallel_cnt,           // 并行数
      const std::string exe_path);                  // exe路径
  bool Exec();

private:
  const std::function<void(std::string)> m_func;
  const std::vector<std::string> m_run_targets;
  const std::size_t m_num_parallel_cnt;
  const std::string m_exe_path;
};
} // namespace frm
