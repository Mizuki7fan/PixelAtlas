#include "process_parallel_executor.h"
#include <boost/process.hpp>

namespace frm {
namespace bp = boost::process;

ProcessParallelExecutor::ProcessParallelExecutor(
    const std::function<void(std::string)> &func, // 需要并行处理的函数
    const std::vector<std::string> &run_targets,  // 需要并行处理的目标
    const std::size_t num_parallel_cnt,           // 并行数
    const std::string exe_path)
    : m_func(func),                         //
      m_run_targets(run_targets),           //
      m_num_parallel_cnt(num_parallel_cnt), //
      m_exe_path(exe_path) {}

bool ProcessParallelExecutor::Exec() {
  //
  return true;
}

} // namespace frm
