#include "process_parallel_executor.h"
#include "global_args.h"
#include <boost/process.hpp>
#include <iostream>

namespace frm {
using GA = GlobalArguments;
namespace bp = boost::process;

ProcessParallelExecutor::ProcessParallelExecutor(
    const std::function<void()> &func,        // 需要并行处理的函数
    const std::vector<fs::path> &run_targets, // 需要并行处理的目标
    const std::size_t num_parallel_cnt,       // 并行数
    const std::string exe_path)
    : m_func(func),                         //
      m_run_targets(run_targets),           //
      m_num_parallel_cnt(num_parallel_cnt), //
      m_exe_path(exe_path) {}

bool ProcessParallelExecutor::Exec() {
  std::size_t running_processes = 0, curr_target_idx = 0;
  std::deque<std::pair<bp::child, std::chrono::steady_clock::time_point>>
      processes;
  // 设置允许进程运行的最大时间
  const auto max_duration = std::chrono::seconds(global::MaxTimeElapsed());

  while (curr_target_idx < m_run_targets.size() || !processes.empty()) {
    // 如果当前运行进程数小于并行数, 且当前例子idx小于总的例子数
    while (running_processes < m_num_parallel_cnt &&
           curr_target_idx < m_run_targets.size()) {
      // 启动一个新的进程
      try {
        const fs::path &target = m_run_targets[curr_target_idx];
        // 传递变量
        std::vector<std::string> args{"-d",                                 //
                                      std::to_string(global::DebugLevel()), //
                                      "--dataset",
                                      global::DatasetName(),
                                      "--single", //
                                      target.filename().string(),
                                      "--work_name",
                                      global::WorkName()};
        if (global::UseIndividualInstanceDir())
          args.emplace_back("--use_individual_instance_dir");

        processes.emplace_back(
            bp::child(bp::exe = m_exe_path,    //
                      bp::args = args),        //
            std::chrono::steady_clock::now()); // 记录启动时间
        ++running_processes;                   // 当前运行的进程数+1
        ++curr_target_idx;                     // 当前处理的例子idx+1
      } catch (const bp::process_error &e) {
        std::cerr << "Fail to launch process: " << e.what() << std::endl;
        return false;
      }
    }
    // 等待所有进程结束
    for (auto it = processes.begin(); it != processes.end();) {
      auto &[child_proc, start_time] = *it;
      // 判断进程是否超时
      if (std::chrono::steady_clock::now() - start_time > max_duration) {
        std::cerr << "Timeout: " << m_run_targets[curr_target_idx - 1]
                  << std::endl;
        child_proc.terminate();
      }

      if (child_proc.running()) {
        ++it;
      } else {
        child_proc.wait();
        it = processes.erase(it);
        --running_processes;
      }
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }
  return true;
}
} // namespace frm
