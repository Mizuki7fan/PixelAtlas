#pragma once
#include <boost/thread.hpp>
#include <deque>
#include <filesystem>
#include <functional>

// 用来运行线程级并行
namespace fs = std::filesystem;
namespace frm {
class ThreadParallelExecutor {
public:
  explicit ThreadParallelExecutor(
      const std::function<void(fs::path)> &func, // 需要并行处理的函数
      const std::vector<fs::path> &run_targets,  // 需要并行处理的目标
      const std::size_t num_parallel_cnt);
  bool Exec();

private:
  void ThreadWorker();

private:
  const std::function<void(fs::path)> m_func;
  const std::vector<fs::path> m_run_targets;
  const std::size_t m_num_targets;
  const std::size_t m_num_parallel_cnt;

private:
  boost::mutex queue_mutex;
  boost::condition_variable queue_cv;
  std::deque<fs::path> task_queue;
  std::atomic<bool> stop_flag{false};
  // 结果统计
  std::atomic<size_t> success_count{0}, processed_count{0};
};
} // namespace frm
