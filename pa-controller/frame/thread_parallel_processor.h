#pragma once
#include <boost/asio.hpp>
#include <boost/thread.hpp>
#include <deque>
#include <functional>
#include <string>

// 用来运行线程级并行

namespace frm {
class ThreadParallelProcessor {
public:
  explicit ThreadParallelProcessor(
      const std::function<void(std::string)> &func, // 需要并行处理的函数
      const std::vector<std::string> &run_targets,  // 需要并行处理的目标
      const std::size_t num_parallel_cnt);
  void Exec();

private:
  void ThreadWorker();

private:
  const std::function<void(std::string)> m_func;
  const std::vector<std::string> m_run_targets;
  const std::size_t m_num_parallel_cnt;

private:
  boost::mutex queue_mutex;
  boost::condition_variable queue_cv;
  std::deque<std::string> task_queue;
  std::atomic<bool> stop_flag{false};
  // 结果统计
  std::atomic<size_t> success_count{0}, processed_count{0};
};
} // namespace frm
