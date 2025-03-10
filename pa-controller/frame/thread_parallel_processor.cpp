#include "thread_parallel_processor.h"
#include <iostream>

namespace frm {
ThreadParallelProcessor::ThreadParallelProcessor(
    const std::function<void(std::string)> &func, // 需要并行处理的函数
    const std::vector<std::string> &run_targets,  // 需要并行处理的目标
    const std::size_t num_parallel_cnt)
    : m_func(func), m_run_targets(run_targets),
      m_num_parallel_cnt(num_parallel_cnt) {
  const std::size_t num_targets = run_targets.size();
  std::cout << std::format("共{}个例子, 通过{}个线程并行执行:", num_targets,
                           m_num_parallel_cnt)
            << std::endl;
  task_queue =
      std::deque<std::string>(m_run_targets.begin(), m_run_targets.end());
}

void ThreadParallelProcessor::Exec() {

  std::atomic<bool> stop_flag{false};
  // 结果统计
  std::atomic<size_t> success_count{0}, processed_count{0};
  // 进度监控
  boost::asio::io_context io;
  boost::asio::deadline_timer timer(io);
}

void ThreadParallelProcessor::ThreadWorker() {
  while (true) {
    std::string model_name;
    {
      boost::unique_lock<boost::mutex> lock(queue_mutex);
      queue_cv.wait(lock, [&] { return !task_queue.empty() || stop_flag; });
    }
  }
}
} // namespace frm
