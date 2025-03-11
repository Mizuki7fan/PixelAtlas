#include "thread_parallel_executor.h"
#include <boost/asio.hpp>
#include <iostream>

namespace frm {
ThreadParallelExecutor::ThreadParallelExecutor(
    const std::function<void(std::string)> &func, // 需要并行处理的函数
    const std::vector<std::string> &run_targets,  // 需要并行处理的目标
    const std::size_t num_parallel_cnt)
    : m_func(func),                        //
      m_run_targets(run_targets),          //
      m_num_targets(m_run_targets.size()), //
      m_num_parallel_cnt(num_parallel_cnt) //
{
  std::cout << std::format("共{}个例子, 通过{}个线程并行执行:", m_num_targets,
                           m_num_parallel_cnt)
            << std::endl;
  task_queue =
      std::deque<std::string>(m_run_targets.begin(), m_run_targets.end());
}

bool ThreadParallelExecutor::Exec() {
  // 进度监控
  boost::asio::io_context io;
  boost::asio::deadline_timer timer(io);

  // 使用std::function包装回调以实现递归调用
  std::function<void(const boost::system::error_code &)> progress_monitor;

  // 统计执行状态
  progress_monitor = [&](const boost::system::error_code &ec) {
    if (ec)
      return;
    const size_t done = processed_count.load();
    std::cout << "\r进度: " << done << "/" << m_num_targets << " ("
              << (done * 100 / m_num_targets) << "%)" << std::flush;
    timer.expires_from_now(boost::posix_time::milliseconds(500));
    timer.async_wait(progress_monitor);
  };

  timer.expires_from_now(boost::posix_time::milliseconds(500));
  timer.async_wait(progress_monitor);
  boost::thread io_thread([&io]() { io.run(); });

  boost::thread_group workers;
  for (std::size_t i = 0; i < m_num_targets; ++i) {
    // 也可以写成匿名函数
    workers.create_thread(
        boost::bind(&ThreadParallelExecutor::ThreadWorker, this));
  }

  // 清理资源
  stop_flag.store(true);
  queue_cv.notify_all();
  workers.join_all();

  io.stop();
  io_thread.join();

  // 输出最终结果
  std::cout << "\n完成: " << success_count << "/" << m_num_targets << " 成功"
            << std::endl;

  return success_count == m_num_targets ? 0 : 1;
}

void ThreadParallelExecutor::ThreadWorker() {
  while (true) {
    std::string model_name;
    {
      boost::unique_lock<boost::mutex> lock(queue_mutex);
      queue_cv.wait(lock, [&] { return !task_queue.empty() || stop_flag; });
      if (stop_flag && task_queue.empty())
        break;
      if (task_queue.empty())
        continue;
      model_name = task_queue.front(); // 取出队列头
      task_queue.pop_front();          // 移除队列头
    }

    try {
      std::cout << "model_name: " << model_name << std::endl;
      m_func(model_name);
      // std::memory_order_relaxed: 无序, 有锁
      success_count.fetch_add(1, std::memory_order_relaxed);
    } catch (const std::exception &e) {
      std::cerr << "\n发现错误: " << model_name << ": " << e.what()
                << std::endl;
    }
    processed_count.fetch_add(1, std::memory_order_relaxed);
  }
}
} // namespace frm
