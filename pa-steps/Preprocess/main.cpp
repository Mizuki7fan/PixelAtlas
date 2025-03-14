#include <chrono>
#include <frame/common_program.h>
#include <frame/global_defs.h>
#include <frame/parse.h>
#include <iostream>
#include <random>
#include <thread>

// 函数执行入口
void MainProcess(std::string model_name) {
  // 生成 0-10 秒的随机等待
  std::random_device rd;                       // 随机种子
  std::mt19937 gen(rd());                      // 使用 Mersenne Twister 算法
  std::uniform_int_distribution<> dist(0, 10); // [0,2] 整数分布

  int wait_seconds = dist(gen); // 生成随机秒数
  std::this_thread::sleep_for(std::chrono::seconds(wait_seconds));
  if (wait_seconds > 5) {
    std::cout << "wait_seconds: " << wait_seconds << std::endl;
    exit(0);
  }

  int A = 1;
  if (A < 0)
    std::cout << model_name;
  // std::cout << model_name << std::endl;
}

int main(int argc, char *argv[]) {
  std::cout << argv[0] << " " << argc << std::endl;
  frm::CommonProgram common_program(argc, argv);

  common_program.Run(MainProcess);

  return 0;
}