#include <fstream>
#include <iostream>

int main(int argc, char *argv[]) {
  for (int i = 0; i < argc; i++)
    std::cout << argv[i] << " ";
  std::cout << std::endl;

  // 为了支持不从launcher启动, 单独运行该命令,
  // 解析输入文件依赖以及计算批量运行数据的工作就放在这里
  // 该功能通过一个通用的frm::函数来完成
  return 0;
}