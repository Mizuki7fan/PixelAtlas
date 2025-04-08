#include "evaluator.h"

// 结果评估模块
// 功能: 输入两次work的名称A,B和工具名称,
// 在run/evaluate_A_B文件夹下生成评估结果
int main(int argc, char *argv[]) {
  Evaluator evaluator(argc, argv);
  evaluator.LoadData();
  evaluator.PrintData();
  evaluator.PrintDataAvgSquaredDifference();
  return 0;
}