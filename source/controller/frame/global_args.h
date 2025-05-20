#pragma once
#include "action_args.h"
#include "common_program.h"
#include <filesystem>

// 使用单例类来达成全局变量的功能, 通过友元函数控制访问权限
namespace fs = std::filesystem;
// 前置声明全局访问函数
namespace global {
// 外部可见的函数, 前置声明
fs::path WorkDir();
fs::path InstanceFullPath();
// 考虑设置ActionDir()和CurrActionDir()
fs::path ActionDir();
fs::path ActionDebugDir();
fs::path ActionLogDir();
fs::path ActionResultDir();
fs::path ActionResultDir(std::size_t action_idx);
bool UseIndividualInstanceDir();
int DebugLevel();
bool CleanActionCache();
std::string SingleInstanceName();
std::string BatchInstanceRegex();
fs::path DatasetDir();
int NumParallelCnt();
std::size_t CurrActionIdx();
fs::path CurrActionPath();
int MaxTimeElapsed();
std::string DatasetName();
std::string WorkName();
const std::unordered_map<std::string, std::size_t> &ActionInputs();

} // namespace global

namespace frm {

// 全局单例类
class GlobalArguments {
public:
  GlobalArguments();           // 显式声明构造函数
  static GlobalArguments &I(); // 取单例实例
  class Token {
  private:
    Token() = default;
    friend class CommonProgram; // 授权整个类为友元
  };

public:
  // 设置为非静态的方法
  void Initialize(int argc, char **argv, const Token &);
  GlobalArguments(const GlobalArguments &) = delete; // 禁用拷贝
  void operator=(const GlobalArguments &) = delete;  // 禁用赋值

private:
  // 类构造时候就可以确定的值
  std::vector<ActionArguments> all_action_list_;
  // 类的成员变量, 由外部输入确定的值
  fs::path curr_action_path_; // 当前action的路径
  friend fs::path global::CurrActionPath();
  int debug_level_; // 调试等级
  friend int global::DebugLevel();
  std::string batch_instance_regex_;
  friend std::string global::BatchInstanceRegex();
  std::string single_instance_name_;
  friend std::string global::SingleInstanceName();
  std::string dataset_name_;
  friend std::string global::DatasetName();
  int num_parallel_cnt_;
  friend int global::NumParallelCnt();
  bool use_individual_instance_dir_;
  friend bool global::UseIndividualInstanceDir();
  int max_time_elapsed_;
  friend int global::MaxTimeElapsed();
  bool clean_action_cache_;
  friend bool global::CleanActionCache();
  std::string work_name_;
  friend std::string global::WorkName();

  // 根据外部输入的值计算得到的衍生值, 计算衍生值的函数应该标记Impl后缀
  mutable fs::path work_dir_;
  mutable fs::path action_dir_;
  friend fs::path global::ActionDir();
  mutable std::size_t curr_action_idx_ =
      std::numeric_limits<std::size_t>::max();
  friend std::size_t global::CurrActionIdx();
  mutable std::string curr_action_name_;
  mutable fs::path action_result_dir_;
  mutable fs::path action_debug_dir_;
  mutable fs::path action_log_dir_;
  mutable fs::path dataset_dir_;
  friend fs::path global::DatasetDir();
  mutable fs::path instance_full_path_;

  // 记录action的所有输入依赖action_idx
  std::vector<std::unordered_map<std::string, std::size_t>>
      map_input_to_action_;
  friend const std::unordered_map<std::string, std::size_t> &
  global::ActionInputs();

private: // 成员函数读写
  fs::path WorkDirImpl() const;
  friend fs::path global::WorkDir();
  fs::path ActionDirImpl() const;
  std::size_t CurrActionIdxImpl() const;
  std::string CurrActionNameImpl() const;
  fs::path ActionDebugDirImpl() const;
  friend fs::path global::ActionDebugDir();
  fs::path ActionResultDirImpl() const;
  friend fs::path global::ActionResultDir();
  friend fs::path global::ActionResultDir(std::size_t);
  fs::path ActionLogDirImpl() const;
  friend fs::path global::ActionLogDir();
  fs::path DatasetDirImpl() const;
  fs::path InstanceFullPathImpl() const;
  friend fs::path global::InstanceFullPath();
};

} // namespace frm
