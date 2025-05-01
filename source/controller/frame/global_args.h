#pragma once
#include "action_args.h"
#include "common_program.h"
#include <filesystem>

// 使用单例类来达成全局变量的功能,
// 并通过友元类保证该单例类的初始化仅限于common_program类中

namespace fs = std::filesystem;

namespace frm {
class GlobalArguments {
public:
  GlobalArguments() {}; // 显式声明构造函数
  // 对外开放函数
  static GlobalArguments &I(); // 取单例实例
  fs::path CurrActionPath() const { return curr_action_path_; }
  bool UseIndividualInstanceDir() const { return use_individual_instance_dir_; }
  int DebugLevel() const { return debug_level_; }
  int NumParallelCnt() const { return num_parallel_cnt_; }
  int MaxTimeElapsed() const { return max_time_elapsed_; }
  std::string DataSetName() const { return dataset_name_; }
  std::string WorkName() const { return work_name_; }
  bool CleanActionCache() const { return clean_action_cache_; }
  std::string SingleInstanceName() const { return single_instance_name_; }
  std::string BatchInstanceRegex() const { return batch_instance_regex_; }
  // 推算的值
  fs::path WorkDir() const;
  fs::path ActionDir() const;
  fs::path ActionResultDir() const;
  fs::path ActionDebugDir() const;
  fs::path ActionLogDir() const;
  fs::path DatasetDir() const;
  fs::path InstanceFullPath() const;
  std::size_t CurrActionIdx() const;

  const auto &MapInputNameToFullPath() const {
    return map_input_name_to_full_path_;
  }

private:
  std::string CurrActionName() const;

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
  // 原生值
  fs::path curr_action_path_;
  int debug_level_;
  std::string batch_instance_regex_;
  std::string single_instance_name_;
  std::string dataset_name_;
  int num_parallel_cnt_;
  bool use_individual_instance_dir_;
  int max_time_elapsed_;
  bool clean_action_cache_;
  std::string work_name_;
  std::vector<ActionArguments> all_action_list_;
  // 衍生值
  mutable fs::path work_dir_;
  mutable fs::path action_dir_;
  mutable std::size_t curr_action_idx_ =
      std::numeric_limits<std::size_t>::max();
  mutable std::string curr_action_name_;
  mutable fs::path action_result_dir_;
  mutable fs::path action_debug_dir_;
  mutable fs::path action_log_dir_;
  mutable fs::path dataset_dir_;
  mutable fs::path instance_full_path_;

  // 标记当前步骤的输入所依赖的前置步骤的完整路径
  std::unordered_map<std::string, fs::path> map_input_name_to_full_path_;
};
} // namespace frm