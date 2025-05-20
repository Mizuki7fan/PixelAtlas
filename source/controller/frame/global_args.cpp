#include "global_args.h"
#include "pa_assert.hpp"
#include <boost/program_options.hpp>
#include <iostream>

namespace po = boost::program_options;

namespace global {
fs::path WorkDir() { return frm::GlobalArguments::I().WorkDirImpl(); }
fs::path InstanceFullPath() {
  return frm::GlobalArguments::I().InstanceFullPathImpl();
}

fs::path ActionDir() { return frm::GlobalArguments::I().ActionDirImpl(); }

fs::path ActionDebugDir() {
  return frm::GlobalArguments::I().ActionDebugDirImpl();
}

fs::path ActionLogDir() { return frm::GlobalArguments::I().ActionLogDirImpl(); }

fs::path ActionResultDir() {
  return frm::GlobalArguments::I().ActionResultDirImpl();
}

fs::path ActionResultDir(std::size_t idx) {
  std::string action_name =
      frm::GlobalArguments::I().all_action_list_[idx].name;
  fs::path work_dir = frm::GlobalArguments::I().WorkDirImpl();
  return work_dir / std::format("{}_{}", idx, action_name) / "result";
}

bool UseIndividualInstanceDir() {
  return frm::GlobalArguments::I().use_individual_instance_dir_;
}
int DebugLevel() { return frm::GlobalArguments::I().debug_level_; }

bool CleanActionCache() {
  return frm::GlobalArguments::I().clean_action_cache_;
}

std::string SingleInstanceName() {
  return frm::GlobalArguments::I().single_instance_name_;
}

std::string BatchInstanceRegex() {
  return frm::GlobalArguments::I().batch_instance_regex_;
}

fs::path DatasetDir() { return frm::GlobalArguments::I().DatasetDirImpl(); }

int NumParallelCnt() { return frm::GlobalArguments::I().num_parallel_cnt_; }

std::size_t CurrActionIdx() {
  return frm::GlobalArguments::I().CurrActionIdxImpl();
}

fs::path CurrActionPath() {
  return frm::GlobalArguments::I().curr_action_path_;
}

int MaxTimeElapsed() { return frm::GlobalArguments::I().max_time_elapsed_; }
std::string DatasetName() { return frm::GlobalArguments::I().dataset_name_; }
std::string WorkName() { return frm::GlobalArguments::I().work_name_; }

const std::unordered_map<std::string, std::size_t> &ActionInputs() {
  std::size_t curr_action_idx = frm::GlobalArguments::I().CurrActionIdxImpl();
  return frm::GlobalArguments::I().map_input_to_action_[curr_action_idx];
}

const std::unordered_set<std::string> &ActionOutputs() {
  std::size_t curr_action_idx = frm::GlobalArguments::I().CurrActionIdxImpl();
  return frm::GlobalArguments::I().all_action_list_[curr_action_idx].outputs;
}

const std::unordered_map<std::string, frm::ValueType> &ActionHyperParameters() {
  std::size_t curr_action_idx = frm::GlobalArguments::I().CurrActionIdxImpl();
  return frm::GlobalArguments::I()
      .all_action_list_[curr_action_idx]
      .hpyer_parameters;
}

} // namespace global

namespace frm {
GlobalArguments::GlobalArguments() {
  // 读入all_action_list
  all_action_list_ = LoadAllActionList();
  map_input_to_action_.resize(all_action_list_.size());

  for (std::size_t action_idx = 0; action_idx < all_action_list_.size();
       ++action_idx) {
    for (std::string input_name : all_action_list_[action_idx].inputs) {
      for (std::size_t i = 0; i < action_idx; ++i) {
        for (std::string output_name : all_action_list_[i].outputs) {
          if (input_name == output_name) {
            map_input_to_action_[action_idx][input_name] = i;
          }
        }
      }
    }
  }
}

GlobalArguments &GlobalArguments::I() {
  static GlobalArguments instance;
  return instance;
}

void GlobalArguments::Initialize(int argc, char **argv, const Token &) {
  curr_action_path_ = argv[0];

  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "帮助信息") //
      ("debug,d", po::value<int>(&debug_level_)->default_value(0),
       "调试等级") //
      ("batch",
       po::value<std::string>(&batch_instance_regex_)->default_value(""),
       "批量执行") //
      ("single",
       po::value<std::string>(&single_instance_name_)->default_value(""),
       "单一执行") //
      ("dataset", po::value<std::string>(&dataset_name_)->default_value(""),
       "数据集") //
      ("parallel,p", po::value<int>(&num_parallel_cnt_)->default_value(1),
       "并行数") //
      ("use_individual_instance_dir",
       po::bool_switch(&use_individual_instance_dir_)->default_value(false),
       "是否每个instance独立建立文件夹") //
      ("max_time_elapsed,t",
       po::value<int>(&max_time_elapsed_)->default_value(1800),
       "最大耗时") //
      ("clean", po::bool_switch(&clean_action_cache_)->default_value(false),
       "是否清理单步缓存") //
      ("work_name", po::value<std::string>(&work_name_)->default_value(""),
       "设置运行批次名");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
      std::cout << desc << std::endl;
      exit(0);
    }
    po::notify(vm);
  } catch (const po::error &e) {
    std::cerr << "error: " << e.what() << std::endl;
    std::cerr << desc << std::endl;
    PA_ASSERT_WITH_MSG(0, "输入异常");
  }

  PA_FALSE_ASSERT_WITH_MSG(batch_instance_regex_.empty() &&
                               single_instance_name_.empty(),
                           "single和batch都未指定需要处理的文件");
  PA_FALSE_ASSERT_WITH_MSG(!batch_instance_regex_.empty() &&
                               !single_instance_name_.empty(),
                           "single和batch都指定了需要处理的文件");
}

fs::path GlobalArguments::WorkDirImpl() const {
  PA_ASSERT_WITH_MSG(!work_name_.empty(), "work_name为空");
  if (work_dir_.empty()) // 如果值为空则设置
    work_dir_ = fs::current_path() / std::format("../work_{}", work_name_);
  return work_dir_;
}

fs::path GlobalArguments::ActionDirImpl() const {
  if (action_dir_.empty())
    action_dir_ = WorkDirImpl() / std::format("{}_{}", CurrActionIdxImpl(),
                                              CurrActionNameImpl());
  return action_dir_;
}

std::size_t GlobalArguments::CurrActionIdxImpl() const {
  if (curr_action_idx_ == std::numeric_limits<std::size_t>::max()) {
    for (std::size_t i = 0; i < all_action_list_.size(); ++i) {
      if (all_action_list_[i].name == CurrActionNameImpl()) {
        curr_action_idx_ = i;
        break;
      }
    }
  }
  return curr_action_idx_;
}

std::string GlobalArguments::CurrActionNameImpl() const {
  if (curr_action_name_.empty())
    curr_action_name_ = curr_action_path_.stem().string();
  return curr_action_name_;
}

fs::path GlobalArguments::ActionDebugDirImpl() const {
  if (action_debug_dir_.empty()) {
    if (use_individual_instance_dir_)
      action_debug_dir_ = ActionDirImpl() / "debug" / single_instance_name_;
    else
      action_debug_dir_ = ActionDirImpl() / "debug";
  }
  return action_debug_dir_;
}

fs::path GlobalArguments::ActionResultDirImpl() const {
  if (action_result_dir_.empty()) {
    if (use_individual_instance_dir_)
      action_result_dir_ = ActionDirImpl() / "result" / single_instance_name_;
    else
      action_result_dir_ = ActionDirImpl() / "result";
  }
  return action_result_dir_;
}

fs::path GlobalArguments::ActionLogDirImpl() const {
  if (action_log_dir_.empty()) {
    if (use_individual_instance_dir_)
      action_log_dir_ = ActionDirImpl() / "log" / single_instance_name_;
    else
      action_log_dir_ = ActionDirImpl() / "log";
  }
  return action_log_dir_;
}

fs::path GlobalArguments::DatasetDirImpl() const {
  if (dataset_dir_.empty())
    dataset_dir_ = fs::current_path() / ".." / ".." / "asset" / dataset_name_;
  return dataset_dir_;
}

fs::path GlobalArguments::InstanceFullPathImpl() const {
  if (instance_full_path_.empty())
    instance_full_path_ = DatasetDirImpl() / single_instance_name_;
  return instance_full_path_;
}

} // namespace frm
