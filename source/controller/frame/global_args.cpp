#include "global_args.h"
#include "assert.hpp"
#include <boost/program_options.hpp>
#include <iostream>

namespace po = boost::program_options;
namespace frm {
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

  all_action_list_ = LoadAllActionList();

  std::unordered_map<std::string, std::size_t> map_output_name_to_action_idx;

  for (std::size_t act_idx = 0; act_idx < CurrActionIdx(); ++act_idx)
    for (const std::string &file_name : all_action_list_[act_idx].outputs)
      map_output_name_to_action_idx[file_name] = act_idx;

  for (const std::string &file_name :
       all_action_list_[CurrActionIdx()].inputs) {
    PA_ASSERT_WITH_MSG(map_output_name_to_action_idx.count(file_name) != 0,
                       std::format("输入文件{}在输出文件中未找到", file_name));
    std::size_t act_idx = map_output_name_to_action_idx[file_name];
    map_input_name_to_full_path_[file_name] =
        WorkDir() /
        std::format("{}_{}", act_idx, all_action_list_[act_idx].name) /
        "result" /
        std::format("{}-{}", InstanceFullPath().filename().string(), file_name);
  }
}

fs::path GlobalArguments::WorkDir() const {
  PA_ASSERT_WITH_MSG(!work_name_.empty(), "work_name为空");
  if (work_dir_.empty()) { // 如果值为空则设置
    work_dir_ = fs::current_path() / std::format("../work_{}", work_name_);
  }
  return work_dir_;
}

fs::path GlobalArguments::ActionDir() const {
  if (action_dir_.empty()) {
    action_dir_ =
        WorkDir() / std::format("{}_{}", CurrActionIdx(), CurrActionName());
  }
  return action_dir_;
}

std::string GlobalArguments::CurrActionName() const {
  if (curr_action_name_.empty())
    curr_action_name_ = curr_action_path_.stem().string();
  return curr_action_name_;
}

std::size_t GlobalArguments::CurrActionIdx() const {
  if (curr_action_idx_ == std::numeric_limits<std::size_t>::max()) {
    for (std::size_t i = 0; i < all_action_list_.size(); ++i) {
      if (all_action_list_[i].name == CurrActionName()) {
        curr_action_idx_ = i;
        break;
      }
    }
  }
  PA_ASSERT_WITH_MSG(curr_action_idx_ !=
                         std::numeric_limits<std::size_t>::max(),
                     std::format("Action Name {}无效", CurrActionName()));
  return curr_action_idx_;
}

fs::path GlobalArguments::ActionDebugDir() const {
  if (action_debug_dir_.empty()) {
    if (UseIndividualInstanceDir()) {
      action_debug_dir_ = ActionDir() / "debug" / single_instance_name_;
    } else {
      action_debug_dir_ = ActionDir() / "debug";
    }
  }
  return action_debug_dir_;
}

fs::path GlobalArguments::ActionResultDir() const {
  if (action_result_dir_.empty()) {
    if (UseIndividualInstanceDir()) {
      action_result_dir_ = ActionDir() / "result" / single_instance_name_;
    } else {
      action_result_dir_ = ActionDir() / "result";
    }
  }
  return action_result_dir_;
}

fs::path GlobalArguments::ActionLogDir() const {
  if (action_log_dir_.empty()) {
    if (UseIndividualInstanceDir()) {
      action_log_dir_ = ActionDir() / "log" / single_instance_name_;
    } else {
      action_log_dir_ = ActionDir() / "log";
    }
  }
  return action_log_dir_;
}

fs::path GlobalArguments::DatasetDir() const {
  if (dataset_dir_.empty()) {
    dataset_dir_ = fs::current_path() / ".." / ".." / "asset" / DataSetName();
  }
  return dataset_dir_;
}

fs::path GlobalArguments::InstanceFullPath() const {
  if (instance_full_path_.empty()) {
    instance_full_path_ = DatasetDir() / SingleInstanceName();
  }
  return instance_full_path_;
}

} // namespace frm
