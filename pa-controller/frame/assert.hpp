#pragma once
#ifdef _WIN32
#ifndef NOMINMAX // 仅在未定义时添加
#define NOMINMAX // 阻止 Windows 的 min/max 宏
#endif
#include <windows.h>
#else
#include <iostream>
#endif

#include <string>

// 辅助函数处理消息显示
inline void AssertOutput(const std::string &message) {
#ifdef _WIN32
  MessageBoxA(nullptr, message.c_str(), "Assertion Failed",
              MB_ICONERROR | MB_OK);
#else
  std::cerr << message << std::endl;
#endif
}

// 主断言宏
#define PA_ASSERT(expr)                                                        \
  do {                                                                         \
    if (!(expr)) {                                                             \
      std::ostringstream oss;                                                  \
      oss << "Assertion failed: " #expr "\n"                                   \
          << "Function: " << __FUNCTION__ << "\n"                              \
          << "File: " << __FILE__ << "\n"                                      \
          << "Line: " << __LINE__;                                             \
      AssertOutput(oss.str());                                                 \
      std::abort();                                                            \
    }                                                                          \
  } while (0)

// 带自定义消息的断言宏
#define PA_ASSERT_WITH_MSG(expr, msg)                                          \
  do {                                                                         \
    if (!(expr)) {                                                             \
      std::ostringstream oss;                                                  \
      oss << "Assertion failed: " #expr "\n"                                   \
          << "Function: " << __FUNCTION__ << "\n"                              \
          << "File: " << __FILE__ << "\n"                                      \
          << "Line: " << __LINE__ << "\n"                                      \
          << "Message: " << (msg);                                             \
      AssertOutput(oss.str());                                                 \
      std::abort();                                                            \
    }                                                                          \
  } while (0)