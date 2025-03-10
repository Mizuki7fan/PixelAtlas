#pragma once
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN // 阻止windows.h包含WinSock.h
#include <windows.h>
#else
#include <iostream>
#endif

#include <cstdlib>
#include <sstream>
#include <string>

// 核心断言处理函数（内联）
inline void AssertInternal(bool condition, const std::string &expr,
                           const std::string &message, const char *file,
                           int line, const char *function) {
  if (!condition) {
    std::ostringstream oss;
    oss << "Assertion failed: " << expr << "\n"
        << "Function: " << function << "\n"
        << "File: " << file << "\n"
        << "Line: " << line << "\n";
    if (!message.empty()) {
      oss << "Message: " << message << "\n";
    }

#ifdef _WIN32
    MessageBoxA(nullptr, oss.str().c_str(), "Assertion Failed",
                MB_ICONERROR | MB_OK);
#else
    std::cerr << oss.str() << std::endl;
#endif
    std::abort();
  }
}

// 包装宏（用于捕获 __FILE__、__LINE__ 等编译时信息）
#define PA_ASSERT(expr)                                                        \
  AssertInternal(static_cast<bool>(expr), #expr, "", __FILE__, __LINE__,       \
                 __FUNCTION__)

#define PA_ASSERT_WITH_MSG(expr, msg)                                          \
  AssertInternal(static_cast<bool>(expr), #expr, (msg), __FILE__, __LINE__,    \
                 __FUNCTION__)