#pragma once

#define PA_ASSERT(expr)                                                        \
  do {                                                                         \
    if (!(expr)) {                                                             \
      std::cerr << "Assertion failed: (" #expr "), function " << __FUNCTION__  \
                << ", file " << __FILE__ << ", line " << __LINE__ << "."       \
                << std::endl;                                                  \
      std::abort();                                                            \
    }                                                                          \
  } while (false)

#define PA_ASSERT_WITH_MSG(expr, msg)                                          \
  do {                                                                         \
    if (!(expr)) {                                                             \
      std::cerr << "Assertion failed: (" #expr "), function " << __FUNCTION__  \
                << ", file " << __FILE__ << ", line " << __LINE__ << "."       \
                << std::endl;                                                  \
      std::cerr << "Message: " << msg << std::endl;                            \
      std::abort();                                                            \
    }                                                                          \
  } while (false)