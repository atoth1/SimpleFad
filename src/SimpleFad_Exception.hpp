#ifndef SIMPLEFAD_EXCEPTION_HPP
#define SIMPLEFAD_EXCEPTION_HPP

#include "SimpleFad_config.hpp"

#include <exception>
#include <string>
#include <utility>

namespace SimpleFad {

class ExceptionBase: public std::exception {
public:
  ExceptionBase(const std::string& msg)
      : msg_(msg) {
  }

  ExceptionBase(std::string&& msg)
      : msg_(std::move(msg)) {
  }

  const char* what() const noexcept override {
    return msg_.c_str();
  }

private:
  std::string msg_;
};

// Exception indicating an attempt to access a derivative component outside
// the bounds of the FAD type's derivative container, or that two FAD types
// with incompatible sizes are being combined.
class BoundsError: public ExceptionBase {
public:
  BoundsError(const std::string& msg) : ExceptionBase(msg) { }
  BoundsError(std::string&& msg) : ExceptionBase(std::move(msg)) { }
};

// Exception indicating that a mathematical function is being asked to evaulate
// at a point where its value or derivative is not defined.
class DomainError: public ExceptionBase {
public:
  DomainError(const std::string& msg) : ExceptionBase(msg) { }
  DomainError(std::string&& msg) : ExceptionBase(std::move(msg)) { }
};

}

#ifdef ENABLE_EXCEPTIONS

#define CHECK_CONFORMING_SIZE(exprSize1, exprSize2) \
  if ((exprSize1) > 0 && (exprSize2) > 0 && (exprSize1) != (exprSize2)) { \
    const char* msg = "ERROR: Construction/assignment from source expression with " \
      "size which does not conform with destination expression."; \
    throw SimpleFad::BoundsError(msg); \
  }

#define CHECK_WITHIN_BOUNDS(sizeExpr, index) \
  if ((index) < 0 || (index) >= (sizeExpr)) { \
    const char* msg = "ERROR: Accessing index out of bound for derivative size."; \
    throw SimpleFad::BoundsError(msg); \
  }

#define CHECK_WITHIN_DOMAIN(condition) \
    if (!(condition)) { \
      const char* msg = "ERROR: Invalid input argument for function evaluation."; \
      throw SimpleFad::DomainError(msg); \
    }

#else

#define CHECK_CONFORMING_SIZE(sizeDstExpr, sizeSrcExpr)

#define CHECK_WITHIN_BOUNDS(sizeExpr, index)

#define CHECK_WITHIN_DOMAIN(condition)

#endif

#endif
