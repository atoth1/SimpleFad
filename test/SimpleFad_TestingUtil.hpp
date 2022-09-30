#ifndef SIMPLEFAD_TESTINGUTIL_HPP
#define SIMPLEFAD_TESTINGUTIL_HPP

namespace SimpleFad {
  template<class T>
  constexpr bool checkFloatingEquality(const T a, const T b, const T eps) {
    return a > b ? a - b < eps : b - a < eps;
  }
}

#endif
