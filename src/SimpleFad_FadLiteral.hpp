#ifndef SIMPLEFAD_FADLITERAL_HPP
#define SIMPLEFAD_FADLITERAL_HPP

#include "SimpleFad_ExprBase.hpp"

namespace SimpleFad
{

// Expression wrapper for a scalar literal of type T.
// Simply wraps a constant value, and returns 0 for
// any requested derivative index.
template <class T>
class FadLiteral : public ExprBase<FadLiteral<T>>
{
public:
  using traits_type = ExprTraits<FadLiteral<T>>;
  using value_type = typename traits_type::value_type;
  using index_type = typename traits_type::index_type;

  constexpr explicit FadLiteral(const T& v) noexcept
      : val_(v) {
  }

  // Value accessor
  constexpr value_type valImpl() const noexcept { return val_; }

  // Derivative accessor
  constexpr value_type dvalImpl(index_type) const noexcept {
    return 0.;
  }

  // Derivative array size
  constexpr index_type dsizeImpl() const noexcept { return 0; }

private:
  value_type val_;
};

}

#endif
