#ifndef SIMPLEFAD_EXPRBASE_HPP
#define SIMPLEFAD_EXPRBASE_HPP

#include "SimpleFad_config.hpp"
#include "SimpleFad_ExprTraits.hpp"

#include <type_traits>

namespace SimpleFad {

template<class Expr>
class ExprBase {
public:
  using traits_type = ExprTraits<Expr>;
  using value_type = typename traits_type::value_type;
  using index_type = typename traits_type::index_type;

  static_assert(std::is_floating_point_v<value_type>);

  SIMPLEFAD_CONSTEXPR value_type val() const {
    return self().valImpl();
  }

  SIMPLEFAD_CONSTEXPR value_type dval(index_type id) const {
    return self().dvalImpl(id);
  }

  constexpr index_type dsize() const {
    return self().dsizeImpl();
  }

  constexpr const Expr& self() const {
    return static_cast<const Expr&>(*this);
  }
};

}

#endif
