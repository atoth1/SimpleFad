#ifndef SIMPLEFAD_COMPARISONOPS_HPP
#define SIMPLEFAD_COMPARISONOPS_HPP

#include "SimpleFad_ExprBase.hpp"

#include <type_traits>

namespace SimpleFad {

#define GENERATE_COMPARISON_OVERLOADS(OPNAME, OP) \
  template<class Expr1, class Expr2> \
  constexpr bool OPNAME(const ExprBase<Expr1>& e1, const ExprBase<Expr2>& e2) { \
    return e1.val() OP e2.val(); \
  } \
  template<class Expr> \
  constexpr bool \
  OPNAME(const ExprBase<Expr>& e, const typename Expr::value_type& t) { \
    return e.val() OP t; \
  } \
  template<class Expr> \
  constexpr bool \
  OPNAME(const typename Expr::value_type& t, const ExprBase<Expr>& e) { \
    return t OP e.val(); \
  }

GENERATE_COMPARISON_OVERLOADS(operator==, ==)
GENERATE_COMPARISON_OVERLOADS(operator!=, !=)
GENERATE_COMPARISON_OVERLOADS(operator<, <)
GENERATE_COMPARISON_OVERLOADS(operator>, >)
GENERATE_COMPARISON_OVERLOADS(operator<=, <=)
GENERATE_COMPARISON_OVERLOADS(operator>=, >=)

#undef GENERATE_COMPARISON_OVERLOADS

}

#endif
