#ifndef SIMPLEFAD_BIARYOPS_HPP
#define SIMPLEFAD_BIARYOPS_HPP

#include "SimpleFad_config.hpp"
#include "SimpleFad_Exception.hpp"
#include "SimpleFad_ExprBase.hpp"
#include "SimpleFad_BinaryExpr.hpp"
#include "SimpleFad_FadLiteral.hpp"

#include <cmath>
#include <type_traits>

namespace SimpleFad {

template<class FOp, class DFXOp, class DFYOp>
class BinaryOp {
public:
  constexpr explicit BinaryOp(const FOp& f, const DFXOp& dfx, const DFYOp& dfy)
      : f_(f), dfx_(dfx), dfy_(dfy) {
  }

  template<class T>
  SIMPLEFAD_CONSTEXPR T f(const T& x, const T& y) const {
    return f_(x, y);
  }

  template<class T>
  SIMPLEFAD_CONSTEXPR T dfx(const T& x, const T& y) const {
    return dfx_(x, y);
  }

  template<class T>
  SIMPLEFAD_CONSTEXPR T dfy(const T& x, const T& y) const {
    return dfy_(x, y);
  }

private:
  FOp f_;
  DFXOp dfx_;
  DFYOp dfy_;
};

template<class FOp, class DFXOp, class DFYOp>
constexpr auto makeBinaryOp(const FOp& f, const DFXOp& dfx, const DFYOp& dfy) {
  return BinaryOp<FOp, DFXOp, DFYOp>(f, dfx, dfy);
}

#define GENERATE_BINARY_OVERLOADS(OPNAME, OP) \
  template<class Expr1, class Expr2> \
  constexpr auto \
  OPNAME (const ExprBase<Expr1>& e1, const ExprBase<Expr2>& e2) { \
     return BinaryExpr<decltype(OP), Expr1, Expr2>(OP, e1.self(), e2.self()); \
  } \
  template<class Expr> \
  constexpr auto \
  OPNAME (const ExprBase<Expr>& e, const typename Expr::value_type& t) { \
    using v_t = typename Expr::value_type; \
    auto l = FadLiteral<v_t>(t); \
    return BinaryExpr<decltype(OP), Expr, FadLiteral<v_t>>(OP, e.self(), l); \
  } \
  template<class Expr> \
  constexpr auto \
  OPNAME (const typename Expr::value_type& t, const ExprBase<Expr>& e) { \
    using v_t = typename Expr::value_type; \
    auto l = FadLiteral<v_t>(t); \
    return BinaryExpr<decltype(OP), FadLiteral<v_t>, Expr>(OP, l, e.self()); \
  }

inline constexpr auto ADD_OP = makeBinaryOp(
    [](const auto& x, const auto& y) {return x+y;},
    [](const auto& x, const auto& y) {return static_cast<decltype(x)>(1);},
    [](const auto& x, const auto& y) {return static_cast<decltype(x)>(1);});
GENERATE_BINARY_OVERLOADS(operator+, ADD_OP)

inline constexpr auto SUBTRACT_OP = makeBinaryOp(
    [](const auto& x, const auto& y) {return x-y;},
    [](const auto& x, const auto& y) {return static_cast<decltype(x)>(1);},
    [](const auto& x, const auto& y) {return static_cast<decltype(x)>(-1);});
GENERATE_BINARY_OVERLOADS(operator-, SUBTRACT_OP)

inline constexpr auto MULTIPLY_OP = makeBinaryOp(
    [](const auto& x, const auto& y) {return x*y;},
    [](const auto& x, const auto& y) {return y;},
    [](const auto& x, const auto& y) {return x;});
GENERATE_BINARY_OVERLOADS(operator*, MULTIPLY_OP)

inline constexpr auto DIVIDE_OP = makeBinaryOp(
    [](const auto& x, const auto& y) {
      CHECK_WITHIN_DOMAIN(y != 0)
      return x/y;
    },
    [](const auto& x, const auto& y) {
      return static_cast<decltype(y)>(1)/y;
    },
    [](const auto& x, const auto& y) {
      return -x/(y*y);}
    );
GENERATE_BINARY_OVERLOADS(operator/, DIVIDE_OP)

inline constexpr auto POW_OP = makeBinaryOp(
    [](const auto& x, const auto& y) {return std::pow(x, y);},
    [](const auto& x, const auto& y) {return y * std::pow(x, y-1);},
    [](const auto& x, const auto& y) {return std::pow(x, y) * std::log(x);});
GENERATE_BINARY_OVERLOADS(pow, POW_OP)

/* NOTE: If std::max/min are visible, for unqualified max/min calls the overload:
   template<class T> std::max/min(const T&, const T&);
   will be preferred to the versions generated below in cases where the arguments
   have identical types (i.e. max(const DynamicFadVariable<T>&, const DynamicFadVaraible<T>&)).
   ADL sees namespace std since vector and array are the default deriv storge template types.
   This is fine as long as "SimpleFad_ComparisonOps.hpp" is included. There is a symantic
   difference though, as the std versions evaluate eagerly and return a const reference to
   the type of the inputs, and the following overloads return a lazily evaluated BinaryExpr.*/
inline constexpr auto MAX_OP = makeBinaryOp(
    [](const auto& x, const auto& y) {return x > y ? x : y;},
    [](const auto& x, const auto& y) {
      using value_type = decltype(x);
      return x > y ? static_cast<value_type>(1) : static_cast<value_type>(0);
    },
    [](const auto& x, const auto& y) {
      using value_type = decltype(x);
      return x > y ? static_cast<value_type>(0) : static_cast<value_type>(1);
    });
GENERATE_BINARY_OVERLOADS(max, MAX_OP)

inline constexpr auto MIN_OP = makeBinaryOp(
    [](const auto& x, const auto& y) {return x < y ? x : y;},
    [](const auto& x, const auto& y) {
      using value_type = decltype(x);
      return x < y ? static_cast<value_type>(1) : static_cast<value_type>(0);
    },
    [](const auto& x, const auto& y) {
      using value_type = decltype(x);
      return x < y ? static_cast<value_type>(0) : static_cast<value_type>(1);
    });
GENERATE_BINARY_OVERLOADS(min, MIN_OP)

inline constexpr auto ATAN2_OP = makeBinaryOp(
    [](const auto& x, const auto& y) {
      CHECK_WITHIN_DOMAIN(x != 0 || y != 0)
      return std::atan2(x, y);
    },
    [](const auto& x, const auto& y) {
      return y/(x*x + y*y);
    },
    [](const auto& x, const auto& y) {
      return -x/(x*x + y*y);}
    );
GENERATE_BINARY_OVERLOADS(atan2, ATAN2_OP)

#undef GENERATE_BINARY_OVERLOADS

}

#endif
