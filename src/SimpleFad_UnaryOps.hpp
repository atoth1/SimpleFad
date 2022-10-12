#ifndef SIMPLEFAD_UNARYOPS_HPP
#define SIMPLEFAD_UNARYOPS_HPP

#include "SimpleFad_config.hpp"
#include "SimpleFad_Exception.hpp"
#include "SimpleFad_ExprBase.hpp"
#include "SimpleFad_UnaryExpr.hpp"

#include <cmath>

namespace SimpleFad {

template<class FOp, class DFOp>
class UnaryOp {
public:
  constexpr explicit UnaryOp(const FOp& f, const DFOp& df)
      : f_(f), df_(df) {

  }

  template<class T>
  SIMPLEFAD_CONSTEXPR T f(const T& x) const {
    return f_(x);
  }

  template<class T>
  SIMPLEFAD_CONSTEXPR T df(const T& x) const {
    return df_(x);
  }

private:
  FOp f_;
  DFOp df_;
};

template<class FOp, class DFOp>
constexpr auto makeUnaryOp(const FOp& f, const DFOp& df) {
  return UnaryOp<FOp, DFOp>(f, df);
}

#define GENERATE_UNARY_OVERLOADS(OP_NAME, OP) \
  template<class Expr> \
  constexpr auto \
  OP_NAME (const ExprBase<Expr>& e) { \
    return UnaryExpr<decltype(OP), Expr>(OP, e.self()); \
  }

inline constexpr auto UNARY_PLUS_OP = makeUnaryOp(
    [](const auto& x) {return +x;},
    [](const auto& x) {return static_cast<decltype(x)>(1);});
GENERATE_UNARY_OVERLOADS(operator+, UNARY_PLUS_OP)

inline constexpr auto UNARY_MINUS_OP = makeUnaryOp(
    [](const auto& x) {return -x;},
    [](const auto& x) {return static_cast<decltype(x)>(-1);});
GENERATE_UNARY_OVERLOADS(operator-, UNARY_MINUS_OP)

inline constexpr auto ABS_OP = makeUnaryOp(
    [](const auto& x) {return std::abs(x);},
    [](const auto& x) {
      using value_type = decltype(x);
      return x < static_cast<value_type>(0) ?
          static_cast<value_type>(-1) : static_cast<value_type>(1);
    });
GENERATE_UNARY_OVERLOADS(abs, ABS_OP)

inline constexpr auto EXP_OP = makeUnaryOp(
    [](const auto& x) {return std::exp(x);},
    [](const auto& x) {return std::exp(x);});
GENERATE_UNARY_OVERLOADS(exp, EXP_OP)

inline constexpr auto EXP2_OP = makeUnaryOp(
    [](const auto& x) {return std::exp2(x);},
    [](const auto& x) {
      using value_type = decltype(x);
      return std::exp2(x) * std::log(static_cast<value_type>(2));
    });
GENERATE_UNARY_OVERLOADS(exp2, EXP2_OP)

inline constexpr auto LOG_OP = makeUnaryOp(
    [](const auto& x) {
      CHECK_WITHIN_DOMAIN(x > static_cast<decltype(x)>(0))
      return std::log(x);
    },
    [](const auto& x) {
      using value_type = decltype(x);
      return static_cast<value_type>(1)/x;
    });
GENERATE_UNARY_OVERLOADS(log, LOG_OP)

inline constexpr auto LOG2_OP = makeUnaryOp(
    [](const auto& x) {
      CHECK_WITHIN_DOMAIN(x > static_cast<decltype(x)>(0))
      return std::log2(x);
    },
    [](const auto& x) {
      using value_type = decltype(x);
      auto coef = std::log(static_cast<value_type>(2));
      return static_cast<value_type>(1)/(coef*x);
    });
GENERATE_UNARY_OVERLOADS(log2, LOG2_OP)

inline constexpr auto LOG10_OP = makeUnaryOp(
    [](const auto& x) {
      CHECK_WITHIN_DOMAIN(x > static_cast<decltype(x)>(0))
      return std::log10(x);
    },
    [](const auto& x) {
      using value_type = decltype(x);
      auto coef = std::log(static_cast<value_type>(10));
      return static_cast<value_type>(1)/(coef*x);
    });
GENERATE_UNARY_OVERLOADS(log10, LOG10_OP)

inline constexpr auto SQRT_OP = makeUnaryOp(
    [](const auto& x) {
      // Value obviously fine at 0, but infinite derivative.
      using value_type = decltype(x);
      CHECK_WITHIN_DOMAIN(x > static_cast<value_type>(0))
      return std::sqrt(x);
    },
    [](const auto& x) {
      using value_type = decltype(x);
      auto coef = static_cast<value_type>(2);
      return static_cast<value_type>(1)/(coef*std::sqrt(x));
    });
GENERATE_UNARY_OVERLOADS(sqrt, SQRT_OP)

inline constexpr auto CBRT_OP = makeUnaryOp(
    [](const auto& x) {
      // Value obviously fine, but infinite derivative at 0.
      using value_type = decltype(x);
      CHECK_WITHIN_DOMAIN(x != static_cast<value_type>(0))
      return std::cbrt(x);
    },
    [](const auto& x) {
      using value_type = decltype(x);
      auto coef = static_cast<value_type>(3);
      auto cbrt = std::cbrt(x);
      return static_cast<value_type>(1)/(coef*cbrt*cbrt);
    });
GENERATE_UNARY_OVERLOADS(cbrt, CBRT_OP)

inline constexpr auto COS_OP = makeUnaryOp(
    [](const auto& x) {return std::cos(x);},
    [](const auto& x) {return -std::sin(x);});
GENERATE_UNARY_OVERLOADS(cos, COS_OP)

inline constexpr auto SIN_OP = makeUnaryOp(
    [](const auto& x) {return std::sin(x);},
    [](const auto& x) {return std::cos(x);});
GENERATE_UNARY_OVERLOADS(sin, SIN_OP)

inline constexpr auto TAN_OP = makeUnaryOp(
    [](const auto& x) {return std::tan(x);},
    [](const auto& x) {
      using value_type = decltype(x);
      auto sec = static_cast<value_type>(1)/std::cos(x);
      return sec*sec;
    });
GENERATE_UNARY_OVERLOADS(tan, TAN_OP)

inline constexpr auto ACOS_OP = makeUnaryOp(
    [](const auto& x) {
      // Infinite derivatives at +/-1.
      using value_type = decltype(x);
      CHECK_WITHIN_DOMAIN(x > static_cast<value_type>(-1) && x < static_cast<value_type>(1))
      return std::acos(x);
    },
    [](const auto& x) {
      auto one = static_cast<decltype(x)>(1);
      return -one/std::sqrt(one - x*x);
    });
GENERATE_UNARY_OVERLOADS(acos, ACOS_OP)

inline constexpr auto ASIN_OP = makeUnaryOp(
    [](const auto& x) {
      // Infinite derivatives at +/-1.
      using value_type = decltype(x);
      CHECK_WITHIN_DOMAIN(x > static_cast<value_type>(-1) && x < static_cast<value_type>(1))
      return std::asin(x);
    },
    [](const auto& x) {
      auto one = static_cast<decltype(x)>(1);
      return one/std::sqrt(one - x*x);
    });
GENERATE_UNARY_OVERLOADS(asin, ASIN_OP)

inline constexpr auto ATAN_OP = makeUnaryOp(
    [](const auto& x) {return std::atan(x);},
    [](const auto& x) {
      auto one = static_cast<decltype(x)>(1);
      return one/(one + x*x);
    });
GENERATE_UNARY_OVERLOADS(atan, ATAN_OP)

inline constexpr auto COSH_OP = makeUnaryOp(
    [](const auto& x) {return std::cosh(x);},
    [](const auto& x) {return std::sinh(x);});
GENERATE_UNARY_OVERLOADS(cosh, COSH_OP)

inline constexpr auto SINH_OP = makeUnaryOp(
    [](const auto& x) {return std::sinh(x);},
    [](const auto& x) {return std::cosh(x);});
GENERATE_UNARY_OVERLOADS(sinh, SINH_OP)

inline constexpr auto TANH_OP = makeUnaryOp(
    [](const auto& x) {return std::tanh(x);},
    [](const auto& x) {
      using value_type = decltype(x);
      auto sech = static_cast<value_type>(1)/std::cosh(x);
      return sech*sech;
    });
GENERATE_UNARY_OVERLOADS(tanh, TANH_OP)

inline constexpr auto ACOSH_OP = makeUnaryOp(
    [](const auto& x) {
      // Value fine at 1, but derivative is infinite.
      using value_type = decltype(x);
      CHECK_WITHIN_DOMAIN(x > static_cast<value_type>(1))
      return std::acosh(x);
    },
    [](const auto& x) {
      auto one = static_cast<decltype(x)>(1);
      return one/std::sqrt(x*x - one);
    });
GENERATE_UNARY_OVERLOADS(acosh, ACOSH_OP)

inline constexpr auto ASINH_OP = makeUnaryOp(
    [](const auto& x) {return std::asinh(x);},
    [](const auto& x) {
      auto one = static_cast<decltype(x)>(1);
      return one/std::sqrt(x*x + one);
    });
GENERATE_UNARY_OVERLOADS(asinh, ASINH_OP)

inline constexpr auto ATANH_OP = makeUnaryOp(
    [](const auto& x) {
      using value_type = decltype(x);
      CHECK_WITHIN_DOMAIN(x > static_cast<value_type>(-1) && x < static_cast<value_type>(1))
      return std::atanh(x);
    },
    [](const auto& x) {
      auto one = static_cast<decltype(x)>(1);
      return one/(one - x*x);
    });
GENERATE_UNARY_OVERLOADS(atanh, ATANH_OP)

#undef GENERATE_UNARY_OVERLOADS
}

#endif
