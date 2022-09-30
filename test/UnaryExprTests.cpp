#include "SimpleFad_config.hpp"
#include "SimpleFad_UnaryExpr.hpp"
#include "SimpleFad_UnaryOps.hpp"
#include "SimpleFad_FadVariable.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <limits>
#include <type_traits>

namespace {
  const long double E = std::exp(1.l);
  const long double LN2 = std::log(2.l);
  const long double LN10 = std::log(10.l);
  const long double PI = std::acos(-1.l);
  const long double SQRT2 = std::sqrt(2.l);
}

template <class T>
class UnaryExprTest : public ::testing::Test {
protected:
  static constexpr std::size_t derivSize = 2;
  using ValueType = T;
  using DynamicFad = SimpleFad::DynamicFadVariable<T>;
  using StaticFad = SimpleFad::StaticFadVariable<T, derivSize>;
  static constexpr T eps = 100 * std::numeric_limits<T>::epsilon();
};
using MyTypes = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(UnaryExprTest, MyTypes);

#define UNARY_OPS_TEST(TEST_NAME, OP_NAME, OP, FAD_TYPE, \
    IN_X, IN_DX0, IN_DX1, \
    OUT_X, OUT_DX0, OUT_DX1) \
TYPED_TEST(UnaryExprTest, TEST_NAME) { \
  using T = typename TestFixture::ValueType; \
  using DynamicFad = typename TestFixture::DynamicFad; \
  using StaticFad = typename TestFixture::StaticFad; \
  using FadType = FAD_TYPE; \
  static constexpr auto eps = TestFixture::eps; \
 \
  FadType x(IN_X, {IN_DX0, IN_DX1}); \
  auto a = OP(x); \
  static_assert(std::is_base_of_v< \
      SimpleFad::ExprBase<SimpleFad::UnaryExpr<decltype(SimpleFad::OP_NAME), FadType>>, \
      decltype(a)>); \
  EXPECT_NEAR(a.val(), OUT_X, eps); \
  EXPECT_EQ(a.dsize(), TestFixture::derivSize); \
  EXPECT_NEAR(a.dval(0), OUT_DX0, eps); \
  EXPECT_NEAR(a.dval(1), OUT_DX1, eps); \
 \
  T zero = 0.; \
  FadType b(zero, {zero, zero}); \
  b = a; \
  EXPECT_NEAR(b.val(), a.val(), eps); \
  EXPECT_EQ(a.dsize(), TestFixture::derivSize); \
  EXPECT_NEAR(b.dval(0), a.dval(0), eps); \
  EXPECT_NEAR(b.dval(1), a.dval(1), eps); \
}
UNARY_OPS_TEST(StaticFadOperatorPlus, UNARY_PLUS_OP, operator+, StaticFad,
    1., 2., 3.,
    1., 2., 3.)
UNARY_OPS_TEST(DynamicFadOperatorPlus, UNARY_PLUS_OP, operator+, DynamicFad,
    1., 2., 3.,
    1., 2., 3.)

UNARY_OPS_TEST(StaticFadOperatorMinus, UNARY_MINUS_OP, operator-, StaticFad,
    1., 2., 3.,
    -1., -2., -3.)
UNARY_OPS_TEST(DynamicFadOperatorMinus, UNARY_MINUS_OP, operator-, DynamicFad,
    1., 2., 3.,
    -1., -2., -3.)

UNARY_OPS_TEST(StaticFadAbs, ABS_OP, abs, StaticFad,
    1., 2., -3.,
    1., 2., -3.)
UNARY_OPS_TEST(DynamicFadAbs, ABS_OP, abs, DynamicFad,
    -1., 2., -3.,
    1., -2., 3.)

UNARY_OPS_TEST(StaticFadExp, EXP_OP, exp, StaticFad,
    1., 2., -3.,
    E, 2 * E, -3 * E)
UNARY_OPS_TEST(DynamicFadExp, EXP_OP, exp, DynamicFad,
    1., 2., -3.,
    E, 2 * E, -3 * E)

UNARY_OPS_TEST(StaticFadExp2, EXP2_OP, exp2, StaticFad,
    1., 2., -3.,
    2., 4 * LN2, -6 * LN2)
UNARY_OPS_TEST(DynamicFadExp2, EXP2_OP, exp2, DynamicFad,
    1., 2., -3.,
    2., 4 * LN2, -6 * LN2)

UNARY_OPS_TEST(StaticFadLog, LOG_OP, log, StaticFad,
    1., 2., -3.,
    0., 2., -3.)
UNARY_OPS_TEST(DynamicFadLog, LOG_OP, log, DynamicFad,
    1., 2., -3.,
    0., 2., -3.)

UNARY_OPS_TEST(StaticFadLog2, LOG2_OP, log2, StaticFad,
    4., 4., -8.,
    2., 1./LN2, -2./LN2)
UNARY_OPS_TEST(DynamicFadLog2, LOG2_OP, log2, DynamicFad,
    4., 4., -8.,
    2., 1./LN2, -2./LN2)

UNARY_OPS_TEST(StaticFadLog10, LOG10_OP, log10, StaticFad,
    100., 2., -1.,
    2., 1./(50. * LN10), -1./(100. * LN10))
UNARY_OPS_TEST(DynamicFadLog10, LOG10_OP, log10, DynamicFad,
    100., 2., -1.,
    2., 1./(50. * LN10), -1./(100. * LN10))

UNARY_OPS_TEST(StaticFadSqrt, SQRT_OP, sqrt, StaticFad,
    4., 4., -8.,
    2., 1., -2.)
UNARY_OPS_TEST(DynamicFadSqrt, SQRT_OP, sqrt, DynamicFad,
    4., 4., -8.,
    2., 1., -2.)

UNARY_OPS_TEST(StaticFadCbrt, CBRT_OP, cbrt, StaticFad,
    8., 12., -24.,
    2., 1., -2.)
UNARY_OPS_TEST(DynamicFadCbrt, CBRT_OP, cbrt, DynamicFad,
    8., 12., -24.,
    2., 1., -2.)

UNARY_OPS_TEST(StaticFadCos, COS_OP, cos, StaticFad,
    PI/2., 1., 2.,
    0., -1., -2.)
UNARY_OPS_TEST(DynamicFadCos, COS_OP, cos, DynamicFad,
    PI, 1., 2.,
    -1., 0., 0.)

UNARY_OPS_TEST(StaticFadSin, SIN_OP, sin, StaticFad,
    PI/2, 1., 2.,
    1., 0., 0.)
UNARY_OPS_TEST(DynamicFadSin, SIN_OP, sin, DynamicFad,
    PI, 1., 2.,
    0, -1., -2.)

UNARY_OPS_TEST(StaticFadTan, TAN_OP, tan, StaticFad,
    PI/4., 1., -2.,
    1., 2., -4.)
UNARY_OPS_TEST(DynamicFadTan, TAN_OP, tan, DynamicFad,
    -PI/4, 1., -2.,
    -1., 2., -4.)

UNARY_OPS_TEST(StaticFadACos, ACOS_OP, acos, StaticFad,
    1/SQRT2, 1., -2.,
    PI/4, -SQRT2, 2*SQRT2)
UNARY_OPS_TEST(DynamicFadACos, ACOS_OP, acos, DynamicFad,
    -1/SQRT2, 1, -2,
    3*PI/4, -SQRT2, 2*SQRT2)

UNARY_OPS_TEST(StaticFadASin, ASIN_OP, asin, StaticFad,
    1/SQRT2, 1., -2.,
    PI/4, SQRT2, -2*SQRT2)
UNARY_OPS_TEST(DynamicFadASin, ASIN_OP, asin, DynamicFad,
    -1/SQRT2, 1., -2.,
    -PI/4, SQRT2, -2*SQRT2)

UNARY_OPS_TEST(StaticFadATan, ATAN_OP, atan, StaticFad,
    1., 2., -4.,
    PI/4, 1., -2.)
UNARY_OPS_TEST(DynamicFadATan, ATAN_OP, atan, DynamicFad,
    -1., 2., -4.,
    -PI/4, 1., -2.)

UNARY_OPS_TEST(StaticFadCosh, COSH_OP, cosh, StaticFad,
    1., 1., -2.,
    std::cosh(1.l), std::sinh(1.l), -2*std::sinh(1.l))
UNARY_OPS_TEST(DynamicFadCosh, COSH_OP, cosh, DynamicFad,
    -1., 1., -2.,
    std::cosh(-1.l), std::sinh(-1.l), -2*std::sinh(-1.l))

UNARY_OPS_TEST(StaticFadSinh, SINH_OP, sinh, StaticFad,
    1., 1., -2.,
    std::sinh(1.l), std::cosh(1.l), -2*std::cosh(1.l))
UNARY_OPS_TEST(DynamicFadSinh, SINH_OP, sinh, DynamicFad,
    -1., 1., -2.,
    std::sinh(-1.l), std::cosh(-1.l), -2*std::cosh(-1.l))

UNARY_OPS_TEST(StaticFadTanh, TANH_OP, tanh, StaticFad,
    1., 1., -2.,
    std::tanh(1.l), 1/(std::cosh(1.l)*std::cosh(1.l)), -2/(std::cosh(1.l)*std::cosh(1.l)))
UNARY_OPS_TEST(DynamicFadTanh, TANH_OP, tanh, DynamicFad,
    -1., 1., -2.,
    std::tanh(-1.l), 1/(std::cosh(-1.l)*std::cosh(-1.l)), -2/(std::cosh(-1.l)*std::cosh(-1.l)))

UNARY_OPS_TEST(StaticFadACosh, ACOSH_OP, acosh, StaticFad,
    std::sqrt(5.l), 2., -2.,
    std::acosh(std::sqrt(5.l)), 1., -1.)
UNARY_OPS_TEST(DynamicFadACosh, ACOSH_OP, acosh, DynamicFad,
    std::sqrt(17.l), 4., -4.,
    std::acosh(std::sqrt(17.l)), 1., -1.)

UNARY_OPS_TEST(StaticFadASinh, ASINH_OP, asinh, StaticFad,
    std::sqrt(3.l), 2., -2.,
    std::asinh(std::sqrt(3.l)), 1., -1.)
UNARY_OPS_TEST(DynamicFadASinh, ASINH_OP, asinh, DynamicFad,
    std::sqrt(15.l), 4., -4.,
    std::asinh(std::sqrt(15.l)), 1., -1.)

UNARY_OPS_TEST(StaticFadATanh, ATANH_OP, atanh, StaticFad,
    1/SQRT2, 1., -2.,
    std::atanh(1/SQRT2), 2., -4.)
UNARY_OPS_TEST(DynamicFadATanh, ATANH_OP, atanh, DynamicFad,
    1/SQRT2, 1., -2.,
    std::atanh(1/SQRT2), 2., -4.)

#ifdef CMATH_SUPPORTS_CONSTEXPR

namespace {

// Let above tests check correctness, just check callability in conxtexpr context.
template<class Scalar>
constexpr SimpleFad::StaticFadVariable<Scalar, 1> constexprTest() {
  constexpr Scalar zero = 0.;
  constexpr Scalar two = 2.;
  SimpleFad::StaticFadVariable<Scalar, 1> x(zero, {zero});
  SimpleFad::StaticFadVariable<Scalar, 1> y(x);
  SimpleFad::StaticFadVariable<Scalar, 1> z(two, {two});

  y = +x;
  y = -x;
  y = abs(x);
  y = exp(x);
  y = exp2(x);
  y = log(z);
  y = log2(z);
  y = log10(z);
  y = sqrt(z);
  y = cbrt(z);
  y = cos(x);
  y = sin(x);
  y = tan(x);
  y = acos(x);
  y = asin(x);
  y = atan(x);
  y = cosh(x);
  y = sinh(x);
  y = tanh(x);
  y = acosh(z);
  y = asinh(x);
  y = atanh(x);
  return y;
}

TYPED_TEST(UnaryExprTest, ConstexprTest) {
  static constexpr auto x = constexprTest<TypeParam>();
}

}

#endif
