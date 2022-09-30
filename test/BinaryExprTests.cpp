#include "SimpleFad_config.hpp"
#include "SimpleFad_BinaryExpr.hpp"
#include "SimpleFad_BinaryOps.hpp"
#include "SimpleFad_ComparisonOps.hpp"
#include "SimpleFad_FadLiteral.hpp"
#include "SimpleFad_FadVariable.hpp"
#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include "SimpleFad_TestingUtil.hpp"

namespace {
  const long double LN2 = std::log(2.l);
  const long double LN3 = std::log(3.l);
  const long double PI = std::acos(-1.l);
}

template <class T>
class BinaryExprTest : public ::testing::Test {
protected:
  static constexpr std::size_t derivSize = 2;
  using ValueType = T;
  using DynamicFad = SimpleFad::DynamicFadVariable<T>;
  using StaticFad = SimpleFad::StaticFadVariable<T, derivSize>;
  using FadLiteral = SimpleFad::FadLiteral<T>;
  static constexpr T eps = 100 * std::numeric_limits<T>::epsilon();
};
using MyTypes = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(BinaryExprTest, MyTypes);

#define FAD_BINARY_OPS_TEST(TEST_NAME, OP_NAME, OP, FAD_TYPE, \
    IN_X0, IN_DX00, IN_DX01, \
    IN_X1, IN_DX10, IN_DX11, \
    IN_C0, IN_C1, \
    OUT_X0, OUT_DX00, OUT_DX01, \
    OUT_X1, OUT_DX10, OUT_DX11, \
    OUT_X2, OUT_DX20, OUT_DX21) \
TYPED_TEST(BinaryExprTest, TEST_NAME) { \
  using T = typename TestFixture::ValueType; \
  using DynamicFad = typename TestFixture::DynamicFad; \
  using StaticFad = typename TestFixture::StaticFad; \
  using FadType = FAD_TYPE; \
  using FadLiteral = typename TestFixture::FadLiteral; \
  static constexpr auto eps = TestFixture::eps; \
 \
  FadType x0(IN_X0, {IN_DX00, IN_DX01}); \
  FadType x1(IN_X1, {IN_DX10, IN_DX11}); \
  auto a = OP(x0, x1); \
  static_assert(std::is_base_of_v< \
      SimpleFad::ExprBase<SimpleFad::BinaryExpr<decltype(SimpleFad::OP_NAME), FadType, FadType>>, \
      decltype(a)>); \
  EXPECT_NEAR(a.val(), static_cast<T>(OUT_X0), eps); \
  EXPECT_EQ(a.dsize(), TestFixture::derivSize); \
  EXPECT_NEAR(a.dval(0), static_cast<T>(OUT_DX00), eps); \
  EXPECT_NEAR(a.dval(1), static_cast<T>(OUT_DX01), eps); \
 \
  T zero = 0.; \
  FadType b(zero, {zero, zero}); \
  b = a; \
  EXPECT_NEAR(b.val(), a.val(), eps); \
  EXPECT_EQ(b.dsize(), TestFixture::derivSize); \
  EXPECT_NEAR(b.dval(0), a.dval(0), eps); \
  EXPECT_NEAR(b.dval(1), a.dval(1), eps); \
 \
  T c0 = IN_C0; \
  auto c = OP(c0, x0); \
  static_assert(std::is_base_of_v< \
      SimpleFad::ExprBase<SimpleFad::BinaryExpr<decltype(SimpleFad::OP_NAME), FadLiteral, FadType>>, \
      decltype(c)>); \
  EXPECT_NEAR(c.val(), static_cast<T>(OUT_X1), eps); \
  EXPECT_EQ(c.dsize(), TestFixture::derivSize); \
  EXPECT_NEAR(c.dval(0), static_cast<T>(OUT_DX10), eps); \
  EXPECT_NEAR(c.dval(1), static_cast<T>(OUT_DX11), eps); \
 \
  b = c; \
  EXPECT_NEAR(b.val(), c.val(), eps); \
  EXPECT_EQ(b.dsize(), TestFixture::derivSize); \
  EXPECT_NEAR(b.dval(0), c.dval(0), eps); \
  EXPECT_NEAR(b.dval(1), c.dval(1), eps); \
 \
  T c1 = IN_C1; \
  auto d = OP(x1, c1); \
  static_assert(std::is_base_of_v< \
      SimpleFad::ExprBase<SimpleFad::BinaryExpr<decltype(SimpleFad::OP_NAME), FadType, FadLiteral>>, \
      decltype(d)>); \
  EXPECT_NEAR(d.val(), static_cast<T>(OUT_X2), eps); \
  EXPECT_EQ(d.dsize(), TestFixture::derivSize); \
  EXPECT_NEAR(d.dval(0), static_cast<T>(OUT_DX20), eps); \
  EXPECT_NEAR(d.dval(1), static_cast<T>(OUT_DX21), eps); \
 \
  b = d; \
  EXPECT_NEAR(b.val(), d.val(), eps); \
  EXPECT_EQ(b.dsize(), TestFixture::derivSize); \
  EXPECT_NEAR(b.dval(0), d.dval(0), eps); \
  EXPECT_NEAR(b.dval(1), d.dval(1), eps); \
}

FAD_BINARY_OPS_TEST(StaticFadOperatorPlus, ADD_OP, operator+, StaticFad,
    2., 2., 1.,
    1., 1., 0.,
    3., 3.,
    3., 3., 1.,
    5., 2., 1.,
    4., 1., 0.)
FAD_BINARY_OPS_TEST(DynamicFadOperatorPlus, ADD_OP, operator+, DynamicFad,
    2., 2., 1.,
    1., 1., 0.,
    3., 3.,
    3., 3., 1.,
    5., 2., 1.,
    4., 1., 0.)

FAD_BINARY_OPS_TEST(StaticFadOperatorMinus, SUBTRACT_OP, operator-, StaticFad,
    2., 2., 1.,
    1., 1., 0.,
    3., 3.,
    1., 1., 1.,
    1., -2., -1.,
    -2., 1., 0.)
FAD_BINARY_OPS_TEST(DynamicFadOperatorMinus, SUBTRACT_OP, operator-, DynamicFad,
    2., 2., 1.,
    1., 1., 0.,
    3., 3.,
    1., 1., 1.,
    1., -2., -1.,
    -2., 1., 0.)

FAD_BINARY_OPS_TEST(StaticFadOperatorTimes, MULTIPLY_OP, operator*, StaticFad,
    2., 2., 1.,
    1., 1., 0.,
    3., 3.,
    2., 4., 1.,
    6., 6., 3.,
    3., 3., 0.)
FAD_BINARY_OPS_TEST(DynamicFadOperatorTimes, MULTIPLY_OP, operator*, DynamicFad,
    2., 2., 1.,
    1., 1., 0.,
    3., 3.,
    2., 4., 1.,
    6., 6., 3.,
    3., 3., 0.)

FAD_BINARY_OPS_TEST(StaticFadOperatorDivide, DIVIDE_OP, operator/, StaticFad,
    2., 2., 1.,
    1., 1., 0.,
    4., 2.,
    2., 0., 1.,
    2., -2., -1.,
    0.5, 0.5, 0.)
FAD_BINARY_OPS_TEST(DynamicFadOperatorDivide, DIVIDE_OP, operator/, DynamicFad,
    2., 2., 1.,
    1., 1., 0.,
    4., 2.,
    2., 0., 1.,
    2., -2., -1.,
    0.5, 0.5, 0.)

FAD_BINARY_OPS_TEST(StaticFadPow, POW_OP, pow, StaticFad,
    2., 2., 1.,
    1., 1., 0.,
    3., 3.,
    2., 2. + 2 * LN2, 1.,
    9., 18 * LN3, 9 * LN3,
    1., 3., 0.)
FAD_BINARY_OPS_TEST(DynamicFadPow, POW_OP, pow, DynamicFad,
    2., 2., 1.,
    1., 1., 0.,
    3., 3.,
    2., 2. + 2 * LN2, 1.,
    9., 18 * LN3, 9 * LN3,
    1., 3., 0.)

// Qualify max/min calls, otherwise versions in std namespace are called
// in first application with two FadType arguments.
FAD_BINARY_OPS_TEST(StaticFadMax, MAX_OP, SimpleFad::max, StaticFad,
    2., 2., 1.,
    1., 1., 0.,
    3., 0.,
    2., 2., 1.,
    3., 0., 0.,
    1., 1., 0.)
FAD_BINARY_OPS_TEST(DynamicFadMax, MAX_OP, SimpleFad::max, DynamicFad,
    2., 2., 1.,
    1., 1., 0.,
    3., 0.,
    2., 2., 1.,
    3., 0., 0.,
    1., 1., 0.)

FAD_BINARY_OPS_TEST(StaticFadMin, MIN_OP, SimpleFad::min, StaticFad,
    2., 2., 1.,
    1., 1., 0.,
    3., 0.,
    1., 1., 0.,
    2., 2., 1.,
    0., 0., 0.)
FAD_BINARY_OPS_TEST(DynamicFadMin, MIN_OP, SimpleFad::min, DynamicFad,
    2., 2., 1.,
    1., 1., 0.,
    3., 0.,
    1., 1., 0.,
    2., 2., 1.,
    0., 0., 0.)

FAD_BINARY_OPS_TEST(StaticFadATan2, ATAN2_OP, atan2, StaticFad,
    2., 2., 1.,
    1., 1., 0.,
    -2., 1.,
    std::atan2(2.l, 1.l), 0., 0.2,
    -PI/4, 0.5, 0.25,
    PI/4, 0.5, 0.)
FAD_BINARY_OPS_TEST(DynamicFadATan2, ATAN2_OP, atan2, DynamicFad,
    2., 2., 1.,
    1., 1., 0.,
    -2., 1.,
    std::atan2(2.l, 1.l), 0., 0.2,
    -PI/4, 0.5, 0.25,
    PI/4, 0.5, 0.)

#ifdef CMATH_SUPPORTS_CONSTEXPR

namespace {

// Let above tests check correctness, just check callability in conxtexpr context.
template<class Scalar>
constexpr SimpleFad::StaticFadVariable<Scalar, 1> constexprTest() {
  constexpr Scalar one = 1.;
  constexpr Scalar two = 2.;
  SimpleFad::StaticFadVariable<Scalar, 1> x(one, {one});
  SimpleFad::StaticFadVariable<Scalar, 1> y(two, {two});
  SimpleFad::StaticFadVariable<Scalar, 1> z(x);

  z = x + y;
  z = x - y;
  z = x * y;
  z = x / y;
  z = pow(x, y);
  // Qualify to avoid std version.
  z = SimpleFad::max(x, y);
  z = SimpleFad::min(x, y);
  // With std::max/std::min.
  z = max(x, y);
  z = min(x, y);
  z = atan2(y, x);
  return z;
}

TYPED_TEST(BinaryExprTest, ConstexprTests) {
  static constexpr auto x = constexprTest<TypeParam>();
}

}

#endif
