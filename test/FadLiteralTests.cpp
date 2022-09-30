#include "SimpleFad_config.hpp"
#include "SimpleFad_FadLiteral.hpp"
#include "gtest/gtest.h"

#include <limits>
#include "SimpleFad_TestingUtil.hpp"

template <class T>
class FadLiteralTest : public ::testing::Test {
protected:
  using FadType = SimpleFad::FadLiteral<T>;
  static constexpr T zero = 0.;
  static constexpr T one = 1.;
  static constexpr T eps = 10 * std::numeric_limits<T>::epsilon();
};
using MyTypes = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(FadLiteralTest, MyTypes);

TYPED_TEST(FadLiteralTest, Basic) {
  using FadType = typename TestFixture::FadType;
  static constexpr auto zero = TestFixture::zero;
  static constexpr auto one = TestFixture::one;
  static constexpr auto eps = TestFixture::eps;

  FadType fad(one);
  EXPECT_NEAR(fad.val(), one, eps);
  for (typename FadType::index_type i = 0; i < 100; ++i)
  {
    EXPECT_NEAR(fad.dval(i), zero, eps);
  }
  EXPECT_EQ(fad.dsize(), zero);
}

#ifdef CMATH_SUPPORTS_CONSTEXPR

TYPED_TEST(FadLiteralTest, Constexpr) {
  using FadType = typename TestFixture::FadType;
  static constexpr auto zero = TestFixture::zero;
  static constexpr auto one = TestFixture::one;
  static constexpr auto eps = TestFixture::eps;

  constexpr FadType fad(one);
  static_assert(SimpleFad::checkFloatingEquality(fad.val(), one, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad.dval(0), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad.dval(1), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad.dval(100), zero, eps));
  static_assert(fad.dsize() == 0);
}

#endif
