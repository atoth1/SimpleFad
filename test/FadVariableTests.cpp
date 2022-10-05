#include "SimpleFad_config.hpp"
#include "SimpleFad_Exception.hpp"
#include "SimpleFad_FadLiteral.hpp"
#include "SimpleFad_FadVariable.hpp"
#include "gtest/gtest.h"

#include <cstddef>
#include <limits>
#include "SimpleFad_TestingUtil.hpp"

template <class T>
class FadVariableTest : public ::testing::Test {
protected:
  using FadLiteral = SimpleFad::FadLiteral<T>;
  static constexpr std::size_t dynamicSize = 5;
  using DynamicFad = SimpleFad::DynamicFadVariable<T>;
  static constexpr std::size_t staticSize = 5;
  using StaticFad = SimpleFad::StaticFadVariable<T, staticSize>;
  static constexpr T zero = 0.;
  static constexpr T one = 1.;
  static constexpr T two = 2.;
  static constexpr T eps = 10 * std::numeric_limits<T>::epsilon();
};
using MyTypes = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(FadVariableTest, MyTypes);

TYPED_TEST(FadVariableTest, DynamicFadConstructors) {
  using DynamicFad = typename TestFixture::DynamicFad;
  using FadLiteral = typename TestFixture::FadLiteral;
  static constexpr auto zero = TestFixture::zero;
  static constexpr auto one = TestFixture::one;
  static constexpr auto two = TestFixture::two;
  static constexpr auto size = TestFixture::dynamicSize;
  static constexpr auto eps = TestFixture::eps;

  // Construct with 0 derivatives
  DynamicFad fad1(size, one);
  EXPECT_NEAR(fad1.val(), one, eps);
  EXPECT_EQ(fad1.dsize(), size);
  EXPECT_NEAR(fad1.dval(0), zero, eps);
  EXPECT_NEAR(fad1.dval(1), zero, eps);
  EXPECT_NEAR(fad1.dval(2), zero, eps);
  EXPECT_NEAR(fad1.dval(3), zero, eps);
  EXPECT_NEAR(fad1.dval(4), zero, eps);

  // Construct with single component derivative initialization
  DynamicFad fad2(size, one, 2, one);
  EXPECT_NEAR(fad2.val(), one, eps);
  EXPECT_EQ(fad2.dsize(), size);
  EXPECT_NEAR(fad2.dval(0), zero, eps);
  EXPECT_NEAR(fad2.dval(1), zero, eps);
  EXPECT_NEAR(fad2.dval(2), one, eps);
  EXPECT_NEAR(fad2.dval(3), zero, eps);
  EXPECT_NEAR(fad2.dval(4), zero, eps);

  // Construct from initializer list
  DynamicFad fad3(one, {one, two, one, two, one});
  EXPECT_NEAR(fad3.val(), one, eps);
  EXPECT_EQ(fad3.dsize(), size);
  EXPECT_NEAR(fad3.dval(0), one, eps);
  EXPECT_NEAR(fad3.dval(1), two, eps);
  EXPECT_NEAR(fad3.dval(2), one, eps);
  EXPECT_NEAR(fad3.dval(3), two, eps);
  EXPECT_NEAR(fad3.dval(4), one, eps);

  // Generated copy constructor
  DynamicFad fad4(fad3);
  EXPECT_NEAR(fad4.val(), fad4.val(), eps);
  EXPECT_EQ(fad4.dsize(), size);
  EXPECT_NEAR(fad4.dval(0), fad3.dval(0), eps);
  EXPECT_NEAR(fad4.dval(1), fad3.dval(1), eps);
  EXPECT_NEAR(fad4.dval(2), fad3.dval(2), eps);
  EXPECT_NEAR(fad4.dval(3), fad3.dval(3), eps);
  EXPECT_NEAR(fad4.dval(4), fad3.dval(4), eps);

  // Copy constructor from other expression type
  FadLiteral literal(one);
  DynamicFad fad5(literal);
  EXPECT_NEAR(fad5.val(), one, eps);
  EXPECT_EQ(fad5.dsize(), 0);
}

TYPED_TEST(FadVariableTest, DynamicFadAssignmentOperator) {
  using DynamicFad = typename TestFixture::DynamicFad;
  using FadLiteral = typename TestFixture::FadLiteral;
  static constexpr auto zero = TestFixture::zero;
  static constexpr auto one = TestFixture::one;
  static constexpr auto two = TestFixture::two;
  static constexpr auto size = TestFixture::dynamicSize;
  static constexpr auto eps = TestFixture::eps;

  // Generated copy assignment
  DynamicFad fad1(one, {one, two, one, two, one});
  DynamicFad fad2(size);
  fad2 = fad1;
  EXPECT_NEAR(fad2.val(), fad1.val(), eps);
  EXPECT_EQ(fad2.dsize(), size);
  EXPECT_NEAR(fad2.dval(0), fad1.dval(0), eps);
  EXPECT_NEAR(fad2.dval(1), fad1.dval(1), eps);
  EXPECT_NEAR(fad2.dval(2), fad1.dval(2), eps);
  EXPECT_NEAR(fad2.dval(3), fad1.dval(3), eps);
  EXPECT_NEAR(fad2.dval(4), fad1.dval(4), eps);

  // Copy assign from other expression type
  FadLiteral literal(one);
  fad2 = literal;
  EXPECT_NEAR(fad2.val(), one, eps);
  EXPECT_EQ(fad2.dsize(), 0);
}

TYPED_TEST(FadVariableTest, DynamicFadUpdateOperators) {
  using DynamicFad = typename TestFixture::DynamicFad;
  static constexpr auto one = TestFixture::one;
  static constexpr auto two = TestFixture::two;
  static constexpr auto eps = TestFixture::eps;

  DynamicFad x(one, {one});
  x += one;
  EXPECT_NEAR(x.val(), two, eps);
  EXPECT_NEAR(x.dval(0), one, eps);
  x -= one;
  EXPECT_NEAR(x.val(), one, eps);
  EXPECT_NEAR(x.dval(0), one, eps);
  x *= two;
  EXPECT_NEAR(x.val(), two, eps);
  EXPECT_NEAR(x.dval(0), two, eps);
  x /= two;
  EXPECT_NEAR(x.val(), one, eps);
  EXPECT_NEAR(x.dval(0), one, eps);

  DynamicFad y(x);
  x += y;
  EXPECT_NEAR(x.val(), two, eps);
  EXPECT_NEAR(x.dval(0), two, eps);
  x -= y;
  EXPECT_NEAR(x.val(), one, eps);
  EXPECT_NEAR(x.dval(0), one, eps);
  x *= y;
  EXPECT_NEAR(x.val(), one, eps);
  EXPECT_NEAR(x.dval(0), two, eps);
  x /= y;
  EXPECT_NEAR(x.val(), one, eps);
  EXPECT_NEAR(x.dval(0), one, eps);
}

#ifdef ENABLE_EXCEPTIONS
TYPED_TEST(FadVariableTest, DynamicFadExpectedExceptions) {
  using DynamicFad = typename TestFixture::DynamicFad;
  static constexpr auto one = TestFixture::one;
  static constexpr auto size = TestFixture::dynamicSize;
  DynamicFad fad(size);
  EXPECT_THROW(fad.dvalNonConst(size) = 0., SimpleFad::BoundsError);
}
#endif

TYPED_TEST(FadVariableTest, StaticFadConstructors) {
  using StaticFad = typename TestFixture::StaticFad;
  using FadLiteral = typename TestFixture::FadLiteral;
  static constexpr auto zero = TestFixture::zero;
  static constexpr auto one = TestFixture::one;
  static constexpr auto two = TestFixture::two;
  static constexpr auto size = TestFixture::staticSize;
  static constexpr auto eps = TestFixture::eps;

  // Construct with 0 derivatives
  StaticFad fad1(size, one);
  EXPECT_NEAR(fad1.val(), one, eps);
  EXPECT_EQ(fad1.dsize(), size);
  EXPECT_NEAR(fad1.dval(0), zero, eps);
  EXPECT_NEAR(fad1.dval(1), zero, eps);
  EXPECT_NEAR(fad1.dval(2), zero, eps);
  EXPECT_NEAR(fad1.dval(3), zero, eps);
  EXPECT_NEAR(fad1.dval(4), zero, eps);

  // Construct with single component derivative initialization
  StaticFad fad2(size, one, 2, one);
  EXPECT_NEAR(fad2.val(), one, eps);
  EXPECT_EQ(fad2.dsize(), size);
  EXPECT_NEAR(fad2.dval(0), zero, eps);
  EXPECT_NEAR(fad2.dval(1), zero, eps);
  EXPECT_NEAR(fad2.dval(2), one, eps);
  EXPECT_NEAR(fad2.dval(3), zero, eps);
  EXPECT_NEAR(fad2.dval(4), zero, eps);

  // Construct from initializer list
  StaticFad fad3(one, {one, two, one, two, one});
  EXPECT_NEAR(fad3.val(), one, eps);
  EXPECT_EQ(fad3.dsize(), size);
  EXPECT_NEAR(fad3.dval(0), one, eps);
  EXPECT_NEAR(fad3.dval(1), two, eps);
  EXPECT_NEAR(fad3.dval(2), one, eps);
  EXPECT_NEAR(fad3.dval(3), two, eps);
  EXPECT_NEAR(fad3.dval(4), one, eps);

  // Generated copy constructor
  StaticFad fad4(fad3);
  EXPECT_NEAR(fad4.val(), fad3.val(), eps);
  EXPECT_EQ(fad4.dsize(), size);
  EXPECT_NEAR(fad4.dval(0), fad3.dval(0), eps);
  EXPECT_NEAR(fad4.dval(1), fad3.dval(1), eps);
  EXPECT_NEAR(fad4.dval(2), fad3.dval(2), eps);
  EXPECT_NEAR(fad4.dval(3), fad3.dval(3), eps);
  EXPECT_NEAR(fad4.dval(4), fad3.dval(4), eps);

  // Copy constructor from other expression type
  FadLiteral literal(one);
  StaticFad fad5(literal);
  EXPECT_NEAR(fad5.val(), one, eps);
  EXPECT_EQ(fad5.dsize(), size);
  EXPECT_NEAR(fad5.dval(0), zero, eps);
  EXPECT_NEAR(fad5.dval(1), zero, eps);
  EXPECT_NEAR(fad5.dval(2), zero, eps);
  EXPECT_NEAR(fad5.dval(3), zero, eps);
  EXPECT_NEAR(fad5.dval(4), zero, eps);
}

TYPED_TEST(FadVariableTest, StaticFadAssignmentOperator) {
  using StaticFad = typename TestFixture::StaticFad;
  using FadLiteral = typename TestFixture::FadLiteral;
  static constexpr auto zero = TestFixture::zero;
  static constexpr auto one = TestFixture::one;
  static constexpr auto two = TestFixture::two;
  static constexpr auto size = TestFixture::staticSize;
  static constexpr auto eps = TestFixture::eps;

  // Generated copy assignment
  StaticFad fad1(one, {one, two, one, two, one});
  StaticFad fad2(size);
  fad2 = fad1;
  EXPECT_NEAR(fad2.val(), fad2.val(), eps);
  EXPECT_EQ(fad2.dsize(), size);
  EXPECT_NEAR(fad2.dval(0), fad2.dval(0), eps);
  EXPECT_NEAR(fad2.dval(1), fad2.dval(1), eps);
  EXPECT_NEAR(fad2.dval(2), fad2.dval(2), eps);
  EXPECT_NEAR(fad2.dval(3), fad2.dval(3), eps);
  EXPECT_NEAR(fad2.dval(4), fad2.dval(4), eps);

  // Copy assign from other expression type
  FadLiteral literal(one);
  fad2 = literal;
  EXPECT_NEAR(fad2.val(), one, eps);
  EXPECT_EQ(fad2.dsize(), size);
  EXPECT_NEAR(fad2.dval(0), zero, eps);
  EXPECT_NEAR(fad2.dval(1), zero, eps);
  EXPECT_NEAR(fad2.dval(2), zero, eps);
  EXPECT_NEAR(fad2.dval(3), zero, eps);
  EXPECT_NEAR(fad2.dval(4), zero, eps);
}

TYPED_TEST(FadVariableTest, StaticFadUpdateOperators) {
  using StaticFad = SimpleFad::StaticFadVariable<TypeParam, 1>;
  static constexpr auto one = TestFixture::one;
  static constexpr auto two = TestFixture::two;
  static constexpr auto eps = TestFixture::eps;

  StaticFad x(one, {one});
  x += one;
  EXPECT_NEAR(x.val(), two, eps);
  EXPECT_NEAR(x.dval(0), one, eps);
  x -= one;
  EXPECT_NEAR(x.val(), one, eps);
  EXPECT_NEAR(x.dval(0), one, eps);
  x *= two;
  EXPECT_NEAR(x.val(), two, eps);
  EXPECT_NEAR(x.dval(0), two, eps);
  x /= two;
  EXPECT_NEAR(x.val(), one, eps);
  EXPECT_NEAR(x.dval(0), one, eps);

  StaticFad y(x);
  x += y;
  EXPECT_NEAR(x.val(), two, eps);
  EXPECT_NEAR(x.dval(0), two, eps);
  x -= y;
  EXPECT_NEAR(x.val(), one, eps);
  EXPECT_NEAR(x.dval(0), one, eps);
  x *= y;
  EXPECT_NEAR(x.val(), one, eps);
  EXPECT_NEAR(x.dval(0), two, eps);
  x /= y;
  EXPECT_NEAR(x.val(), one, eps);
  EXPECT_NEAR(x.dval(0), one, eps);
}

#ifdef ENABLE_EXCEPTIONS
TYPED_TEST(FadVariableTest, StaticFadExpectedExceptions) {
  using StaticFad = typename TestFixture::StaticFad;
  using MismatchedStaticFad = SimpleFad::StaticFadVariable<TypeParam, 4>;
  static constexpr auto one = TestFixture::one;
  static constexpr auto size = TestFixture::staticSize;

  StaticFad fad1(5);
  EXPECT_THROW(fad1.dvalNonConst(size) = one, SimpleFad::BoundsError);
  MismatchedStaticFad fad2(4);
  EXPECT_THROW(StaticFad fad3(fad2), SimpleFad::BoundsError);
  EXPECT_THROW(fad1 = fad2, SimpleFad::BoundsError);
}
#endif

#ifdef CMATH_SUPPORTS_CONSTEXPR

TYPED_TEST(FadVariableTest, StaticFadConstexprConstructors) {
  using StaticFad = typename TestFixture::StaticFad;
  using FadLiteral = typename TestFixture::FadLiteral;
  static constexpr auto zero = TestFixture::zero;
  static constexpr auto one = TestFixture::one;
  static constexpr auto two = TestFixture::two;
  static constexpr auto size = TestFixture::staticSize;
  static constexpr auto eps = TestFixture::eps;

  // Construct with 0 derivatives
  constexpr StaticFad fad1(size, one);
  static_assert(SimpleFad::checkFloatingEquality(fad1.val(), one, eps));
  static_assert(fad1.dsize() == size);
  static_assert(SimpleFad::checkFloatingEquality(fad1.dval(0), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad1.dval(1), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad1.dval(2), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad1.dval(3), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad1.dval(4), zero, eps));

  // Construct with single component derivative initialization
  constexpr StaticFad fad2(size, one, 2, one);
  static_assert(SimpleFad::checkFloatingEquality(fad2.val(), one, eps));
  static_assert(fad2.dsize() == size);
  static_assert(SimpleFad::checkFloatingEquality(fad2.dval(0), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad2.dval(1), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad2.dval(2), one, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad2.dval(3), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad2.dval(4), zero, eps));

  // Construct from initializer list
  constexpr StaticFad fad3(one, {one, two, one, two, one});
  static_assert(SimpleFad::checkFloatingEquality(fad3.val(), one, eps));
  static_assert(fad3.dsize() == size);
  static_assert(SimpleFad::checkFloatingEquality(fad3.dval(0), one, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad3.dval(1), two, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad3.dval(2), one, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad3.dval(3), two, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad3.dval(4), one, eps));

  // Generated copy constructor
  constexpr StaticFad fad4(fad3);
  static_assert(SimpleFad::checkFloatingEquality(fad4.val(), fad3.val(), eps));
  static_assert(fad4.dsize() == size);
  static_assert(SimpleFad::checkFloatingEquality(fad4.dval(0), fad3.dval(0), eps));
  static_assert(SimpleFad::checkFloatingEquality(fad4.dval(1), fad3.dval(1), eps));
  static_assert(SimpleFad::checkFloatingEquality(fad4.dval(2), fad3.dval(2), eps));
  static_assert(SimpleFad::checkFloatingEquality(fad4.dval(3), fad3.dval(3), eps));
  static_assert(SimpleFad::checkFloatingEquality(fad4.dval(4), fad3.dval(4), eps));

  // Copy constructor from other expression type
  constexpr FadLiteral literal(one);
  constexpr StaticFad fad5(literal);
  static_assert(SimpleFad::checkFloatingEquality(fad5.val(), one, eps));
  static_assert(fad5.dsize() == size);
  static_assert(SimpleFad::checkFloatingEquality(fad5.dval(0), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad5.dval(1), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad5.dval(2), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad5.dval(3), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(fad5.dval(4), zero, eps));
}

namespace {
template<class Scalar>
constexpr SimpleFad::StaticFadVariable<Scalar, 5> constexprGeneratedAssignment() {
  using FadType = SimpleFad::StaticFadVariable<Scalar, 5>;
  FadType fad1(5, 1.);
  FadType fad2(5, 2., 2, 3.);
  fad2 = fad1;
  return fad2;
}

template<class Scalar>
constexpr SimpleFad::StaticFadVariable<Scalar, 5> constexprTemplatedAssignment() {
  SimpleFad::FadLiteral<Scalar> fad1(1.);
  SimpleFad::StaticFadVariable<Scalar, 5> fad2(5, 2., 2, 3.);
  fad2 = fad1;
  return fad2;
}
}

TYPED_TEST(FadVariableTest, StaticFadConstexprAssignmentOperator) {
  static constexpr auto zero = TestFixture::zero;
  static constexpr auto one = TestFixture::one;
  static constexpr auto two = TestFixture::two;
  static constexpr auto size = TestFixture::staticSize;
  static constexpr auto eps = TestFixture::eps;

  constexpr auto result1 = constexprGeneratedAssignment<TypeParam>();
  static_assert(SimpleFad::checkFloatingEquality(result1.val(), one, eps));
  static_assert(result1.dsize() == size);
  static_assert(SimpleFad::checkFloatingEquality(result1.dval(0), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(result1.dval(1), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(result1.dval(2), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(result1.dval(3), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(result1.dval(4), zero, eps));

  constexpr auto result2 = constexprTemplatedAssignment<TypeParam>();
  static_assert(SimpleFad::checkFloatingEquality(result2.val(), one, eps));
  static_assert(result2.dsize() == 5);
  static_assert(SimpleFad::checkFloatingEquality(result2.dval(0), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(result2.dval(1), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(result2.dval(2), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(result2.dval(3), zero, eps));
  static_assert(SimpleFad::checkFloatingEquality(result2.dval(4), zero, eps));
}

#endif
