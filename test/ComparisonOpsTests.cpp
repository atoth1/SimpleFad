#include "SimpleFad_ComparisonOps.hpp"
#include "SimpleFad_FadVariable.hpp"

#include "gtest/gtest.h"

template <class T>
class ComparisonOpsTest : public ::testing::Test {
protected:
  static constexpr std::size_t derivSize = 1;
  using ValueType = T;
  using DynamicFad = SimpleFad::DynamicFadVariable<T>;
  using StaticFad = SimpleFad::StaticFadVariable<T, derivSize>;
};
using MyTypes = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(ComparisonOpsTest, MyTypes);

#define FAD_COMPARISON_OPS_TEST(TEST_NAME, FAD_TYPE) \
TYPED_TEST(ComparisonOpsTest, TEST_NAME) { \
  using T = typename TestFixture::ValueType; \
  using DynamicFad = typename TestFixture::DynamicFad; \
  using StaticFad = typename TestFixture::StaticFad; \
  using FadType = FAD_TYPE; \
   \
  static constexpr T zero = 0.; \
  static constexpr T one = 1.; \
  static constexpr T two = 2.; \
  FadType x(one, {one}); \
  EXPECT_TRUE(x == one); \
  EXPECT_TRUE(one == x); \
  EXPECT_TRUE(x != zero); \
  EXPECT_TRUE(zero != x); \
  EXPECT_TRUE(x < two); \
  EXPECT_TRUE(zero < x); \
  EXPECT_TRUE(x > zero); \
  EXPECT_TRUE(two > x); \
  EXPECT_TRUE(x <= one); \
  EXPECT_TRUE(one <= x); \
  EXPECT_TRUE(x >= one); \
  EXPECT_TRUE(one >= x); \
   \
  FadType y(zero, {two}); \
  EXPECT_TRUE(x == x); \
  EXPECT_TRUE(x != y); \
  EXPECT_TRUE(y < x); \
  EXPECT_TRUE(x > y); \
  EXPECT_TRUE(y <= x); \
  EXPECT_TRUE(x >= y); \
}

FAD_COMPARISON_OPS_TEST(StaticFadComparisonOps, StaticFad)
FAD_COMPARISON_OPS_TEST(DynamicFadComparisonOps, DynamicFad)
