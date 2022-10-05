#include "SimpleFad.hpp"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <vector>

template <class T>
class EvaluationTest : public ::testing::Test {
protected:
  using StaticFad = SimpleFad::StaticFadVariable<T, 3>;
  using DynamicFad = SimpleFad::DynamicFadVariable<T>;
  static constexpr T eps = 10 * std::numeric_limits<T>::epsilon();
};
using MyTypes = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(EvaluationTest, MyTypes);

TYPED_TEST(EvaluationTest, LinearCombination) {
  using FadType = typename TestFixture::StaticFad;
  static constexpr auto eps = TestFixture::eps;

  FadType x{3, 1., 0, 1.};
  FadType y{3, 2., 1, 1.};
  FadType z{3, 3., 2, 1.};

  FadType result = 3*x + 2*y + z;
  EXPECT_NEAR(result.val(), 10., eps);
  EXPECT_NEAR(result.dval(0), 3., eps);
  EXPECT_NEAR(result.dval(1), 2., eps);
  EXPECT_NEAR(result.dval(2), 1., eps);
}

TYPED_TEST(EvaluationTest, MatVec) {
  using FadType = typename TestFixture::StaticFad;
  using VectorType = std::array<FadType, 3>;
  using MatrixType = std::array<std::array<TypeParam, 3>, 3>;
  static constexpr auto eps = TestFixture::eps;

  MatrixType A{std::array<TypeParam, 3>{1., 2., 3.},
               std::array<TypeParam, 3>{4., 5., 6.},
               std::array<TypeParam, 3>{7., 8., 9.}};
  VectorType x{FadType{3, 1., 0, 1.}, FadType{3, 2., 1, 1.}, FadType{3, 3., 2, 1.}};
  VectorType y{FadType{}, FadType{}, FadType{}};

  std::array<TypeParam, 3> expected{14., 32., 50.};
  for (int row = 0; row < A.size(); ++row) {
    for (int col = 0; col < A[row].size(); ++col) {
      y[row] += A[row][col] * x[col];
    }
    EXPECT_NEAR(expected[row], y[row].val(), eps);
    for (int col = 0; col < A[row].size(); ++col) {
      EXPECT_NEAR(A[row][col], y[row].dval(col), eps);
    }
  }
}

TYPED_TEST(EvaluationTest, Polynomials) {
  using FadType = typename TestFixture::StaticFad;
  static constexpr auto eps = TestFixture::eps;

  FadType x{3, 1., 0, 1.};
  FadType y{3, 2., 1, 1.};
  FadType z{3, 3., 2, 1.};
  FadType w;

  w = +x*y*z;
  EXPECT_NEAR(w.val(), 6., eps);
  EXPECT_NEAR(w.dval(0), 6., eps);
  EXPECT_NEAR(w.dval(1), 3., eps);
  EXPECT_NEAR(w.dval(2), 2., eps);

  w = x*x + y*y;
  EXPECT_NEAR(w.val(), 5., eps);
  EXPECT_NEAR(w.dval(0), 2., eps);
  EXPECT_NEAR(w.dval(1), 4., eps);
  EXPECT_NEAR(w.dval(2), 0., eps);

  w = 3*x*x*z - y*z + 4*x*z*z*z;
  EXPECT_NEAR(w.val(), 111., eps);
  EXPECT_NEAR(w.dval(0), 126., eps);
  EXPECT_NEAR(w.dval(1), -3., eps);
  EXPECT_NEAR(w.dval(2), 109., eps);
}

TYPED_TEST(EvaluationTest, PowChainRule) {
  using FadType = typename TestFixture::StaticFad;
  static constexpr auto eps = TestFixture::eps;

  FadType x{3, 1., 0, 1.};
  FadType y{3, 2., 1, 1.};
  FadType z{3, 3., 2, 1.};
  FadType w;

  w = pow(x+y, z/2.);
  EXPECT_NEAR(w.val(), std::pow(x.val()+y.val(), 1.5), eps);
  EXPECT_NEAR(w.dval(0), 1.5*std::pow(x.val()+y.val(), 0.5), eps);
  EXPECT_NEAR(w.dval(1), 1.5*std::pow(x.val()+y.val(), 0.5), eps);
  EXPECT_NEAR(w.dval(2), std::log(x.val()+y.val())*std::pow(x.val()+y.val(), 1.5)/2., eps);
}

TYPED_TEST(EvaluationTest, ExpChainRule) {
  using FadType = typename TestFixture::StaticFad;
  static constexpr auto eps = TestFixture::eps;
  FadType x{3, 1., 0, 1.};
  FadType y{3, 2., 1, 1.};
  FadType z;

  z = exp(x*y);
  EXPECT_NEAR(z.val(), std::exp(x.val()*y.val()), eps);
  EXPECT_NEAR(z.dval(0), 2.*std::exp(x.val()*y.val()), eps);
  EXPECT_NEAR(z.dval(1), std::exp(x.val()*y.val()), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);

  z = exp2(x+y);
  EXPECT_NEAR(z.val(), std::exp2(x.val()+y.val()), eps);
  EXPECT_NEAR(z.dval(0), std::log(2.)*std::exp2(x.val()+y.val()), eps);
  EXPECT_NEAR(z.dval(1), std::log(2.)*std::exp2(x.val()+y.val()), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);
}

TYPED_TEST(EvaluationTest, LogChainRule) {
  using FadType = typename TestFixture::StaticFad;
  static constexpr auto eps = TestFixture::eps;
  FadType x{3, 2., 0, 1.};
  FadType y{3, 3., 1, 1.};
  FadType z;

  z = log(x*y);
  EXPECT_NEAR(z.val(), std::log(6.l), eps);
  EXPECT_NEAR(z.dval(0), 1.l/2, eps);
  EXPECT_NEAR(z.dval(1), 1.l/3, eps);
  EXPECT_NEAR(z.dval(2), 0., eps);

  z = log2(x+y);
  EXPECT_NEAR(z.val(), std::log2(5.l), eps);
  EXPECT_NEAR(z.dval(0), 1./(5*std::log(2.l)), eps);
  EXPECT_NEAR(z.dval(1), 1./(5*std::log(2.l)), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);

  z = log10(y-x);
  EXPECT_NEAR(z.val(), std::log10(1.l), eps);
  EXPECT_NEAR(z.dval(0), -1./(std::log(10.l)), eps);
  EXPECT_NEAR(z.dval(1), 1./(std::log(10.l)), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);
}

TYPED_TEST(EvaluationTest, RootChainRule) {
  using FadType = typename TestFixture::StaticFad;
  static constexpr auto eps = TestFixture::eps;

  FadType x{3, 4., 0, 1.};
  FadType y{3, 9., 1, 1.};
  FadType z;

  z = sqrt(x*y);
  EXPECT_NEAR(z.val(), 6., eps);
  EXPECT_NEAR(z.dval(0), 3.l/4, eps);
  EXPECT_NEAR(z.dval(1), 1.l/3, eps);
  EXPECT_NEAR(z.dval(2), 0., eps);

  x = FadType{3, 8., 0, 1};
  y = FadType{3, 27., 1, 1};
  z = cbrt(x*y);
  EXPECT_NEAR(z.val(), 6., eps);
  EXPECT_NEAR(z.dval(0), 1.l/4, eps);
  EXPECT_NEAR(z.dval(1), 2.l/27, eps);
  EXPECT_NEAR(z.dval(2), 0., eps);
}

TYPED_TEST(EvaluationTest, TrigChainRule) {
  using FadType = typename TestFixture::StaticFad;
  static constexpr auto eps = TestFixture::eps;

  FadType x{3, 1., 0, 1.};
  FadType y{3, 2., 1, 1.};
  FadType z;

  z = cos(x*y);
  EXPECT_NEAR(z.val(), std::cos(2.l), eps);
  EXPECT_NEAR(z.dval(0), -2*std::sin(2.l), eps);
  EXPECT_NEAR(z.dval(1), -std::sin(2.l), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);

  z = sin(x*y);
  EXPECT_NEAR(z.val(), std::sin(2.l), eps);
  EXPECT_NEAR(z.dval(0), 2*std::cos(2.l), eps);
  EXPECT_NEAR(z.dval(1), std::cos(2.l), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);

  z = tan(x*y);
  EXPECT_NEAR(z.val(), std::tan(2.l), eps);
  EXPECT_NEAR(z.dval(0), 2/(std::cos(2.l)*std::cos(2.l)), eps);
  EXPECT_NEAR(z.dval(1), 1/(std::cos(2.l)*std::cos(2.l)), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);
}

TYPED_TEST(EvaluationTest, InverseTrigChainRule) {
  using FadType = typename TestFixture::StaticFad;
  static constexpr auto eps = TestFixture::eps;

  FadType x{3, 1./2, 0, 1.};
  FadType y{3, -1./2, 1, 1.};
  FadType z;

  z = acos(x*y);
  EXPECT_NEAR(z.val(), std::acos(-0.25l), eps);
  EXPECT_NEAR(z.dval(0), 1/(2*std::sqrt(0.9375l)), eps);
  EXPECT_NEAR(z.dval(1), -1/(2*std::sqrt(0.9375l)), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);

  z = asin(x*y);
  EXPECT_NEAR(z.val(), std::asin(-0.25l), eps);
  EXPECT_NEAR(z.dval(0), -1/(2*std::sqrt(0.9375l)), eps);
  EXPECT_NEAR(z.dval(1), 1/(2*std::sqrt(0.9375l)), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);

  z = atan(x*y);
  EXPECT_NEAR(z.val(), std::atan(-0.25l), eps);
  EXPECT_NEAR(z.dval(0), -1/2.125l, eps);
  EXPECT_NEAR(z.dval(1), 1/2.125l, eps);
  EXPECT_NEAR(z.dval(2), 0., eps);
}

TYPED_TEST(EvaluationTest, HyperbolicTrigChainRule) {
  using FadType = typename TestFixture::StaticFad;
  static constexpr auto eps = TestFixture::eps;

  FadType x{3, 1., 0, 1.};
  FadType y{3, 2., 1, 1.};
  FadType z;

  z = cosh(x*y);
  EXPECT_NEAR(z.val(), std::cosh(2.l), eps);
  EXPECT_NEAR(z.dval(0), 2*std::sinh(2.l), eps);
  EXPECT_NEAR(z.dval(1), std::sinh(2.l), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);

  z = sinh(x*y);
  EXPECT_NEAR(z.val(), std::sinh(2.l), eps);
  EXPECT_NEAR(z.dval(0), 2*std::cosh(2.l), eps);
  EXPECT_NEAR(z.dval(1), std::cosh(2.l), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);

  z = tanh(x*y);
  EXPECT_NEAR(z.val(), std::tanh(2.l), eps);
  EXPECT_NEAR(z.dval(0), 2/(std::cosh(2.l)*std::cosh(2.l)), eps);
  EXPECT_NEAR(z.dval(1), 1/(std::cosh(2.l)*std::cosh(2.l)), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);
}

TYPED_TEST(EvaluationTest, InverseHyperbolicTrigChainRule) {
  using FadType = typename TestFixture::StaticFad;
  static constexpr auto eps = TestFixture::eps;

  FadType x{3, 1., 0, 1.};
  FadType y{3, 2., 1, 1.};
  FadType z;

  z = acosh(x*y);
  EXPECT_NEAR(z.val(), std::acosh(2.l), eps);
  EXPECT_NEAR(z.dval(0), 2/(std::sqrt(3.l)), eps);
  EXPECT_NEAR(z.dval(1), 1/(std::sqrt(3.l)), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);

  z = asinh(x*y);
  EXPECT_NEAR(z.val(), std::asinh(2.l), eps);
  EXPECT_NEAR(z.dval(0), 2/(std::sqrt(5.l)), eps);
  EXPECT_NEAR(z.dval(1), 1/(std::sqrt(5.l)), eps);
  EXPECT_NEAR(z.dval(2), 0., eps);

  y.valNonConst() = 0.5;
  z = atanh(x*y);
  EXPECT_NEAR(z.val(), std::atanh(0.5l), eps);
  EXPECT_NEAR(z.dval(0), 2.l/3, eps);
  EXPECT_NEAR(z.dval(1), 4.l/3, eps);
  EXPECT_NEAR(z.dval(2), 0., eps);
}

TYPED_TEST(EvaluationTest, HEquation) {
  using FadType = typename TestFixture::DynamicFad;
  static constexpr auto eps = TestFixture::eps;

  static constexpr std::size_t N = 20;
  static constexpr TypeParam c = 0.5;
  std::vector<FadType> x(N);
  std::vector<TypeParam> mu(N);
  for (std::size_t i = 0; i < N; ++i) {
    x[i] = FadType{N, 1., i, 1.};
    mu[i] = 1./(2*N) + static_cast<TypeParam>(i)/N;
  }

  std::vector<FadType> f(x);
  std::vector<TypeParam> coefs(N, 1.);
  for (std::size_t i = 0; i < N; ++i) {
    FadType delta{N, 1.};
    for (std::size_t j = 0; j < N; ++j) {
      auto tmp = c * mu[i] * x[j]/(2 * N * (mu[i] + mu[j]));
      delta -= tmp;
      coefs[i] -= tmp.val();
    }
    f[i] -= 1/delta;
  }

  auto identity = [] (const std::size_t i, const std::size_t j) {
    if (i == j) return static_cast<TypeParam>(1);
    else return static_cast<TypeParam>(0);
  };

  for (std::size_t i = 0; i < N; ++i) {
    EXPECT_NEAR(f[i].val(), 1. - 1./coefs[i], eps);
    for (std::size_t j = 0; j < N; ++j) {
      auto expected = identity(i, j) - c * mu[i]/(2 * N * (mu[i] + mu[j]) * coefs[i] * coefs[i]);
      EXPECT_NEAR(f[i].dval(j), expected, eps);
    }
  }
}
