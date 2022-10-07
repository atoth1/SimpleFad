#include "SimpleFad.hpp"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <vector>

// LAPACK declarations. May not be portable due to Fortran name mangling.
// If LAPACK library symbols don't match this, disable compilation of this
// file with the CMake option ENABLE_NEWTON_TESTS (enabled by default).
extern "C" {
  void sgesv_(int* n, int* nrhs,  float* a,  int* lda,
    int* ipivot, float* b, int* ldb, int* info);

  void dgesv_(int* n, int* nrhs,  double* a,  int* lda,
      int* ipivot, double* b, int* ldb, int* info);

  void sgtsv_(int* n, int* nrhs,  float* dl, float* d, float* du,
      float* b, int* ldb, int* info);

  void dgtsv_(int* n, int* nrhs,  double* dl, double* d, double* du,
      double* b, int* ldb, int* info);
}

namespace {
template<class T>
void solveDense(std::vector<T>& A, std::vector<T>& x);

template<>
void solveDense(std::vector<float>& A, std::vector<float>& x) {
  int n = x.size();
  int nrhs = 1;
  int info = 0;
  std::vector<int> ipivot(n);
  sgesv_(&n, &nrhs, A.data(), &n, ipivot.data(), x.data(), &n, &info);
}

template<>
void solveDense(std::vector<double>& A, std::vector<double>& x) {
  int n = x.size();
  int nrhs = 1;
  int info = 0;
  std::vector<int> ipivot(n);
  dgesv_(&n, &nrhs, A.data(), &n, ipivot.data(), x.data(), &n, &info);
}

template<class T>
void solveTridiag(std::vector<T>& lower, std::vector<T>& diag,
    std::vector<T>& upper, std::vector<T>& x);

template<>
void solveTridiag(std::vector<float>& lower, std::vector<float>& diag,
    std::vector<float>& upper, std::vector<float>& x) {
  int n = diag.size();
  int nrhs = 1;
  int info = 0;
  sgtsv_(&n, &nrhs, lower.data(), diag.data(), upper.data(), x.data(), &n, &info);
}

template<>
void solveTridiag(std::vector<double>& lower, std::vector<double>& diag,
    std::vector<double>& upper, std::vector<double>& x) {
  int n = diag.size();
  int nrhs = 1;
  int info = 0;
  dgtsv_(&n, &nrhs, lower.data(), diag.data(), upper.data(), x.data(), &n, &info);
}
}

template<class T>
class NewtonsMethodTest : public ::testing::Test {
};
using MyTypes = ::testing::Types<float, double>;
TYPED_TEST_SUITE(NewtonsMethodTest, MyTypes);

// Solving simple 2D nonlinear equation:
// F_0(x_0, x_1) = -x_0 + x_1^2 - k_0
// F_1(x_0, x_1) = c * (x_0^2 - x_1 - k_1)
// Constants chosen to give solution x_0 = 3, x_1 = 4.
TYPED_TEST(NewtonsMethodTest, Simple2D) {
  static constexpr std::size_t N = 2;
  using FadType = SimpleFad::StaticFadVariable<TypeParam, N>;
  const auto eps = std::sqrt(std::numeric_limits<TypeParam>::epsilon());

  std::vector<FadType> x{ FadType{N, 1., 0, 1.},  FadType{N, 1., 1, 1.} };
  std::vector<FadType> f(N);
  std::vector<TypeParam> jac(N*N);
  std::vector<TypeParam> update(N);

  static constexpr TypeParam k0 = 13;
  static constexpr TypeParam k1 = 5;
  static constexpr TypeParam c = 3;
  auto eval = [&] () {
    f[0] = -x[0] + x[1]*x[1] - k0;
    f[1] = c * (x[0]*x[0] - x[1] - k1);
  };

  auto solveLinSys = [&] () {
    // Jacobian expected in column-major format
    jac[0] = f[0].dval(0);
    jac[1] = f[1].dval(0);
    jac[2] = f[0].dval(1);
    jac[3] = f[1].dval(1);
    update[0] = -f[0].val();
    update[1] = -f[1].val();
    solveDense(jac, update);
  };

  auto updateSol = [&] () {
    x[0] += update[0];
    x[1] += update[1];
    return std::max(std::abs(update[0]), std::abs(update[1]));
  };

  TypeParam del = 1;
  int itCount = 0;
  static constexpr int maxIts = 10;
  while (itCount < maxIts && del > eps) {
    ++itCount;
    eval();
    solveLinSys();
    del = updateSol();
  }
  EXPECT_NEAR(x[0].val(), 3., eps);
  EXPECT_NEAR(x[1].val(), 4., eps);
}

// Solving  -u'' = sin(u) + 2 - sin(x(1-x)) 0 <= x <= 1, subject to u(0) = u(1) = 0
// using finite element method with linear elements. This has exact solution u(x) = x(1-x).
TYPED_TEST(NewtonsMethodTest, 1DFEM) {
  static constexpr std::size_t N = 5001;
  using FadType = SimpleFad::StaticFadVariable<TypeParam, 2>;
  const auto eps = std::sqrt(std::numeric_limits<TypeParam>::epsilon());

  std::vector<TypeParam> x(N);
  std::vector<TypeParam> f(N);
  std::vector<TypeParam> mesh(N);
  for (std::size_t i = 0; i < N; ++i) mesh[i] = static_cast<TypeParam>(i)/(N-1);
  std::vector<TypeParam> lower(N-1);
  std::vector<TypeParam> diag(N);
  std::vector<TypeParam> upper(N-1);
  std::vector<TypeParam> update(N);

  auto basisFn = [&] (const TypeParam pos, const int elem, const int basisFnId) {
    if (basisFnId == elem) {
      return (mesh[elem+1] - pos)/(mesh[elem+1] - mesh[elem]);
    } else if (basisFnId == elem+1) {
      return (pos - mesh[elem])/(mesh[elem+1] - mesh[elem]);
    } else {
      return static_cast<TypeParam>(0);
    }
  };

  auto dbasisFn = [&] (const int elem, const int basisFnId) {
    if (basisFnId == elem) {
      return -1/(mesh[elem+1] - mesh[elem]);
    } else if (basisFnId == elem+1) {
      return 1/(mesh[elem+1] - mesh[elem]);
    } else {
      return static_cast<TypeParam>(0);
    }
  };

  auto reconstruction = [&] (const TypeParam pos, const int elem,
      const FadType& l, const FadType& r) {
    return l*basisFn(pos, elem, elem) + r*basisFn(pos, elem, elem+1);
  };

  auto dreconstruction = [&] (const int elem, const FadType& l, const FadType& r) {
    return l*dbasisFn(elem, elem) + r*dbasisFn(elem, elem+1);
  };

  auto translate = [&] (const TypeParam pt, const TypeParam l, const TypeParam r) {
    return 0.5 * (r - l) * pt + 0.5 * (r + l);
  };

  auto contribution = [&] (const int elem, const int basisFnId) {
    using std::sin;
    FadType l{2, x[elem], 0, 1.};
    FadType r{2, x[elem+1], 1, 1.};

    static const std::array<TypeParam, 2>
      gaussPts{-1/std::sqrt(static_cast<TypeParam>(3)), 1/std::sqrt(static_cast<TypeParam>(3))};

    FadType f{2};
    for (const auto pt : gaussPts) {
      auto trans = translate(pt, mesh[elem], mesh[elem+1]);
      f += dreconstruction(elem, l, r) * dbasisFn(elem, basisFnId);
      f += sin(trans*(1-trans)) * basisFn(trans, elem, basisFnId);
      f -= sin(reconstruction(trans, elem, l, r)) * basisFn(trans, elem, basisFnId);
      f -= 2 * basisFn(trans, elem, basisFnId);
      f *= 0.5 * (mesh[elem+1] - mesh[elem]);
    }
    return f;
  };

  auto eval = [&] () {
    std::fill(f.begin(), f.end(), static_cast<TypeParam>(0));
    std::fill(diag.begin(), diag.end(), static_cast<TypeParam>(0));
    for (std::size_t elem = 0; elem < N-1; ++elem) {
      auto fi = contribution(elem, elem);
      f[elem] += fi.val();
      diag[elem] += fi.dval(0);
      upper[elem] = fi.dval(1);
      auto fip1 = contribution(elem, elem+1);
      f[elem+1] += fip1.val();
      lower[elem] = fip1.dval(0);
      diag[elem+1] += fip1.dval(1);
    }
    // Adjust for 0 Dirichlett BCs
    f[0] = x[0];
    diag[0] = 1.;
    upper[0] = 0.;
    f[N-1] = x[N-1];
    lower[N-2] = 0.;
    diag[N-1] = 1.;
  };

  auto solveLinSys = [&] () {
    std::transform(f.cbegin(), f.cend(), update.begin(),
        [] (const TypeParam resid) { return -resid; });
    solveTridiag(lower, diag, upper, update);
  };

  auto updateSol = [&] () {
    TypeParam del = 0;
    for (std::size_t row = 0; row < N; ++row) {
      x[row] += update[row];
      del = std::max(del, std::abs(update[row]));
    }
    return del;
  };

  TypeParam del = 1;
  int itCount = 0;
  static constexpr int maxIts = 10;
  while (itCount < maxIts && del > eps) {
    ++itCount;
    eval();
    solveLinSys();
    del = updateSol();
  }

  auto exact = [] (const TypeParam x) {
    return x * (1-x);
  };

  for (std::size_t i = 0; i < N; ++i) {
    EXPECT_NEAR(x[i], exact(mesh[i]), eps);
  }
}

// Solving Chandrasekhar H-equation with Newton's method.
TYPED_TEST(NewtonsMethodTest, HEquation) {
  using FadType = SimpleFad::DynamicFadVariable<TypeParam>;
  const auto eps = std::sqrt(std::numeric_limits<TypeParam>::epsilon());

  static constexpr std::size_t N = 50;
  static constexpr TypeParam c = 0.5;
  std::vector<FadType> x(N);
  std::vector<TypeParam> mu(N);
  for (std::size_t i = 0; i < N; ++i) {
    x[i] = FadType{N, 1., i, 1.};
    mu[i] = 1./(2*N) + static_cast<TypeParam>(i)/N;
  }
  std::vector<FadType> f(N, FadType{N});
  std::vector<TypeParam> jac(N*N);
  std::vector<TypeParam> update(N);

  auto eval = [&] () {
    for (std::size_t i = 0; i < N; ++i) {
      FadType delta{N, 1.};
      for (std::size_t j = 0; j < N; ++j) {
        auto tmp = c * mu[i] * x[j]/(2 * N * (mu[i] + mu[j]));
        delta -= tmp;
      }
      f[i] = x[i] - 1./delta;
    }
  };

  auto solveLinSys = [&] () {
    for (std::size_t row = 0; row < N; ++row) {
      update[row] = -f[row].val();
      // Jacobian expected in column-major format
      for (std::size_t col = 0; col < N; ++col) {
        jac[row + N*col] = f[row].dval(col);
      }
    }
    solveDense(jac, update);
  };

  auto updateSol = [&] () {
    TypeParam del = 0;
    for (std::size_t row = 0; row < N; ++row) {
      x[row] += update[row];
      del = std::max(del, std::abs(update[row]));
    }
    return del;
  };

  TypeParam del = 1;
  int itCount = 0;
  static constexpr int maxIts = 10;
  while (itCount < maxIts && del > eps) {
    eval();
    solveLinSys();
    del = updateSol();
  }
  // Different termination threshold for float and double.
  // float terminates in 3 and double in 4 on my machine.
  EXPECT_TRUE(itCount <= 4);
}
