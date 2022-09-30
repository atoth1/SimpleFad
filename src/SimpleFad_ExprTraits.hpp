#ifndef SIMPLEFAD_EXPRTRAITS_HPP
#define SIMPLEFAD_EXPRTRAITS_HPP

#include "SimpleFad_DerivStorageTraits.hpp"

#include <cstddef>
#include <type_traits>

namespace SimpleFad {

template<class Expr>
struct ExprTraits;

// Traits for FadLiteral
template<class T>
class FadLiteral;

template<class T>
struct ExprTraits<FadLiteral<T>> {
  using value_type = T;
  using index_type = std::size_t;
};

// Traits for FadVariable
template<class T, class DerivStorageType>
class FadVariable;

template<class T, std::size_t N>
struct ExprTraits<FadVariable<T, DefaultStaticStorageType<T, N>>> {
  using deriv_storage_type = DefaultStaticStorageType<T, N>;
  using deriv_storage_traits = DerivStorageTraits<deriv_storage_type>;
  using value_type = typename deriv_storage_traits::value_type;
  using index_type = typename deriv_storage_traits::index_type;
};

template<class T>
struct ExprTraits<FadVariable<T, DefaultDynamicStorageType<T>>> {
  using deriv_storage_type = DefaultDynamicStorageType<T>;
  using deriv_storage_traits = DerivStorageTraits<deriv_storage_type>;
  using value_type = typename deriv_storage_traits::value_type;
  using index_type = typename deriv_storage_traits::index_type;
};

// Traits for UnaryExpr
template<class UnaryOp, class Expr>
class UnaryExpr;

template<class UnaryOp, class Expr>
struct ExprTraits<UnaryExpr<UnaryOp, Expr>> {
  using expr_traits = ExprTraits<Expr>;
  using value_type = typename expr_traits::value_type;
  using index_type = typename expr_traits::index_type;
};

// Traits for BinaryExpr
template<class BinaryOp, class Expr1, class Expr2>
class BinaryExpr;

template<class BinaryOp, class Expr1, class Expr2>
struct ExprTraits<BinaryExpr<BinaryOp, Expr1, Expr2>> {
  using expr1_traits = ExprTraits<Expr1>;
  using expr2_traits = ExprTraits<Expr2>;
  using value_type = std::common_type_t<typename expr1_traits::value_type,
      typename expr2_traits::value_type>;
  using index_type = std::common_type_t<typename expr1_traits::index_type,
      typename expr2_traits::index_type>;
};

}

#endif
