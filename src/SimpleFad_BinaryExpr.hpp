#ifndef SIMPLEFAD_BINARYEXPR_HPP
#define SIMPLEFAD_BINARYEXPR_HPP

#include "SimpleFad_config.hpp"
#include "SimpleFad_Exception.hpp"
#include "SimpleFad_ExprBase.hpp"

namespace SimpleFad {

template<class BinaryOp, class Expr1, class Expr2>
class BinaryExpr: public ExprBase<BinaryExpr<BinaryOp, Expr1, Expr2>> {
public:
  using traits_type = ExprTraits<BinaryExpr<BinaryOp, Expr1, Expr2>>;
  using value_type = typename traits_type::value_type;
  using index_type = typename traits_type::index_type;

  constexpr explicit BinaryExpr(const BinaryOp& op, const Expr1& expr1,
      const Expr2& expr2)
      : op_(op)
      , expr1_(expr1)
      , expr2_(expr2)
      , expr1Val_(expr1_.val())
      , expr2Val_(expr2_.val()) {
    if (expr1_.dsize() != 0 && expr2_.dsize() != 0) {
      CHECK_CONFORMING_SIZE(expr1_.dsize(), expr2_.dsize())
    }
  }

  SIMPLEFAD_CONSTEXPR value_type valImpl() const {
    return op_.f(expr1Val_, expr2Val_);
  }

  SIMPLEFAD_CONSTEXPR value_type dvalImpl(index_type id) const {
    value_type ret = static_cast<value_type>(0.);
    if (expr1_.dsize()) ret += op_.dfx(expr1Val_, expr2Val_) * expr1_.dval(id);
    if (expr2_.dsize()) ret += op_.dfy(expr1Val_, expr2Val_) * expr2_.dval(id);
    return ret;
  }

  constexpr index_type dsizeImpl() const {
    return expr1_.dsize() ? expr1_.dsize() : expr2_.dsize();
  }

private:
  BinaryOp op_;
  Expr1 expr1_;
  Expr2 expr2_;
  value_type expr1Val_;
  value_type expr2Val_;
};

}

#endif
