#ifndef SIMPLEFAD_UNARYEXPR_HPP
#define SIMPLEFAD_UNARYEXPR_HPP

#include "SimpleFad_config.hpp"
#include "SimpleFad_ExprBase.hpp"

namespace SimpleFad {

template<class UnaryOp, class Expr>
class UnaryExpr: public ExprBase<UnaryExpr<UnaryOp, Expr>> {
public:
  using traits_type = ExprTraits<UnaryExpr<UnaryOp, Expr>>;
  using value_type = typename traits_type::value_type;
  using index_type = typename traits_type::index_type;

  constexpr explicit UnaryExpr(const UnaryOp& op, const Expr& expr)
      : op_(op), expr_(expr) {
  }

  SIMPLEFAD_CONSTEXPR value_type valImpl() const {
    return op_.f(expr_.val());
  }

  SIMPLEFAD_CONSTEXPR value_type dvalImpl(index_type id) const {
    return op_.df(expr_.val()) * expr_.dval(id);
  }

  constexpr index_type dsizeImpl() const {
    return expr_.dsize();
  }

private:
  UnaryOp op_;
  Expr expr_;
};

}

#endif
