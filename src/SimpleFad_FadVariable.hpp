#ifndef SIMPLEFAD_FADVARIABLE_HPP
#define SIMPLEFAD_FADVARIABLE_HPP

#include "SimpleFad_config.hpp"
#include "SimpleFad_ExprBase.hpp"
#include "SimpleFad_ExprTraits.hpp"
#include "SimpleFad_DerivStorageTraits.hpp"

#include <cstddef>
#include <initializer_list>

namespace SimpleFad
{

// Class template representing a variable forward-mode automatic differentiation type.
template<class T, class DerivStorageType>
class FadVariable: public ExprBase<FadVariable<T, DerivStorageType>>
{
public:
  using storage_traits_type = DerivStorageTraits<DerivStorageType>;
  using value_type = typename storage_traits_type::value_type;
  using deriv_storage_type = typename storage_traits_type::container_type;
  using index_type = typename storage_traits_type::index_type;

  // Default constructor. Value initializes value and storage, so value is initialized
  // to 0. Statically sized derivative container is initialized to array of 0's, and
  // dynamically sized derivative container is initialized with size 0.
  constexpr FadVariable() = default;

  // Construct a FadVariable with value set to v and zero-initialized derivatives
  constexpr FadVariable(const index_type n, const T& v = T())
      : val_(v), dval_(storage_traits_type::create(n)) {
  }

  // Construct a FadVariable with value set to v and derivative component id value set to dv
  constexpr FadVariable(const index_type n, const T& v, const index_type id, const T& dv)
      : val_(v), dval_(storage_traits_type::create(n, id, dv)) {
  }

  // Construct a FadVariable with value set to v and derivative values given by initializer_list
  constexpr FadVariable(const T& v, const std::initializer_list<T>& dv)
      : val_(v), dval_(storage_traits_type::create(dv)) {
  }

  // Construct a FadVariable initialized with value/derivative from another expression
  template<class Expr>
  constexpr FadVariable(const ExprBase<Expr>& other)
      : val_(other.val()), dval_(storage_traits_type::create(other.dsize())) {
    for (index_type id = 0; id < other.dsize(); ++id) {
      dvalNonConst(id) = other.dval(id);
    }
  }

  // Assign from another expression
  template<class Expr>
  constexpr FadVariable& operator=(const ExprBase<Expr>& other) {
    val_ = other.val();
    storage_traits_type::reinit(dval_, other.dsize());
    for (index_type id = 0; id < other.dsize(); ++id) {
      dvalNonConst(id) = other.dval(id);
    }
    return *this;
  }

  // Update assignment operators
  constexpr FadVariable<T, DerivStorageType>&
  operator+=(const T& update) {
    val_ += update;
    return *this;
  }

  template<class Expr>
  constexpr FadVariable<T, DerivStorageType>&
  operator+=(const ExprBase<Expr>& update) {
    CHECK_CONFORMING_SIZE(this->dsize(), update.dsize())
    val_ += update.val();
    for (index_type id = 0; id < update.dsize(); ++id) {
      dvalNonConst(id) += update.dval(id);
    }
    return *this;
  }

  constexpr FadVariable<T, DerivStorageType>&
  operator-=(const T& update) {
    val_ -= update;
    return *this;
  }

  template<class Expr>
  constexpr FadVariable<T, DerivStorageType>&
  operator-=(const ExprBase<Expr>& update) {
    CHECK_CONFORMING_SIZE(this->dsize(), update.dsize())
    val_ -= update.val();
    for (index_type id = 0; id < update.dsize(); ++id) {
      dvalNonConst(id) -= update.dval(id);
    }
    return *this;
  }

  constexpr FadVariable<T, DerivStorageType>&
  operator*=(const T& update) {
    val_ *= update;
    for (index_type id = 0; id < this->dsize(); ++id) {
      dvalNonConst(id) *= update;
    }
    return *this;
  }

  template<class Expr>
  constexpr FadVariable<T, DerivStorageType>&
  operator*=(const ExprBase<Expr>& update) {
    CHECK_CONFORMING_SIZE(this->dsize(), update.dsize())
    const auto x = val_;
    const auto y = update.val();
    val_ *= y;
    for (index_type id = 0; id < update.dsize(); ++id) {
      dvalNonConst(id) = y * this->dval(id) + x * update.dval(id);
    }
    return *this;
  }

  constexpr FadVariable<T, DerivStorageType>&
  operator/=(const T& update) {
    val_ /= update;
    for (index_type id = 0; id < this->dsize(); ++id) {
      dvalNonConst(id) /= update;
    }
    return *this;
  }

  template<class Expr>
  constexpr FadVariable<T, DerivStorageType>&
  operator/=(const ExprBase<Expr>& update) {
    CHECK_CONFORMING_SIZE(this->dsize(), update.dsize())
    const auto x = val_;
    const auto y = update.val();
    val_ /= y;
    for (index_type id = 0; id < update.dsize(); ++id) {
      dvalNonConst(id) = this->dval(id) / y - x * update.dval(id) / (y * y);
    }
    return *this;
  }

  // Value accessors
  constexpr value_type& valNonConst() {
    return val_;
  }
  
  constexpr value_type valImpl() const {
    return val_;
  }

  // Derivative accessors
  constexpr value_type& dvalNonConst(index_type id) {
    return storage_traits_type::getNonConst(dval_, id);
  }

  constexpr value_type dvalImpl(index_type id) const {
    return storage_traits_type::get(dval_, id);
  }

  // Derivative array size
  constexpr index_type dsizeImpl() const {
    return storage_traits_type::derivSize(dval_);
  }

private:
  value_type val_;
  deriv_storage_type dval_;
};

// Alias for default statically sized FAD type.
template<class T, std::size_t N>
using StaticFadVariable = FadVariable<T, DefaultStaticStorageType<T, N>>;

//Alis for default dynamically sized FAD type.
template<class T>
using DynamicFadVariable = FadVariable<T, DefaultDynamicStorageType<T>>;

}

#endif
