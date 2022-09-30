#ifndef SIMPLEFAD_DERIVSTORAGETRAITS_HPP
#define SIMPLEFAD_DERIVSTORAGETRAITS_HPP

#include "SimpleFad_Exception.hpp"

#include <array>
#include <initializer_list>
#include <iterator>
#include <vector>

template<class StorageType>
struct DerivStorageTraits;

template<class Scalar, std::size_t N>
using DefaultStaticStorageType = std::array<Scalar, N>;

template<class Scalar, std::size_t N>
struct DerivStorageTraits<DefaultStaticStorageType<Scalar, N>> {
  using container_type = DefaultStaticStorageType<Scalar, N>;
  using value_type = typename container_type::value_type;
  using index_type = typename container_type::size_type;
  static constexpr index_type derivSize(const container_type&) {
    return N;
  }

  static constexpr const value_type& get(const container_type& container,
      const index_type index) {
    CHECK_WITHIN_BOUNDS(container.size(), index)
    return container[index];
  }

  static constexpr value_type& getNonConst(container_type& container,
      const index_type index) {
    CHECK_WITHIN_BOUNDS(container.size(), index)
    return container[index];
  }

  static constexpr container_type create(const index_type size) {
    CHECK_CONFORMING_SIZE(N, size)
    return container_type{ };
  }

  static constexpr container_type create(const index_type size, const index_type id,
      const value_type val) {
    CHECK_CONFORMING_SIZE(N, size)
    container_type tmp{ };
    tmp[id] = val;
    return tmp;
  }

  static constexpr container_type create(
      const std::initializer_list<value_type>& init) {
    CHECK_CONFORMING_SIZE(N, init.size())
    container_type tmp{ };
    for (index_type i = 0; i < init.size(); ++i)
      tmp[i] = *std::next(std::begin(init), i);
    return tmp;
  }

  static constexpr void reinit(container_type& container, const index_type size) {
    CHECK_CONFORMING_SIZE(N, size)
    for (index_type i = 0; i < N; ++i) container[i] = static_cast<value_type>(0.);
  }
};

template<class Scalar>
using DefaultDynamicStorageType = std::vector<Scalar>;

template<class Scalar>
struct DerivStorageTraits<DefaultDynamicStorageType<Scalar>> {
  using container_type = DefaultDynamicStorageType<Scalar>;
  using value_type = typename container_type::value_type;
  using index_type = typename container_type::size_type;
  static index_type derivSize(const container_type& container) {
    return container.size();
  }

  static const value_type& get(const container_type& container,
      const index_type index) {
    CHECK_WITHIN_BOUNDS(container.size(), index)
    return container[index];
  }

  static value_type& getNonConst(container_type& container,
      const index_type index) {
    CHECK_WITHIN_BOUNDS(container.size(), index)
    return container[index];
  }

  static container_type create(const index_type size) {
    return container_type(size);
  }

  static container_type create(const index_type size, const index_type id,
      const value_type val) {
    container_type tmp(size);
    tmp[id] = val;
    return tmp;
  }

  static container_type create(const std::initializer_list<value_type>& init) {
    return container_type(init);
  }

  static void reinit(container_type& container, const index_type size) {
    container.resize(size);
  }
};

#endif
