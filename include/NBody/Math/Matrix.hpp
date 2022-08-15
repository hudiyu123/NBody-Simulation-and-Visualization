#ifndef CPP_PROJECT_INCLUDE_NBODY_MATH_MATRIX_HPP
#define CPP_PROJECT_INCLUDE_NBODY_MATH_MATRIX_HPP

#include <NBody/Concept.hpp>
#include <NBody/Math/Vector.hpp>

namespace NBody::Math {

/// \brief Matrix.
/// \tparam Row The number of rows
/// \tparam Column The number of columns
template<Arithmetic T, std::size_t Row, std::size_t Column>
class Matrix {
 public:
  /// \brief The type of elements
  using ValueType = T;
  /// \brief The type of rows
  using RowType = Vector<T, Column>;
  /// \brief The type of array
  using ArrayType = std::array<RowType, Row>;
  /// \brief The reference to the row
  using Reference = typename ArrayType::reference;
  /// \brief The const reference to the row
  using ConstReference = typename ArrayType::const_reference ;
  /// \brief The type of iterator
  using Iterator = typename ArrayType::iterator;
  /// \brief The type of const iterator
  using ConstIterator = typename ArrayType::const_iterator;

 private:
  ArrayType arr_;

 public:
  /// \brief Default constructor.
  constexpr Matrix() {
    std::ranges::fill(arr_, RowType());
  }

  /// \brief Constructs the matrix with the value of elements at diagonal.
  /// \param value The value of elements at diagonal
  constexpr explicit Matrix(ValueType value) {
    for (auto it = begin(); it != end(); ++it) {
      const auto row_index = it - cbegin();
      for (auto vector_it = it->begin(); vector_it != it->end(); ++vector_it) {
        const auto element_index = vector_it - it->begin();
        *vector_it = row_index == element_index ? value : ValueType{};
      }
    }
  }

  /// \brief Returns an iterator to the first row of the matrix.
  /// \return Iterator to the first row.
  constexpr Iterator begin() noexcept {
    return arr_.begin();
  }

  /// \brief Returns an iterator to the first row of the matrix.
  /// \return Iterator to the first row.
  constexpr ConstIterator begin() const noexcept {
    return arr_.begin();
  }

  /// \brief Returns an iterator to the row following the last row of the
  /// matrix.
  /// \return Iterator to the row following the last row.
  constexpr Iterator end() noexcept {
    return arr_.end();
  }

  /// \brief Returns an iterator to the row following the last row of the
  /// matrix.
  /// \return Iterator to the row following the last row.
  constexpr ConstIterator end() const noexcept {
    return arr_.end();
  }

  /// \brief Returns an iterator to the first row of the matrix.
  /// \return Iterator to the first row.
  constexpr ConstIterator cbegin() noexcept {
    return arr_.cbegin();
  }

  /// \brief Returns an iterator to the first row of the matrix.
  /// \return Iterator to the first row.
  constexpr ConstIterator cbegin() const noexcept {
    return arr_.cbegin();
  }

  /// \brief Returns an iterator to the row following the last row of the
  /// matrix.
  /// \return Iterator to the row following the last row.
  constexpr ConstIterator cend() noexcept {
    return arr_.cend();
  }

  /// \brief Returns an iterator to the row following the last row of the
  /// matrix.
  /// \return Iterator to the row following the last row.
  constexpr ConstIterator cend() const noexcept {
    return arr_.cend();
  }

  /// \brief Returns a reference to the row at location index.
  /// \param index Position of the row
  /// \return Reference to the requested row
  constexpr Reference operator[](std::size_t index) {
    return arr_[index];
  }

  /// \brief Returns a reference to the row at location index.
  /// \param index Position of the row to return
  /// \return Reference to the requested row
  constexpr ConstReference operator[](std::size_t index) const {
    return arr_[index];
  }
};

} // namespace NBody::Math

#endif //CPP_PROJECT_INCLUDE_NBODY_MATH_MATRIX_HPP
