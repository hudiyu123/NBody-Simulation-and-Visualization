#ifndef CPP_PROJECT_INCLUDE_NBODY_MATH_VECTOR_HPP
#define CPP_PROJECT_INCLUDE_NBODY_MATH_VECTOR_HPP

#include <NBody/Concept.hpp>
#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <iostream>

namespace NBody::Math {

/// \brief Vector.
/// \tparam Length The length of the vector
template<Arithmetic T, std::size_t Length>
class Vector {
 public:
  /// \brief The type of elements
  using ValueType = T;
  /// \brief The type of Array
  using ArrayType = std::array<ValueType, Length>;
  /// \brief The reference to the element
  using Reference = typename ArrayType::reference;
  /// \brief The const reference to the element
  using ConstReference = typename ArrayType::const_reference;
  /// \brief The type of iterator
  using Iterator = typename ArrayType::iterator;
  /// \brief The type of const iterator
  using ConstIterator = typename ArrayType::const_iterator;

 private:
  ArrayType arr_;

 public:
  /// \brief Constructs a vector with the value of elements.
  /// \param value The value of elements
  constexpr explicit Vector(ValueType value = ValueType{}) {
    std::ranges::fill(arr_, value);
  }

  /// \brief Constructs a vector with initializer list.
  /// \details Precondition: The length of initializer list equals the length
  /// of vector.
  /// \param list The initializer list
  constexpr Vector(std::initializer_list<T> list) {
    std::ranges::move(list, arr_.begin());
  }

  /// \brief Returns an iterator to the first element of the vector.
  /// \return Iterator to the first element.
  constexpr Iterator begin() noexcept {
    return arr_.begin();
  }

  /// \brief Returns an iterator to the first element of the vector.
  /// \return Iterator to the first element.
  constexpr ConstIterator begin() const noexcept {
    return arr_.begin();
  }

  /// \brief Returns an iterator to the element following the last element of
  /// the vector.
  /// \return Iterator to the element following the last element.
  constexpr Iterator end() noexcept {
    return arr_.end();
  }

  /// \brief Returns an iterator to the element following the last element of
  /// the vector.
  /// \return Iterator to the element following the last element.
  constexpr ConstIterator end() const noexcept {
    return arr_.end();
  }

  /// \brief Returns an iterator to the first element of the vector.
  /// \return Iterator to the first element.
  constexpr ConstIterator cbegin() noexcept {
    return arr_.cbegin();
  }

  /// \brief Returns an iterator to the first element of the vector.
  /// \return Iterator to the first element.
  constexpr ConstIterator cbegin() const noexcept {
    return arr_.cbegin();
  }

  /// \brief Returns an iterator to the element following the last element of
  /// the vector.
  /// \return Iterator to the element following the last element.
  constexpr ConstIterator cend() noexcept {
    return arr_.cend();
  }

  /// \brief Returns an iterator to the element following the last element of
  /// the vector.
  /// \return Iterator to the element following the last element.
  constexpr ConstIterator cend() const noexcept {
    return arr_.cend();
  }

  /// \brief Returns a reference to the element at location index.
  /// \param index Position of the element to return
  /// \return Reference to the requested element
  constexpr Reference operator[](std::size_t index) {
    return arr_[index];
  }

  /// \brief Returns a reference to the element at location index.
  /// \param index Position of the element to return
  /// \return Reference to the requested element
  constexpr ConstReference operator[](std::size_t index) const {
    return arr_[index];
  }

  /// \brief Adds other the vector.
  /// \param other Vector other
  /// \return Reference to the vector
  constexpr Vector& operator+=(const Vector& other) {
    std::ranges::transform(arr_, other.arr_, arr_.begin(), std::plus{});
    return *this;
  }

  /// \brief Subtracts other from the vector.
  /// \param other Vector other
  /// \return Reference to the vector
  constexpr Vector& operator-=(const Vector& other) {
    std::ranges::transform(arr_, other.arr_, arr_.begin(), std::minus{});
    return *this;
  }

  /// \brief Multiplies the vector by scalar other.
  /// \param other Scalar other
  /// \return Reference to the vector
  constexpr Vector& operator*=(ValueType other) {
    std::ranges::transform(arr_, arr_.begin(), [=](auto x) {
      return x * other;
    });
    return *this;
  }

  /// \brief Divides the vector by scalar other.
  /// \note The divisor should not be zero.
  /// \param other Scalar other
  /// \return Reference to the vector
  constexpr Vector& operator/=(ValueType other) {
    std::ranges::transform(arr_, arr_.begin(), [=](auto x) {
      return x / other;
    });
    return *this;
  }
};

/// \brief Returns the negation of the vector.
/// \tparam Length The length of the vector
/// \param vector Vector
/// \return Negation of the vector
template<Arithmetic T, std::size_t Length>
constexpr Vector<T, Length> operator-(const Vector<T, Length>& vector) {
  auto result = vector;
  std::ranges::transform(vector, result.begin(), std::negate{});
  return result;
}

/// \brief Returns the sum of two vectors.
/// \tparam Length The length of the vector
/// \param lhs Vector lhs
/// \param rhs Vector rhs
/// \return The result of addition
template<Arithmetic T, std::size_t Length>
constexpr Vector<T, Length> operator+(
    const Vector<T, Length>& lhs,
    const Vector<T, Length>& rhs
) {
  return Vector<T, Length>{lhs} += rhs;
}

/// \brief Returns the result of subtracting vector rhs from vector lhs.
/// \tparam Length The length of the vector
/// \param lhs Vector lhs
/// \param rhs Vector rhs
/// \return The result of subtraction
template<Arithmetic T, std::size_t Length>
constexpr Vector<T, Length> operator-(
    const Vector<T, Length>& lhs,
    const Vector<T, Length>& rhs
) {
  return Vector<T, Length>{lhs} -= rhs;
}

/// \brief Returns the product of vector lhs and scalar rhs.
/// \tparam Length The length of the vector
/// \param lhs Vector lhs
/// \param rhs Scalar rhs
/// \return The product.
template<Arithmetic T, std::size_t Length>
constexpr Vector<T, Length> operator*(const Vector<T, Length>& lhs, T rhs) {
  return Vector<T, Length>{lhs} *= rhs;
}

/// \brief Returns the product of scalar lhs and vector rhs.
/// \tparam Length The length of the vector
/// \param lhs Scalar lhs
/// \param rhs Vector rhs
/// \return The product.
template<Arithmetic T, std::size_t Length>
constexpr Vector<T, Length> operator*(T lhs, const Vector<T, Length>& rhs) {
  return Vector<T, Length>{rhs} *= lhs;
}

/// \brief Returns the quotient of vector rhs and scalar lhs.
/// \tparam Length The length of the vector
/// \param lhs Vector lhs
/// \param rhs Scalar rhs
/// \return The quotient
template<Arithmetic T, std::size_t Length>
constexpr Vector<T, Length> operator/(const Vector<T, Length>& lhs, T rhs) {
  return Vector<T, Length>{lhs} /= rhs;
}

/// \brief Writes the vector to a character output stream.
/// \details Elements are separated by white space.
/// \tparam Length The length of the vector
/// \param ostream The character output stream
/// \param vector The vector to be extracted
/// \return The ostream
template<Arithmetic T, std::size_t Length>
std::ostream& operator<<(
    std::ostream& ostream,
    const Vector<T, Length>& vector
) {
  for (auto it = vector.cbegin(); it != vector.cend(); ++it) {
    ostream << *it;
    if (std::next(it) != vector.cend()) {
      ostream << " ";
    }
  }
  return ostream;
}

/// \brief Reads a vector from a character input stream.
/// \tparam Length The length of the vector
/// \param istream The character input stream
/// \param vector The vector to be inserted
/// \return The istream
template<Arithmetic T, std::size_t Length>
std::istream& operator>>(std::istream& istream, Vector<T, Length>& vector) {
  for (auto& value: vector) {
    if (!(istream >> value)) { break; }
  }
  return istream;
}

} // namespace NBody::Math

#endif // CPP_PROJECT_INCLUDE_NBODY_MATH_VECTOR_HPP
