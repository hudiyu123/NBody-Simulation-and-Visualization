#ifndef CPP_PROJECT_INCLUDE_NBODY_MATH_COMMON_HPP
#define CPP_PROJECT_INCLUDE_NBODY_MATH_COMMON_HPP

#include <NBody/Concept.hpp>
#include <NBody/Math/Vector.hpp>
#include <NBody/Math/Matrix.hpp>
#include <numbers>
#include <cmath>
#include <concepts>

namespace NBody::Math {

/// \brief Converts degrees to radians and returns the result.
/// \param degrees Degrees to be converted
/// \return Radians
template<Arithmetic T>
constexpr T radians(T degrees) {
  return degrees * std::numbers::pi_v<T> / static_cast<T>(180);
}

/// \brief Returns the length of vector x.
/// \tparam Length Length of the vector
/// \param x Vector
/// \return Length
template<Arithmetic T, std::size_t Length>
constexpr T length(const Vector<T, Length>& x) {
  T result = 0;
  for (auto it = x.cbegin(); it != x.cend(); ++it) {
    result += *it * *it;
  }
  return std::sqrt(result);
}

/// \brief Returns the length of vector x with softening.
/// \tparam Length Length of the vector
/// \param x Vector
/// \param softening Softening
/// \return Length
template<Arithmetic T, std::size_t Length>
constexpr T length(const Vector<T, Length>& x, T softening) {
  T result = 0;
  for (auto it = x.cbegin(); it != x.cend(); ++it) {
    result += *it * *it;
  }
  result += softening * softening;
  return std::sqrt(result);
}

/// \brief Returns the dot product of vector x and vector y.
/// \tparam Length Length of the vector
/// \param x Vector x
/// \param y Vector y
/// \return Dot product
template<Arithmetic T, std::size_t Length>
constexpr T dot(const Vector<T, Length>& x, const Vector<T, Length>& y) {
  T result = 0;
  auto it_x = x.cbegin();
  auto it_y = y.cbegin();
  for (; it_x != x.cend() && it_y != y.cend(); ++it_x, ++it_y) {
    result += *it_x * *it_y;
  }
  return result;
}

/// \brief Returns the inverse sqrt of v.
/// \param v Floating-point scalar v
/// \return Inverse sqrt
template<Arithmetic T>
T inverseSqrt(T v) {
  return static_cast<T>(1) / std::sqrt(v);
}

/// \brief Returns a vector in the same direction as x but with length of 1.
/// \tparam Length Length of the vector
/// \param x Vector x
/// \return Normalized vector
template<Arithmetic T, std::size_t Length>
Vector<T, Length> normalize(const Vector<T, Length>& x) {
  return x * inverseSqrt(dot(x, x));
}

/// \brief Returns the cross product of vector x and vector y.
/// \param x Vector x
/// \param y Vector y
/// \return Cross product
template<Arithmetic T>
constexpr Vector<T, 3> cross(const Vector<T, 3>& x, const Vector<T, 3>& y) {
  Vector<T, 3> result;
  result[0] = x[1] * y[2] - y[1] * x[2];
  result[1] = x[2] * y[0] - y[2] * x[0];
  result[2] = x[0] * y[1] - y[0] * x[1];
  return result;
}

/// \brief Build a look at view matrix
/// \param eye Position of the camera
/// \param center Position where the camera is looking at
/// \param up Normalized up vector, how the camera is oriented
/// \return Loot at view matrix
template<Arithmetic T>
Matrix<T, 4, 4> lookAt(
    const Vector<T, 3>& eye,
    const Vector<T, 3>& center,
    const Vector<T, 3>& up
) {
  const auto f = normalize(center - eye);
  const auto s = normalize(cross(f, up));
  const auto u = cross(s, f);

  Matrix<T, 4, 4> result;
  result[0][0] = s[0];
  result[0][1] = s[1];
  result[0][2] = s[2];
  result[1][0] = u[0];
  result[1][1] = u[1];
  result[1][2] = u[2];
  result[2][0] = -f[0];
  result[2][1] = -f[1];
  result[2][2] = -f[2];
  result[0][3] = -dot(s, eye);
  result[1][3] = -dot(u, eye);
  result[2][3] = dot(f, eye);
  result[3][3] = T{1};
  return result;
}

/// \brief Creates a matrix for a symmetric perspective-view frustum with far
/// plane at infinite.
/// \param fovy Angle between the upper and lower sides of the viewing frustum
/// \param aspect Aspect ratio of the viewing window. (width / height)
/// \param zNear Distance to the near clipping plane along the -Z axis.
/// \return Perspective transformation matrix.
template<Arithmetic T>
Matrix<T, 4, 4> infinitePerspective(
    const T& fovy,
    const T& aspect,
    const T& zNear
) {
  auto range = std::tan(fovy / T{2}) * zNear;
  auto left = -range * aspect;
  auto right = range * aspect;
  auto bottom = -range;
  auto top = range;

  Matrix<T, 4, 4> result;
  result[0][0] = (T{2} * zNear) / (right - left);
  result[1][1] = (T{2} * zNear) / (top - bottom);
  result[2][2] = -T{1};
  result[3][2] = -T{1};
  result[2][3] = -T{2} * zNear;
  return result;
}

} // namespace NBody::Math

#endif // CPP_PROJECT_INCLUDE_NBODY_MATH_COMMON_HPP
