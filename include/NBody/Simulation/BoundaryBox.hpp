#ifndef CPP_PROJECT_INCLUDE_NBODY_SIMULATION_BOUNDARYBOX_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_SIMULATION_BOUNDARYBOX_HPP_

#include <NBody/Math/Vector.hpp>
#include <concepts>

namespace NBody::Simulation {

/// \brief The type of the boundary plane.
/// \details Six enums represent the six faces of a cube.
enum class BoundaryPlaneType {
  NegativeX,
  PositiveX,
  NegativeY,
  PositiveY,
  NegativeZ,
  PositiveZ
};

/// \brief Represents the boundary cube.
template<std::floating_point T>
struct BoundaryBox {
  using ValueType = T;
  using VectorType = Math::Vector<ValueType, 3>;

  /// \brief Low positions on x, y and z axis of the vertices of the cube.
  VectorType low;
  /// \brief High positions on x, y and z axis of the vertices of the cube.
  VectorType high;

  /// Returns whether the boundary box is valid
  /// \return true if the boundary box is valid, false otherwise
  [[nodiscard]] bool isValid() const {
    return low[0] < high[0] && low[1] < high[1] && low[2] < high[2];
  }

  /// \brief Returns the maximum length of edges.
  /// \return Maximum length of edges
  [[nodiscard]] ValueType maxLength() const {
    auto diff = high - low;
    return std::max({diff[0], diff[1], diff[2]});
  }

  /// \brief Returns the center of the boundary box.
  [[nodiscard]] VectorType center() const {
    return low + (high - low) / ValueType{2};
  }
};

/// \brief Writes the boundary box to a character output stream.
/// \param os The character output stream
/// \param boundaryBox The boundary box
/// \return The character output stream
template<std::floating_point T>
std::ostream& operator<<(std::ostream& os, const BoundaryBox<T>& boundaryBox) {
  os << boundaryBox.low << " " << boundaryBox.high;
  return os;
}

/// \brief Reads the boundary box from a character input stream.
/// \param is The character input stream
/// \param boundaryBox The boundary box
/// \return The character input stream
template<std::floating_point T>
std::istream& operator>>(std::istream& is, BoundaryBox<T>& boundaryBox) {
  if (!(is >> boundaryBox.low)) { return is; }
  if (!(is >> boundaryBox.high)) { return is; }
  return is;
}

} // namespace NBody::Simulation

#endif // CPP_PROJECT_INCLUDE_NBODY_SIMULATION_BOUNDARYBOX_HPP_
