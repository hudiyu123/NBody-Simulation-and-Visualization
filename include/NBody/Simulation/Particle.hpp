#ifndef CPP_PROJECT_INCLUDE_NBODY_PARTICLE_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_PARTICLE_HPP_

#include <NBody/Math/Vector.hpp>
#include <NBody/Math/Common.hpp>
#include <NBody/Simulation/BoundaryBox.hpp>
#include <concepts>
#include <vector>
#include <iostream>

namespace NBody::Simulation {

/// \brief An instance of data type Particle represents a particle in nbody
/// simulation.
template<std::floating_point T>
struct Particle {
  using ValueType = T;
  using VectorType = Math::Vector<ValueType, 3>;

  /// \brief The index of the particle.
  /// \details The index should be unique in a simulation process.
  std::size_t index;
  /// \brief The mass of the particle.
  ValueType mass;
  /// \brief The radius of the particle.
  /// \details The particles in the simulation are considered as sphere.
  ValueType radius;
  /// \brief The position of the particle in the simulation space.
  VectorType position;
  /// \brief The velocity of the particle.
  VectorType velocity;
  /// \brief The acceleration of the particle.
  VectorType acceleration;
  /// \brief The RGB color of the particle.
  VectorType color;

  /// Checks whether the particle is collided with the other particle.
  /// \param other The other particle.
  /// \return true if the two particles are collided, false otherwise
  bool isCollidedWith(const Particle& other) const {
    if (this == &other) {
      return false;
    }

    const auto r = radius + other.radius;
    const auto d = position - other.position;
    const auto dv = velocity - other.velocity;
    if (Math::length(d) > r) {
      return false;
    } else {
      return Math::dot(d, dv) <= ValueType{0};
    }
  }

  /// \brief Checks whether the particle is passed through the given boundary
  /// box.
  /// \param boundaryBox The boundary box
  /// \return true if the particle passed the boundary box, false otherwise
  bool isPassedThroughBoundary(const BoundaryBox<T>& boundaryBox) const {
    return position[0] - radius < boundaryBox.low[0]
        || position[0] + radius > boundaryBox.high[0]
        || position[1] - radius < boundaryBox.low[1]
        || position[1] + radius > boundaryBox.high[1]
        || position[2] - radius < boundaryBox.low[2]
        || position[2] + radius > boundaryBox.high[2];
  }

  /// \brief Returns the boundary plane types collided with the particle from
  /// the given boundary box.
  /// \details The return vector will be empty if there is no collision.
  /// \param boundaryBox The boundary box
  /// \return The vector of boundary types collided with the particle
  std::vector<BoundaryPlaneType> getCollidedBoundaryPlaneTypes(
      const BoundaryBox<T>& boundaryBox
  ) const {
    auto v = std::vector<BoundaryPlaneType>{};
    if (position[0] - radius < boundaryBox.low[0] && velocity[0] < T{0}) {
      v.emplace_back(BoundaryPlaneType::NegativeX);
    }
    if (position[0] + radius > boundaryBox.high[0] && velocity[0] > T{0}) {
      v.emplace_back(BoundaryPlaneType::PositiveX);
    }
    if (position[1] - radius < boundaryBox.low[1] && velocity[1] < T{0}) {
      v.emplace_back(BoundaryPlaneType::NegativeY);
    }
    if (position[1] + radius > boundaryBox.high[1] && velocity[1] > T{0}) {
      v.emplace_back(BoundaryPlaneType::PositiveY);
    }
    if (position[2] - radius < boundaryBox.low[2] && velocity[2] < T{0}) {
      v.emplace_back(BoundaryPlaneType::NegativeZ);
    }
    if (position[2] + radius > boundaryBox.high[2] && velocity[2] > T{0}) {
      v.emplace_back(BoundaryPlaneType::PositiveZ);
    }
    return v;
  }

};

/// \brief Writes the particle to a character output stream.
/// \details Format:\n
/// \n
/// Index Mass Radius Position Velocity
/// \param ostream The character output stream
/// \param particle The particle to be extracted
/// \return The character output stream
template<std::floating_point T>
std::ostream& operator<<(std::ostream& ostream, const Particle<T>& particle) {
  ostream << particle.index << " "
          << particle.mass << " "
          << particle.radius << " "
          << particle.position << " "
          << particle.velocity;
  return ostream;
}

/// \brief Reads the particle from a character input stream
/// \details Format:\n
/// \n
/// Index Mass Radius Position Velocity
/// \param istream The character input stream
/// \param particle The particle to be inserted
/// \return The character input stream
template<std::floating_point T>
std::istream& operator>>(std::istream& istream, Particle<T>& particle) {
  istream >> particle.index
          >> particle.mass
          >> particle.radius
          >> particle.position
          >> particle.velocity;
  return istream;
}

} // namespace NBody::Simulation

#endif // CPP_PROJECT_INCLUDE_NBODY_PARTICLE_HPP_
