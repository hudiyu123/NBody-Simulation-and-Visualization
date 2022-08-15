#ifndef CPP_PROJECT_INCLUDE_NBODY_SIMULATION_INITIALDATA_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_SIMULATION_INITIALDATA_HPP_

#include <NBody/Simulation/BoundaryBox.hpp>
#include <NBody/Simulation/Particle.hpp>
#include <concepts>
#include <vector>
#include <iostream>

namespace NBody::Simulation {

/// \brief This data class stores the initial data required for nbody simulation.
template<std::floating_point T>
struct InitialData {
  using ParticleType = Particle<T>;
  using BoundaryBoxType = BoundaryBox<T>;

  /// \brief Default constructor.
  InitialData() = default;

  /// \brief The simulation space of the simulation.
  BoundaryBoxType boundaryBox;
  /// \brief The particles participate in the simulation.
  std::vector<ParticleType> particles;
};

/// \brief Reads the initial data from a character input stream.
/// \param is The character input stream
/// \param initialData The initial data to be inserted
/// \return The character input stream
template<std::floating_point T>
std::istream& operator>>(std::istream& is, InitialData<T>& initialData) {
  if (!(is >> initialData.boundaryBox)) { return is; }
  std::uint64_t size;
  if (!(is >> size)) { return is; }
  Particle<T> particle;
  for (std::uint64_t i = 0; i < size; ++i) {
    if (!(is >> particle)) { return is; }
    initialData.particles.emplace_back(std::move(particle));
  }
  return is;
}

} // namespace NBody::Simulation

#endif //CPP_PROJECT_INCLUDE_NBODY_SIMULATION_INITIALDATA_HPP_
