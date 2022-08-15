#ifndef CPP_PROJECT_INCLUDE_NBODY_SIMULATION_INTEGRATOR_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_SIMULATION_INTEGRATOR_HPP_

#include <NBody/Simulation/Particle.hpp>
#include <omp.h>
#include <concepts>
#include <vector>

namespace NBody::Simulation {

/// \brief Base class of integrators.
template<std::floating_point T>
class Integrator {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;

  virtual ~Integrator() = default;

  /// \brief Computes the integration of given particles in place.
  /// \note Needs to be implemented in the derived class.
  /// \param particles The particles to be forwarded
  /// \param activeParticleIndices The indices of active particles
  virtual void forward(
      std::vector<ParticleType>& particles,
      const std::vector<std::size_t>& activeParticleIndices
  ) = 0;
};

/// \brief Integrator uses Leapfrog method.
template<std::floating_point T>
class LeapfrogIntegrator : public Integrator<T> {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;

 private:
  const ValueType timeStep_;

 public:
  /// \brief Constructs a Leapfrog integrator.
  /// \param timeStep The time step
  explicit LeapfrogIntegrator(ValueType timeStep) : timeStep_{timeStep} {}

  /// \brief Forwards particles in place using Leapfrog method.
  /// \param particles The particles to be forwarded
  /// \param activeParticleIndices The indices of active particles
  void forward(
      std::vector<ParticleType>& particles,
      const std::vector<std::size_t>& activeParticleIndices
  ) override {
#pragma omp parallel for default(none) \
  shared(particles, activeParticleIndices, timeStep_)
    for (std::size_t i = 0; i < activeParticleIndices.size(); ++i) {
      auto& particle = particles[activeParticleIndices[i]];
      // First stage.
      particle.position += ValueType{0.5} * timeStep_ * particle.velocity;
    }
#pragma omp parallel for default(none) \
  shared(particles, activeParticleIndices, timeStep_)
    for (std::size_t i = 0; i < activeParticleIndices.size(); ++i) {
      auto& particle = particles[activeParticleIndices[i]];
      // Second stage.
      particle.velocity += particle.acceleration * timeStep_;
      // Third stage.
      particle.position += ValueType{0.5} * timeStep_ * particle.velocity;
    }
  }
};

} // namespace NBody::Simulation

#endif // CPP_PROJECT_INCLUDE_NBODY_SIMULATION_INTEGRATOR_HPP_
