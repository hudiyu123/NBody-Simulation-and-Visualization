#ifndef CPP_PROJECT_INCLUDE_NBODY_SIMULATION_GRAVITYSOLVER_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_SIMULATION_GRAVITYSOLVER_HPP_

#include <NBody/Simulation/Particle.hpp>
#include <NBody/Simulation/Octree.hpp>
#include <omp.h>
#include <concepts>
#include <memory>
#include <vector>

namespace NBody::Simulation {

/// \brief Base class of gravity solvers.
template<std::floating_point T>
class BaseGravitySolver {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;

  virtual ~BaseGravitySolver() = default;

  /// \brief Computes the acceleration of given particles in place.
  /// \note Needs to be implemented in the derived class.
  /// \param particles The particles to be computed
  /// \param activeParticleIndices The indices of active particles
  virtual void computeAcceleration(
      std::vector<ParticleType>& particles,
      const std::vector<std::size_t>& activeParticleIndices
  ) const = 0;
};

/// \brief The gravity solver uses direct algorithm.
template<std::floating_point T>
class DirectGravitySolver : public BaseGravitySolver<T> {
 public:
  using ValueType = T;
  using VectorType = Math::Vector<ValueType, 3>;
  using ParticleType = Particle<ValueType>;

 private:
  const ValueType G_;
  const ValueType softening_;

 public:
  /// \brief Constructs a direct gravity solver.
  /// \param G The gravitational constant
  /// \param softening  The softening
  DirectGravitySolver(ValueType G, ValueType softening)
      : G_{G}, softening_{softening} {}

  /// \brief Computes the acceleration of given particles in place using direct
  /// algorithm.
  /// \param particles The particles to be computed
  /// \param activeParticleIndices The indices of active particles
  void computeAcceleration(
      std::vector<ParticleType>& particles,
      const std::vector<std::size_t>& activeParticleIndices
  ) const override {
    // Clear accelerations.
#pragma omp parallel for default(none) shared(particles, activeParticleIndices)
    for (std::size_t i = 0; i < activeParticleIndices.size(); ++i) {
      auto& particle = particles[activeParticleIndices[i]];
      particle.acceleration = VectorType{};
    }

#pragma omp parallel for default(none) shared(particles, activeParticleIndices)
    for (std::size_t i = 0; i < activeParticleIndices.size(); ++i) {
      const auto numThreads = omp_get_num_threads();
      if (numThreads == 1) {
        // Optimize algorithm under single thread context. The gravitational
        // force between two particles is only computed once.
        if (i >= (activeParticleIndices.size() + 1) / 2) { continue; }
        auto& particle = particles[activeParticleIndices[i]];
        for (std::size_t j = i + 1; j < activeParticleIndices.size(); ++j) {
          auto& otherParticle = particles[activeParticleIndices[j]];
          const auto difference = otherParticle.position - particle.position;
          const auto distance = Math::length(difference, softening_);
          const auto factor =
              G_ / std::pow(distance, ValueType{3}) * difference;
          particle.acceleration += factor * otherParticle.mass;
          otherParticle.acceleration -= factor * particle.mass;
        }
      } else {
        auto& particle = particles[activeParticleIndices[i]];
        for (std::size_t j = 0; j < activeParticleIndices.size(); ++j) {
          if (i == j) { continue; }
          const auto& otherParticle = particles[activeParticleIndices[j]];
          const auto difference = otherParticle.position - particle.position;
          const auto distance = Math::length(difference, softening_);
          const auto factor =
              G_ / std::pow(distance, ValueType{3}) * difference;
          particle.acceleration += factor * otherParticle.mass;
        }
      }
    }
  }
};

/// \brief Gravity solver using Barnes-Hut algorithm.
template<std::floating_point T>
class BarnesHutGravitySolver : public BaseGravitySolver<T> {
 public:
  using ValueType = T;
  using VectorType = Math::Vector<ValueType, 3>;
  using ParticleType = Particle<ValueType>;
  using OctreeType = Octree<ValueType>;

 private:
  const std::shared_ptr<OctreeType> octreePtr_;
  const ValueType G_;
  const ValueType softening_;
  const ValueType openingAngle_;

 public:
  /// \brief Constructs a Barnes-Hut gravity solver.
  /// \param octree A pointer to an octree
  /// \param G The gravitational constant
  /// \param softening The softening
  /// \param openingAngle The opening angle
  BarnesHutGravitySolver(
      const std::shared_ptr<Octree<ValueType>>& octree,
      ValueType G,
      ValueType softening,
      ValueType openingAngle
  ) : octreePtr_{octree},
      G_{G},
      softening_{softening},
      openingAngle_{openingAngle} {}

  /// \brief Copy constructor of BarnesHutGravitySolver is deleted.
  BarnesHutGravitySolver(const BarnesHutGravitySolver& other) = delete;
  /// \brief Copy assignment operator of BarnesHutGravitySolver is deleted.
  BarnesHutGravitySolver& operator==(
      const BarnesHutGravitySolver& other
  ) = delete;
  /// \brief Move constructor of BarnesHutGravitySolver is deleted.
  BarnesHutGravitySolver(BarnesHutGravitySolver&& other) = delete;
  /// \brief Move assignment operator of BarnesHutGravitySolver is deleted.
  BarnesHutGravitySolver& operator==(BarnesHutGravitySolver&& other) = delete;

  /// \brief Computes the acceleration of particles in place using Barnes-Hut
  /// algorithm.
  /// \param particles The particles to be computed
  /// \param activeParticleIndices The indices of active particles
  void computeAcceleration(
      std::vector<ParticleType>& particles,
      const std::vector<std::size_t>& activeParticleIndices
  ) const override {
#pragma omp parallel for default(none) \
    shared(particles, activeParticleIndices) schedule(guided)
    for (std::size_t i = 0; i < activeParticleIndices.size(); ++i) {
      auto approximateParticles = octreePtr_->getApproximateParticles(
          activeParticleIndices[i],
          particles,
          openingAngle_
      );
      auto& particle = particles[activeParticleIndices[i]];
      particle.acceleration = VectorType{};
      for (const auto& approximateParticle: approximateParticles) {
        const auto difference =
            particle.position - approximateParticle.position;
        const auto distance = Math::length(difference, softening_);
        const auto factor = -G_ / std::pow(distance, ValueType{3}) * difference;
        particle.acceleration += factor * approximateParticle.mass;
      }
    }
  }
};

} // namespace NBody::Simulation

#endif // CPP_PROJECT_INCLUDE_NBODY_SIMULATION_GRAVITYSOLVER_HPP_
