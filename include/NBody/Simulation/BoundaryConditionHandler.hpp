#ifndef CPP_PROJECT_INCLUDE_NBODY_SIMULATION_BOUNDARYCONDITIONHANDLER_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_SIMULATION_BOUNDARYCONDITIONHANDLER_HPP_

#include <NBody/Simulation/Particle.hpp>
#include <algorithm>
#include <concepts>
#include <vector>
#include <omp.h>

namespace NBody::Simulation {

/// \brief Base class of boundary condition handlers.
template<std::floating_point T>
class BaseBoundaryConditionHandler {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;

  virtual ~BaseBoundaryConditionHandler() = default;

  /// \brief Handles the boundary condition of given particles.
  /// \note Needs to be implemented in the derived class.
  /// \param particles The particles to be handled
  /// \param activeParticleIndices The indices of active particles
  virtual void handle(
      std::vector<ParticleType>& particles,
      std::vector<std::size_t>& activeParticleIndices
  ) const = 0;
};

/// \brief The boundary handler assumes the open boundary condition.
template<std::floating_point T>
class OpenBoundaryConditionHandler : public BaseBoundaryConditionHandler<T> {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;
  using BoundaryBoxType = BoundaryBox<ValueType>;

 private:
  const BoundaryBoxType boundaryBox_;

 public:
  /// \brief Returns the indices of particles passed boundary.
  /// \param particles The particles to be handled
  /// \param activeParticleIndices The indices of active particles
  /// \param boundaryBox The boundary box
  /// \return The indices of particles passed boundary
  static std::vector<std::size_t> computePassedParticleIndices(
      const std::vector<ParticleType>& particles,
      const std::vector<std::size_t>& activeParticleIndices,
      const BoundaryBoxType& boundaryBox
  ) {
    auto passedParticleIndices = std::vector<std::size_t>{};
#pragma omp parallel for default(none) \
    shared( \
    particles, \
    activeParticleIndices, \
    passedParticleIndices, \
    boundaryBox \
    ) \
    schedule(guided)
    for (std::size_t i = 0; i < activeParticleIndices.size(); ++i) {
      const auto& particle = particles[activeParticleIndices[i]];
      if (particle.isPassedThroughBoundary(boundaryBox)) {
#pragma omp critical
        passedParticleIndices.emplace_back(activeParticleIndices[i]);
      }
    }
    return passedParticleIndices;
  }

  /// \brief Constructs an open boundary condition handler.
  /// \param boundaryBox The boundary box
  explicit OpenBoundaryConditionHandler(const BoundaryBoxType& boundaryBox)
      : boundaryBox_{boundaryBox} {}

  /// \brief Handles the boundary condition by assuming the open boundary
  /// condition.
  /// \param particles The particles to be handled
  /// \param activeParticleIndices The indices of active particles
  /// \param boundaryBox The boundary box
  static void handleOpenBoundaryCondition(
      std::vector<ParticleType>& particles,
      std::vector<std::size_t>& activeParticleIndices,
      const BoundaryBoxType& boundaryBox
  ) {
    auto passedParticleIndices = computePassedParticleIndices(
        particles,
        activeParticleIndices,
        boundaryBox
    );
    // Removes the particle from the active particle indices vector if it passes
    // the boundary.
    for (auto index: passedParticleIndices) {
      activeParticleIndices.erase(std::find(
          activeParticleIndices.cbegin(),
          activeParticleIndices.cend(),
          index
      ));
    }
  }

  /// \brief Handles the boundary condition by assuming the open boundary
  /// condition.
  /// \note This function simply calls the #handleOpenBoundaryCondition
  /// function.
  /// \param particles The particles to be handled
  /// \param activeParticleIndices The indices of active particles
  void handle(
      std::vector<ParticleType>& particles,
      std::vector<std::size_t>& activeParticleIndices
  ) const override {
    handleOpenBoundaryCondition(particles, activeParticleIndices, boundaryBox_);
  }
};

/// \brief The boundary handler assumes the periodic boundary condition.
template<std::floating_point T>
class PeriodicBoundaryConditionHandler
    : public BaseBoundaryConditionHandler<T> {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;
  using BoundaryBoxType = BoundaryBox<ValueType>;

 private:
  const BoundaryBoxType boundaryBox_;

 public:
  /// \brief Constructs a periodic boundary condition handler.
  /// \param boundaryBox The boundary box
  explicit PeriodicBoundaryConditionHandler(const BoundaryBoxType& boundaryBox)
      : boundaryBox_{boundaryBox} {}

  /// \brief Handles the boundary condition by assuming the periodic boundary
  /// condition.
  /// \param particles The particles to be handled
  /// \param activeParticleIndices The indices of active particles
  void handle(
      std::vector<ParticleType>& particles,
      std::vector<std::size_t>& activeParticleIndices
  ) const override {
    constexpr auto handler = [](
        ValueType& position,
        ValueType low,
        ValueType high
    ) {
      const auto width = high - low;
      if (position < low) {
        position = high - std::fmod(low - position, width);
      } else if (position > high) {
        position = low + std::fmod(position - high, width);
      }
    };

#pragma omp parallel for default(none) \
  shared(particles, activeParticleIndices, handler) schedule(guided)
    for (std::size_t i = 0; i < activeParticleIndices.size(); ++i) {
      auto& particle = particles[activeParticleIndices[i]];
      handler(particle.position[0], boundaryBox_.low[0], boundaryBox_.high[0]);
      handler(particle.position[1], boundaryBox_.low[1], boundaryBox_.high[1]);
      handler(particle.position[2], boundaryBox_.low[2], boundaryBox_.high[2]);
    }
  }
};

/// \brief The boundary handler assumes the rebound boundary condition.
template<std::floating_point T>
class ReboundBoundaryConditionHandler : public BaseBoundaryConditionHandler<T> {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;
  using BoundaryBoxType = BoundaryBox<ValueType>;

 private:
  const BoundaryBoxType boundaryBox_;
  const ValueType e_;

 public:
  /// \brief Constructs a rebound boundary condition handler.
  /// \param boundaryBox The boundary box
  /// \param e The coefficient of restitution of the boundaries
  ReboundBoundaryConditionHandler(BoundaryBoxType boundaryBox, ValueType e)
      : boundaryBox_{boundaryBox}, e_{e} {}

  /// \brief Handles the boundary condition by assuming the rebound boundary
  /// condition.
  /// \note Treats rebound as open if the coefficient of restitution is 0.
  /// \param particles The particles to be handled
  /// \param activeParticleIndices The indices of active particles
  void handle(
      std::vector<ParticleType>& particles,
      std::vector<std::size_t>& activeParticleIndices
  ) const override {
    if (e_ == ValueType{0}) {
      // Treats rebound as open if coefficient of restitution is 0.
      OpenBoundaryConditionHandler<ValueType>::handleOpenBoundaryCondition(
          particles,
          activeParticleIndices,
          boundaryBox_
      );
    } else {
#pragma omp parallel for default(none) \
  shared(particles, activeParticleIndices) schedule(guided)
      for (std::size_t i = 0; i < activeParticleIndices.size(); ++i) {
        auto& particle = particles[activeParticleIndices[i]];
        const auto collidedBoundaryTypes =
            particle.getCollidedBoundaryPlaneTypes(boundaryBox_);
        for (auto type: collidedBoundaryTypes) {
          switch (type) {
            case BoundaryPlaneType::NegativeX:
              [[fallthrough]];
            case BoundaryPlaneType::PositiveX:
              // Reverts the x coordinate if the particle cross the boundaries
              // on the x-axis.
              particle.velocity[0] = -particle.velocity[0];
              break;
            case BoundaryPlaneType::NegativeY:
              [[fallthrough]];
            case BoundaryPlaneType::PositiveY:
              // Reverts the y coordinate if the particle cross the boundaries
              // on the y-axis.
              particle.velocity[1] = -particle.velocity[1];
              break;
            case BoundaryPlaneType::NegativeZ:
              [[fallthrough]];
            case BoundaryPlaneType::PositiveZ:
              // Reverts the z coordinate if the particle cross the boundaries
              // on the z-axis.
              particle.velocity[2] = -particle.velocity[2];
              break;
          }
        }
        // Handles the restitution.
        if (!collidedBoundaryTypes.empty()) {
          particle.velocity *= e_;
        }
      }
    }
  }
};

} // namespace NBody::Simulation

#endif // CPP_PROJECT_INCLUDE_NBODY_SIMULATION_BOUNDARYCONDITIONHANDLER_HPP_
