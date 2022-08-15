#ifndef CPP_PROJECT_INCLUDE_NBODY_SIMULATION_SIMULATOR_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_SIMULATION_SIMULATOR_HPP_

#include <NBody/Simulation/Configuration.hpp>
#include <NBody/Simulation/BoundaryBox.hpp>
#include <NBody/Simulation/Particle.hpp>
#include <NBody/Simulation/InitialData.hpp>
#include <NBody/Simulation/GravitySolver.hpp>
#include <NBody/Simulation/Integrator.hpp>
#include <NBody/Simulation/BoundaryConditionHandler.hpp>
#include <NBody/Simulation/Octree.hpp>
#include <NBody/Simulation/CollisionDetector.hpp>
#include <NBody/Simulation/CollisionResolver.hpp>

#include <concepts>
#include <vector>
#include <memory>

namespace NBody::Simulation {

/// \brief The nbody problem simulator.
template<std::floating_point T>
class Simulator {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;
  using OctreeType = Octree<ValueType>;

 private:
  /// \brief The delta time for each step.
  const ValueType timeStep_;
  /// \brief The current step of simulation.
  std::uint64_t currentStep_;
  /// \brief The current time elapsed of simulation.
  ValueType currentTime_;
  /// \brief The maximum finishStep will be forwarded by the simulation.
  const std::uint64_t finishStep_;
  /// \brief The boundary of the simulation space.
  const BoundaryBox<ValueType> boundaryBox_;
  /// \brief The pointer to the gravity solver.
  std::unique_ptr<BaseGravitySolver<ValueType>> gravitySolverPtr_;
  /// \brief The pointer to the integrator.
  std::unique_ptr<Integrator<ValueType>> integratorPtr_;
  /// \brief The pointer to the boundary condition handler.
  std::unique_ptr<BaseBoundaryConditionHandler<ValueType>>
      boundaryConditionHandlerPtr_;
  /// \brief The pointer to the collision detector.
  std::unique_ptr<BaseCollisionDetector<ValueType>> collisionDetectorPtr_;
  /// \brief The pointer to the collision resolver.
  std::unique_ptr<BaseCollisionResolver<ValueType>> collisionResolverPtr_;
  /// \brief The pointer to the octree used in the gravity solver and the
  /// collision detector.
  std::shared_ptr<OctreeType> octreePtr_;
  /// \brief The particles participate in the simulation.
  /// \details Particles will not be removed from the particles vector after
  /// insertion.
  std::vector<ParticleType> particles_;
  /// \brief The vector of the active particle indices.
  /// \details This vector determines whether a particle is still active in
  /// simulation. If a particle becomes inactive due to mergence or passing the
  /// boundary, its index will be removed from this vector instead of removing
  /// the particle from the particles vector.\n
  /// \n
  /// This vector is sorted so binary search is used to find the activeness of
  /// arbitrary particle.
  /// \note The elements stored in this vector are particles' index in particles
  /// vector (particles_), not the index read from input.
  std::vector<std::size_t> activeParticleIndices_;

 public:
  /// \brief Constructs a simulator with configuration and initial data.
  /// \param configuration The configuration
  /// \param initialData  The initial data
  Simulator(
      const Configuration<ValueType>& configuration,
      const InitialData<ValueType>& initialData
  ) : timeStep_{configuration.timeStep},
      finishStep_{configuration.finishStep},
      boundaryBox_{initialData.boundaryBox},
      currentStep_{0},
      currentTime_{0},
      particles_{initialData.particles},
      activeParticleIndices_{} {
    // Initialize the active particle indices vector.
    for (std::size_t i = 0; i < particles_.size(); ++i) {
      activeParticleIndices_.emplace_back(i);
    }

    // Constructs the gravity solver.
    switch (configuration.gravitySolverType) {
      case GravitySolverType::Direct:
        gravitySolverPtr_ = std::make_unique<DirectGravitySolver<ValueType>>(
            configuration.G,
            configuration.softening
        );
        break;
      case GravitySolverType::BarnesHut:
        if (!octreePtr_) {
          octreePtr_ = std::make_shared<Octree<ValueType>>(
              initialData.boundaryBox
          );
          octreePtr_->update(particles_, activeParticleIndices_);
        }

        gravitySolverPtr_ = std::make_unique<BarnesHutGravitySolver<ValueType>>(
            octreePtr_,
            configuration.G,
            configuration.softening,
            configuration.openingAngle
        );
        break;
    }

    // Constructs the integrator.
    switch (configuration.integratorType) {
      case IntegratorType::Leapfrog:
        integratorPtr_ = std::make_unique<LeapfrogIntegrator<ValueType>>(
            configuration.timeStep
        );
        break;
    }

    // Constructs the boundary condition handler.
    switch (configuration.boundaryConditionHandlerType) {
      case BoundaryConditionHandlerType::None:
        break;
      case BoundaryConditionHandlerType::Open:
        boundaryConditionHandlerPtr_ =
            std::make_unique<OpenBoundaryConditionHandler<ValueType>>(
                initialData.boundaryBox
            );
        break;
      case BoundaryConditionHandlerType::Rebound:
        boundaryConditionHandlerPtr_ =
            std::make_unique<ReboundBoundaryConditionHandler<ValueType>>(
                initialData.boundaryBox,
                configuration.boundaryCOR
            );
        break;
      case BoundaryConditionHandlerType::Periodic:
        boundaryConditionHandlerPtr_ =
            std::make_unique<PeriodicBoundaryConditionHandler<ValueType>>(
                initialData.boundaryBox
            );
        break;
    }

    // Constructs the collision detector.
    switch (configuration.collisionDetectorType) {
      case CollisionDetectorType::None:
        break;
      case CollisionDetectorType::Direct:
        collisionDetectorPtr_ =
            std::make_unique<DirectCollisionDetector<ValueType>>();
        break;
      case CollisionDetectorType::Octree:
        if (!octreePtr_) {
          octreePtr_ = std::make_shared<Octree<ValueType>>(
              initialData.boundaryBox
          );
          octreePtr_->update(particles_, activeParticleIndices_);
        }

        collisionDetectorPtr_ =
            std::make_unique<OctreeCollisionDetector<ValueType>>(octreePtr_);
        break;
    }

    // Constructs the collision resolver if the type of collision detector is
    // not None.
    if (configuration.collisionDetectorType != CollisionDetectorType::None) {
      switch (configuration.collisionResolverType) {
        case CollisionResolverType::Hard:
          collisionResolverPtr_ =
              std::make_unique<HardCollisionResolver<ValueType>>(
                  configuration.COR,
                  configuration.mergeThreshold
              );
          break;
        case CollisionResolverType::Merge:
          collisionResolverPtr_ =
              std::make_unique<MergeCollisionResolver<ValueType>>();
          break;
      }
    }
  }

  /// \brief The copy constructor is deleted.
  Simulator(const Simulator& other) = delete;
  /// \brief The copy assignment operator is deleted.
  Simulator& operator==(const Simulator& other) = delete;
  /// \brief The move constructor is deleted.
  Simulator(Simulator&& other) = delete;
  /// \brief The move assignment operator is deleted.
  Simulator& operator==(Simulator&& other) = delete;

  /// \brief Forwards the simulation by delta steps and returns the particles.
  /// \param deltaStep The delta steps
  /// \return The particles
  std::vector<ParticleType> forwardStep(std::uint64_t deltaStep) {
    if (isFinished()) { return particles_; }
    const auto remainingSteps = finishStep_ - currentStep_;
    if (remainingSteps >= deltaStep) {
      currentStep_ += deltaStep;
      return forward(deltaStep);
    } else {
      currentStep_ += remainingSteps;
      return forward(remainingSteps);
    }
  }

  /// \brief Forwards the simulation by delta time and returns the particles.
  /// \param deltaTime The delta time
  /// \return The particles
  std::vector<ParticleType> forwardTime(ValueType deltaTime) {
    if (isFinished()) { return particles_; }
    currentTime_ += deltaTime;
    constexpr auto max = std::numeric_limits<std::uint64_t>::max();
    std::uint64_t deltaStep;
    // Check overflow during conversion.
    if (currentTime_ / timeStep_ > static_cast<ValueType>(max)) {
      deltaStep = max - currentStep_;
    } else {
      deltaStep =
          static_cast<std::uint64_t>(currentTime_ / timeStep_) - currentStep_;
    }
    return forwardStep(deltaStep);
  }

  /// \brief Checks whether the simulation is finished.
  /// \return true if the simulation is finished, false otherwise
  [[nodiscard]] bool isFinished() const {
    return currentStep_ == finishStep_;
  }

  /// \brief Returns the boundary box of the simulation.
  /// \return The boundary box
  BoundaryBox<ValueType> getBoundaryBox() const {
    return boundaryBox_;
  }

  /// \brief Returns the particles at current step.
  /// \return The particles
  std::vector<ParticleType> getParticles() const {
    return particles_;
  }

  /// \brief Returns the current step.
  /// \return The current step
  [[nodiscard]] std::uint64_t getCurrentStep() const {
    return currentStep_;
  }

  /// \brief Returns the time step of the simulation.
  /// \return The time step
  [[nodiscard]] ValueType getTimeStep() const {
    return timeStep_;
  }

 private:
  std::vector<ParticleType> forward(std::uint64_t deltaStep) {
    for (std::uint64_t i = 0; i < deltaStep; ++i) {
      // Compute interaction (acceleration).
      gravitySolverPtr_->computeAcceleration(
          particles_,
          activeParticleIndices_
      );
      // Compute integration.
      integratorPtr_->forward(particles_, activeParticleIndices_);
      // Reconstruct Octree after integration.
      if (octreePtr_) {
        octreePtr_->update(particles_, activeParticleIndices_);
      }
      // Handle boundary condition.
      if (boundaryConditionHandlerPtr_) {
        boundaryConditionHandlerPtr_->handle(
            particles_,
            activeParticleIndices_
        );
      }
      // Detect and resolve collisions.
      if (collisionDetectorPtr_) {
        auto pairs = collisionDetectorPtr_->getCollidedParticleIndexPairs(
            particles_,
            activeParticleIndices_
        );
        collisionResolverPtr_->resolveCollisions(
            particles_,
            activeParticleIndices_,
            pairs
        );
        // Reconstruct Octree after collision.
        if (octreePtr_) {
          octreePtr_->update(particles_, activeParticleIndices_);
        }
      }
    }

    // Only return active particles.
    auto activeParticles = std::vector<ParticleType>{};
    activeParticles.reserve(activeParticleIndices_.size());
    for (auto index: activeParticleIndices_) {
      activeParticles.emplace_back(particles_[index]);
    }
    return activeParticles;
  }
};

} // namespace NBody::Simulation

#endif // CPP_PROJECT_INCLUDE_NBODY_SIMULATION_SIMULATOR_HPP_
