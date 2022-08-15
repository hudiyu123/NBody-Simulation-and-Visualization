#ifndef CPP_PROJECT_INCLUDE_NBODY_SIMULATION_COLLISIONDETECTOR_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_SIMULATION_COLLISIONDETECTOR_HPP_

#include <NBody/Simulation/Particle.hpp>
#include <NBody/Simulation/Octree.hpp>
#include <concepts>
#include <memory>
#include <set>
#include <utility>
#include <vector>

namespace NBody::Simulation {

/// \brief Stores a pair of particle indices.
/// \details The first index is always less than the second index.
struct ParticleIndexPair {
 public:
  /// \brief The first index (always the smaller one).
  const std::size_t first;
  /// \brief The second index (always the greater one).
  const std::size_t second;

  /// \brief Constructs a particle index pair with given indices.
  /// \param index The index
  /// \param otherIndex The other index
  ParticleIndexPair(std::size_t index, std::size_t otherIndex)
      : first{std::min(index, otherIndex)},
        second{std::max(index, otherIndex)} {}

  /// \brief Returns the result of three-way comparison.
  /// \param other The other particle
  /// \return The result of three-way comparison
  auto operator<=>(const ParticleIndexPair& other) const = default;
};

/// \brief Base class of collision detectors.
template<std::floating_point T>
class BaseCollisionDetector {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;

  virtual ~BaseCollisionDetector() = default;

  /// \brief Computes and returns the collided particle index pairs of given
  /// particles.
  /// \note Needs to be implemented in the derived class.
  /// \param particles The particles to be detected
  /// \param activeParticleIndices The indices of active particles
  /// \return The set of collided particle index pairs
  virtual std::set<ParticleIndexPair> getCollidedParticleIndexPairs(
      const std::vector<ParticleType>& particles,
      const std::vector<std::size_t>& activeParticleIndices
  ) const = 0;
};

/// \brief Collision detector uses direct algorithm.
template<std::floating_point T>
class DirectCollisionDetector : public BaseCollisionDetector<T> {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;

  /// \brief Computes and returns the collided particle index pairs using the
  /// direct algorithm.
  /// \param particles The particles to be detected
  /// \param activeParticleIndices The indices of active particles
  /// \return The set of collided particle index pairs
  static std::set<ParticleIndexPair> computeCollidedParticleIndexPairs(
      const std::vector<ParticleType>& particles,
      const std::vector<std::size_t>& activeParticleIndices
  ) {
    auto pairs = std::set<ParticleIndexPair>{};
#pragma omp parallel for default(none) collapse(2) \
  shared(particles, activeParticleIndices, pairs) schedule(guided)
    for (std::size_t i = 0; i < activeParticleIndices.size(); ++i) {
      for (std::size_t j = 0; j < activeParticleIndices.size(); ++j) {
        if (i == j) { continue; }
        const auto& particle = particles[activeParticleIndices[i]];
        const auto& theOtherParticle = particles[activeParticleIndices[j]];
        if (particle.isCollidedWith(theOtherParticle)) {
#pragma omp critical
          pairs.emplace(activeParticleIndices[i], activeParticleIndices[j]);
        }
      }
    }
    return pairs;
  }

  /// \brief Computes and returns the collided particle index pairs using the
  /// direct algorithm.
  /// \note This function simply calls the computeCollidedParticleIndexPairs
  /// function.
  /// \param particles The particles to be detected
  /// \param activeParticleIndices The indices of active particles
  /// \return The set of collided particle index pairs
  std::set<ParticleIndexPair> getCollidedParticleIndexPairs(
      const std::vector<ParticleType>& particles,
      const std::vector<std::size_t>& activeParticleIndices
  ) const override {
    return computeCollidedParticleIndexPairs(particles, activeParticleIndices);
  }
};

/// \brief Collision detector uses the octree method.
template<std::floating_point T>
class OctreeCollisionDetector : public BaseCollisionDetector<T> {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;
  using OctreeType = Octree<ValueType>;

 private:
  const std::shared_ptr<OctreeType> octreePtr_;

 public:
  /// \brief Constructs an octree collision detector.
  /// \param octreePtr The pointer to an octree
  explicit OctreeCollisionDetector(
      const std::shared_ptr<OctreeType>& octreePtr
  ) : octreePtr_{octreePtr} {}

  /// \brief Copy constructor of OctreeCollisionDetector is deleted.
  OctreeCollisionDetector(const OctreeCollisionDetector& other) = delete;
  /// \brief Copy assignment operator of OctreeCollisionDetector is deleted.
  OctreeCollisionDetector& operator==(
      const OctreeCollisionDetector& other
  ) = delete;
  /// \brief Move constructor of OctreeCollisionDetector is deleted.
  OctreeCollisionDetector(OctreeCollisionDetector&& other) = delete;
  /// \brief Move assignment operator of OctreeCollisionDetector is deleted.
  OctreeCollisionDetector& operator==(OctreeCollisionDetector&& other) = delete;

  /// \brief Computes and returns the collided particle index pairs using the
  /// octree method.
  /// \param particles The particles to be detected
  /// \param activeParticleIndices The indices of active particles
  /// \return The set of collided particle index pairs
  std::set<ParticleIndexPair> getCollidedParticleIndexPairs(
      const std::vector<ParticleType>& particles,
      const std::vector<std::size_t>& activeParticleIndices
  ) const override {
    auto pairs = std::set<ParticleIndexPair>{};
#pragma omp parallel for default(none) \
  shared(particles, activeParticleIndices, pairs) schedule(guided)
    for (std::size_t i = 0; i < activeParticleIndices.size(); ++i) {
      const auto collidedParticleIndices =
          octreePtr_->getCollidedParticleIndices(
              activeParticleIndices[i],
              particles
          );
      if (!collidedParticleIndices.empty()) {
#pragma omp critical
        {
          for (auto index: collidedParticleIndices) {
            pairs.emplace(activeParticleIndices[i], index);
          }
        }
      }
    }
    return pairs;
  }
};

} // namespace NBody::Simulation

#endif //CPP_PROJECT_INCLUDE_NBODY_SIMULATION_COLLISIONDETECTOR_HPP_
