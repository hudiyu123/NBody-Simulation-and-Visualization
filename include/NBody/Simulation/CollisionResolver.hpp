#ifndef CPP_PROJECT_INCLUDE_NBODY_SIMULATION_COLLISIONRESOLVER_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_SIMULATION_COLLISIONRESOLVER_HPP_

#include <NBody/Simulation/CollisionDetector.hpp>
#include <algorithm>
#include <cmath>
#include <concepts>
#include <map>
#include <numbers>
#include <set>
#include <vector>

namespace NBody::Simulation {

/// \brief Base class of collision resolvers.
template<std::floating_point T>
class BaseCollisionResolver {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;

  virtual ~BaseCollisionResolver() = default;

  /// \brief Resolves the collisions of particles in place with given collided
  /// particle index pairs.
  /// \note Needs to be implemented in the derived class.
  /// \param particles The particles to be resolved
  /// \param activeParticleIndices The indices of active particles
  /// \param collidedParticleIndexPairs The collided particle index pairs
  virtual void resolveCollisions(
      std::vector<ParticleType>& particles,
      std::vector<std::size_t>& activeParticleIndices,
      const std::set<ParticleIndexPair>& collidedParticleIndexPairs
  ) const = 0;
};

/// \brief Collision resolver assumes merge after collision.
template<std::floating_point T>
class MergeCollisionResolver : public BaseCollisionResolver<T> {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;

  /// \brief Resolves the collisions by merging the two colliding particles.
  /// \note This function simply calls the static resolveMergeCollisions
  /// function.
  /// \param particles The particles to be resolved
  /// \param activeParticleIndices The indices of active particles
  /// \param collidedParticleIndexPairs The collided particle index pairs
  void resolveCollisions(
      std::vector<ParticleType>& particles,
      std::vector<std::size_t>& activeParticleIndices,
      const std::set<ParticleIndexPair>& collidedParticleIndexPairs
  ) const override {
    resolveMergeCollisions(
        particles,
        activeParticleIndices,
        collidedParticleIndexPairs
    );
  }

  /// \brief Resolves the collisions by merging the two colliding particles.
  /// \param particles The particles to be resolved
  /// \param activeParticleIndices The indices of active particles
  /// \param collidedParticleIndexPairs The collided particle index pairs
  static void resolveMergeCollisions(
      std::vector<ParticleType>& particles,
      std::vector<std::size_t>& activeParticleIndices,
      const std::set<ParticleIndexPair>& collidedParticleIndexPairs
  ) {
    for (const auto& pair: collidedParticleIndexPairs) {
      // Since
      const auto isParticleActive = std::ranges::binary_search(
          activeParticleIndices,
          pair.first
      );
      if (isParticleActive) {
        const auto& p1 = particles[pair.first];
        auto& p2 = particles[pair.second];
        p2.radius = std::cbrt(std::pow(p1.radius, 3) + std::pow(p2.radius, 3));
        p2.velocity = (p1.mass * p1.velocity + p2.mass * p2.velocity)
            / (p1.mass + p2.mass);
        p2.position = (p1.mass * p1.position + p2.mass * p2.position)
            / (p1.mass + p2.mass);
        p2.mass += p1.mass;
        activeParticleIndices.erase(std::find(
            activeParticleIndices.cbegin(),
            activeParticleIndices.cend(),
            pair.first)
        );
      }
    }
  }
};

/// \brief Collision resolver assumes hard-sphere collision.
template<std::floating_point T>
class HardCollisionResolver : public BaseCollisionResolver<T> {
 public:
  using ValueType = T;
  using ParticleType = Particle<ValueType>;
  using VectorType = Math::Vector<ValueType, 3>;

 private:
  const ValueType e_;
  const std::uint64_t mergeThreshold_;
  /// \brief A recoder stores the count of consecutive collisions of two
  /// particles.
  mutable std::map<ParticleIndexPair, std::uint64_t> collisionRecorder_;

 public:
  /// \brief Constructs a hard collision resolver that assumes hard-sphere
  /// collision.
  /// \param e The coefficient of restitution
  /// \param mergeThreshold The merge threshold
  explicit HardCollisionResolver(ValueType e, std::uint64_t mergeThreshold)
      : e_{e}, mergeThreshold_{mergeThreshold}, collisionRecorder_{} {}

  /// \brief Resolves the collisions of particles by assuming hard-sphere
  /// collision.
  /// \note The two colliding particles will be merged if they consecutively
  /// collide with each other and the count of the consecutive collision reaches
  /// the merge threshold.
  /// \param particles The particles to be resolved
  /// \param activeParticleIndices The indices of active particles
  /// \param collidedParticleIndexPairs The collided particle index pairs
  void resolveCollisions(
      std::vector<ParticleType>& particles,
      std::vector<std::size_t>& activeParticleIndices,
      const std::set<ParticleIndexPair>& collidedParticleIndexPairs
  ) const override {
    resolveReboundCollisions(
        particles,
        activeParticleIndices,
        collidedParticleIndexPairs,
        e_
    );

    // Merge two particles if they consecutively collide with each other over
    // mergeThreshold times.
    if (mergeThreshold_ > 0) {
      // Clear the record if the tow particle do not collide.
      for (auto it = collisionRecorder_.begin();
           it != collisionRecorder_.end();) {
        if (!collidedParticleIndexPairs.contains(it->first)) {
          it = collisionRecorder_.erase(it);
        } else {
          ++it;
        }
      }

      // Increase the count for each pair of collided particles. Add the pair
      // into merge list if its count is over the threshold.
      auto mergedParticleIndexPairs = std::set<ParticleIndexPair>{};
      for (const auto& pair: collidedParticleIndexPairs) {
        if (++collisionRecorder_[pair] >= mergeThreshold_) {
          mergedParticleIndexPairs.emplace(pair);
        }
      }

      MergeCollisionResolver<ValueType>::resolveMergeCollisions(
          particles,
          activeParticleIndices,
          mergedParticleIndexPairs
      );
    }
  }

  /// \brief Resolves the
  /// <a href="https://www.plasmaphysics.org.uk/collision3d.htm"> elastic and
  /// inelastic collision </a> in place.
  /// \param particles The particles to be resolved
  /// \param activeParticleIndices The indices of active particles
  /// \param collidedParticleIndexPairs The collided particle index pairs
  /// \param e The coefficient of restitution
  static void resolveReboundCollisions(
      std::vector<ParticleType>& particles,
      std::vector<std::size_t>& activeParticleIndices,
      const std::set<ParticleIndexPair>& collidedParticleIndexPairs,
      ValueType e
  ) {
    for (const auto& pair: collidedParticleIndexPairs) {
      auto& particle = particles[pair.first];
      auto& otherParticle = particles[pair.second];
      auto m1 = particle.mass;
      auto m2 = otherParticle.mass;
      auto& v1 = particle.velocity;
      auto& v2 = otherParticle.velocity;
      auto v_cm = (m1 * v1 + m2 * v2) / (m1 + m2);

      auto p1 = particle.position;
      auto p2 = otherParticle.position;

      auto r12 = particle.radius + otherParticle.radius;
      auto m21 = otherParticle.mass / particle.mass;
      auto p21 = p2 - p1;
      auto v21 = v2 - v1;
      auto d = Math::length(p21);
      auto v = Math::length(v21);
      p2 = p21;
      v1 = -v21;

      auto theta2 = std::acos(p2[2] / d);
      auto phi2 = ValueType{};
      if (p2[0] == ValueType{0} && p2[1] == ValueType{0}) {
        phi2 = ValueType{};
      } else {
        phi2 = std::atan2(p2[1], p2[0]);
      }
      auto st = std::sin(theta2);
      auto ct = std::cos(theta2);
      auto sp = std::sin(phi2);
      auto cp = std::cos(phi2);

      auto v1r = VectorType{};
      v1r[0] = ct * cp * v1[0] + ct * sp * v1[1] - st * v1[2];
      v1r[1] = cp * v1[1] - sp * v1[0];
      v1r[2] = st * cp * v1[0] + st * sp * v1[1] + ct * v1[2];

      auto fvz1r = v1r[2] / v;
      if (fvz1r > ValueType{1}) {
        fvz1r = ValueType{1};
      } else if (fvz1r < ValueType{-1}) {
        fvz1r = ValueType{-1};
      }

      auto thetav = std::acos(fvz1r);

      auto phiv = ValueType{};
      if (v1r[0] == ValueType{0} && v1r[1] == ValueType{0}) {
        phiv = ValueType{0};
      } else {
        phiv = std::atan2(v1r[1], v1r[0]);
      }

      auto dr = d * std::sin(thetav) / r12;

      auto alpha = std::asin(-dr);
      auto beta = phiv;
      auto sbeta = std::sin(beta);
      auto cbeta = std::cos(beta);

      auto a = std::tan(thetav + alpha);
      auto dvz2 =
          ValueType{2} * (v1r[2] + a * (cbeta * v1r[0] + sbeta * v1r[1]))
              / ((ValueType{1} + a * a) * (ValueType{1} + m21));
      auto v2r = VectorType{};
      v2r[2] = dvz2;
      v2r[0] = a * cbeta * dvz2;
      v2r[1] = a * sbeta * dvz2;
      v1r -= m21 * v2r;

      v1[0] = ct * cp * v1r[0] - sp * v1r[1] + st * cp * v1r[2] + v2[0];
      v1[1] = ct * sp * v1r[0] + cp * v1r[1] + st * sp * v1r[2] + v2[1];
      v1[2] = ct * v1r[2] - st * v1r[0] + v2[2];

      v2[0] = ct * cp * v2r[0] - sp * v2r[1] + st * cp * v2r[2] + v2[0];
      v2[1] = ct * sp * v2r[0] + cp * v2r[1] + st * sp * v2r[2] + v2[1];
      v2[2] = ct * v2r[2] - st * v2r[0] + v2[2];

      v1 = (v1 - v_cm) * e + v_cm;
      v2 = (v2 - v_cm) * e + v_cm;
    }
  }
};

} // namespace NBody::Simulation

#endif // CPP_PROJECT_INCLUDE_NBODY_SIMULATION_COLLISIONRESOLVER_HPP_
