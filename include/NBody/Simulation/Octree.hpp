#ifndef CPP_PROJECT_INCLUDE_NBODY_SIMULATION_OCTREE_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_SIMULATION_OCTREE_HPP_

#include <NBody/Simulation/Particle.hpp>
#include <NBody/Simulation/BoundaryBox.hpp>
#include <NBody/Math/Vector.hpp>
#include <cmath>
#include <concepts>
#include <numbers>
#include <set>
#include <vector>

namespace NBody::Simulation {

/// \brief The base data class of the octree node.
/// \details Each octree node represents a cube space in simulation space. This
/// base class stores the width and the center of a cube space.
template<std::floating_point T>
struct BaseOctreeNode {
  using ValueType = T;
  using VectorType = Math::Vector<ValueType, 3>;

  /// \brief The width of the cube space
  ValueType width;
  /// \brief The center of the cube space
  VectorType center;

  virtual ~BaseOctreeNode() = default;

  /// \brief Resets all data members.
  virtual void reset() = 0;
};

/// \brief Resets the width and the center。
template<std::floating_point T>
void BaseOctreeNode<T>::reset() {
  width = {};
  center = VectorType{};
}

/// \brief The branch node of the octree.
/// \details One branch node contains more than one particle. It has eight
/// octants (sub-nodes) and each of them could be a branch node, a leaf node or
/// an empty node.
template<std::floating_point T>
struct BranchOctreeNode : public BaseOctreeNode<T> {
 public:
  using ValueType = T;
  using VectorType = Math::Vector<ValueType, 3>;
  using BaseNodeType = BaseOctreeNode<ValueType>;

  /// \brief The total mass of all particles in the branch node.
  ValueType mass;
  /// \brief The mass center of all particles in the branch node.
  VectorType massCenter;
  /// \brief The maximum radius of all particles in the branch node.
  ValueType maxRadius;
  /// \brief The array of eight pointers point to octant.
  BaseNodeType* octants[8];

  /// \brief Constructs a branch octree node with default values.
  BranchOctreeNode() : mass{}, massCenter{}, maxRadius{}, octants{} {
    clearOctants();
  }

  /// \brief Resets all data members.
  void reset() override {
    BaseNodeType::reset();
    mass = {};
    massCenter = VectorType{};
    maxRadius = {};
    clearOctants();
  }

 private:
  void clearOctants() {
    for (auto& value: octants) {
      value = nullptr;
    }
  }
};

/// \brief The leaf node of the octree.
/// \details One leaf node only contains one particle and stores the index of
/// that particle.
template<std::floating_point T>
struct LeafOctreeNode : public BaseOctreeNode<T> {
  using BaseNodeType = BaseOctreeNode<T>;

  /// \brief The index of the particle in the leaf node.
  std::size_t particleIndex;

  /// \brief Resets the particle index.
  void reset() override {
    BaseNodeType::reset();
    particleIndex = {};
  }
};

/// \brief The concept of octree node.
/// \tparam T The type of the node
template<typename T>
concept OctreeNode =
std::is_base_of_v<BaseOctreeNode<typename T::ValueType>, T>;

/// \brief A pool stores the octree nodes to be reused.
/// \details This class reduce memory allocation, construction and destruction
/// of octree nodes when reconstructing octrees at a high rate.
/// \tparam T The octree node type
template<OctreeNode T>
class OctreeNodePool {
 public:
  using OctreeNodeType = T;
  using Pointer = OctreeNodeType*;

 private:
  std::vector<Pointer> pool_;

 public:
  /// \brief Constructs an empty node pool.
  OctreeNodePool() : pool_{} {}

  /// \brief Destructs the node pool, free all nodes in pool.
  ~OctreeNodePool() {
    clear();
  }

  /// \brief Pushes the node into the node pool.
  /// \param nodePtr The pointer to the node to be pushed
  void push(Pointer nodePtr) {
    pool_.emplace_back(nodePtr);
  }

  /// \brief Pops and returns a node.
  /// \details The node pool will update a node on heap if it is empty.
  /// \return The pointer to the popped node.
  T* pop() {
    if (pool_.empty()) {
      return new T();
    } else {
      auto nodePtr = pool_.back();
      pool_.pop_back();
      nodePtr->reset();
      return nodePtr;
    }
  }

  /// \brief Pops and frees all nodes in the node pool.
  void clear() {
    while (!pool_.empty()) {
      auto ptr = pool_.back();
      delete ptr;
      pool_.pop_back();
    }
  }
};

/// \brief The approximate particle of an octree node.
/// \details Treats all particles within a node as a single node.
/// \tparam T The floating-point scalar type
template<std::floating_point T>
struct ApproximateParticle {
  using ValueType = T;
  using VectorType = Math::Vector<ValueType, 3>;

  /// \brief The mass of the approximate particle.
  /// \details The sum of all particles' mass within the node.
  ValueType mass;
  /// \brief The position of the approximate particle.
  /// \details The center of mass of the node.
  VectorType position;

  ApproximateParticle(ValueType mass_, const VectorType& position_)
      : mass{mass_}, position{position_} {}
};

/// \brief Octree container used in Barnes-Hut gravity solver and Octree
/// collision detector.
/// \note This class is not well implemented currently and some common functions
/// are missed (i.e., iterator).
template<std::floating_point T>
class Octree {
 public:
  using ValueType = T;
  using VecterType = Math::Vector<ValueType, 3>;
  using BoundaryBoxType = BoundaryBox<ValueType>;
  using ParticleType = Particle<ValueType>;
  using BaseNodeType = BaseOctreeNode<ValueType>;
  using BranchNodeType = BranchOctreeNode<ValueType>;
  using LeafNodeType = LeafOctreeNode<ValueType>;
  using ApproximateParticleType = ApproximateParticle<ValueType>;

 private:
  const BoundaryBoxType boundaryBox_;
  OctreeNodePool<BranchNodeType> branchNodePool_;
  OctreeNodePool<LeafNodeType> leafNodePool_;
  BaseNodeType* rootNode_;

 public:
  /// \brief Constructs an octree with boundary box.
  /// \details The boundary box is used to construct the root node.
  /// \param boundaryBox The boundary box
  explicit Octree(const BoundaryBoxType& boundaryBox)
      : boundaryBox_{boundaryBox},
        branchNodePool_{},
        leafNodePool_{},
        rootNode_{} {}

  /// \brief The copy constructor of octree is deleted.
  Octree(const Octree&) = delete;
  /// \brief The copy assignment operator of octree is deleted.
  Octree& operator==(const Octree&) = delete;
  /// \brief The move constructor of octree is deleted.
  Octree(Octree&&) = delete;
  /// \brief The move assignment operator of octree is deleted.
  Octree& operator==(Octree&&) = delete;

  /// \brief Destructs the octree.
  ~Octree() {
    clear();
  }

  /// \brief Updates the octree with particles and active particle indices.
  /// \param particles The particles
  /// \param activeParticleIndices The indices of active particles
  void update(
      const std::vector<ParticleType>& particles,
      const std::vector<std::size_t>& activeParticleIndices
  ) {
    // Assume all particles are unsettled at beginning.
    auto unsettledParticleIndices = std::set<std::size_t>{
        activeParticleIndices.cbegin(),
        activeParticleIndices.cend()
    };
    updateNodeStructure(
        rootNode_,
        particles,
        activeParticleIndices,
        unsettledParticleIndices
    );
    // Insert unsettled particles into octree.
    for (auto index: unsettledParticleIndices) {
      insert(index, particles);
    }
    updateNodeData(rootNode_, particles);
  }

  /// \brief Clears the octree.
  void clear() {
    clear(rootNode_);
  }

  /// \brief Returns the approximate particles of the given particle.
  /// \param particleIndex The index of particle
  /// \param particles The particles
  /// \param openingAngle The critical opening angle
  /// \return The approximate particles of the particle with give index.
  std::vector<ApproximateParticleType> getApproximateParticles(
      const std::size_t particleIndex,
      const std::vector<ParticleType>& particles,
      ValueType openingAngle
  ) const {
    auto approximateParticles = std::vector<ApproximateParticleType>{};
    getApproximateParticles(
        particleIndex,
        particles,
        openingAngle,
        rootNode_,
        approximateParticles
    );
    return approximateParticles;
  }

  /// \brief Returns the indices of particles that are collided with the given
  /// particle.
  /// \param particleIndex The index of particle
  /// \param particles The particles
  /// \return The indices of collided particles
  std::vector<std::size_t> getCollidedParticleIndices(
      const std::size_t particleIndex,
      const std::vector<ParticleType>& particles
  ) const {
    auto collidedParticleIndices = std::vector<std::size_t>{};
    computeCollidedParticleIndices(
        particleIndex,
        particles,
        rootNode_,
        collidedParticleIndices
    );
    return collidedParticleIndices;
  }

 private:
  static void clear(BaseNodeType* node) {
    // If the node is a branch node, recursively clear its octants first.
    if (auto branchNode = dynamic_cast<BranchNodeType*>(node)) {
      for (auto octant: branchNode->octants) {
        clear(octant);
      }
    }
    // Then, free the node.
    delete node;
  }

  /// \brief Inserts a particle into the octree.
  /// \param particleIndex The index of the particle
  /// \param particles The particles
  void insert(
      const std::size_t particleIndex,
      const std::vector<ParticleType>& particles
  ) {
    const auto& particle = particles[particleIndex];
    if (!rootNode_) {
      // Create the root node as a leaf node if the root node is not exist.
      auto leafNode = leafNodePool_.pop();
      leafNode->particleIndex = particleIndex;
      leafNode->width = boundaryBox_.maxLength();
      leafNode->center = boundaryBox_.center();
      rootNode_ = leafNode;
    } else {
      // Enlarge the root node until the particle to be inserted is inside the
      // root node.
      while (isOutsideOfNode(particle, rootNode_)) {
        auto octantIndex = getOctantIndex(particle, rootNode_);
        auto newRootNode = branchNodePool_.pop();
        newRootNode->width = rootNode_->width * ValueType{2};
        newRootNode->center = rootNode_->center;
        newRootNode->center[0] +=
            getOctantSign(octantIndex, 0) * rootNode_->width / ValueType{2};
        newRootNode->center[1] +=
            getOctantSign(octantIndex, 1) * rootNode_->width / ValueType{2};
        newRootNode->center[2] +=
            getOctantSign(octantIndex, 2) * rootNode_->width / ValueType{2};
        newRootNode->octants[7 - octantIndex] = rootNode_;
        rootNode_ = newRootNode;
      }
      insert(particleIndex, particles, rootNode_);
    }
  }

  /// \brief Inserts a particle into the given node.
  /// \param particleIndex The index of the particle
  /// \param particles The particles
  /// \param nodePtr The reference to the node pointer
  void insert(
      const std::size_t particleIndex,
      const std::vector<ParticleType>& particles,
      BaseNodeType*& nodePtr
  ) {
    auto& particle = particles[particleIndex];
    if (auto branchNode = dynamic_cast<BranchNodeType*>(nodePtr)) {
      // If the node is a branch node, get the octant index of the particle in
      // this branch node.
      auto octantIndex = getOctantIndex(particle, branchNode);
      auto& octant = branchNode->octants[octantIndex];
      if (octant) {
        // If the octant exist, insert the particle into this octant.
        insert(particleIndex, particles, octant);
      } else {
        // Otherwise, get a new leaf node from leaf node pool to hold the
        // particle.
        auto leafNode = leafNodePool_.pop();
        leafNode->particleIndex = particleIndex;
        updateOctant(branchNode, leafNode, octantIndex);
        // Insert this new leaf node into the branch node.
        branchNode->octants[octantIndex] = leafNode;
      }
    } else if (auto leafNode = dynamic_cast<LeafNodeType*>(nodePtr)) {
      // If the node is a leaf node, replace this leaf node by a branch node.
      auto newBranchNode = branchNodePool_.pop();
      newBranchNode->width = leafNode->width;
      newBranchNode->center = leafNode->center;
      nodePtr = newBranchNode;
      // Reinsert the particle in this leaf node to the new branch node.
      insert(leafNode->particleIndex, particles, nodePtr);
      // Insert the particle to be inserted to the new branch node too.
      insert(particleIndex, particles, nodePtr);
      // Push this leaf node to the leaf node pool for reuse.
      leafNodePool_.push(leafNode);
    }
  }

  /// \brief Updates the structure of given node.
  /// \param nodePtr The node pointer
  /// \param particles The particles
  /// \param activeParticleIndices The active particle indices
  /// \param unsettledParticleIndices The unsettled particle indices
  void updateNodeStructure(
      BaseNodeType*& nodePtr,
      const std::vector<ParticleType>& particles,
      const std::vector<std::size_t>& activeParticleIndices,
      std::set<std::size_t>& unsettledParticleIndices
  ) {
    if (auto branchNode = dynamic_cast<BranchNodeType*>(nodePtr)) {
      // Update the structure of a branch node recursively.
      for (std::size_t i = 0; i < 8; ++i) {
        updateNodeStructure(
            branchNode->octants[i],
            particles,
            activeParticleIndices,
            unsettledParticleIndices
        );
      }

      // Compute the number of active octant and get the index of the last
      // active octant.
      std::size_t numOctants = 0;
      std::size_t octantIndex = 0;
      for (std::size_t i = 0; i < 8; ++i) {
        if (branchNode->octants[i]) {
          ++numOctants;
          octantIndex = i;
        }
      }

      if (numOctants == 0) {
        // If this node has no active octant, remove it from the octree and push
        // it to the branch node pool for reuse.
        nodePtr = nullptr;
        branchNodePool_.push(branchNode);
      } else if (numOctants == 1) {
        auto& octant = branchNode->octants[octantIndex];
        // If this branch node only has one active octant and this octant is a
        // leaf node, replace this branch node by this leaf node and push the
        // branch node to branch node pool for reuse.
        if (auto leafOctant = dynamic_cast<LeafNodeType*>(octant)) {
          leafOctant->width = branchNode->width;
          leafOctant->center = branchNode->center;
          nodePtr = octant;
          branchNodePool_.push(branchNode);
        }
      }
    } else if (auto leafNode = dynamic_cast<LeafNodeType*>(nodePtr)) {
      if (!isParticleActive(leafNode->particleIndex, activeParticleIndices)) {
        // If the particle in this leaf node is no longer active (e.g., due to
        // collision), remove the leaf node from the octree and push the leaf
        // node to leaf node poll for reuse.
        nodePtr = nullptr;
        leafNodePool_.push(leafNode);
      } else {
        // Otherwise, the particle in this leaf node is still active, further
        // determines whether it is outside the node.

        // The particle is going to be handled, remove it from unsettled
        // particle indices vector.
        unsettledParticleIndices.erase(leafNode->particleIndex);
        // If the particle in this leaf node is outside, reinsert the particle
        // and push the leaf node into leaf node pool for reuse.
        if (isOutsideOfNode(particles[leafNode->particleIndex], leafNode)) {
          nodePtr = nullptr;
          insert(leafNode->particleIndex, particles);
          leafNodePool_.push(leafNode);
        }
      }
    }
  }

  void updateNodeData(
      BaseNodeType* nodePtr,
      const std::vector<ParticleType>& particles
  ) {
    if (auto branchNode = dynamic_cast<BranchNodeType*>(nodePtr)) {
      branchNode->mass = {};
      branchNode->massCenter = VecterType{};
      branchNode->maxRadius = {};
      for (auto octant: branchNode->octants) {
        // Recursively update the octants.
        updateNodeData(octant, particles);
        if (auto branchOctant = dynamic_cast<BranchNodeType*>(octant)) {
          branchNode->mass += branchOctant->mass;
          branchNode->massCenter +=
              branchOctant->massCenter * branchOctant->mass;
          branchNode->maxRadius = std::max(
              branchNode->maxRadius,
              branchOctant->maxRadius
          );
        } else if (auto leafOctant = dynamic_cast<LeafNodeType*>(octant)) {
          const auto& particle = particles[leafOctant->particleIndex];
          branchNode->mass += particle.mass;
          branchNode->massCenter += particle.position * particle.mass;
          branchNode->maxRadius = std::max(
              branchNode->maxRadius,
              particle.radius
          );
        }
      }
      if (branchNode->mass > ValueType{0}) {
        branchNode->massCenter /= branchNode->mass;
      }
    }
  }

  void getApproximateParticles(
      const std::size_t particleIndex,
      const std::vector<ParticleType>& particles,
      ValueType openingAngle,
      BaseNodeType* nodePtr,
      std::vector<ApproximateParticleType>& approximateParticles
  ) const {
    const auto& particle = particles[particleIndex];
    if (auto branchNode = dynamic_cast<BranchNodeType*>(nodePtr)) {
      auto difference = particle.position - branchNode->massCenter;
      if (nodePtr->width / Math::length(difference) < openingAngle) {
        // The particle is far enough to the node, treat the node as a single
        // approximate particle.
        approximateParticles.emplace_back(
            branchNode->mass,
            branchNode->massCenter
        );
      } else {
        // If the opening angle of the given particle is not less than the
        // critical opening angle, going into the node further.
        for (auto octant: branchNode->octants) {
          getApproximateParticles(
              particleIndex,
              particles,
              openingAngle,
              octant,
              approximateParticles
          );
        }
      }
    } else if (auto leafNode = dynamic_cast<LeafNodeType*>(nodePtr)) {
      if (particleIndex != leafNode->particleIndex) {
        const auto& otherParticle = particles[leafNode->particleIndex];
        approximateParticles.emplace_back(
            otherParticle.mass,
            otherParticle.position
        );
      }
    }
  }

  void computeCollidedParticleIndices(
      const std::size_t particleIndex,
      const std::vector<ParticleType>& particles,
      BaseNodeType* nodePtr,
      std::vector<std::size_t>& collidedParticleIndices
  ) const {
    const auto& particle = particles[particleIndex];
    if (auto branchNode = dynamic_cast<BranchNodeType*>(nodePtr)) {
      // The distance of the particle to the center of node.
      const auto distance =
          Math::length(particle.position - branchNode->center);
      // The distance of the center of node to the vertices of the node.
      const auto nodeRadius =
          branchNode->width * std::numbers::sqrt3_v<ValueType> / ValueType{2};
      const auto critical =
          particle.radius + branchNode->maxRadius + nodeRadius;
      if (distance <= critical) {
        for (auto octant: branchNode->octants) {
          computeCollidedParticleIndices(
              particleIndex,
              particles,
              octant,
              collidedParticleIndices
          );
        }
      }
    } else if (auto leafNode = dynamic_cast<LeafNodeType*>(nodePtr)) {
      if (particleIndex != leafNode->particleIndex) {
        const auto& otherParticle = particles[leafNode->particleIndex];
        if (particle.isCollidedWith(otherParticle)) {
          collidedParticleIndices.emplace_back(leafNode->particleIndex);
        }
      }
    }
  }

  /// \brief Returns the index of the octant that the given particle should be
  /// placed in the given node.
  /// \param particle The particle
  /// \param nodePtr The node pointer
  /// \return The octant index
  std::size_t getOctantIndex(
      const ParticleType& particle,
      const BaseNodeType* nodePtr
  ) const {
    std::size_t octantIndex = 0;
    if (particle.position[0] < nodePtr->center[0]) { octantIndex += 1; }
    if (particle.position[1] < nodePtr->center[1]) { octantIndex += 2; }
    if (particle.position[2] < nodePtr->center[2]) { octantIndex += 4; }
    return octantIndex;
  }

  /// \brief Determines whether the particle is outside the node.
  /// \param particle The particle
  /// \param nodePtr The node pointer
  /// \return true if the particle is outside the node, false otherwise
  bool isOutsideOfNode(
      const ParticleType& particle,
      const BaseNodeType* nodePtr
  ) const {
    const auto difference = particle.position - nodePtr->center;
    const auto nodeHalfWidth = nodePtr->width / ValueType{2};
    return std::fabs(difference[0]) > nodeHalfWidth
        || std::fabs(difference[1]) > nodeHalfWidth
        || std::fabs(difference[2]) > nodeHalfWidth;
  }

  /// \brief Returns the sign of the given octant index in the given dimension.
  /// \param octantIndex The octant index
  /// \param dimension The dimension. 0 for x, 1 for y, and 2 for z
  /// \return The sign of the octant. 1 for positive and -1 for negative
  ValueType getOctantSign(
      std::size_t octantIndex,
      std::size_t dimension
  ) const {
    return (octantIndex >> dimension) % 2 == 0 ? ValueType{1} : ValueType{-1};
  }

  /// \brief Updates the width and center of the newly created octant.
  /// \details This function is call after a new octant is created and inserted
  /// into a branch node.
  /// \param parentNodePtr The parent node of the given octant
  /// \param octantPtr The octant pointer
  /// \param octantIndex The octant index of the given octant
  void updateOctant(
      const BaseNodeType* parentNodePtr,
      BaseNodeType* octantPtr,
      const std::size_t octantIndex
  ) const {
    octantPtr->width = parentNodePtr->width / ValueType{2};
    auto octantHalfWidth = octantPtr->width / ValueType{2};
    octantPtr->center = parentNodePtr->center;
    octantPtr->center[0] += getOctantSign(octantIndex, 0) * octantHalfWidth;
    octantPtr->center[1] += getOctantSign(octantIndex, 1) * octantHalfWidth;
    octantPtr->center[2] += getOctantSign(octantIndex, 2) * octantHalfWidth;
  }

  /// \details Determines whether a particle is active.
  /// \note Precondition: The active particle indices vector is sorted.
  /// \param particleIndex The index of the particle
  /// \param activeParticleIndices The active particle indices
  /// \return true is the particle is active, false otherwise
  [[nodiscard]] bool isParticleActive(
      const std::size_t particleIndex,
      const std::vector<std::size_t>& activeParticleIndices
  ) const {
    return std::ranges::binary_search(
        activeParticleIndices,
        particleIndex
    );
  }
};

} // namespace NBody::Simulation

#endif // CPP_PROJECT_INCLUDE_NBODY_SIMULATION_OCTREE_HPP_
