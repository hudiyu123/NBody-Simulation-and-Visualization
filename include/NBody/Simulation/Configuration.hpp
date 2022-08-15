#ifndef CPP_PROJECT_INCLUDE_NBODY_SIMULATION_CONFIGURATION_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_SIMULATION_CONFIGURATION_HPP_

#include <concepts>
#include <iostream>

namespace NBody::Simulation {

/// \brief The type of gravity solver.
enum class GravitySolverType {
  /// \brief Uses Direct algorithm for interaction computation.
  /// \details The Direct algorithm involves computing the gravitational forces
  /// influenced by other particles for each particle. The complexity of this
  /// algorithm is Θ(N^2).
  /// \note Use this method only when the number of particles is small or the
  /// precision of simulation matters. Otherwise, uses Barnes-Hut algorithm
  /// instead to reduce the simulation time.
  Direct,
  /// \brief Uses Barnes-Hut algorithm for interaction computation.
  /// \details The Barnes-Hut algorithm is an approximate method. This algorithm
  /// reduces the time complexity in comparison with the  Direct algorithm.\n
  /// In Barnes–Hut algorithm, the 3D space is recursively subdivided into cubic
  /// nodes via an octree. Only particles from nearby cells need to be treated
  /// individually. Particles in distant cells can be treated as a single large
  /// particle centered at the cell's center of mass. The complexity of
  /// Barnes-Hut algorithm is Θ(NlogN).
  /// \note Uses the Direct algorithm instead if the precision of simulation
  /// matters.
  BarnesHut
};

/// \brief The type of integrator.
enum class IntegratorType {
  /// \brief Uses the Leapfrog method in computing integration.
  /// \details The Leapfrog integration algorithm involves three stages. In the
  /// first stage, the particle is advanced for half a time step while keeping
  /// the velocities fixed:\n
  /// \n
  /// x(i+1/2) = x(i) + v(i) * Δt / 2\n
  /// \n
  /// In the second stage, the velocity is advanced for one time step while
  /// keeping the position fixed:\n
  /// \n
  /// v(i+1) = v(i) + a(i+1) * Δt\n
  /// \n
  /// In the third stage, the velocity is again advanced for half a time step:\n
  /// \n
  /// x(i+1) = x(i+1/2) + v(i+1) * Δt / 2\n
  /// \n
  /// \note Currently, only Leapfrog method is supported for integration.
  Leapfrog
};

/// \brief The type of boundary condition handler.
enum class BoundaryConditionHandlerType {
  /// \brief Boundaries are ignored. The particles can have arbitrary
  /// coordinates in the space.
  None,
  /// \brief The particles will be removed if they cross the boundary.
  Open,
  /// \brief The particles are considered as colliding with the boundary box and
  /// rebound if they reach the boundary. Boundaries are considered as static
  /// solid bodies with infinite mass.
  Rebound,
  /// \brief The particles are reinserted on the opposite side of the boundary
  /// box if they cross the boundary.
  Periodic
};

/// \brief The type of collision detector.
enum class CollisionDetectorType {
  /// \brief The collisions of particles are ignored.
  /// \note This might cause strong force gradients on small scales. Use
  /// non-zero softening to reduce this effect.
  None,
  /// \brief Uses a brute force algorithm to detect the collision.
  /// \details This algorithm simply checks the overlaps between every particle
  /// pair at the end of the time step. The time complexity of this algorithm is
  /// Θ(N2).
  /// \note Use Direct algorithm only when the number of particles is small.
  Direct,
  /// \brief Uses an octree to detect the collision.
  /// \details To find overlapping particles, one starts at the root node and
  /// recursively searches into sub-nodes if the distance between the particle
  /// and the cell center is smaller than a critical value:\n
  /// \n
  /// r(ic) < R(i) + R(max) + R(c)\n
  /// \n
  /// where r(ic) is distance of particle i to the center of node c, R(i) is the
  /// radius of the particle i, R(max) is the maximum radius of the particles in
  /// node c, R(c) is the distance of the center of cell c to its vertices.
  /// \note The direct algorithm would be more efficient if the number of
  /// particles is small.
  Octree
};

/// \brief The type of collision resolver.
enum class CollisionResolverType {
  /// \brief This method assumes hard-sphere collision.
  /// \details This method conserves momentum, mass, and kinetic energy
  /// depending on the COR (Coefficient of Restitution). The COR must be in
  /// range [0, 1]. If COR equals 0, the collision is perfectly inelastic,
  /// in which two colliding particles have the same velocity after collision.
  Hard,
  /// \brief This method merges the two colliding particles.
  /// \details This method conserves momentum, mass, volume, but not kinetic
  /// energy. The particle has the lower index will be removed after collision.
  Merge
};

/// \brief The configuration class of simulator.
template<std::floating_point T>
struct Configuration {
  using ValueType = T;

  /// \brief The delta step. The default value is 0.05.
  /// \note The value of must be greater than or equal to 0.
  ValueType timeStep{0.05};
  /// \brief The nbody program will integrate exactly up to the request step.
  /// The default value is 0.
  std::uint64_t finishStep{0};
  /// \brief The gravitational constant. The default value is 1.
  ValueType G{1};
  /// \brief The gravitational softening parameter. The default value is 0.
  /// \details The gravitational force of a particle in the x direction is
  /// calculated as\n
  /// \n
  /// F(x) = -x * G * m(1) * m(2) / (x^2 + y^2 + z^2 + b^2)^(3/2)\n
  /// \n
  /// where b is the gravitational softening parameter. This can be used to
  /// avoid numerical divergence when a particle comes too close to another, in
  /// which the force goes to infinity.
  /// \note The value of this parameter must be greater than or equal to 0.
  ValueType softening{0};
  /// \brief Then opening angle parameter for the Barnes-Hut algorithm. The
  /// default value is 1.
  /// \note The value of this parameter must be greater than or equal to 0.
  /// \details This parameter determines the accuracy of the gravity calculation
  /// when the Barnes-Hut algorithm is used. Let θ = w/R denotes the opening
  /// angle of a particle, where w is the width of the target cell and R is the
  /// distance from the cell’s center of mass to the particle. If the opening
  /// angle is smaller than a critical angle θ(critical) > θ, the total mass and
  /// center of mass of the cell are used to calculate the contribution to the
  /// gravitational force. Otherwise, the sub-cells are opened until the
  /// criterion is met. One has to choose θ(critical) appropriately to achieve a
  /// balance between accuracy and speed.
  /// \note The simulation result of using Barnes-Hut algorithm will be
  /// identical to the result of using direct algorithm if the opening angle
  /// equals 0.
  ValueType openingAngle{1};
  /// \brief The coefficient of restitution (COR) for the hard collision
  /// resolver. The default value is 1.
  /// \note The COR must be in range [0, 1].
  /// \details If the COR equals 0, the collision is perfectly inelastic. If the
  /// COR equals 1, the collision is perfectly elastic.
  /// \note The computation might yield consecutive collision of two particles
  /// when COR is less than 1, which leads two particles keep getting close to
  /// each other due to the Leapfrog integration method. As a result, the two
  /// particles will bounce off in a high velocity. The parameter merge
  /// threshold is introduced to avoid this situation.
  ValueType COR{1};
  /// \brief The coefficient of restitution (COR) for the rebound boundary
  /// condition. The default value is 1.
  /// \note The boundary COR must be in range [0, 1].
  /// \details If the COR of boundary equals 0, the particle will be removed as
  /// same as open boundary condition.
  ValueType boundaryCOR{1};
  /// \brief The merge threshold parameter for handling consecutive collision
  /// when hard collision resolver is used. Consecutive collision handling will
  /// be disabled if the value of this parameter equals 0. The default value is
  /// 0.
  /// \details The two colliding particles will be merged if they consecutively
  /// collide with each other and the count of the consecutive collision reaches
  /// the threshold.
  std::uint64_t mergeThreshold{0};

  /// \brief The type of the gravity solver will be used in the simulation. The
  /// default value is Direct.
  GravitySolverType gravitySolverType{};
  /// \brief The type of the integrator will be used in the simulation. The
  /// default value is Leapfrog.
  IntegratorType integratorType{};
  /// \brief The type of the boundary condition handler will be used in the
  /// simulation. The default value is None.
  BoundaryConditionHandlerType boundaryConditionHandlerType{};
  /// \brief The type of the collision detector will be used in the simulation.
  /// The default value is None.
  CollisionDetectorType collisionDetectorType{};
  /// \brief The type of the collision resolver will be used in the simulation.
  /// The default value is Hard.
  CollisionResolverType collisionResolverType{};
};

/// \brief Reads the type of gravity solver from a character input stream.
/// \param istream The character input stream
/// \param type The type of gravity solver
/// \return The character input stream
std::istream& operator>>(std::istream& istream, GravitySolverType& type) {
  std::string typeStr;
  if (!(istream >> typeStr)) { return istream; }
  if (typeStr == "Direct") {
    type = GravitySolverType::Direct;
  } else if (typeStr == "BarnesHut") {
    type = GravitySolverType::BarnesHut;
  } else {
    istream.setstate(std::ios_base::failbit);
  }
  return istream;
}

/// \brief Writes the gravity solver type to a character output stream.
/// \param ostream The character output stream
/// \param type The gravity solver type
/// \return The character output stream
std::ostream& operator<<(std::ostream& ostream, const GravitySolverType& type) {
  switch (type) {
    case GravitySolverType::Direct:
      ostream << "Direct";
      break;
    case GravitySolverType::BarnesHut:
      ostream << "BarnesHut";
      break;
  }
  return ostream;
}

/// \brief Reads the type of integrator from a character input stream.
/// \param istream The character input stream
/// \param type The type of integrator
/// \return The character input stream
std::istream& operator>>(std::istream& istream, IntegratorType& type) {
  std::string typeStr;
  if (!(istream >> typeStr)) { return istream; }
  if (typeStr == "Leapfrog") {
    type = IntegratorType::Leapfrog;
  } else {
    istream.setstate(std::ios_base::failbit);
  }
  return istream;
}

/// \brief Writes the integrator type to a character output stream.
/// \param ostream The character output stream
/// \param type The integrator type
/// \return The character output stream
std::ostream& operator<<(std::ostream& ostream, const IntegratorType& type) {
  switch (type) {
    case IntegratorType::Leapfrog:
      ostream << "Leapfrog";
      break;
  }
  return ostream;
}

/// \brief Reads the type of boundary condition handler from a character input stream.
/// \param istream The character input stream
/// \param type The type of boundary condition handler
/// \return The character input stream
std::istream& operator>>(
    std::istream& istream,
    BoundaryConditionHandlerType& type
) {
  std::string typeStr;
  if (!(istream >> typeStr)) { return istream; }
  if (typeStr == "None") {
    type = BoundaryConditionHandlerType::None;
  } else if (typeStr == "Open") {
    type = BoundaryConditionHandlerType::Open;
  } else if (typeStr == "Rebound") {
    type = BoundaryConditionHandlerType::Rebound;
  } else if (typeStr == "Periodic") {
    type = BoundaryConditionHandlerType::Periodic;
  } else {
    istream.setstate(std::ios_base::failbit);
  }
  return istream;
}

/// \brief Writes the boundary condition handler type to a character output
/// stream.
/// \param ostream The character output stream
/// \param type The boundary condition handler type
/// \return The character output stream
std::ostream& operator<<(
    std::ostream& ostream,
    const BoundaryConditionHandlerType& type
) {
  switch (type) {
    case BoundaryConditionHandlerType::None:
      ostream << "None";
      break;
    case BoundaryConditionHandlerType::Open:
      ostream << "Open";
      break;
    case BoundaryConditionHandlerType::Rebound:
      ostream << "Rebound";
      break;
    case BoundaryConditionHandlerType::Periodic:
      ostream << "Periodic";
      break;
  }
  return ostream;
}

/// \brief Reads the type of collision detector from a character input stream.
/// \param istream The character input stream
/// \param type The type of collision detector
/// \return The character input stream
std::istream& operator>>(std::istream& istream, CollisionDetectorType& type) {
  std::string typeStr;
  if (!(istream >> typeStr)) { return istream; }
  if (typeStr == "None") {
    type = CollisionDetectorType::None;
  } else if (typeStr == "Direct") {
    type = CollisionDetectorType::Direct;
  } else if (typeStr == "Octree") {
    type = CollisionDetectorType::Octree;
  } else {
    istream.setstate(std::ios_base::failbit);
  }
  return istream;
}

/// \brief Writes the collision detector type to a character output stream.
/// \param ostream The character output stream
/// \param type The collision detector type
/// \return The character output stream
std::ostream& operator<<(
    std::ostream& ostream,
    const CollisionDetectorType& type
) {
  switch (type) {
    case CollisionDetectorType::None:
      ostream << "None";
      break;
    case CollisionDetectorType::Direct:
      ostream << "Direct";
      break;
    case CollisionDetectorType::Octree:
      ostream << "Octree";
      break;
  }
  return ostream;
}

/// \brief Reads the type of collision resolver from a character input stream.
/// \param istream The character input stream
/// \param type The type of collision resolver
/// \return The character input stream
std::istream& operator>>(std::istream& istream, CollisionResolverType& type) {
  std::string typeStr;
  if (!(istream >> typeStr)) { return istream; }
  if (typeStr == "Merge") {
    type = CollisionResolverType::Merge;
  } else if (typeStr == "Hard") {
    type = CollisionResolverType::Hard;
  } else {
    istream.setstate(std::ios_base::failbit);
  }
  return istream;
}

/// \brief Writes the collision resolver type to a character output stream.
/// \param ostream The character output stream
/// \param type The collision resolver type
/// \return The character output stream
std::ostream& operator<<(
    std::ostream& ostream,
    const CollisionResolverType& type
) {
  switch (type) {
    case CollisionResolverType::Merge:
      ostream << "Merge";
      break;
    case CollisionResolverType::Hard:
      ostream << "Hard";
      break;
  }
  return ostream;
}

} // namespace NBody::Simulation

#endif // CPP_PROJECT_INCLUDE_NBODY_SIMULATION_CONFIGURATION_HPP_
