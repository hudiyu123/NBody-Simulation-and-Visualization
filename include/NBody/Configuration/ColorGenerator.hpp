#ifndef CPP_PROJECT_INCLUDE_NBODY_CONFIGURATION_COLORGENERATOR_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_CONFIGURATION_COLORGENERATOR_HPP_

#include <NBody/Math/Vector.hpp>
#include <concepts>
#include <random>

namespace NBody::Configuration {

/// \brief This class is used to generator random color for particles.
template<std::floating_point T>
class ColorGenerator {
 public:
  using ValueType = T;
  using VectorType = Math::Vector<ValueType, 3>;

 private:
  std::mt19937 engine_;
  std::uniform_real_distribution<ValueType> distribution_;

 public:
  /// \brief Construct a color generator with minimum and maximum value of RGB
  /// channel.
  /// \note Precondition: Min and max are in range [0, 1] and min is less than
  /// or equal to max.
  /// \param min Minimum value of RGB channels
  /// \param max Maximum value of RGB channels
  ColorGenerator(ValueType min, ValueType max)
      : engine_{std::random_device{}()},
        distribution_{min, max} {}

  /// \brief Returns a random RGB color.
  /// \return The RGB color
  VectorType operator()() {
    return {
        distribution_(engine_),
        distribution_(engine_),
        distribution_(engine_)
    };
  }
};

} // namespace NBody::Configuration

#endif //CPP_PROJECT_INCLUDE_NBODY_CONFIGURATION_COLORGENERATOR_HPP_
