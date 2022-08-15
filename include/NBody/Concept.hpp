#ifndef CPP_PROJECT_INCLUDE_NBODY_CONCEPTS_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_CONCEPTS_HPP_

#include <concepts>

namespace NBody {

template<typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

}

#endif //CPP_PROJECT_INCLUDE_NBODY_CONCEPTS_HPP_
