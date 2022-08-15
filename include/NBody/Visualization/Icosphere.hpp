#ifndef CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_ICOSPHERE_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_ICOSPHERE_HPP_

#include <NBody/Math/Vector.hpp>
#include <cmath>
#include <numbers>
#include <vector>

namespace NBody::Visualization {

/// \brief The data class stores the information of an icosphere.
/// \tparam ValueT Floating-point scalar type
/// \tparam IndexT The type of index
/// \tparam SizeT The type of size
/// \tparam CountT The type of count
template<typename ValueT, typename IndexT, typename SizeT, typename CountT>
class Icosphere {
 public:
  using ValueType = ValueT;
  using IndexType = IndexT;
  using SizeType = SizeT;
  using CountType = CountT;
  using VertexType = Math::Vector<ValueType, 3>;
  using TriangleIndexType = Math::Vector<IndexType, 3>;

 private:
  const ValueType radius_;
  const std::uint64_t subdivision_;
  std::vector<VertexType> vertices_;
  std::vector<TriangleIndexType> triangleIndices_;

 public:
  /// \brief Constructs an icosphere with radius and the level of subdivision.
  /// \param radius The radius
  /// \param subdivision The level of subdivision
  Icosphere(ValueType radius, std::uint64_t subdivision)
      : radius_{radius},
        subdivision_{subdivision},
        vertices_{},
        triangleIndices_{} {
    buildVertices();
  }

  /// \brief Returns the raw vertices.
  /// \details Vertices are flattened. For example, [[x1, y1, z1], [x2, y2, z2]]
  /// will be flattened to [x1, y1, z1, x2, y2, z2] and returned.
  /// \return The raw vertices
  std::vector<ValueType> getRawVertices() const {
    auto rawVertices = std::vector<ValueType>{};
    for (const auto& vertex: vertices_) {
      rawVertices.emplace_back(vertex[0]);
      rawVertices.emplace_back(vertex[1]);
      rawVertices.emplace_back(vertex[2]);
    }
    return rawVertices;
  }

  /// \brief Returns the count of raw vertices. Equals getRawVertices().size().
  /// \return The count of raw vertices
  [[nodiscard]] CountType getRawVerticesCount() const {
    return vertices_.size() * 3;
  }

  /// \brief Returns size of raw vertices in bytes.
  /// \return The size of raw vertices.
  [[nodiscard]] SizeType getRawVerticesSize() const {
    return getRawVerticesCount() * sizeof(ValueType);
  }

  /// \brief Returns the vertex indices.
  /// \return The vertex indices
  std::vector<IndexType> getVertexIndices() const {
    auto rawIndices = std::vector<IndexType>{};
    for (const auto& index: triangleIndices_) {
      rawIndices.emplace_back(index[0]);
      rawIndices.emplace_back(index[1]);
      rawIndices.emplace_back(index[2]);
    }
    return rawIndices;
  }

  /// \brief Returns the count of vertex indices.
  /// \return The count of vertex indices
  [[nodiscard]] CountType getVertexIndicesCount() const {
    return triangleIndices_.size() * 3;
  }

  /// \brief Returns the size of vertex indices in bytes.
  /// \return The size of vertex indices
  [[nodiscard]] SizeType getVertexIndicesSize() const {
    return getVertexIndicesCount() * sizeof(IndexType);
  }

 private:
  static ValueType computeScaleForLength(
      const VertexType& v,
      ValueType length
  ) {
    return length / Math::length(v);
  }

  static VertexType computeHalfVertex(
      const VertexType& v1,
      const VertexType& v2,
      ValueType radius
  ) {
    auto newVertex = VertexType{v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]};
    const auto scale = computeScaleForLength(newVertex, radius);
    newVertex *= scale;
    return newVertex;
  }

  std::vector<VertexType> computeIcosahedronVertices() const {
    const auto pi = std::numbers::pi_v<ValueType>;
    const auto hAngle = pi / ValueType{180} * ValueType{72};
    const auto vAngle = std::atan(ValueType{1} / ValueType{2});
    const auto z = radius_ * std::sin(vAngle);
    const auto xy = radius_ * std::cos(vAngle);

    auto vertices = std::vector<VertexType>{};
    vertices.emplace_back(std::initializer_list<ValueType>{
        ValueType{0},
        ValueType{0},
        radius_
    });

    auto hAngle1 = -pi / ValueType{2} - hAngle / ValueType{2};
    for (std::size_t i = 0; i < 5; ++i) {
      vertices.emplace_back(std::initializer_list<ValueType>{
          xy * std::cos(hAngle1),
          xy * std::sin(hAngle1),
          z
      });
      hAngle1 += hAngle;
    }

    auto hAngle2 = -pi / ValueType{2};
    for (std::size_t i = 0; i < 5; ++i) {
      vertices.emplace_back(std::initializer_list<ValueType>{
          xy * std::cos(hAngle2),
          xy * std::sin(hAngle2),
          -z
      });
      hAngle2 += hAngle;
    }

    vertices.emplace_back(std::initializer_list<ValueType>{
        ValueType{0},
        ValueType{0},
        -radius_
    });

    return vertices;
  }

  void buildVertices() {
    vertices_ = computeIcosahedronVertices();
    triangleIndices_.clear();
    for (IndexType i = 1; i <= 5; ++i) {
      const auto index1 = i == 5 ? 1 : i + 1;
      const auto index2 = i == 5 ? 6 : i + 6;
      triangleIndices_.emplace_back(
          std::initializer_list<IndexT>{0, i, index1}
      );
      triangleIndices_.emplace_back(
          std::initializer_list<IndexT>{index1, i, i + 5}
      );
      triangleIndices_.emplace_back(
          std::initializer_list<IndexT>{i + 5, index2, index1}
      );
      triangleIndices_.emplace_back(
          std::initializer_list<IndexT>{index2, i + 5, 11}
      );
    }

    for (std::size_t i = 0; i < subdivision_; ++i) {
      auto tempIndices = std::vector<TriangleIndexType>{};
      for (const auto& triangleIndex: triangleIndices_) {
        const auto vertexA = vertices_[triangleIndex[0]];
        const auto vertexB = vertices_[triangleIndex[1]];
        const auto vertexC = vertices_[triangleIndex[2]];
        const auto indexAB = static_cast<IndexT>(vertices_.size());
        const auto indexBC = indexAB + 1;
        const auto indexCA = indexAB + 2;
        vertices_.emplace_back(computeHalfVertex(vertexA, vertexB, radius_));
        vertices_.emplace_back(computeHalfVertex(vertexB, vertexC, radius_));
        vertices_.emplace_back(computeHalfVertex(vertexC, vertexA, radius_));
        tempIndices.emplace_back(
            std::initializer_list<IndexT>{triangleIndex[0], indexAB, indexCA}
        );
        tempIndices.emplace_back(
            std::initializer_list<IndexT>{indexAB, triangleIndex[1], indexBC}
        );
        tempIndices.emplace_back(
            std::initializer_list<IndexT>{indexCA, indexAB, indexBC}
        );
        tempIndices.emplace_back(
            std::initializer_list<IndexT>{indexCA, indexBC, triangleIndex[2]}
        );
      }
      triangleIndices_ = std::move(tempIndices);
    }
  };
};

} // namespace NBody::Visualization

#endif // CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_ICOSPHERE_HPP_
