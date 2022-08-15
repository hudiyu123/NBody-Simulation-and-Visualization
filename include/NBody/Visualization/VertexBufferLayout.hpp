#ifndef CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_VERTEXBUFFERLAYOUT_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_VERTEXBUFFERLAYOUT_HPP_

#include <NBody/Visualization/GLLoader.hpp>
#include <vector>

namespace NBody::Visualization {

/// \brief The abstraction of the vertex buffer element.
struct VertexBufferElement {
  /// \brief The type of the element.
  GLenum type;
  /// \brief The count of the element.
  GLint count;
  /// \brief The size in bytes of the element.
  std::size_t size;
  /// \brief Determines whether fixed-point data values should be normalized.
  GLboolean normalized;
};

class VertexBufferLayout {
 private:
  std::vector<VertexBufferElement> elements_;
  GLsizei stride_;

 public:
  /// \brief Constructs an empty vertex buffer layout.
  VertexBufferLayout();

  /// \brief Pushes a float type element into the layout with count.
  /// \param count The count of the element
  void pushGLfloat(GLint count);

  /// \brief Returns the elements in the layout.
  /// \return The elements
  [[nodiscard]] const std::vector<VertexBufferElement>& getElements() const;

  /// \brief Returns the stride of this layout
  /// \return The stride
  [[nodiscard]] GLsizei getStride() const;
};

} // namespace NBody::Visualization

#endif // CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_VERTEXBUFFERLAYOUT_HPP_
