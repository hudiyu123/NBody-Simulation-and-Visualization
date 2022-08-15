#include <NBody/Visualization/VertexBufferLayout.hpp>

namespace NBody::Visualization {

VertexBufferLayout::VertexBufferLayout() : elements_{}, stride_{} {}

void VertexBufferLayout::pushGLfloat(GLint count) {
  auto size = count * sizeof(GLfloat);
  elements_.push_back({GL_FLOAT, count, size, GL_FALSE});
  stride_ += static_cast<GLsizei>(size);
}

[[nodiscard]] const std::vector<VertexBufferElement>&
VertexBufferLayout::getElements() const {
  return elements_;
}

[[nodiscard]] GLsizei VertexBufferLayout::getStride() const {
  return stride_;
}

} // namespace NBody::Visualization