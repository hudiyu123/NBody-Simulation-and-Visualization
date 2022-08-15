#include <NBody/Visualization/VertexArray.hpp>

namespace NBody::Visualization {

VertexArray::VertexArray() : id_{} {
  glGenVertexArrays(1, &id_);
}

VertexArray::~VertexArray() {
  glDeleteVertexArrays(1, &id_);
}

void VertexArray::addVertexBuffer(
    const VertexBuffer& vb,
    const VertexBufferLayout& layout
) const {
  bind();
  vb.bind();
  const auto& elements = layout.getElements();
  std::size_t offset = 0;
  for (std::size_t i = 0; i < elements.size(); ++i) {
    const auto& element = elements[i];
    glEnableVertexAttribArray(i);
    glVertexAttribPointer(
        i,
        element.count,
        element.type,
        element.normalized,
        layout.getStride(),
        (const void*) offset
    );
    offset += element.size;
  }
}

void VertexArray::bind() const {
  glBindVertexArray(id_);
}

void VertexArray::unbind() {
  glBindVertexArray(0);
}

} // namespace NBody::Visualization