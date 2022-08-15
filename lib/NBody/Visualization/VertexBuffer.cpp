#include <NBody/Visualization/VertexBuffer.hpp>

namespace NBody::Visualization {

VertexBuffer::VertexBuffer(const void* data, GLsizeiptr size, GLenum usage)
    : id_{} {
  glGenBuffers(1, &id_);
  glBindBuffer(GL_ARRAY_BUFFER, id_);
  glBufferData(GL_ARRAY_BUFFER, size, data, usage);
}

VertexBuffer::~VertexBuffer() {
  glDeleteBuffers(1, &id_);
}

void VertexBuffer::bind() const {
  glBindBuffer(GL_ARRAY_BUFFER, id_);
}

void VertexBuffer::unbind() {
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

} // namespace NBody::Visualization