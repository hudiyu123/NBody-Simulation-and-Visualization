#include <NBody/Visualization/IndexBuffer.hpp>

namespace NBody::Visualization {

IndexBuffer::IndexBuffer(
    const void* data,
    GLsizeiptr size,
    GLsizei count,
    GLenum type,
    GLenum usage
) : id_{}, count_{count}, type_{type} {
  glGenBuffers(1, &id_);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, id_);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, data, usage);
}

IndexBuffer::~IndexBuffer() {
  glDeleteBuffers(1, &id_);
}

void IndexBuffer::bind() const {
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, id_);
}

void IndexBuffer::unbind() {
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

GLsizei IndexBuffer::getCount() const {
  return count_;
}

GLenum IndexBuffer::getType() const {
  return type_;
}

} // namespace NBody::Visualization
