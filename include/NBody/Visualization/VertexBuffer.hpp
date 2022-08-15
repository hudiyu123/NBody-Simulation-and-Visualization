#ifndef CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_VERTEXBUFFER_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_VERTEXBUFFER_HPP_

#include <NBody/Visualization/GLLoader.hpp>

namespace NBody::Visualization {

/// \brief The abstraction of vertex buffer in OpenGL.
class VertexBuffer {
 private:
  GLuint id_;

 public:
  /// \brief Constructs a vertex buffer with data, size and usage. The vertex
  /// buffer will be automatically bound to the curren OpenGL context.
  /// \param data The pointer to the vertex data that will be copied into the
  /// data store for initialization.
  /// \param size he size in bytes of the vertex buffer
  /// \param usage The expected usage pattern of the data store, i.e.,
  /// GL_STATIC_DRAW
  VertexBuffer(const void* data, GLsizeiptr size, GLenum usage);

  /// \brief Destructs the vertex buffer.
  ~VertexBuffer();

  /// \brief Binds the vertex buffer to current OpenGL context.
  void bind() const;

  /// \brief Unbinds the current vertex buffer on the OpenGL context.
  static void unbind();
};

} // namespace NBody::Visualization

#endif // CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_VERTEXBUFFER_HPP_
