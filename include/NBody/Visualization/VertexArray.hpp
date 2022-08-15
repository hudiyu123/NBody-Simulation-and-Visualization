#ifndef CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_VERTEXARRAY_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_VERTEXARRAY_HPP_

#include <NBody/Visualization/GLLoader.hpp>
#include <NBody/Visualization/VertexBuffer.hpp>
#include <NBody/Visualization/VertexBufferLayout.hpp>

namespace NBody::Visualization {

/// \brief The abstraction of vertex array in OpenGL.
class VertexArray {
 private:
  GLuint id_;

 public:
  /// \brief Constructs a vertex array and binds it to current OpenGL context.
  VertexArray();

  /// \brief Destructs the vertex array.
  ~VertexArray();

  /// \brief Adds a vertex buffer to the vertex array with layout.
  /// \param vb The vertex buffer
  /// \param layout The layout of the vertex buffer
  void addVertexBuffer(
      const VertexBuffer& vb,
      const VertexBufferLayout& layout
  ) const;

  /// \brief Binds the vertex array to current OpenGL context.
  void bind() const;

  /// \brief Unbinds the vertex array on the current OpenGL context.
  static void unbind();
};

} // namespace NBody::Simulation

#endif // CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_VERTEXARRAY_HPP_
