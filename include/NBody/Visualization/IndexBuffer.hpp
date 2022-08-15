#ifndef CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_INDEXBUFFER_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_INDEXBUFFER_HPP_

#include <NBody/Visualization/GLLoader.hpp>

namespace NBody::Visualization {

/// \brief This class the abstraction of OpenGL index buffer.
class IndexBuffer {
 private:
  GLuint id_;
  GLsizei count_;
  GLenum type_;

 public:
  /// \brief Constructs an index buffer and binds it to current OpenGL context.
  /// \param data The pointer to the index data that will be copied into the
  /// data store for initialization.
  /// \param size The size in bytes of the index buffer
  /// \param count The count  of the inde
  /// x buffer
  /// \param type The type of the index, is usually GL_UNSIGNED_INT
  /// \param usage The expected usage pattern of the data store, i.e.,
  /// GL_STATIC_DRAW
  IndexBuffer(
      const void* data,
      GLsizeiptr size,
      GLsizei count,
      GLenum type,
      GLenum usage
  );

  /// \brief Destructs the index buffer. Automatically unbinds the index
  /// buffer from OpenGL context.
  ~IndexBuffer();

  /// \brief Binds the index buffer.
  void bind() const;

  /// \brief unbind the current index buffer on OpenGL context.
  static void unbind();

  /// \brief Returns the count of indices.
  /// \return The count of indices.
  [[nodiscard]] GLsizei getCount() const;

  /// \brief Returns the type of index.
  /// \return The type of index
  [[nodiscard]] GLenum getType() const;
};

} // namespace NBody::Visualization

#endif // CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_INDEXBUFFER_HPP_
